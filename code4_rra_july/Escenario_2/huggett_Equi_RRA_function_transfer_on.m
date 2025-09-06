function [r_opt, ir, pop1_vector, a, g_opt, c_opt] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vector, sI_vector1, sF_vector1, cfg)
% HUGGETT_EQUI_RRA_FUNCTION_TRANSFER_ON
%  Versión con soporte para COVID shock y transfers que aumentan:
%   Transfer_covid = (transfer_multiplier) * phi * z1_base
%
%  Uso:
%    [r_opt, ir, pop, a, g_opt, c_opt] = huggett_Equi_RRA_function_transfer_on(eta, sI, sF)
%    [ ... ] = huggett_Equi_RRA_function_transfer_on(eta, sI, sF, cfg)
%
%  cfg (struct) principales:
%    scenario             : "baseline" | "covid_uptransfer"
%    psi_I, psi_F         : multiplicadores de ingreso (default 1,1)
%    transfer_multiplier  : κ>=1 (default 1)  -> solo se usa si scenario ~= baseline
%    keep_transfers_level : true (mantiene nivel basado en z1_base)
%    amin_mode            : "baseline" | "shocked" (cómo se fija a_min)
%    gov_nonneg           : true -> trunca bien público a >= 0
%    phi                  : 0.13 (share sobre z1_base)
%    taxF, taxI, tauc     : tasas fiscales (tauc aún no se usa en el presupuesto)
%    theta                : 0.02 (prima endeudamiento informal)
%    r0, rmin, rmax       : bracket para r*
%
%  Retorna por agente jj:
%    r_opt(jj), g_opt{jj} (I x 2), c_opt{jj} (I x 2), malla a (I x 1)

% -----------------------------
% Defaults de configuración
% -----------------------------
if nargin < 4 || isempty(cfg), cfg = struct(); end
def = struct( ...
    'scenario',            "baseline", ...
    'psi_I',               1.0, ...
    'psi_F',               1.0, ...
    'transfer_multiplier', 1.0, ...
    'keep_transfers_level', true, ...
    'amin_mode',           "baseline", ...
    'gov_nonneg',          true, ...
    'phi',                 0.13, ...
    'taxF',                0.10, ...
    'taxI',                0.00, ...
    'tauc',                0.00, ...
    'theta',               0.02, ...
    'r0',                  0.03, ...
    'rmin',                0.01, ...
    'rmax',                0.04 ...
);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end

% -----------------------------
% Parámetros del modelo
% -----------------------------
rho     = 0.05;      % descuento
z1_base = 0.33;      % ingreso informal base
z2_base = 1.00;      % ingreso formal base
I       = 700;       % puntos de malla
amax    = 5.0;

% Transiciones ocupacionales
p22 = 0.75;
la2 = -log(p22);

% Salidas
pop1_vector = eta_vector;
J  = numel(pop1_vector);
r_opt = zeros(1, J);
g_opt = cell(1, J);
c_opt = cell(1, J);

% -----------------------------
% Bucle en agentes jj
% -----------------------------
for jj = 1:J
    pop1 = pop1_vector(jj);

    % Intensidades ajustadas por tamaño del sector
    la1 = (1 - pop1) * la2 / pop1;   % informal -> formal
    p11 = exp(-la1); %#ok<NASGU>
    la  = [la1, la2];

    taxF = cfg.taxF; taxI = cfg.taxI; tauc = cfg.tauc; %#ok<NASGU>
    phi  = cfg.phi;  theta = cfg.theta;

    % Ingresos con/ sin shock
    if cfg.scenario == "covid_uptransfer"
        z1 = cfg.psi_I * z1_base;
        z2 = cfg.psi_F * z2_base;
        kappa = cfg.transfer_multiplier;   % κ > 1
    else
        z1 = z1_base; z2 = z2_base;
        kappa = 1.0;
    end
    z  = [z1, z2];

    % Transfer: NIVEL basado en z1_base, multiplicado por κ si hay shock
    if cfg.keep_transfers_level
        Transfer = kappa * phi * z1_base;   % contracíclico y reforzado
    else
        Transfer = kappa * phi * z1;        % proporcional al ingreso actual
    end

    % Malla de activos
    if cfg.amin_mode == "shocked"
        amin = -0.3 * z1;        % límite cae con el ingreso actual
    else
        amin = -0.3 * z1_base;   % límite del benchmark (comparativa pura)
    end
    a  = linspace(amin, amax, I)'; 
    da = (amax - amin) / (I - 1);
    aa = [a, a];
    zz = ones(I,1) * [z1, z2];

    % Parámetros de solución
    r     = cfg.r0;
    rmin  = cfg.rmin;  rmax = cfg.rmax;
    maxit = 100;       crit = 1e-6;
    Delta = 1000;
    Ir    = 100;       crit_S = 1e-5;

    % RRA
    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    % Matriz de transición ocupacional
    Aswitch = [-speye(I)*la(1),  speye(I)*la(1);
                speye(I)*la(2), -speye(I)*la(2)];

    % Valor inicial V0
    rr = zeros(I,1);
    for i = 1:I
        if a(i) >= 0, rr(i) = r; else, rr(i) = r + theta; end
    end
    V0 = zeros(I,2);
    V0(:,1) = (((1 - taxI) * zz(:,1) + rr .* a + Transfer).^(1 - sI) / (1 - sI)) / rho;
    V0(:,2) = (((1 - taxF) * zz(:,2) + r  .* a).^(1 - sF)   / (1 - sF)) / rho;

    % ===== Encontrar r* por bisección =====
    for ir = 1:Ir
        if ir > 1, V0 = V_r(:,:,ir-1); end
        V = V0;

        dVf = zeros(I,2); dVb = zeros(I,2);

        for n = 1:maxit
            % Derivadas forward/backward
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVf(I,1) = ((1 - taxI)*z(1) + r*amax)^(-sI);
            dVf(I,2) = ((1 - taxF)*z(2) + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVb(1,1) = ((1 - taxI)*z(1) + (r+theta)*amin)^(-sI);
            dVb(1,2) = ((1 - taxF)*z(2) + r*amin)^(-sF);

            % rr(a) para informales
            for i = 1:I
                if a(i) >= 0, rr(i) = r; else, rr(i) = r + theta; end
            end

            % Consumo/ahorro forward & backward
            cf  = [max(dVf(:,1),1e-10).^(-1/sI), max(dVf(:,2),1e-10).^(-1/sF)];
            ssf = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cf;

            cb  = [max(dVb(:,1),1e-10).^(-1/sI), max(dVb(:,2),1e-10).^(-1/sF)];
            ssb = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cb;

            % Consumo en punto fijo (drift ~ 0)
            c0  = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa;

            % Upwind
            If = ssf > 0;
            Ib = ssb < 0;
            I0 = (1 - If - Ib);

            c = cf.*If + cb.*Ib + c0.*I0;

            % Utilidad (informal incluye Transfer fijo reforzado)
            U1 = ((c(:,1) + Transfer).^(1 - sI)) / (1 - sI);
            U2 =  (c(:,2).^(1 - sF)) / (1 - sF);

            % Presupuesto del Gobierno (bien público + transfer)
            popI = la2 / (la1 + la2);
            popF = la1 / (la1 + la2);
            TaxRevenue = taxI*z(1)*popI + taxF*z(2)*popF;  % + tauc*E[c] si implementamos
            Gov_pub    = TaxRevenue - Transfer * popI;
            if cfg.gov_nonneg, Gov_pub = max(0, Gov_pub); end
            Gov = Gov_pub + Transfer;

            u = [U1, U2] + Gov * ones(I,2);

            % Matriz A para FP
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z =  max(ssf,0)/da;

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            if max(abs(sum(A,2))) > 1e-9
                error('Improper Transition Matrix');
            end

            % Resolver HJB implícito
            B = (1/Delta + rho)*speye(2*I) - A;
            u_st = [u(:,1); u(:,2)];
            V_st = [V(:,1); V(:,2)];
            b = u_st + V_st/Delta;
            V_st = B \ b;
            V_new = [V_st(1:I), V_st(I+1:2*I)];

            if max(max(abs(V_new - V))) < crit
                V = V_new; 
                break
            end
            V = V_new;
        end % end HJB

        % ===== Fokker-Planck =====
        AT = A';
        bb = zeros(2*I,1);
        i_fix = 1;  bb(i_fix) = .1;
        row = zeros(1, 2*I); row(i_fix) = 1;
        AT(i_fix,:) = row;

        gg = AT \ bb;
        g_sum = gg' * ones(2*I,1) * da;
        gg = gg ./ g_sum;

        g = [gg(1:I), gg(I+1:2*I)];

        % Clearing: ahorro agregado = 0
        S = g(:,1)'*a*da + g(:,2)'*a*da;

        % Actualiza r (bisección)
        if S >  crit_S
            rmax = r;  r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r;  r = 0.5*(r + rmax);
        else
            break
        end

        % Guarda para la próxima iteración de r
        V_r(:,:,ir)  = V;     %#ok<AGROW>
        g_r(:,:,ir)  = g;     %#ok<AGROW>
    end % end loop r*

    % Salidas por agente
    r_opt(jj) = r;
    g_opt{jj} = g;
    c_opt{jj} = c;
end

end
