function [r_opt, ir, pop1_vector, a, g_opt, c_opt] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vector, sI_vector1, sF_vector1, cfg)
% HUGGETT_EQUI_RRA_FUNCTION_TRANSFER_ON
%  Versión corregida:
%   - Transferencia entra al PRESUPUESTO del informal (no a utilidad)
%   - Impuesto al consumo en la FOC: c = [(1+tauc)*V_a]^(-1/sigma)
%   - Sin "Gov" sumado al flujo de utilidad del HJB
%   - rr(a) vectorizado y mejoras de estabilidad
%
%  cfg (struct) principales:
%    scenario             : "baseline" | "covid_uptransfer"
%    psi_I, psi_F         : multiplicadores de ingreso (default 1,1)
%    transfer_multiplier  : κ>=1 (default 1)  -> solo si scenario ~= baseline
%    keep_transfers_level : true (nivel basado en z1_base) | false (proporcional a z1)
%    amin_mode            : "baseline" | "shocked"
%    phi, taxF, taxI, tauc, theta
%    r0, rmin, rmax
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

% Transiciones ocupacionales (Poisson)
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
    la1 = (1 - pop1) * la2 / pop1;        % informal -> formal
    % p11 = exp(-la1);                    % sin usar directamente
    la  = [la1, la2];

    taxF = cfg.taxF; taxI = cfg.taxI; tauc = cfg.tauc;
    phi  = cfg.phi;  theta = cfg.theta;

    % Ingresos con/ sin shock
    if cfg.scenario == "covid_uptransfer"
        z1 = cfg.psi_I * z1_base;
        z2 = cfg.psi_F * z2_base;
        kappa = cfg.transfer_multiplier;            % κ > 1
    else
        z1 = z1_base; z2 = z2_base;
        kappa = 1.0;
    end
    z  = [z1, z2];

    % Transfer: nivel basado en z1_base (si keep_transfers_level) o proporcional a z1
    if cfg.keep_transfers_level
        Transfer = kappa * phi * z1_base;          % contracíclico y reforzado
    else
        Transfer = kappa * phi * z1;               % proporcional al ingreso actual
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

    % Parámetros de solución
    r     = cfg.r0;
    rmin  = cfg.rmin;  rmax = cfg.rmax;
    maxit = 100;       crit = 1e-6;
    Delta = 50; %1000;
    Ir    = 100;       crit_S = 1e-5;

    % RRA
    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    % Matriz de transición ocupacional
    Aswitch = [-speye(I)*la(1),  speye(I)*la(1);
                speye(I)*la(2), -speye(I)*la(2)];

    % Valor inicial V0 (regla de pulgar)
    rr_vec = r + theta*(a<0); % vectorizado
    V0 = zeros(I,2);
    % Consumo de "no ahorro" inicial, con Transfer solo en informales
    cI0 = max((1 - taxI)*z(1) + rr_vec.*a + Transfer, 1e-10);
    cF0 = max((1 - taxF)*z(2) + r*a,            1e-10);
    V0(:,1) = (cI0.^(1 - sI)) / (1 - sI) / rho;
    V0(:,2) = (cF0.^(1 - sF)) / (1 - sF) / rho;

    % ===== Encontrar r* por bisección =====
    for ir = 1:Ir
        if ir > 1, V0 = V_r(:,:,ir-1); end
        V = V0;

        dVf = zeros(I,2); dVb = zeros(I,2);

        for n = 1:maxit
            % Derivadas forward/backward (con condiciones de frontera consistentes)
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            % Frontera superior: usar envelope V_a = u'(c)/(1+tauc)
            cI_sup = max((1 - taxI)*z(1) + r*amax + Transfer, 1e-10);
            cF_sup = max((1 - taxF)*z(2) + r*amax,            1e-10);
            dVf(I,1) = cI_sup^(-sI) / (1 + tauc);
            dVf(I,2) = cF_sup^(-sF) / (1 + tauc);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            % Frontera inferior
            cI_inf = max((1 - taxI)*z(1) + (r+theta)*amin + Transfer, 1e-10);
            cF_inf = max((1 - taxF)*z(2) + r*amin,                      1e-10);
            dVb(1,1) = cI_inf^(-sI) / (1 + tauc);
            dVb(1,2) = cF_inf^(-sF) / (1 + tauc);

            % rr(a) vectorizado
            rr_vec = r + theta*(a<0);

            % Recursos (incluye Transfer SOLO en informales)
            res_I = (1 - taxI)*z(1) + rr_vec.*aa(:,1) + Transfer;   % informales
            res_F = (1 - taxF)*z(2) + r*aa(:,2);                    % formales

            % FOC con tauc: u'(c) = (1+tauc)*V_a  => c = [(1+tauc)V_a]^(-1/s)
            dvf_I = max(dVf(:,1), 1e-10);  dvf_F = max(dVf(:,2), 1e-10);
            dvb_I = max(dVb(:,1), 1e-10);  dvb_F = max(dVb(:,2), 1e-10);

            cf  = [ ((1+tauc)*dvf_I).^(-1/sI) , ((1+tauc)*dvf_F).^(-1/sF) ];
            cb  = [ ((1+tauc)*dvb_I).^(-1/sI) , ((1+tauc)*dvb_F).^(-1/sF) ];

            % Drifts
            ssf = [res_I, res_F] - cf;
            ssb = [res_I, res_F] - cb;

            % Consumo en punto fijo (drift ~ 0)
            c0  = [res_I, res_F];

            % Upwind
            If = ssf > 0;
            Ib = ssb < 0;
            I0 = ~(If | Ib);

            c = cf.*If + cb.*Ib + c0.*I0;

            % Utilidades (sin Transfer en gusto, sin Gov)
            U1 = (max(c(:,1),1e-12).^(1 - sI)) / (1 - sI);
            U2 = (max(c(:,2),1e-12).^(1 - sF)) / (1 - sF);
            u  = [U1, U2];

            % Matriz A para FP (con chequeo de filas ~ 0)
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z =  max(ssf,0)/da;

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            if max(abs(sum(A,2))) > 1e-7
                error('Improper Transition Matrix (rows do not sum to zero)');
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

    % Aviso si r* quedó al borde
    if ir == Ir
        warning('r* at bracket edge: rmin=%.4f r=%.4f rmax=%.4f (S=%.3e)', rmin, r, rmax, S);
    end

    % Salidas por agente
    r_opt(jj) = r;
    g_opt{jj} = g;
    c_opt{jj} = c;
end

end
