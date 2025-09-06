function [r_opt, ir, pop1_vector, a, g_opt, c_opt] = ...
    huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1, cfg)
%HUGGETT_EQUI_RRA_FUNCTION_TRANSFER
%  Versión con soporte para shock COVID y reglas de transferencias.
%
%  Uso retro‑compatible:
%    [r_opt, ir, pop1, a, g_opt, c_opt] = huggett_Equi_RRA_function_transfer(eta, sI, sF)
%  Uso con opciones (recomendado):
%    [r_opt, ir, pop1, a, g_opt, c_opt] = huggett_Equi_RRA_function_transfer(eta, sI, sF, cfg)
%
%  cfg (struct) campos principales:
%    scenario              : "baseline" (default) | "covid"
%    psi_I, psi_F          : multiplicadores de ingreso (default 1,1)
%    keep_transfers_level  : true mantiene NIVEL pre‑shock (default: true en covid)
%    amin_mode             : "baseline" | "shocked" (límite inferior de activos)
%    gov_nonneg            : true → trunca bien público a >= 0 (default true)
%    phi                   : 0.13 (transfer rule as % of informal BASE income)
%    taxF, taxI, tauc      : tasas fiscales
%    theta                 : 0.02 (prima endeudamiento informal)
%    r0, rmin, rmax        : bracket para r*
%
%  Retorna por agente (jj):
%    r_opt(jj), a (malla), g_opt{jj}(Ix2), c_opt{jj}(Ix2), ir (iter alcanzada)
%
%  Autor: (tu nombre)
%  Fecha: (hoy)

% -----------------------------
% Defaults (si cfg no se pasa)
% -----------------------------
if nargin < 4 || isempty(cfg), cfg = struct(); end
def = struct( ...
    'scenario',             "baseline", ...
    'psi_I',                1.0, ...
    'psi_F',                1.0, ...
    'keep_transfers_level', true, ...
    'amin_mode',            "baseline", ...   % "baseline" o "shocked"
    'gov_nonneg',           true, ...
    'phi',                  0.13, ...
    'taxF',                 0.10, ...
    'taxI',                 0.00, ...
    'tauc',                 0.00, ...
    'theta',                0.02, ...
    'r0',                   0.03, ...
    'rmin',                 0.01, ...
    'rmax',                 0.04 ...
);
fields = fieldnames(def);
for k=1:numel(fields)
    f = fields{k};
    if ~isfield(cfg, f), cfg.(f) = def.(f); end
end

% -----------------------------
% Parámetros “globales” modelo
% -----------------------------
rho = 0.05;                % tasa de descuento
z1_base = 0.33;            % ingreso informal base
z2_base = 1.00;            % ingreso formal base
I  = 700;                  % puntos de la malla de activos
amax = 5.0;                % cota superior activos

% Saltos ocupacionales (ajuste p11 con tamaño del sector informal)
p22 = 0.75;                 % Prob[y2|y2] = exp(-lambda_2)
la2 = -log(p22);            % intensidad Formal -> Informal

% Resultados por agente
pop1_vector = eta_vector;
r_opt = zeros(1, numel(pop1_vector));
g_opt = cell(1, numel(pop1_vector));
c_opt = cell(1, numel(pop1_vector));

% -----------------------------
% Bucle por agente jj
% -----------------------------
for jj = 1:numel(pop1_vector)
    pop1 = pop1_vector(jj);

    % Intensidad de salto para que p1 dependa de pop1 (estacionario)
    la1 = (1 - pop1) * la2 / pop1;   % informal -> formal
    p11 = exp(-la1);
    la  = [la1, la2];

    % Parámetros fiscales / precios
    taxF = cfg.taxF; taxI = cfg.taxI; tauc = cfg.tauc;
    phi  = cfg.phi;  theta = cfg.theta;

    % Ingresos con/ sin shock
    if cfg.scenario == "covid"
        z1 = cfg.psi_I * z1_base;
        z2 = cfg.psi_F * z2_base;
    else
        z1 = z1_base; z2 = z2_base;
    end
    z  = [z1, z2];

    % Transfer: NIVEL fijo (phi * ingreso BASE) o proporcional al ingreso actual
    if cfg.keep_transfers_level
        Transfer = phi * z1_base;  % contracíclico: mantiene monto base
    else
        Transfer = phi * z1;       % proporcional: cae con ingreso
    end

    % Malla de activos
    if cfg.amin_mode == "shocked"
        amin = -0.3 * z1;          % límite cae con ingreso actual
    else
        amin = -0.3 * z1_base;     % límite del benchmark
    end
    a  = linspace(amin, amax, I)'; 
    da = (amax - amin)/(I - 1);
    aa = [a, a];
    zz = ones(I,1) * z;

    % Parámetros de solución HJB/FP
    r   = cfg.r0;                 % guess inicial
    rmin = cfg.rmin; rmax = cfg.rmax;
    maxit = 100;                  % iter máx. HJB
    crit  = 1e-6;                 % tolerancia HJB
    Delta = 1000;                 % paso temporal implícito
    Ir    = 100;                  % iter máx. para hallar r*
    crit_S = 1e-5;                % tolerancia clearing

    % RRA de este agente
    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    % Matriz de transición entre sectores
    Aswitch = [-speye(I)*la(1),  speye(I)*la(1);
                speye(I)*la(2), -speye(I)*la(2)];

    % ----- Valor inicial V0 (con Transfer para informales) -----
    rr = zeros(I,1);
    for i=1:I
        rr(i) = (a(i) >= 0) * r + (a(i) < 0) * (r + theta); % informal lending/borrowing
    end

    V0 = zeros(I,2);
    V0(:,1) = (((1 - taxI) * zz(:,1) + rr .* a + Transfer).^(1 - sI) / (1 - sI)) / rho;
    V0(:,2) = (((1 - taxF) * zz(:,2) + r  .* a).^(1 - sF) / (1 - sF)) / rho;

    % ====== Iteración para r* (mercado de activos) ======
    for ir = 1:Ir
        if ir > 1
            V0 = V_r(:,:,ir-1);
        end

        V = V0;
        dVf = zeros(I,2);
        dVb = zeros(I,2);

        for n = 1:maxit
            % Derivadas forward/backward
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVf(I,1) = ((1 - taxI) * z(1) + r  * amax)^(-sI);
            dVf(I,2) = ((1 - taxF) * z(2) + r  * amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;   % backward en 2..I
            dVb(1,1)   = ((1 - taxI) * z(1) + (r + theta) * amin)^(-sI);
            dVb(1,2)   = ((1 - taxF) * z(2) + r * amin)^(-sF);

            % Recalcular rr (informal) por signo de a
            for i=1:I
                rr(i) = (a(i) >= 0) * r + (a(i) < 0) * (r + theta);
            end

            % Consumo & ahorro (forward/backward)
            cf  = [max(dVf(:,1),1e-10).^(-1/sI), max(dVf(:,2),1e-10).^(-1/sF)];
            ssf = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cf;

            cb  = [max(dVb(:,1),1e-10).^(-1/sI), max(dVb(:,2),1e-10).^(-1/sF)];
            %cb  = [max(dVb(:,1),1e-10).^(-1/sI), max(dVb(:,2),1e-10).^(-1/sF)];
            %ssb = [(1 - taxI)*zz[:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cb; %#ok<NBRAK>
            ssb = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cb;

            % Consumo en estado estacionario
            c0  = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa;

            % Upwind
            If = ssf > 0;      % drift positivo -> forward
            Ib = ssb < 0;      % drift negativo -> backward
            I0 = (1 - If - Ib);

            c = cf.*If + cb.*Ib + c0.*I0;

            % Utilidad (informal con Transfer fijo)
            U1 = ((c(:,1) + Transfer).^(1 - sI)) / (1 - sI);
            U2 = (c(:,2).^(1 - sF)) / (1 - sF);

            % Presupuesto del Gobierno (bien público + Transfer)
            popI = la2 / (la1 + la2);
            popF = la1 / (la1 + la2);

            TaxRevenue = taxI * z(1) * popI + taxF * z(2) * popF + tauc * (0); %#ok<*NASGU>
            Gov_pub    = TaxRevenue - Transfer * popI;
            if cfg.gov_nonneg
                Gov_pub = max(0, Gov_pub);
            end
            Gov = Gov_pub + Transfer;

            u = [U1, U2] + Gov * ones(I,2);

            % Matriz de transición A (Fokker-Planck)
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z =  max(ssf,0)/da;

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            if max(abs(sum(A,2))) > 1e-9
                error('Improper Transition Matrix');
            end

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
        end % end HJB loop

        % ========= Fokker-Planck en estacionario =========
        AT = A';
        bb = zeros(2*I,1);
        i_fix = 1; bb(i_fix) = .1;
        row = zeros(1, 2*I); row(i_fix) = 1;
        AT(i_fix, :) = row;

        gg = AT \ bb;                 % densidad apilada
        g_sum = gg' * ones(2*I,1) * da;
        gg = gg ./ g_sum;

        g = [gg(1:I), gg(I+1:2*I)];   % col1 informal, col2 formal

        % Ahorro agregado (clearing)
        S = g(:,1)' * a * da + g(:,2)' * a * da;

        % Actualizar r (bisección)
        if S >  crit_S
            rmax = r;  r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r;  r = 0.5*(r + rmax);
        else
            % equilibrio encontrado
            break
        end

        % almacenar para próxima iter
        V_r(:,:,ir)   = V;
        g_r(:,:,ir)   = g;
        adot(:,:,ir)  = [(1 - taxI)*zz(:,1), (1 - taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - c;
    end % end loop r*

    % Salidas por agente
    r_opt(jj) = r;
    g_opt{jj} = g;
    c_opt{jj} = c;

end % end loop jj

end
