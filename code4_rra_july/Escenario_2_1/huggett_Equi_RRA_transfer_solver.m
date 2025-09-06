function res = huggett_Equi_RRA_transfer_solver(eta_vector, sI_vector1, sF_vector1, cfg)
% HUGGETT_EQUI_RRA_TRANSFER_SOLVER
% Resuelve el equilibrio para una configuración dada (shock y transfer),
% devolviendo políticas, distribuciones y estadísticas útiles.
%
% res:
%   .r_opt (1xJ)              : r* por agente
%   .a (Ix1)                  : malla de activos
%   .g_opt {1xJ} (Ix2)        : densidades (col1 informal, col2 formal)
%   .c_opt {1xJ} (Ix2)        : consumo óptimo
%   .statsC_mean_inf (1xJ)    : consumo medio informal por agente (media ponderada por g)
%
% Notas:
%   - Transferencia a informales: Transfer = kappa * phi * z1_base (si keep_transfers_level=true)
%   - Si scenario="shock", psi_I y psi_F escalan z1,z2 (ingresos).
%   - Límite inferior de activos (amin) según cfg.amin_mode.

% Defaults si faltan campos
if nargin < 4 || isempty(cfg), cfg = struct(); end
def = struct( ...
    'scenario',             "baseline", ...
    'psi_I',                1.0, ...
    'psi_F',                1.0, ...
    'keep_transfers_level', true, ...
    'transfer_multiplier',  1.0, ...
    'amin_mode',            "baseline", ...
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
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end

% Parámetros globales
rho     = 0.05;
z1_base = 0.33;
z2_base = 1.00;
I       = 700;
amax    = 5.0;

% Otras constantes
p22 = 0.75;            % formal->formal
la2 = -log(p22);       % intensidad formal->informal

J = numel(eta_vector);

% Salidas
res.r_opt            = zeros(1,J);
res.g_opt            = cell(1,J);
res.c_opt            = cell(1,J);
res.statsC_mean_inf  = zeros(1,J);
res.a                = [];        % se setea en el loop (constante entre agentes)

for jj = 1:J
    pop1 = eta_vector(jj);

    % Intensidad informal->formal para sostener pop estacionaria
    la1 = (1 - pop1) * la2 / pop1;
    la  = [la1, la2];

    % Ingresos con/sin shock
    if cfg.scenario == "shock"
        z1 = cfg.psi_I * z1_base;
        z2 = cfg.psi_F * z2_base;
    else
        z1 = z1_base; z2 = z2_base;
    end
    z  = [z1, z2];

    % Transfer nivel vs proporcional
    if cfg.keep_transfers_level
        Transfer = cfg.transfer_multiplier * cfg.phi * z1_base;
    else
        Transfer = cfg.transfer_multiplier * cfg.phi * z1;
    end

    % Malla de activos
    if cfg.amin_mode == "shocked"
        amin = -0.3 * z1;
    else
        amin = -0.3 * z1_base;
    end
    a  = linspace(amin, amax, I)'; 
    da = (amax - amin) / (I - 1);
    aa = [a, a];
    zz = ones(I,1) * z;

    if jj == 1
        res.a = a; % guardar malla
    end

    % Parámetros solución
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
        if a(i) >= 0, rr(i) = r; else, rr(i) = r + cfg.theta; end
    end
    V0 = zeros(I,2);
    V0(:,1) = (((1 - cfg.taxI) * zz(:,1) + rr .* a + Transfer).^(1 - sI) / (1 - sI)) / rho;
    V0(:,2) = (((1 - cfg.taxF) * zz(:,2) + r  .* a).^(1 - sF)   / (1 - sF)) / rho;

    % ===== Buscar r* por bisección =====
    for ir = 1:Ir
        if ir > 1, V0 = V_r(:,:,ir-1); end
        V = V0;

        dVf = zeros(I,2); dVb = zeros(I,2);

        for n = 1:maxit
            % Derivadas
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVf(I,1) = ((1 - cfg.taxI)*z(1) + r*amax)^(-sI);
            dVf(I,2) = ((1 - cfg.taxF)*z(2) + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVb(1,1)   = ((1 - cfg.taxI)*z(1) + (r+cfg.theta)*amin)^(-sI);
            dVb(1,2)   = ((1 - cfg.taxF)*z(2) + r*amin)^(-sF);

            % rr(a)
            for i = 1:I
                if a(i) >= 0, rr(i) = r; else, rr(i) = r + cfg.theta; end
            end

            % Consumo/ahorro
            cf  = [max(dVf(:,1),1e-10).^(-1/sI), max(dVf(:,2),1e-10).^(-1/sF)];
            ssf = [(1 - cfg.taxI)*zz(:,1), (1 - cfg.taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cf;

            cb  = [max(dVb(:,1),1e-10).^(-1/sI), max(dVb(:,2),1e-10).^(-1/sF)];
            ssb = [(1 - cfg.taxI)*zz(:,1), (1 - cfg.taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa - cb;

            c0  = [(1 - cfg.taxI)*zz(:,1), (1 - cfg.taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa;

            If = ssf > 0;  Ib = ssb < 0;  I0 = (1 - If - Ib);
            c  = cf.*If + cb.*Ib + c0.*I0;

            % Utilidad
            U1 = ((c(:,1) + Transfer).^(1 - sI)) / (1 - sI);
            U2 =  (c(:,2).^(1 - sF)) / (1 - sF);

            % Gobierno
            popI = la2 / (la1 + la2);
            popF = la1 / (la1 + la2);
            TaxRevenue = cfg.taxI*z(1)*popI + cfg.taxF*z(2)*popF; % (tauc no implementado)
            Gov_pub    = TaxRevenue - Transfer * popI;
            if cfg.gov_nonneg, Gov_pub = max(0, Gov_pub); end
            Gov = Gov_pub + Transfer;

            u = [U1, U2] + Gov * ones(I,2);

            % Matriz A (FP)
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z =  max(ssf,0)/da;

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            if max(abs(sum(A,2))) > 1e-9
                error('Improper Transition Matrix');
            end

            % HJB implícito
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

        % FP estacionario
        AT = A';
        bb = zeros(2*I,1);
        i_fix = 1; bb(i_fix) = .1;
        row = zeros(1, 2*I); row(i_fix) = 1;
        AT(i_fix,:) = row;

        gg = AT \ bb;
        g_sum = gg' * ones(2*I,1) * da;
        gg = gg ./ g_sum;

        g = [gg(1:I), gg(I+1:2*I)];

        % Clearing (ahorro agregado)
        S = g(:,1)'*a*da + g(:,2)'*a*da;

        % Actualiza r por bisección
        if S >  1e-5
            rmax = r;  r = 0.5*(r + rmin);
        elseif S < -1e-5
            rmin = r;  r = 0.5*(r + rmax);
        else
            break
        end

        V_r(:,:,ir) = V; %#ok<AGROW>
    end % end loop r*

    % Salidas por agente
    res.r_opt(jj) = r;
    res.g_opt{jj} = g;
    res.c_opt{jj} = c;

    % Consumo medio informal ponderado por g_inf(a)
    g1 = g(:,1); g1(g1<0)=0;
    ci = c(:,1);
    res.statsC_mean_inf(jj) = (g1' * ci * da) / (sum(g1)*da);
end

end
