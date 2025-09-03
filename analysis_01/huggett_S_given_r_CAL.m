function [S_by_agent, S_total, g_out, c_out, a_grid] = ...
    huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg, r_fixed)
% Resuelve HJB+FP con r = r_fixed (común a todos los j) y retorna
% S_j(r), S_total(r), densidades y políticas.

% Defaults (subset relevante)
def = struct( ...
    'scenario',            "baseline", ...
    'psi_I',               1.0, ...
    'psi_F',               1.0, ...
    'transfer_multiplier', 1.0, ...
    'keep_transfers_level', true, ...
    'amin_policy',         "absolute", ...
    'amin_abs',           -1.00, ...
    'phi',                 0.05, ...
    'taxF',                0.10, ...
    'taxI',                0.00, ...
    'tauc',                0.00, ...
    'theta',               0.005 ...
);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end

% Parámetros numéricos y del modelo
rho   = 0.05;
z1b   = 0.33;  z2b = 1.00;
I     = 700;   amax = 5.0;
p22   = 0.75;  la2  = -log(p22);
Delta = 50;    maxit = 100; crit = 1e-6;
tauc  = cfg.tauc; taxF = cfg.taxF; taxI = cfg.taxI; theta = cfg.theta;

J = numel(eta_vector);
S_by_agent = zeros(J,1);
g_out = cell(J,1); c_out = cell(J,1);

for jj=1:J
    pop1 = eta_vector(jj);
    la1  = (1 - pop1) * la2 / pop1;
    Asw  = [-speye(I)*la1, speye(I)*la1; speye(I)*la2, -speye(I)*la2];

    % ingresos & transfer
    if cfg.scenario == "covid_uptransfer"
        z1 = cfg.psi_I * z1b; z2 = cfg.psi_F * z2b; kappa = cfg.transfer_multiplier;
    else
        z1 = z1b; z2 = z2b; kappa = 1.0;
    end
    if cfg.keep_transfers_level, Transfer = kappa*cfg.phi*z1b; else, Transfer = kappa*cfg.phi*z1; end

    % malla a
    switch string(cfg.amin_policy)
        case "absolute", amin = cfg.amin_abs;
        case "relative", amin = -1.0*z1b;
        case "shocked",  amin = -0.3*z1;
        otherwise,       amin = -0.3*z1b;
    end
    a_grid = linspace(amin, amax, I)'; da = (amax-amin)/(I-1);

    % inicialización V0
    rr = r_fixed + theta*(a_grid<0);
    cI0 = max((1 - taxI)*z1 + rr.*a_grid + Transfer, 1e-10);
    cF0 = max((1 - taxF)*z2 + r_fixed*a_grid,        1e-10);
    V = zeros(I,2);
    sI = sI_vector1(jj); sF = sF_vector1(jj);
    V(:,1) = (cI0.^(1 - sI)) / (1 - sI) / rho;
    V(:,2) = (cF0.^(1 - sF)) / (1 - sF) / rho;

    for n=1:maxit
        % derivadas
        dVf = zeros(I,2); dVb = zeros(I,2);
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        cI_sup = max((1 - taxI)*z1 + r_fixed*amax + Transfer, 1e-10);
        cF_sup = max((1 - taxF)*z2 + r_fixed*amax,            1e-10);
        dVf(I,1) = cI_sup^(-sI) / (1+tauc);
        dVf(I,2) = cF_sup^(-sF) / (1+tauc);
        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        cI_inf = max((1 - taxI)*z1 + (r_fixed+theta)*amin + Transfer, 1e-10);
        cF_inf = max((1 - taxF)*z2 + r_fixed*amin,                     1e-10);
        dVb(1,1) = cI_inf^(-sI) / (1+tauc);
        dVb(1,2) = cF_inf^(-sF) / (1+tauc);

        rr = r_fixed + theta*(a_grid<0);
        res_I = (1 - taxI)*z1 + rr.*a_grid + Transfer;
        res_F = (1 - taxF)*z2 + r_fixed*a_grid;

        cf = [ ((1+tauc)*max(dVf(:,1),1e-10)).^(-1/sI) , ...
               ((1+tauc)*max(dVf(:,2),1e-10)).^(-1/sF) ];
        cb = [ ((1+tauc)*max(dVb(:,1),1e-10)).^(-1/sI) , ...
               ((1+tauc)*max(dVb(:,2),1e-10)).^(-1/sF) ];

        ssf = [res_I, res_F] - cf;
        ssb = [res_I, res_F] - cb;
        If = ssf>0; Ib = ssb<0; I0 = ~(If|Ib);
        c  = cf.*If + cb.*Ib + [res_I, res_F].*I0;

        U1 = (max(c(:,1),1e-12).^(1 - sI)) / (1 - sI);
        U2 = (max(c(:,2),1e-12).^(1 - sF)) / (1 - sF);
        u  = [U1, U2];

        X = -min(ssb,0)/da;  Y = -max(ssf,0)/da + min(ssb,0)/da;  Z = max(ssf,0)/da;
        A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1, sparse(I,I); sparse(I,I), A2] + Asw;

        B = (1/Delta + rho)*speye(2*I) - A;
        V_st = B \ ([u(:,1);u(:,2)] + [V(:,1);V(:,2)]/Delta);
        V_new = [V_st(1:I), V_st(I+1:2*I)];

        if max(max(abs(V_new - V))) < crit, V = V_new; break; end
        V = V_new;
    end

    % FP con normalización
    AT = A'; bb = zeros(2*I,1);
    AT(1,:) = ones(1, 2*I)*da; bb(1) = 1;
    gg = AT \ bb; gg = max(gg,0); gg = gg / (sum(gg)*da);
    g  = [gg(1:I), gg(I+1:2*I)];

    S_by_agent(jj) = (g(:,1)'*a_grid + g(:,2)'*a_grid)*da;
    g_out{jj} = g; c_out{jj} = c;
end

S_total = sum(S_by_agent);
end
