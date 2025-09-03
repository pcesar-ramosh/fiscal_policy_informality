function [r_opt, ir_out, pop1_vector, a, g_opt, c_opt] = ...
    huggett_Equi_RRA_function_transfer_on_CAL(eta_vector, sI_vector1, sF_vector1, cfg)
% HUGGETT_EQUI_RRA_FUNCTION_TRANSFER_ON_CAL
% - Transferencia entra al PRESUPUESTO del informal (no a utilidad)
% - FOC con tauc: c = [(1+tauc) V_a]^(-1/sigma)
% - Sin "Gov" en utilidad
% - Delta reducido (mejor condición numérica)
% - Fokker–Planck: normalización sum(g)*da = 1
% - r*_j por agente (bisección en cada j)

% Defaults
if nargin < 4 || isempty(cfg), cfg = struct(); end
def = struct( ...
    'scenario',            "baseline", ...
    'psi_I',               1.0, ...
    'psi_F',               1.0, ...
    'transfer_multiplier', 1.0, ...
    'keep_transfers_level', true, ...
    'amin_mode',           "baseline", ...
    'amin_policy',         "absolute", ...
    'amin_abs',           -1.00, ...
    'phi',                 0.05, ...
    'taxF',                0.10, ...
    'taxI',                0.00, ...
    'tauc',                0.00, ...
    'theta',               0.005, ...
    'r0',                  0.02, ...
    'rmin',               -0.01, ...
    'rmax',                0.06 ...
);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end

% Parámetros del modelo
rho     = 0.05;
z1_base = 0.33;
z2_base = 1.00;
I       = 700;
amax    = 5.0;

% Transiciones ocupacionales
p22 = 0.75; la2 = -log(p22);

% Salidas
pop1_vector = eta_vector;
J  = numel(pop1_vector);
r_opt = zeros(1, J);
g_opt = cell(1, J);
c_opt = cell(1, J);
ir_out = zeros(1, J);  % iteraciones usadas en r para cada j

for jj = 1:J
    pop1 = pop1_vector(jj);
    la1  = (1 - pop1) * la2 / pop1;
    Aswitch = [-speye(I)*la1,  speye(I)*la1;
                speye(I)*la2, -speye(I)*la2];

    taxF = cfg.taxF; taxI = cfg.taxI; tauc = cfg.tauc;
    phi  = cfg.phi;  theta = cfg.theta;

    if cfg.scenario == "covid_uptransfer"
        z1 = cfg.psi_I * z1_base;
        z2 = cfg.psi_F * z2_base;
        kappa = cfg.transfer_multiplier;
    else
        z1 = z1_base; z2 = z2_base; kappa = 1.0;
    end
    if cfg.keep_transfers_level
        Transfer = kappa * phi * z1_base;
    else
        Transfer = kappa * phi * z1;
    end

    % Malla de activos
    switch string(cfg.amin_policy)
        case "absolute", amin = cfg.amin_abs;
        case "relative", amin = -1.0*z1_base;
        case "shocked",  amin = -0.3*z1;
        otherwise,       amin = -0.3*z1_base;
    end
    a  = linspace(amin, amax, I)'; 
    da = (amax - amin) / (I - 1);
    aa = [a, a];

    % Numérico
    r     = cfg.r0;
    rmin  = cfg.rmin;  rmax = cfg.rmax;
    maxit = 100;       crit = 1e-6;
    Delta = 50;        % más diagonal
    Ir    = 100;       crit_S = 1e-5;

    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    % Inicial V0 (consumo sin ahorro)
    rr_vec = r + theta*(a<0);
    cI0 = max((1 - taxI)*z1 + rr_vec.*a + Transfer, 1e-10);
    cF0 = max((1 - taxF)*z2 + r*a,                 1e-10);
    V0 = zeros(I,2);
    V0(:,1) = (cI0.^(1 - sI)) / (1 - sI) / rho;
    V0(:,2) = (cF0.^(1 - sF)) / (1 - sF) / rho;

    % Búsqueda de r*
    for ir = 1:Ir
        V = V0;  % warm start con la solución previa

        for n = 1:maxit
            % Derivadas forward/backward (fronteras con envelope)
            dVf = zeros(I,2); dVb = zeros(I,2);
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            cI_sup = max((1 - taxI)*z1 + r*amax + Transfer, 1e-10);
            cF_sup = max((1 - taxF)*z2 + r*amax,            1e-10);
            dVf(I,1) = cI_sup^(-sI) / (1 + tauc);
            dVf(I,2) = cF_sup^(-sF) / (1 + tauc);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            cI_inf = max((1 - taxI)*z1 + (r+theta)*amin + Transfer, 1e-10);
            cF_inf = max((1 - taxF)*z2 + r*amin,                      1e-10);
            dVb(1,1) = cI_inf^(-sI) / (1 + tauc);
            dVb(1,2) = cF_inf^(-sF) / (1 + tauc);

            % rr(a)
            rr_vec = r + theta*(a<0);

            % Recursos
            res_I = (1 - taxI)*z1 + rr_vec.*aa(:,1) + Transfer;
            res_F = (1 - taxF)*z2 + r*aa(:,2);

            % FOC con tauc
            cf  = [ ((1+tauc)*max(dVf(:,1),1e-10)).^(-1/sI) , ...
                    ((1+tauc)*max(dVf(:,2),1e-10)).^(-1/sF) ];
            cb  = [ ((1+tauc)*max(dVb(:,1),1e-10)).^(-1/sI) , ...
                    ((1+tauc)*max(dVb(:,2),1e-10)).^(-1/sF) ];

            % Drifts
            ssf = [res_I, res_F] - cf;
            ssb = [res_I, res_F] - cb;

            % Upwind
            If = ssf > 0;  Ib = ssb < 0;  I0 = ~(If | Ib);
            c  = cf.*If + cb.*Ib + [res_I, res_F].*I0;

            % Utilidad instantánea
            U1 = (max(c(:,1),1e-12).^(1 - sI)) / (1 - sI);
            U2 = (max(c(:,2),1e-12).^(1 - sF)) / (1 - sF);
            u  = [U1, U2];

            % Matriz A
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z =  max(ssf,0)/da;

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            % HJB implícito
            B = (1/Delta + rho)*speye(2*I) - A;
            V_st = B \ ([u(:,1);u(:,2)] + [V(:,1);V(:,2)]/Delta);
            V_new = [V_st(1:I), V_st(I+1:2*I)];

            if max(max(abs(V_new - V))) < crit
                V = V_new; 
                break
            end
            V = V_new;
        end % end HJB

        % FP con normalización directa
        AT = A'; bb = zeros(2*I,1);
        AT(1,:) = ones(1, 2*I) * da; bb(1) = 1;
        gg = AT \ bb;
        gg = max(gg,0); gg = gg / (sum(gg)*da);
        g  = [gg(1:I), gg(I+1:2*I)];

        % Clearing por agente
        S = g(:,1)'*a*da + g(:,2)'*a*da;

        % Bisección en r
        if S >  crit_S
            rmax = r;  r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r;  r = 0.5*(r + rmax);
        else
            break
        end

        % Warm-start para la prox iter de r
        V0 = V;
    end % end loop r

    r_opt(jj) = r;
    g_opt{jj} = g;
    c_opt{jj} = c;
    ir_out(jj) = ir;
end
end
