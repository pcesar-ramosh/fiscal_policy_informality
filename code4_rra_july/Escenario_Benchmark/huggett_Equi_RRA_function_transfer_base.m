function [r_opt, ir, pop1_vector, a, g_opt, c_opt] = ...
    huggett_Equi_RRA_function_transfer_base(eta_vector, sI_vector1, sF_vector1, cfg)
% HUGGETT_EQUI_RRA_FUNCTION_TRANSFER_BASE
% Equilibrio general con informalidad y transferencia focalizada a informales.
% Soporta escenario baseline o "covid" via cfg, pero aqui lo usamos para baseline.
%
% Entradas:
%   eta_vector   : 1xJ, proporcion del sector informal por tipo
%   sI_vector1   : 1xJ, RRA para informales
%   sF_vector1   : 1xJ, RRA para formales
%   cfg struct   :
%     scenario              = 'baseline' | 'covid'
%     psi_I, psi_F          = multiplicadores de ingreso (default 1)
%     keep_transfers_level  = true -> Transfer = phi * z1_base
%     transfer_multiplier   = kappa (default 1)
%     amin_mode             = 'baseline' | 'shocked'
%     gov_nonneg            = true (truncar bien publico a >= 0)
%     phi                   = 0.13
%     taxF, taxI, tauc      = tasas (tauc no usado en presupuesto aqui)
%     theta                 = 0.02 (prima de endeudamiento informal)
%     r0, rmin, rmax        = bracket para r*
%
% Salidas:
%   r_opt(1xJ)   : tasa de interes de equilibrio por tipo
%   ir           : iteraciones usadas en la ultima busqueda de r
%   pop1_vector  : copia de eta_vector
%   a(Ix1)       : malla de activos
%   g_opt{j}(Ix2): distribuciones estacionarias (col1 informal, col2 formal)
%   c_opt{j}(Ix2): politicas de consumo

% -------- defaults --------
if nargin < 4 || isempty(cfg), cfg = struct(); end
def = struct( ...
    'scenario',            'baseline', ...
    'psi_I',               1.0, ...
    'psi_F',               1.0, ...
    'keep_transfers_level', true, ...
    'transfer_multiplier', 1.0, ...
    'amin_mode',           'baseline', ...
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
    if ~isfield(cfg, fn{k})
        cfg.(fn{k}) = def.(fn{k});
    end
end

% -------- parametros base --------
rho     = 0.05;
z1_base = 0.33;
z2_base = 1.00;
I       = 700;
amax    = 5.0;

% Transiciones ocupacionales (CTMC 2 estados)
p22 = 0.75;           % Prob[y2|y2] formal->formal
la2 = -log(p22);      % intensidad formal->informal

J = numel(eta_vector);

% Salidas
pop1_vector = eta_vector;
r_opt = zeros(1, J);
g_opt = cell(1, J);
c_opt = cell(1, J);
ir    = 0;

for jj = 1:J
    pop1 = pop1_vector(jj);

    % Elegir la1 para que masa estacionaria coincida con pop1
    la1 = (1 - pop1) * la2 / pop1;   % informal->formal
    la  = [la1, la2];

    % Ingresos (baseline o covid)
    if strcmp(cfg.scenario, 'covid')
        z1 = cfg.psi_I * z1_base;
        z2 = cfg.psi_F * z2_base;
    else
        z1 = z1_base;
        z2 = z2_base;
    end
    z  = [z1, z2];

    % Transfer focalizada a informales
    if cfg.keep_transfers_level
        Transfer = cfg.transfer_multiplier * cfg.phi * z1_base;
    else
        Transfer = cfg.transfer_multiplier * cfg.phi * z1;
    end

    % Malla de activos
    if strcmp(cfg.amin_mode, 'shocked')
        amin = -0.3 * z1;       % mas restrictivo bajo shock
    else
        amin = -0.3 * z1_base;  % fijar con benchmark
    end
    a  = linspace(amin, amax, I)';
    da = (amax - amin) / (I - 1);
    aa = [a, a];
    zz = ones(I,1) * z;

    % Parametros numericos
    r     = cfg.r0;
    rmin  = cfg.rmin;
    rmax  = cfg.rmax;
    maxit = 100;
    crit  = 1e-6;
    Delta = 1000;
    Ir    = 100;
    crit_S = 1e-5;

    % Preferencias
    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    % Matriz de switching ocupacional
    Aswitch = [-speye(I)*la(1),  speye(I)*la(1); ...
                speye(I)*la(2), -speye(I)*la(2)];

    % Valor inicial V0
    rr = zeros(I,1);
    for i = 1:I
        if a(i) >= 0
            rr(i) = r;                 % prestar al r
        else
            rr(i) = r + cfg.theta;     % endeudarse al r+theta
        end
    end
    V0 = zeros(I,2);
    V0(:,1) = (((1 - cfg.taxI) * zz(:,1) + rr .* a + Transfer).^(1 - sI) / (1 - sI)) / rho;
    V0(:,2) = (((1 - cfg.taxF) * zz(:,2) + r  .* a).^(1 - sF)   / (1 - sF)) / rho;

    % ===== busqueda de r* por biseccion =====
    for ir = 1:Ir
        if ir > 1
            V0 = V_r(:,:,ir-1);
        end
        V = V0;

        dVf = zeros(I,2);
        dVb = zeros(I,2);

        for n = 1:maxit
            % derivadas forward/backward
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVf(I,1) = ((1 - cfg.taxI)*z(1) + r*amax)^(-sI);
            dVf(I,2) = ((1 - cfg.taxF)*z(2) + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVb(1,1)   = ((1 - cfg.taxI)*z(1) + (r+cfg.theta)*amin)^(-sI);
            dVb(1,2)   = ((1 - cfg.taxF)*z(2) + r*amin)^(-sF);

            % rr(a) con r actual
            for i = 1:I
                if a(i) >= 0
                    rr(i) = r;
                else
                    rr(i) = r + cfg.theta;
                end
            end

            % Consumo y drift
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

            % Gobierno (simple)
            popI = la2 / (la1 + la2);
            popF = la1 / (la1 + la2);
            TaxRevenue = cfg.taxI*z(1)*popI + cfg.taxF*z(2)*popF;
            Gov_pub    = TaxRevenue - Transfer * popI;
            if cfg.gov_nonneg
                Gov_pub = max(0, Gov_pub);
            end
            Gov = Gov_pub + Transfer;

            u = [U1, U2] + Gov * ones(I,2);

            % Generador A
            X = -min(ssb,0)/da;
            Y = -max(ssf,0)/da + min(ssb,0)/da;
            Z =  max(ssf,0)/da;

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            if max(abs(sum(A,2))) > 1e-9
                error('Improper Transition Matrix');
            end

            % HJB implicito
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

        % Fokker-Planck estacionario
        AT = A';
        bb = zeros(2*I,1);
        i_fix = 1; bb(i_fix) = .1;
        row = zeros(1, 2*I); row(i_fix) = 1;
        AT(i_fix,:) = row;

        gg = AT \ bb;
        g_sum = gg' * ones(2*I,1) * da;
        gg = gg ./ g_sum;

        g = [gg(1:I), gg(I+1:2*I)];

        % Clearing de mercado de activos
        S = g(:,1)'*a*da + g(:,2)'*a*da;

        % Actualizar r por biseccion
        if S >  crit_S
            rmax = r;
            r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r;
            r = 0.5*(r + rmax);
        else
            break
        end

        V_r(:,:,ir) = V; %#ok<AGROW>
    end % end loop r

    % Guardar resultados por tipo
    r_opt(jj) = r;
    g_opt{jj} = g;
    c_opt{jj} = c;
end

end
