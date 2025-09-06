function [r_opt, ir, pop1_vector, a, g_opt, c_opt, stats] = ...
    huggett_Equi_RRA_function_transfer_taxes(eta_vector, sI_vector1, sF_vector1, cfg)
% Equilibrio HA con informalidad, transferencias a informales,
% impuesto al ingreso (solo formales) e impuesto al consumo (ambos).
%
% Salidas: r_opt, ir, pop1_vector, a, g_opt{j}, c_opt{j}, stats (ver abajo)

% -------- Defaults --------
if nargin < 4 || isempty(cfg), cfg = struct(); end
def = struct( ...
    'psi_I', 1.0, 'psi_F', 1.0, ...
    'taxI', 0.0, 'taxF', 0.0, 'tauc', 0.0, ...
    'phi', 0.0, 'transfer_mode','level_base', 'keep_transfers_level', true, ...
    'gshare_public', 0.0, 'gov_nonneg', true, ...
    'theta', 0.02, ...
    'r0', 0.03, 'rmin', 0.01, 'rmax', 0.04, ...
    'amin_mode', 'baseline' ...
);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end

% -------- Basicos --------
rho     = 0.05;
z1_base = 0.33;  % ingreso informal base
z2_base = 1.00;  % ingreso formal base
I       = 700;   % puntos en malla de activos
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

% Contenedores para estadisticos
cmean_inf = zeros(J,1);
cmean_for = zeros(J,1);
frac_bound = zeros(J,2);
gini_inf   = zeros(J,1);

% Agregados fiscales ex post (promedio sobre tipos)
rev_income_acc = 0; rev_tauc_acc = 0; transfers_acc = 0; gov_pub_acc = 0; gov_surplus_acc = 0;

for jj = 1:J
    pop1 = pop1_vector(jj);
    % ajustar la1 para estacionaria = pop1
    la1 = (1 - pop1) * la2 / pop1;
    la  = [la1, la2];

    % ingresos
    z1 = cfg.psi_I * z1_base;
    z2 = cfg.psi_F * z2_base;
    z  = [z1, z2];

    % transferencia focalizada (informales)
    if strcmp(cfg.transfer_mode,'proportional')
        Transfer = cfg.phi * z1;       % proporcional al ingreso informal
    else
        if cfg.keep_transfers_level
            Transfer = cfg.phi * z1_base;
        else
            Transfer = cfg.phi * z1;
        end
    end

    % malla de activos
    if strcmp(cfg.amin_mode,'shocked')
        amin = -0.3 * z1;
    else
        amin = -0.3 * z1_base;
    end
    a  = linspace(amin, amax, I)';
    da = (amax - amin) / (I - 1);
    aa = [a, a];
    zz = ones(I,1) * z;

    % parametros numericos (robustos)
    r     = cfg.r0; rmin = cfg.rmin; rmax = cfg.rmax;
    maxit = 200;        % + iteraciones HJB
    crit  = 1e-6;       
    Delta = 3000;       % mas implicito => mejor condicionamiento
    Ir    = 120;        
    crit_S= 5e-5;       % clearing menos estricto

    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    Aswitch = [-speye(I)*la(1),  speye(I)*la(1); ...
                speye(I)*la(2), -speye(I)*la(2)];

    % valor inicial V0 (simple, suficiente)
    rr = zeros(I,1);
    for i = 1:I
        if a(i) >= 0, rr(i) = r; else, rr(i) = r + cfg.theta; end
    end
    V0 = zeros(I,2);
    V0(:,1) = (((1 - cfg.taxI) * zz(:,1) + rr .* a + Transfer).^(1 - sI) / (1 - sI)) / rho;
    V0(:,2) = (((1 - cfg.taxF) * zz(:,2) + r  .* a).^(1 - sF)   / (1 - sF)) / rho;

    % ===== Busqueda de r* por biseccion =====
    for ir = 1:Ir
        if ir > 1, V0 = V_r(:,:,ir-1); end
        V = V0;

        dVf = zeros(I,2); dVb = zeros(I,2);

        for n = 1:maxit
            % ===== Derivadas forward/backward =====
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVf(I,1) = ((1 - cfg.taxI)*z(1) + r*amax)^(-sI);
            dVf(I,2) = ((1 - cfg.taxF)*z(2) + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:)) / da;
            dVb(1,1)   = ((1 - cfg.taxI)*z(1) + (r+cfg.theta)*amin)^(-sI);
            dVb(1,2)   = ((1 - cfg.taxF)*z(2) + r*amin)^(-sF);

            % ---- Pisos de derivadas ----
            dVf = max(dVf, 1e-12);
            dVb = max(dVb, 1e-12);

            % ===== r(a) con prima para endeudamiento informal =====
            for i = 1:I
                if a(i) >= 0, rr(i) = r; else, rr(i) = r + cfg.theta; end
            end

            % ===== Consumo y drift con impuesto al consumo =====
            tau1 = 1 + cfg.tauc;

            cf  = [dVf(:,1).^(-1/sI) ./ tau1, dVf(:,2).^(-1/sF) ./ tau1];
            cb  = [dVb(:,1).^(-1/sI) ./ tau1, dVb(:,2).^(-1/sF) ./ tau1];

            % ---- Pisos de consumo upwind ----
            cf = max(cf, 1e-12);
            cb = max(cb, 1e-12);

            resources = [(1 - cfg.taxI)*zz(:,1), (1 - cfg.taxF)*zz(:,2)] + [rr, r*ones(I,1)].*aa;
            resources = max(resources, 1e-12);

            ssf = resources - tau1.*cf;
            ssb = resources - tau1.*cb;

            c0  = max(resources ./ tau1, 1e-12);

            If = ssf > 0;  Ib = ssb < 0;  I0 = (1 - If - Ib);
            c  = cf.*If + cb.*Ib + c0.*I0;
            c(~isfinite(c)) = 1e-12;
            c = max(c, 1e-12);

            % ===== Utilidad =====
            U1 = ((c(:,1) + Transfer).^(1 - sI)) / (1 - sI);
            U2 =  (c(:,2).^(1 - sF)) / (1 - sF);

            % ===== Bien público en utilidad (si se usa) =====
            popI = la2 / (la1 + la2);
            popF = la1 / (la1 + la2);
            TaxRev_income  = cfg.taxF*z(2)*popF + cfg.taxI*z(1)*popI;
            Gov_pub_exante = TaxRev_income - Transfer * popI;
            if cfg.gov_nonneg, Gov_pub_exante = max(0, Gov_pub_exante); end
            Gov = cfg.gshare_public * Gov_pub_exante;

            u = [U1, U2] + Gov * ones(I,2);

            % ---------- Generador A (robusto) ----------
            X = max(-min(ssb,0)/da, 0);
            Z = max( max(ssf,0)/da, 0);

            % Piso de flujos (evita filas casi nulas)
            flow_floor = 1e-12;
            X = max(X, flow_floor);
            Z = max(Z, flow_floor);

            % Vectores sub/super (longitud I-1)
            L1 = X(2:I,1);   U1d = Z(1:I-1,1);
            L2 = X(2:I,2);   U2d = Z(1:I-1,2);

            % Padding a longitud I
            L1p = [0; L1];   U1p = [U1d; 0];
            L2p = [0; L2];   U2p = [U2d; 0];

            D1  = -(L1p + U1p);
            D2  = -(L2p + U2p);

            A1 = spdiags([-L1p, D1, U1p], [-1, 0, 1], I, I);
            A2 = spdiags([-L2p, D2, U2p], [-1, 0, 1], I, I);

            % Conmutación ocupacional
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            % Proyección conservativa: suma por fila = 0
            rowsum = A*ones(2*I,1);
            A = A - spdiags(rowsum, 0, 2*I, 2*I);
            A(abs(A) < 1e-14) = 0;

            % ===== HJB implícito con regularización adaptativa + damping =====
            %reg = 1e-8;                 % regularización base
            %B = (1/Delta + rho + reg)*speye(2*I) - A;

            %rc = rcond(B);  
            %if rc < 1e-10
               % reg = 1e-6;
                %B = (1/Delta + rho + reg)*speye(2*I) - A;
            %end
            %if rcond(B) < 1e-12
               % reg = 1e-4;
               % B = (1/Delta + rho + reg)*speye(2*I) - A;
            %end

            reg = 1e-8;
            B = (1/Delta + rho + reg)*speye(2*I) - A;   % B es sparse
            
            % Usa condest (soporta sparse) en vez de rcond
            rc = condest(B);             % estima 1/cond(B)
            if rc < 1e-10
                reg = 1e-6;
                B = (1/Delta + rho + reg)*speye(2*I) - A;
                rc = condest(B);
            end
            if rc < 1e-12
                reg = 1e-4;
                B = (1/Delta + rho + reg)*speye(2*I) - A;
                % opcional: rc = condest(B);
            end
  

            u_st = [u(:,1); u(:,2)];
            V_st = [V(:,1); V(:,2)];
            bvec = u_st + V_st/Delta;

            V_st = B \ bvec;
            V_new = [V_st(1:I), V_st(I+1:2*I)];

            % ---- Damping en la actualización de V ----
            omega = 0.3;  % 0<omega<=1
            V_damped = (1-omega)*V + omega*V_new;

            % Convergencia
            if max(max(abs(V_damped - V))) < crit
                V = V_damped; break;
            end
            V = V_damped;
        end % HJB loop

        % ===== FP estacionaria (normalización exacta) =====
        AT = A';
        bb = zeros(2*I,1);

        % ∑_states ∫ g(a,s) da = 1  -> reemplazamos una fila por ecuación de masa total
        i_fix = 1;
        row = ones(1, 2*I) * da;
        bb(i_fix) = 1.0;
        AT(i_fix,:) = row;

        gg = AT \ bb;

        % Normaliza por seguridad
        g_sum = gg' * ones(2*I,1) * da;
        gg = gg ./ g_sum;

        g = [gg(1:I), gg(I+1:2*I)];

        % Clearing de activos
        S = g(:,1)'*a*da + g(:,2)'*a*da;

        % Actualiza r
        if S >  crit_S
            rmax = r; r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r; r = 0.5*(r + rmax);
        else
            break
        end

        V_r(:,:,ir) = V; %#ok<AGROW>
    end % loop r*

    % Guardar por tipo
    r_opt(jj) = r;
    g_opt{jj} = g;
    c_opt{jj} = c;

    % -------- Estadisticos por tipo --------
    g1 = g(:,1); g2 = g(:,2); g1(g1<0)=0; g2(g2<0)=0;

    % Pesos internos (normalizar por masa de cada grupo)
    mass1 = sum(g1)*da; mass2 = sum(g2)*da;
    if mass1>0, w1 = g1/mass1; else, w1 = zeros(I,1); end
    if mass2>0, w2 = g2/mass2; else, w2 = zeros(I,1); end

    cmean_inf(jj) = sum(w1 .* c(:,1)) * da;
    cmean_for(jj) = sum(w2 .* c(:,2)) * da;

    % Fraccion en la restriccion (aprox masa en a_min)
    frac_bound(jj,:) = (g(1,:) ./ sum(g(:)))';

    % Gini de riqueza para informales (solo a>=0)
    mask_pos = (a >= 0);
    if any(mask_pos)
        gini_inf(jj) = gini_from_pdf(a(mask_pos), g(mask_pos,1), da);
    else
        gini_inf(jj) = NaN;
    end

    % -------- Agregados fiscales ex post (este tipo) --------
    popI = la2 / (la1 + la2);
    popF = la1 / (la1 + la2);

    % Ingreso por impuesto a la renta:
    rev_income_j = cfg.taxF*z2*popF + cfg.taxI*z1*popI;

    % Ingreso por IVA (integral sobre consumo)
    cons_mass = (c(:,1).*g(:,1) + c(:,2).*g(:,2));   % c(a,s)*pdf(a,s)
    rev_tauc_j = cfg.tauc * sum(cons_mass) * da;

    % Costo de transferencias
    transfers_j = Transfer * popI;

    % Remanente a bien publico (ex post)
    gov_pub_j   = rev_income_j + rev_tauc_j - transfers_j;
    if cfg.gov_nonneg, gov_pub_j = max(0, gov_pub_j); end

    % Superavit (si negativo, deficit)
    gov_surplus_j = rev_income_j + rev_tauc_j - transfers_j - gov_pub_j;

    % Acumular para promedio por escenario
    rev_income_acc  = rev_income_acc  + rev_income_j;
    rev_tauc_acc    = rev_tauc_acc    + rev_tauc_j;
    transfers_acc   = transfers_acc   + transfers_j;
    gov_pub_acc     = gov_pub_acc     + gov_pub_j;
    gov_surplus_acc = gov_surplus_acc + gov_surplus_j;
end

% -------- Gini total (sobre promedio de tipos, a>=0) --------
mask_pos = (a >= 0);
G_total = NaN;
if any(mask_pos)
    % promedio de distribuciones totales (informal+formal) a traves de tipos
    Gavg = zeros(sum(mask_pos),1);
    for jj = 1:J
        Gavg = Gavg + (g_opt{jj}(mask_pos,1) + g_opt{jj}(mask_pos,2));
    end
    Gavg = Gavg / J;
    % normalizar a pdf
    massG = sum(Gavg)*da;
    if massG > 0
        Gavg = Gavg / massG;
        G_total = gini_from_pdf(a(mask_pos), Gavg, da);
    end
end

% empaquetar estadisticos
stats = struct();
stats.cmean_inf  = cmean_inf;
stats.cmean_for  = cmean_for;
stats.frac_bound = frac_bound;
stats.gini_inf   = gini_inf;
stats.gini_total = G_total;

% promedios fiscales por escenario (sobre tipos)
stats.rev_income  = rev_income_acc  / J;
stats.rev_tauc    = rev_tauc_acc    / J;
stats.transfers   = transfers_acc   / J;
stats.gov_pub     = gov_pub_acc     / J;
stats.gov_surplus = gov_surplus_acc / J;

end

%% =========================
%   Funcion auxiliar Gini
%% =========================
function G = gini_from_pdf(x, pdf, da)
% Calcula Gini para variable x>=0 con densidad pdf(x) (curva de Lorenz).
w = pdf * da;
mass = sum(w);
if mass <= 0
    G = NaN; return;
end
w = w / mass;
mu = sum(w .* x);
if mu <= 0
    G = NaN; return;
end
[xs, idx] = sort(x);
ws = w(idx);
ys = xs .* ws;
P = cumsum(ws);
L = cumsum(ys) / sum(ys);
A = trapz(P, L);
G = 1 - 2*A;
end
