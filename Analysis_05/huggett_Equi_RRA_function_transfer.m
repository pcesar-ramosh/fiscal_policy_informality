function [r_opt, ir, pop1_vector, statsMatrix, statsCMatrix, GDistribution, a, Distribution, Fiscal, Cpolicies, Spolicies] = ...
         huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1)
% Baseline con:
%  - IVA tau_c en BC y Euler
%  - Impuesto laboral tau_l solo formales (BC)
%  - Transferencias phi*y_i solo informales (BC)
%  - Bien público Gov en utilidad (exógeno baseline)
%  - r endógeno: clearing de activos con B(r) = (Tl+Tc - G - Tr)/r
%
% Salidas extra:
%  - Cpolicies{jj} : c(a) óptima [I x 2]
%  - Spolicies{jj} : \dot a(a)   [I x 2]
%  - Fiscal{jj}    : [Tc Tl Tr G rB Gap popI popF Y Btarget_eq]

tic;

% ---------- Política fiscal (baseline) ----------
tau_l  = 0.15;     % impuesto laboral (formales)
tau_c  = 0.10;     % IVA (ambos)
phi    = 0.15;     % transferencias a informales

% ---------- Bien público en utilidad ----------
Gov   = 0.07;      % shifter (proxy de psi*log G)

% ---------- Ingresos ----------
z1 = 0.33;   % informal
z2 = 1.00;   % formal
z  = [z1, z2];

% ---------- Preferencias / descuento ----------
rho  = 0.05;       % tasa de descuento (continua, esquema HJB)

% ---------- Grilla de activos y controles ----------
r0   = 0.03;       % guess inicial r
rmin = 0.005;      % evita dividir por ~0 en B(r)
rmax = 0.08;

I    = 700;        
amin = -0.3*z1;    % límite de deuda
amax = 5.0;
a    = linspace(amin, amax, I)';
da   = (amax - amin)/(I - 1);
aa   = [a, a];
zz   = ones(I,1)*z;

maxit = 100;
crit  = 1e-6;
Delta = 1000;

pop1_vector = eta_vector;   % share informal target

% Storage
Distribution  = cell(1, numel(pop1_vector));
GDistribution = cell(1, numel(pop1_vector));
Fiscal        = cell(1, numel(pop1_vector));
Cpolicies     = cell(1, numel(pop1_vector));
Spolicies     = cell(1, numel(pop1_vector));

Ir = 1000;
crit_S = 1e-5;

r_opt = nan(1, numel(pop1_vector));
% statsMatrix cols (10):
% [gmean_inf gmean_for gmed_inf gmed_for gmean_tot gmed_tot gini_inf gini_for gini_tot p11]
statsMatrix  = nan(numel(pop1_vector), 10);
statsCMatrix = nan(numel(pop1_vector), 12);  % (no usamos todas las cols)

for jj = 1:numel(pop1_vector)
    pop1 = pop1_vector(jj);

    % ---- Intensidades para empatar pop estacionaria ----
    p22 = 0.8155;               
    la2 = -log(p22);            % formal->informal
    la1 = (1 - pop1) * la2 / pop1;  
    p11 = exp(-la1);                
    la  = [la1, la2];

    % ---- RRA ----
    sI = sI_vector1(jj);
    sF = sF_vector1(jj);

    % ---- Prima de endeudamiento en informales ----
    theta = 0.02;               

    % ---- Guess V ----
    r = r0;
    rr = zeros(I,1);
    for i=1:I
        rr(i) = (a(i)>=0) * r + (a(i)<0) * (r + theta);
    end

    Transfer = phi * z1;

    % Matriz de cambio entre ocupaciones
    Aswitch = [-speye(I)*la(1),  speye(I)*la(1);
                speye(I)*la(2), -speye(I)*la(2)];

    % Guess para V (consumo estático con IVA)
    v0 = zeros(I,2);
    v0(:,1) = ((((1-0)*zz(:,1) + rr.*a + Transfer)/(1+tau_c)).^(1-sI)/(1-sI) + Gov)/rho;  % informal
    v0(:,2) = ((((1-tau_l)*zz(:,2) + r.*a          )/(1+tau_c)).^(1-sF)/(1-sF) + Gov)/rho; % formal

    dVf = zeros(I,2); dVb = zeros(I,2);

    for ir = 1:Ir
        if ir>1, v0 = V; end
        v = v0;

        for n = 1:maxit
            V = v;

            % --- Diferencias para V_a ---
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVf(I,1) = ((1-0)*z(1) + (r)*amax)^(-sI);
            dVf(I,2) = ((1-tau_l)*z(2) + (r)*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVb(1,1)   = ((1-0)*z(1) + (r+theta)*amin)^(-sI);
            dVb(1,2)   = ((1-tau_l)*z(2) + (r)*amin)^(-sF);

            % --- Actualizar rr(a) ---
            for i=1:I
                if a(i)>=0, rr(i) = r; else, rr(i) = r + theta; end
            end

            % ---------- Euler con IVA ----------
            cf = [max((1+tau_c)*dVf(:,1),1e-10).^( -1/sI), max((1+tau_c)*dVf(:,2),1e-10).^( -1/sF)];
            cb = [max((1+tau_c)*dVb(:,1),1e-10).^( -1/sI), max((1+tau_c)*dVb(:,2),1e-10).^( -1/sF)];

            % Recursos por tipo (phi solo en BC informal)
            res_inf = (1-0)*zz(:,1) + rr.*aa(:,1) + Transfer;   % informal
            res_for = (1-tau_l)*zz(:,2) + r .*aa(:,2);          % formal

            % Drifts (el gasto en caja incluye (1+tau_c))
            ssf = [res_inf, res_for] - (1+tau_c)*cf;
            ssb = [res_inf, res_for] - (1+tau_c)*cb;
            c0  = [res_inf, res_for];   % si I0=1 => c = c0/(1+tau_c)

            % Upwind
            If = ssf > 0;
            Ib = ssb < 0;
            I0 = (1-If-Ib);

            c = cf.*If + cb.*Ib + c0.*(I0/(1+tau_c));

            % Utilidad (Gov solo como shifter)
            u = [c(:,1).^(1-sI)/(1-sI), c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);

            % ---------- Generador A (robusto: filas suman 0) ----------
            X = max(-ssb,0)/da;   % flujo a a-1
            Z = max( ssf,0)/da;   % flujo a a+1
            X(1,:) = 0;  Z(I,:) = 0;
            Y = -(X + Z);

            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

            rowsum = sum(A,2);
            A = A - spdiags(rowsum,0,2*I,2*I);

            B = (1/Delta + rho)*speye(2*I) - A;

            u_stacked = [u(:,1); u(:,2)];
            V_stacked = [V(:,1); V(:,2)];

            b_rhs = u_stacked + V_stacked/Delta;
            V_stacked = B \ b_rhs;
            V = [V_stacked(1:I), V_stacked(I+1:2*I)];

            Vchange = V - v;
            v = V;

            if max(max(abs(Vchange))) < crit
                break;
            end
        end % end V-iteration

        % ---------- Distribución estacionaria (Fokker-Planck) ----------
        AT = A';
        b0 = zeros(2*I,1);
        i_fix = 1; b0(i_fix) = .1;
        row = zeros(1,2*I); row(i_fix) = 1;
        AT(i_fix,:) = row;
        gg = AT \ b0;
        g_sum = (gg' * ones(2*I,1)) * da;
        gg = gg ./ g_sum;
        g  = [gg(1:I), gg(I+1:2*I)];

        % Políticas (flujo de activos)
        adot = [res_inf, res_for] - (1+tau_c)*c;

        % ---------- Presupuesto público y B(r) ----------
        popI = sum(g(:,1))*da;
        popF = sum(g(:,2))*da;
        C_tot = sum(c(:,1) .* g(:,1))*da + sum(c(:,2) .* g(:,2))*da;
        Tc    = tau_c * C_tot;                     % IVA
        Tl    = tau_l * z2 * popF;                 % renta laboral
        TrAgg = phi   * z1 * popI;                 % transferencias
        Y     = z1*popI + z2*popF;                 % output
        G_agg = Gov * (popI + popF);               % gasto público (proxy)

        % Oferta de bonos del gobierno endógena
        Btarget = (Tl + Tc - G_agg - TrAgg) / max(r,1e-6);
        rB      = r * Btarget;
        Gap     = (G_agg + TrAgg + rB) - (Tl + Tc);  % ~0 por construcción

        % ---------- Clearing de activos ----------
        S = g(:,1)'*a*da + g(:,2)'*a*da - Btarget;

        % Búsqueda de r
        if S >  crit_S
            rmax = r; r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r; r = 0.5*(r + rmax);
        else
            break;
        end
    end % end r-iteration

    r_opt(jj) = r;

    % ---------- Guardar distribuciones y políticas ----------
    Distribution{jj}  = g;
    GDistribution{jj} = g(:,1) + g(:,2);
    Cpolicies{jj}     = c;
    Spolicies{jj}     = adot;

    % ---------- Estadísticas (riqueza) ----------
    % Medias/medianas
    gmean_inf = sum(a .* g(:,1)) / max(popI,eps);
    gmean_for = sum(a .* g(:,2)) / max(popF,eps);
    gmean_tot = sum(a .* (g(:,1)+g(:,2))) / max(popI+popF,eps);
    gmed_inf  = weightedMedian(a, g(:,1));
    gmed_for  = weightedMedian(a, g(:,2));
    gmed_tot  = weightedMedian(a, g(:,1)+g(:,2));
    % Ginis por tipo y total (maneja valores negativos con shift interno)
    gini_inf  = giniWeighted(a, g(:,1));
    gini_for  = giniWeighted(a, g(:,2));
    gini_tot  = giniWeighted(a, g(:,1)+g(:,2));

    statsMatrix(jj,1:10) = [gmean_inf, gmean_for, gmed_inf, gmed_for, ...
                            gmean_tot, gmed_tot, gini_inf, gini_for, gini_tot, p11];

    % ---------- Estadísticas (consumo) ----------
    c_inf = c(:,1); c_for = c(:,2);
    Cmean_inf = sum(c_inf .* g(:,1)) / max(popI,eps);
    Cmean_for = sum(c_for .* g(:,2)) / max(popF,eps);
    Cmean_tot = sum((c_inf+c_for) .* (g(:,1)+g(:,2))) / max(popI+popF,eps);
    Cmed_inf  = weightedMedian(c_inf, g(:,1));
    Cmed_for  = weightedMedian(c_for, g(:,2));
    Cmed_tot  = weightedMedian([c_inf; c_for], [g(:,1); g(:,2)]);

    statsCMatrix(jj,1:6) = [Cmean_inf, Cmean_for, Cmed_inf, Cmed_for, Cmean_tot, Cmed_tot];

    % ---------- Resumen fiscal ----------
    Fiscal{jj} = [Tc, Tl, TrAgg, G_agg, rB, Gap, popI, popF, Y, Btarget];
end

toc;

end % end main function


% ===================== Helpers =====================

function m = weightedMedian(x, w)
w = w(:); x = x(:);
W = sum(w);
if W <= 0, m = NaN; return; end
[xx, idx] = sort(x);
ww = w(idx);
cw = cumsum(ww);
k  = find(cw >= 0.5*W, 1, 'first');
m = xx(k);
end

function gini = giniWeighted(x, w)
% Gini ponderado robusto a valores negativos: aplica un "shift" mínimo
% para asegurar no-negatividad y luego usa la fórmula estándar.
x = x(:); w = w(:);
W = sum(w);
if W <= 0, gini = NaN; return; end
xshift = x - min(0,min(x)) + 1e-12;  % asegura xshift>=0
[xx, idx] = sort(xshift);
ww = w(idx)/W;
cumw = cumsum(ww);
cumxw = cumsum(xx .* ww);
if cumxw(end) <= 0
    gini = NaN; return;   % evita división por ~0 (datos degenerados)
end
L = cumxw / cumxw(end);        % curva de Lorenz
gini = 1 - 2 * trapz(cumw, L); % 0..1
end
