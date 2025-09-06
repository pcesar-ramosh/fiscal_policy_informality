function [r_opt, ir, pop1_vector, statsMatrix, statsCMatrix, GDistribution, a, Distribution, Fiscal, Cpolicies, Spolicies, Params] = ...
         huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1)
% - IVA tau_c en BC/Euler
% - Impuesto laboral tau_l solo formales (BC)
% - Transferencias phi*y_i solo informales (BC)
% - Bien público Gov en utilidad (exógeno)
% - r endógeno: clearing con B(r) = (Tl + Tc - G - Tr)/r
% - Si eta_vector(jj) es NaN => lambdas dependen de sI (composición endógena)
%   Si eta_vector(jj) en [0,1] => fijamos composición estacionaria pop1 = eta
% Salidas: distribuciones, políticas, fiscal, Gini por tipo y total, 
%          Params: theta,tau_c,tau_l,phi,Gov,p11_vec,p22_vec.

tic;

% ---------- Política fiscal ----------
tau_l  = 0.15;     % labor tax (formales)
tau_c  = 0.10;     % VAT (ambos)
phi    = 0.15;     % transferencias a informales

% ---------- Bien público en utilidad ----------
Gov   = 0.07;

% ---------- Ingresos ----------
z1 = 0.33; z2 = 1.00; z = [z1, z2];

% ---------- Preferencias / numéricos ----------
rho  = 0.05;
r0   = 0.03; rmin = 0.005; rmax = 0.08;

I    = 700;  amin = -0.3*z1; amax = 5.0;
a    = linspace(amin, amax, I)';  da = (amax-amin)/(I-1);
aa   = [a, a];  zz = ones(I,1)*z;

maxit = 100; crit = 1e-6; Delta = 1000;

pop1_vector = eta_vector;

Distribution  = cell(1, numel(pop1_vector));
GDistribution = cell(1, numel(pop1_vector));
Fiscal        = cell(1, numel(pop1_vector));
Cpolicies     = cell(1, numel(pop1_vector));
Spolicies     = cell(1, numel(pop1_vector));

Ir = 1000; crit_S = 1e-5;

r_opt = nan(1, numel(pop1_vector));
% [gmean_inf gmean_for gmed_inf gmed_for gmean_tot gmed_tot gini_inf gini_for gini_tot p11]
statsMatrix  = nan(numel(pop1_vector), 10);
statsCMatrix = nan(numel(pop1_vector), 12);

p11_vec = nan(1, numel(pop1_vector));
p22_vec = nan(1, numel(pop1_vector));

% --- calibración base de lambdas y dependencia en sI cuando eta es NaN
p22_base = 0.8155;                % base: formal permanece formal
la2_base = -log(p22_base);
s_ref    = median(sI_vector1);
k2 = 0.10; k1 = 0.20;             % elasticidades (tuneables)

for jj = 1:numel(pop1_vector)
    eta_target = pop1_vector(jj);
    sI = sI_vector1(jj); sF = sF_vector1(jj);

    % ---- Lambda(sI)
    if isnan(eta_target)
        % composición endógena: ambas intensidades dependen de sI
        la2 = la2_base * exp(k2*(sI - s_ref));       % F->I
        la1 = la2 * exp(k1*(sI - s_ref));            % I->F
        p22 = exp(-la2);
        p11 = exp(-la1);
    else
        % composición objetivo: fija pop1=eta_target
        la2 = la2_base;                  p22 = exp(-la2);
        la1 = (1 - eta_target) * la2 / max(eta_target,1e-12);
        p11 = exp(-la1);
    end
    la  = [la1, la2];

    % ---- Prima de endeudamiento en informales ----
    theta = 0.02;

    % ---- Guess V ----
    r = r0; rr = zeros(I,1);
    for i=1:I, rr(i) = (a(i)>=0)*r + (a(i)<0)*(r+theta); end

    Transfer = phi * z1;

    Aswitch = [-speye(I)*la(1),  speye(I)*la(1);
                speye(I)*la(2), -speye(I)*la(2)];

    % Guess para V
    v0 = zeros(I,2);
    v0(:,1) = ((((1-0)*zz(:,1) + rr.*a + Transfer)/(1+tau_c)).^(1-sI)/(1-sI) + Gov)/rho;
    v0(:,2) = ((((1-tau_l)*zz(:,2) + r.*a          )/(1+tau_c)).^(1-sF)/(1-sF) + Gov)/rho;

    dVf = zeros(I,2); dVb = zeros(I,2);

    for ir = 1:Ir
        if ir>1, v0 = V; end
        v = v0;

        for n = 1:maxit
            V = v;

            % Derivadas
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVf(I,1) = ((1-0)*z1 + r*amax)^(-sI);
            dVf(I,2) = ((1-tau_l)*z2 + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVb(1,1)   = ((1-0)*z1 + (r+theta)*amin)^(-sI);
            dVb(1,2)   = ((1-tau_l)*z2 + r*amin)^(-sF);

            % rr(a)
            for i=1:I, if a(i)>=0, rr(i)=r; else, rr(i)=r+theta; end; end

            % Euler con IVA
            cf = [max((1+tau_c)*dVf(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVf(:,2),1e-10).^(-1/sF)];
            cb = [max((1+tau_c)*dVb(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVb(:,2),1e-10).^(-1/sF)];

            res_inf = (1-0)*zz(:,1) + rr.*aa(:,1) + Transfer;
            res_for = (1-tau_l)*zz(:,2) + r .*aa(:,2);

            ssf = [res_inf, res_for] - (1+tau_c)*cf;
            ssb = [res_inf, res_for] - (1+tau_c)*cb;
            c0  = [res_inf, res_for];

            If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
            c  = cf.*If + cb.*Ib + c0.*(I0/(1+tau_c));

            u = [c(:,1).^(1-sI)/(1-sI), c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);

            % Generador A (filas suman 0)
            X = max(-ssb,0)/da; Z = max(ssf,0)/da; X(1,:)=0; Z(I,:)=0; Y = -(X+Z);
            A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
            A  = A - spdiags(sum(A,2),0,2*I,2*I);

            Bmat = (1/Delta + rho)*speye(2*I) - A;

            u_stacked = [u(:,1); u(:,2)];
            V_stacked = [V(:,1); V(:,2)];
            V_stacked = Bmat \ (u_stacked + V_stacked/Delta);
            V = [V_stacked(1:I), V_stacked(I+1:2*I)];

            if max(max(abs(V - v))) < crit, break; end
            v = V;
        end

        % Fokker-Planck
        AT = A'; b0 = zeros(2*I,1); i_fix = 1; b0(i_fix)=.1;
        row = zeros(1,2*I); row(i_fix)=1; AT(i_fix,:) = row;
        gg = AT \ b0;
        g  = [gg(1:I), gg(I+1:2*I)];
        g  = g / ((sum(g(:))*da));  % normaliza
        adot = [res_inf, res_for] - (1+tau_c)*c;

        % Presupuesto público y B(r)
        popI = sum(g(:,1))*da; popF = sum(g(:,2))*da;
        C_tot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;
        Tc    = tau_c * C_tot;
        Tl    = tau_l * z2 * popF;
        TrAgg = phi   * z1 * popI;
        Y     = z1*popI + z2*popF;
        G_agg = Gov * (popI + popF);
        Btarget = (Tl + Tc - G_agg - TrAgg) / max(r,1e-6);
        rB      = r * Btarget;
        Gap     = (G_agg + TrAgg + rB) - (Tl + Tc);

        % Clearing de activos
        S = g(:,1)'*a*da + g(:,2)'*a*da - Btarget;
        if S >  1e-5, rmax=r; r=0.5*(r+rmin);
        elseif S < -1e-5, rmin=r; r=0.5*(r+rmax);
        else, break; end
    end

    r_opt(jj) = r;

    % Guardar
    Distribution{jj}  = g;
    GDistribution{jj} = g(:,1)+g(:,2);
    Cpolicies{jj}     = c;
    Spolicies{jj}     = adot;

    % Wealth stats
    gmean_inf = sum(a.*g(:,1))/max(popI,eps);
    gmean_for = sum(a.*g(:,2))/max(popF,eps);
    gmean_tot = sum(a.*(g(:,1)+g(:,2)))/max(popI+popF,eps);
    gmed_inf  = weightedMedian(a,g(:,1));
    gmed_for  = weightedMedian(a,g(:,2));
    gmed_tot  = weightedMedian(a,g(:,1)+g(:,2));
    gini_inf  = giniWeighted(a,g(:,1));
    gini_for  = giniWeighted(a,g(:,2));
    gini_tot  = giniWeighted(a,g(:,1)+g(:,2));
    statsMatrix(jj,:) = [gmean_inf gmean_for gmed_inf gmed_for gmean_tot gmed_tot gini_inf gini_for gini_tot p11];

    % Consumo stats
    c_inf = c(:,1); c_for = c(:,2);
    Cmean_inf = sum(c_inf.*g(:,1))/max(popI,eps);
    Cmean_for = sum(c_for.*g(:,2))/max(popF,eps);
    Cmean_tot = sum((c_inf+c_for).*(g(:,1)+g(:,2)))/max(popI+popF,eps);
    Cmed_inf  = weightedMedian(c_inf,g(:,1));
    Cmed_for  = weightedMedian(c_for,g(:,2));
    Cmed_tot  = weightedMedian([c_inf;c_for],[g(:,1);g(:,2)]);
    statsCMatrix(jj,1:6) = [Cmean_inf Cmean_for Cmed_inf Cmed_for Cmean_tot Cmed_tot];

    Fiscal{jj} = [Tc Tl TrAgg G_agg rB Gap popI popF Y Btarget];

    p11_vec(jj) = p11;
    p22_vec(jj) = p22;
end

Params = struct('theta',theta,'tau_c',tau_c,'tau_l',tau_l,'phi',phi,'Gov',Gov, ...
                'p11_vec',p11_vec,'p22_vec',p22_vec);

toc;
end

% -------- Helpers --------
function m = weightedMedian(x,w)
w = w(:); x = x(:); W = sum(w); if W<=0, m=NaN; return; end
[xx,idx]=sort(x); ww=w(idx); cw=cumsum(ww); k=find(cw>=0.5*W,1,'first'); m=xx(k);
end

function gini = giniWeighted(x,w)
x=x(:); w=w(:); W=sum(w); if W<=0, gini=NaN; return; end
xshift = x - min(0,min(x)) + 1e-12; [xx,idx]=sort(xshift); ww=w(idx)/W;
cumw=cumsum(ww); cumxw=cumsum(xx.*ww);
if cumxw(end)<=0, gini=NaN; return; end
L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end
