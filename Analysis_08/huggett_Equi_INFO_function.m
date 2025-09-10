function [r_opt, ir, pop1_vector, statsMatrix, statsCMatrix, ...
          GDistribution, a, Distribution, informal, formal, ...
          C, Saving, CDistribution, GCdistribution, I, amin] = ...
         huggett_Equi_INFO_function(RRA, eta_vector, lambdaSpec)
% -------------------------------------------------------------------------
% HUGGETT 2 tipos (Informal/Formal) — MODELO BASE + LAMBDA DINÁMICO/ENDÓGENO
% - tau_l (formales), tau_c (ambos), transferencias phi*y1 (informales),
%   bien público Gov en utilidad (aditivo), prima theta para informales con a<0
% - Deuda pública B(r) cierra presupuesto: B = (T_l + T_c - G - Tr)/r
% - (NUEVO) Intensidades de transición:
%     * Por defecto (sin lambdaSpec):  λ2 = const, λ1 = (1-eta)/eta * λ2
%     * Con lambdaSpec.mode='custom':  usa funciones λ1(a,eta), λ2(a,eta)
%
% I/O compatible con tus scripts legacy.
% -------------------------------------------------------------------------

tic;

%% --- Política/calibración base
tau_l  = 0.15;    % impuesto laboral (solo formales)
tau_c  = 0.10;    % IVA (ambos)
phi    = 0.15;    % transferencias (solo informales; proporcionales a y1)
Gov    = 0.07;    % bien público (entra a utilidad)
theta  = 0.02;    % prima de endeudamiento informal (a<0)

% Preferencias
sI = RRA; sF = RRA;
rho = 0.05;

% Ingresos
z1 = 0.33; z2 = 1.00; z = [z1,z2];

%% --- Discretización
r0   = 0.03; rmin = 0.005; rmax = 0.08;
I    = 700;
amin = -0.30*z1;
amax = 5.0;
a    = linspace(amin,amax,I)'; 
da   = (amax-amin)/(I-1);
aa   = [a,a];
zz   = ones(I,1)*z;

maxit = 100; crit = 1e-6; Delta = 1000;
Ir    = 1000; crit_S = 1e-5;

%% --- Lambda base (si no viene lambdaSpec)
useCustomLambda = (nargin>=3) && isstruct(lambdaSpec) && isfield(lambdaSpec,'mode') ...
                  && strcmpi(lambdaSpec.mode,'custom');

% p22 base (solo se usa si no hay custom lambda)
p22_bar = 0.8155;             
la2_bar = -log(p22_bar);      

% Prealocar
Distribution    = cell(1, numel(eta_vector));
GDistribution   = cell(1, numel(eta_vector));
C               = cell(1, numel(eta_vector));
Saving          = cell(1, numel(eta_vector));
CDistribution   = cell(1, numel(eta_vector));
GCdistribution  = cell(1, numel(eta_vector));
informal        = cell(1, numel(eta_vector));
formal          = cell(1, numel(eta_vector));
statsMatrix     = nan(numel(eta_vector), 12);
statsCMatrix    = nan(numel(eta_vector), 6);
r_opt           = nan(1, numel(eta_vector));
pop1_vector     = nan(1, numel(eta_vector)); % (NUEVO) llenamos con popI endógeno

for jj = 1:numel(eta_vector)
    eta_target = eta_vector(jj);

    %% ----- λ1(a), λ2(a) -----
    if useCustomLambda
        % Esperamos handles: la1_fun(a, eta_guess), la2_fun(a, eta_guess)
        if ~isfield(lambdaSpec,'la1_fun') || ~isfield(lambdaSpec,'la2_fun')
            error('lambdaSpec.custom requiere la1_fun y la2_fun');
        end
        L1_fun = @(avec,eta_guess) max(lambdaSpec.la1_fun(avec, eta_guess), 1e-10);
        L2_fun = @(avec,eta_guess) max(lambdaSpec.la2_fun(avec, eta_guess), 1e-10);
        % Inicializa con eta_guess = eta_target (solo como semilla)
        L1_vec = L1_fun(a, eta_target);
        L2_vec = L2_fun(a, eta_target);
        p11_rep = exp(-mean(L1_vec)); % resumen para stats (prob permanecer informal)
    else
        % λ constantes: λ2 fijo y λ1 ajustado a eta_target
        la2 = la2_bar;
        la1 = (1-eta_target)*la2 / max(eta_target,1e-12);
        L1_vec = la1*ones(I,1);
        L2_vec = la2*ones(I,1);
        p11_rep = exp(-la1);
    end

    %% ----- Generador de switches -----
    Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
                 spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];

    %% -------- HJB + búsqueda de r --------
    r = r0;
    dVf = zeros(I,2); dVb = zeros(I,2);

    rr = r*ones(I,1); rr(a<0) = r + theta;

    v0 = zeros(I,2);
    % IVA se implementa en Euler; acá solo un guess razonable:
    v0(:,1) = ((((1-0)*zz(:,1) + rr.*a + phi*z1)/(1+tau_c)).^(1-sI)/(1-sI) + Gov)/rho;
    v0(:,2) = ((((1-tau_l)*zz(:,2) + r .*a          )/(1+tau_c)).^(1-sF)/(1-sF) + Gov)/rho;

    for ir = 1:Ir
        if ir>1, v0 = V; end
        v = v0;

        for n=1:maxit
            V = v;

            % derivadas
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVf(I,1) = ((1-0)*z1 + r*amax)^(-sI);
            dVf(I,2) = ((1-tau_l)*z2 + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVb(1,1)   = ((1-0)*z1 + (r+theta)*amin)^(-sI);
            dVb(1,2)   = ((1-tau_l)*z2 + r*amin)^(-sF);

            rr = r*ones(I,1); rr(a<0) = r + theta;

            res_inf = (1-0)*zz(:,1) + rr.*aa(:,1) + phi*z1;
            res_for = (1-tau_l)*zz(:,2) + r .*aa(:,2);

            cf = [max((1+tau_c)*dVf(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVf(:,2),1e-10).^(-1/sF)];
            cb = [max((1+tau_c)*dVb(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVb(:,2),1e-10).^(-1/sF)];

            ssf = [res_inf, res_for] - (1+tau_c)*cf;
            ssb = [res_inf, res_for] - (1+tau_c)*cb;
            c0  = [res_inf, res_for];

            If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
            c  = cf.*If + cb.*Ib + c0.*(I0/(1+tau_c));

            u = [c(:,1).^(1-sI)/(1-sI), c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);

            % generador con filas suma cero
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

        % ---- Fokker–Planck ----
        AT = A';
        b0 = zeros(2*I,1); i_fix = 1; b0(i_fix)=.1;
        row = zeros(1,2*I); row(i_fix)=1; AT(i_fix,:) = row;

        gg = AT \ b0;
        g  = [gg(1:I), gg(I+1:2*I)];
        g  = g / (sum(g(:))*da);

        % políticas de ahorro
        adot = [res_inf, res_for] - (1+tau_c)*c;

        % ---- Gobierno y B(r) ----
        popI = sum(g(:,1))*da; popF = sum(g(:,2))*da;
        C_tot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;
        Tc    = tau_c * C_tot;
        Tl    = tau_l * z2 * popF;
        TrAgg = phi   * z1 * popI;
        G_agg = Gov * (popI + popF);

        Btarget = (Tl + Tc - G_agg - TrAgg) / max(r,1e-9);

        S = g(:,1)'*a*da + g(:,2)'*a*da - Btarget;

        % si lambda custom depende de eta, actualiza Aswitch con popI actual
        if useCustomLambda
            L1_vec = L1_fun(a, popI);
            L2_vec = L2_fun(a, popI);
            Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
                         spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];
        end

        % bisección en r
        if S >  crit_S
            rmax = r; r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r; r = 0.5*(r + rmax);
        else
            break
        end
    end % ir

    r_opt(jj) = r;
    pop1_vector(jj) = popI; % (NUEVO) devolver composición endógena

    % ---- Salidas para tus scripts ----
    Distribution{jj}   = g;
    GDistribution{jj}  = g(:,1) + g(:,2);
    C{jj}              = c;
    Saving{jj}         = adot;
    CDistribution{jj}  = g;                  % proxy para ejes en c
    GCdistribution{jj} = GDistribution{jj};
    informal{jj}       = g(:,1);
    formal{jj}         = g(:,2);

    % ---- Stats riqueza ----
    gmean_inf = sum(a.*g(:,1))/max(popI,eps);
    gmean_for = sum(a.*g(:,2))/max(popF,eps);
    gmean_tot = sum(a.*(g(:,1)+g(:,2)))/max(popI+popF,eps);
    gmed_inf  = weightedMedian(a,g(:,1));
    gmed_for  = weightedMedian(a,g(:,2));
    gmed_tot  = weightedMedian(a,g(:,1)+g(:,2));
    gini_inf  = giniWeighted(a,g(:,1));
    gini_for  = giniWeighted(a,g(:,2));
    gini_tot  = giniWeighted(a,g(:,1)+g(:,2));
    pp_const  = (g(1,1)+g(1,2)) / max(sum(g(:)),eps);

    statsMatrix(jj,:) = [gmean_inf gmean_for gmed_inf gmed_for ...
                         gmean_tot gmed_tot gini_inf gini_for gini_tot ...
                         p11_rep r pp_const];

    % ---- Stats consumo ----
    c_inf = c(:,1); c_for = c(:,2);
    Cmean_inf = sum(c_inf.*g(:,1))/max(popI,eps);
    Cmean_for = sum(c_for.*g(:,2))/max(popF,eps);
    Cmean_tot = sum((c_inf+c_for).*(g(:,1)+g(:,2)))/max(popI+popF,eps);
    Cmed_inf  = weightedMedian(c_inf,g(:,1));
    Cmed_for  = weightedMedian(c_for,g(:,2));
    Cmed_tot  = weightedMedian([c_inf;c_for],[g(:,1);g(:,2)]);
    statsCMatrix(jj,:) = [Cmean_inf Cmean_for Cmed_inf Cmed_for Cmean_tot Cmed_tot];

end % jj

toc;
end

% ---------- Helpers ----------
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
