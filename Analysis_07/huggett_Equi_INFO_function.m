function [r_opt, ir, pop1_vector, statsMatrix, statsCMatrix, ...
          GDistribution, a, Distribution, informal, formal, ...
          C, Saving, CDistribution, GCdistribution, I, amin] = ...
         huggett_Equi_INFO_function(RRA, eta_vector)
% -------------------------------------------------------------------------
% HUGGETT con 2 tipos (Informal/Formal) — MODELO BASE
% - Impuesto laboral tau_l (solo formales)
% - IVA tau_c (ambos)
% - Transferencias phi*y1 (solo informales)
% - Bien público Gov entra a la utilidad (aditivo)
% - Prima theta: informales con a<0 enfrentan r+theta
% - Deuda pública B(r) = (Tl + Tc - G - Tr) / r  y clearing: ∫ a g(a) da = B(r)
% Salidas compatibles con tus scripts legacy (C, Saving, Distribution, etc.)
% -------------------------------------------------------------------------

tic;

%% --- Política y calibración base (LATAM emergente)
tau_l  = 0.15;    % impuesto laboral (formales)
tau_c  = 0.10;    % IVA (ambos)
phi    = 0.15;    % transferencias (solo informales, proporcionales a y1)
Gov    = 0.07;    % bien público (entra a utilidad como shifter)
theta  = 0.02;    % prima de endeudamiento informal (a<0)

% Preferencias
sI = RRA;         % RRA informal  (puedes diferenciarlos si quieres)
sF = RRA;         % RRA formal

rho = 0.05;       % tasa de descuento

% Ingresos calibrados
z1 = 0.33;        % informal
z2 = 1.00;        % formal
z  = [z1, z2];

%% --- Discretización
r0   = 0.03; rmin = 0.005; rmax = 0.08;
I    = 700;
amin = -0.3*z1; 
amax = 5.0;
a    = linspace(amin,amax,I)'; 
da   = (amax-amin)/(I-1);
aa   = [a,a];
zz   = ones(I,1)*z;

maxit = 100; crit = 1e-6; Delta = 1000;
Ir    = 1000; crit_S = 1e-5;

%% --- Estructura ocupacional (misma lógica legacy: fija p22 y ajusta lambda1 por eta)
p22 = 0.8155;             % Prob[formal->formal]
la2 = -log(p22);          % intensidad F->I

pop1_vector = eta_vector(:)';  % (% informal)

% Prealocar salidas
Distribution    = cell(1, numel(pop1_vector));
GDistribution   = cell(1, numel(pop1_vector));
C               = cell(1, numel(pop1_vector));
Saving          = cell(1, numel(pop1_vector));
CDistribution   = cell(1, numel(pop1_vector));
GCdistribution  = cell(1, numel(pop1_vector));
informal        = cell(1, numel(pop1_vector));
formal          = cell(1, numel(pop1_vector));

% statsMatrix (riqueza): 
% [gmean_inf gmean_for gmed_inf gmed_for gmean_tot gmed_tot gini_inf gini_for gini_tot p11 r pp_const]
statsMatrix  = nan(numel(pop1_vector), 12);
% statsCMatrix (consumo): [Cmean_inf Cmean_for Cmed_inf Cmed_for Cmean_tot Cmed_tot]
statsCMatrix = nan(numel(pop1_vector), 6);

r_opt = nan(1, numel(pop1_vector));
ir    = NaN;

%% --- Loop en eta
for jj = 1:numel(pop1_vector)
    pop1 = pop1_vector(jj);

    % fijar lambda1 para alcanzar pop1 con lambda2 dado
    la1 = (1 - pop1) * la2 / max(pop1,1e-12);
    p11 = exp(-la1);
    la  = [la1, la2];

    % Matriz de cambio ocupacional
    Aswitch = [-speye(I)*la(1),  speye(I)*la(1);
                speye(I)*la(2), -speye(I)*la(2)];

    % -------- HJB + búsqueda de r --------
    r = r0;
    dVf = zeros(I,2); dVb = zeros(I,2);

    % rr(a) efectivo (informal paga theta si a<0)
    rr = r*ones(I,1); rr(a<0) = r + theta;

    % Guess V0 consistente con IVA y transferencias, y Gov en utilidad
    v0 = zeros(I,2);
    v0(:,1) = ((((1-0)*zz(:,1) + rr.*a + phi*z1)/(1+tau_c)).^(1-sI)/(1-sI) + Gov)/rho; % informal
    v0(:,2) = ((((1-tau_l)*zz(:,2) + r.*a          )/(1+tau_c)).^(1-sF)/(1-sF) + Gov)/rho; % formal

    for ir = 1:Ir
        if ir>1, v0 = V; end
        v = v0;

        for n=1:maxit
            V = v;

            % Derivadas
            dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVf(I,1) = ((1-0)*z1 + r*amax)^(-sI);
            dVf(I,2) = ((1-tau_l)*z2 + r*amax)^(-sF);

            dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
            dVb(1,1)   = ((1-0)*z1 + (r+theta)*amin)^(-sI);
            dVb(1,2)   = ((1-tau_l)*z2 + r*amin)^(-sF);

            % rr(a) actualizado
            rr = r*ones(I,1); rr(a<0) = r + theta;

            % Recursos netos (antes de IVA)
            res_inf = (1-0)*zz(:,1) + rr.*aa(:,1) + phi*z1;
            res_for = (1-tau_l)*zz(:,2) + r .*aa(:,2);

            % Euler con IVA
            cf = [max((1+tau_c)*dVf(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVf(:,2),1e-10).^(-1/sF)];
            cb = [max((1+tau_c)*dVb(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVb(:,2),1e-10).^(-1/sF)];

            % Drifts
            ssf = [res_inf, res_for] - (1+tau_c)*cf;
            ssb = [res_inf, res_for] - (1+tau_c)*cb;
            c0  = [res_inf, res_for];

            % Upwind
            If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
            c  = cf.*If + cb.*Ib + c0.*(I0/(1+tau_c));

            % Utilidad con bien público
            u = [c(:,1).^(1-sI)/(1-sI), c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);

            % Generador A con filas que suman cero
            X = max(-ssb,0)/da; Z = max(ssf,0)/da; X(1,:)=0; Z(I,:)=0; Y = -(X+Z);
            A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
            A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
            A  = A - spdiags(sum(A,2),0,2*I,2*I);   % fuerza sumas por fila = 0

            % Sistema lineal
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
        b0 = zeros(2*I,1);
        i_fix = 1; b0(i_fix)=.1;
        row = zeros(1,2*I); row(i_fix)=1; AT(i_fix,:) = row;

        gg = AT \ b0;
        g  = [gg(1:I), gg(I+1:2*I)];
        g  = g / (sum(g(:))*da);              % normaliza masa total = 1

        % Políticas de ahorro (adot)
        adot = [res_inf, res_for] - (1+tau_c)*c;

        % ---- Gobierno y B(r) ----
        popI = sum(g(:,1))*da; popF = sum(g(:,2))*da;
        C_tot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;
        Tc    = tau_c * C_tot;
        Tl    = tau_l * z2 * popF;
        TrAgg = phi   * z1 * popI;
        G_agg = Gov * (popI + popF);

        Btarget = (Tl + Tc - G_agg - TrAgg) / max(r,1e-6);
        % rB    = r * Btarget; % (si lo necesitas externamente)

        % Clearing de activos con deuda pública
        S = g(:,1)'*a*da + g(:,2)'*a*da - Btarget;

        % Búsqueda de r
        if S >  crit_S
            rmax = r; r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r; r = 0.5*(r + rmax);
        else
            break
        end
    end % ir

    r_opt(jj) = r;

    %% ---- Almacenar (compatibles con tus scripts) ----
    Distribution{jj}   = g;                   % densidades por tipo (a)
    GDistribution{jj}  = g(:,1) + g(:,2);     % total
    C{jj}              = c;                   % políticas de consumo por tipo
    Saving{jj}         = adot;                % políticas de ahorro por tipo

    % Proxy para distribuciones en eje c (compatibilidad legacy)
    CDistribution{jj}  = g;
    GCdistribution{jj} = GDistribution{jj};

    informal{jj} = g(:,1);
    formal{jj}   = g(:,2);

    %% ---- Stats riqueza (medias/medianas/gini) ----
    gmean_inf = sum(a.*g(:,1))/max(popI,eps);
    gmean_for = sum(a.*g(:,2))/max(popF,eps);
    gmean_tot = sum(a.*(g(:,1)+g(:,2)))/max(popI+popF,eps);
    gmed_inf  = weightedMedian(a,g(:,1));
    gmed_for  = weightedMedian(a,g(:,2));
    gmed_tot  = weightedMedian(a,g(:,1)+g(:,2));
    gini_inf  = giniWeighted(a,g(:,1));
    gini_for  = giniWeighted(a,g(:,2));
    gini_tot  = giniWeighted(a,g(:,1)+g(:,2));
    pp_const  = (g(1,1)+g(1,2)) / max(sum(g(:)),eps);  % masa en a=amin

    statsMatrix(jj,:) = [gmean_inf gmean_for gmed_inf gmed_for ...
                         gmean_tot gmed_tot gini_inf gini_for gini_tot ...
                         p11 r pp_const];

    %% ---- Stats consumo ----
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
