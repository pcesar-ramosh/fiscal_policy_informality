function out = huggett_base_function(params)
% -------------------------------------------------------------------------
% Modelo Huggett con 2 tipos (Informal/Formal) + impuestos y transferencias
% - tau_l: impuesto laboral (solo formales)
% - tau_c: IVA (ambos)
% - phi:   transferencias a informales (phi * z1, por persona informal)
% - Gov:   bien público (aditivo en utilidad), provisto siempre
% - theta: prima de endeudamiento informal (tasa r+theta si a<0)
% - B(r) endógena cierra el presupuesto: B = (Tl + Tc - G - Tr)/r
% - Transiciones ocupacionales:
%     * Por defecto: λ2 constante a partir de p22_bar, y
%       λ1 = ((1-eta_target)/eta_target) * λ2  (fija informalidad objetivo).
%     * Con params.lambdaSpec.mode='custom': usa la1_fun(a,eta), la2_fun(a,eta).
%
% ENTRADA: struct params con (valores por defecto entre []):
%   RRA_I [2.3], RRA_F [2.3], rho [0.05], theta [0.02]
%   tau_l [0.15], tau_c [0.15], phi [0.10], Gov [0.05]
%   z1 [0.33], z2 [1.0]
%   I [700], amin [-0.30*z1], amax [5], r_guess [0.03], rmin [0.005], rmax [0.08]
%   p22_bar [0.8155], eta_target [0.54]
%   fix_r [0]  (0 = busca r por clearing; 1 = usa r fijo sin bisección)
%   lambdaSpec (opcional): struct .mode='custom', .la1_fun(a,eta), .la2_fun(a,eta)
%
% SALIDA: struct out con:
%   r, a, g (I x 2), c (I x 2), s (I x 2), popI, popF, Ctot, Y
%   fiscal: Tl, Tc, Tr, G, PB, B, rB, BB
%   stats: wealth/cons (medias, medianas, Gini por tipo y total), p11
%   borrowers: fracciones y montos (por tipo) para a<0 y a>0
% -------------------------------------------------------------------------

%% ----- Lectura de parámetros y defaults
arg = @(f,def) (isfield(params,f) && ~isempty(params.(f))) * params.(f) + ...
               ~(isfield(params,f) && ~isempty(params.(f))) * def;

sI    = arg('RRA_I', 2.3);
sF    = arg('RRA_F', 2.3);
rho   = arg('rho',   0.05);
theta = arg('theta', 0.02);

tau_l = arg('tau_l', 0.15);
tau_c = arg('tau_c', 0.15);
phi   = arg('phi',   0.10);
Gov   = arg('Gov',   0.05);

z1    = arg('z1',    0.33);
z2    = arg('z2',    1.00);
z     = [z1, z2];

I     = arg('I',     700);
amin  = arg('amin', -0.30*z1);
amax  = arg('amax',  5.0);
a     = linspace(amin, amax, I)'; 
da    = (amax-amin)/(I-1);
aa    = [a, a];
zz    = ones(I,1) * z;

r     = arg('r_guess', 0.03);
rmin  = arg('rmin',     0.005);
rmax  = arg('rmax',     0.08);
fix_r = arg('fix_r',    0);     % <---- NUEVO

p22_bar    = arg('p22_bar', 0.8155);
eta_target = arg('eta_target', 0.54);

useCustomLambda = (isfield(params,'lambdaSpec') && isstruct(params.lambdaSpec) && ...
                   isfield(params.lambdaSpec,'mode') && strcmpi(params.lambdaSpec.mode,'custom'));
if useCustomLambda
    la1_fun = params.lambdaSpec.la1_fun;  % la1_fun(a, eta_guess)
    la2_fun = params.lambdaSpec.la2_fun;  % la2_fun(a, eta_guess)
end

%% ----- Iteración HJB + búsqueda de r (o r fijo)
maxit = 200; crit = 1e-6; Delta = 1000;
Ir    = 1000; crit_S = 1e-5;

% Lambda inicial
if useCustomLambda
    L1_vec = max(la1_fun(a, eta_target), 1e-10);
    L2_vec = max(la2_fun(a, eta_target), 1e-10);
    p11_rep = exp(-mean(L1_vec));
else
    la2 = -log(p22_bar);
    la1 = (1-eta_target)*la2 / max(eta_target,1e-12);
    L1_vec = la1 * ones(I,1);
    L2_vec = la2 * ones(I,1);
    p11_rep = exp(-la1);
end

Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
             spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];

% Valor inicial (consumo de todo ingreso y Gov en flujo)
rr = r*ones(I,1); rr(a<0) = r + theta;
v0(:,1) = ((((1-0)*zz(:,1) + rr.*a + phi*z1)/(1+tau_c)).^(1-sI)/(1-sI) + Gov)/rho;
v0(:,2) = ((((1-tau_l)*zz(:,2) + r .*a          )/(1+tau_c)).^(1-sF)/(1-sF) + Gov)/rho;

for ir = 1:Ir
    v = v0;
    for n=1:maxit
        V = v;

        % derivadas forward/backward
        dVf = zeros(I,2); dVb = zeros(I,2);
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVf(I,1) = ((1-0)*z1 + r*amax)^(-sI);
        dVf(I,2) = ((1-tau_l)*z2 + r*amax)^(-sF);

        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVb(1,1)   = ((1-0)*z1 + (r+theta)*amin)^(-sI);
        dVb(1,2)   = ((1-tau_l)*z2 + r*amin)^(-sF);

        % tasas efectivas
        rr = r*ones(I,1); rr(a<0) = r + theta;

        % “recursos” (lado derecho de BC sin consumo)
        res_inf = (1-0)*zz(:,1) + rr.*aa(:,1) + phi*z1;
        res_for = (1-tau_l)*zz(:,2) + r .*aa(:,2);

        % consumo con IVA (upwind)
        cf = [max((1+tau_c)*dVf(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVf(:,2),1e-10).^(-1/sF)];
        cb = [max((1+tau_c)*dVb(:,1),1e-10).^(-1/sI), max((1+tau_c)*dVb(:,2),1e-10).^(-1/sF)];
        ssf = [res_inf, res_for] - (1+tau_c)*cf;
        ssb = [res_inf, res_for] - (1+tau_c)*cb;
        c0  = [res_inf, res_for];

        If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
        c  = cf.*If + cb.*Ib + c0.*(I0/(1+tau_c));  % drift≈0 -> c = res/(1+tau_c)

        u = [c(:,1).^(1-sI)/(1-sI), c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);

        % generador (suma de filas = 0)
        X = max(-ssb,0)/da;  Z = max(ssf,0)/da;
        X(1,:) = 0; Z(I,:) = 0;  Y = -(X+Z);
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
        A  = A - spdiags(sum(A,2),0,2*I,2*I);

        % resolver HJB estacionario (implícito)
        Bmat = (1/Delta + rho)*speye(2*I) - A;
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        V_stacked = Bmat \ (u_stacked + V_stacked/Delta);
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];

        if max(max(abs(V - v))) < crit
            break
        end
        v = V;
    end

    % ---- Fokker–Planck ----
    AT = A';
    b0 = zeros(2*I,1); i_fix = 1; b0(i_fix)=.1;
    row = zeros(1,2*I); row(i_fix)=1; AT(i_fix,:) = row;
    gg = AT \ b0;
    g  = [gg(1:I), gg(I+1:2*I)];
    g  = g / (sum(g(:))*da);

    % políticas de ahorro s(a)
    s = [res_inf, res_for] - (1+tau_c)*c;

    % ---- Gobierno: cuentas y deuda B(r)
    popI = sum(g(:,1))*da; popF = sum(g(:,2))*da;
    Ctot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;
    Tl   = tau_l * z2 * popF;
    Tc   = tau_c * Ctot;
    Tr   = phi   * z1 * popI;
    Gx   = Gov * (popI + popF);

    Btarget = (Tl + Tc - (Gx + Tr)) / max(r,1e-9);
    S = g(:,1)'*a*da + g(:,2)'*a*da - Btarget;

    % Si λ depende de composición, actualiza y repetir dentro del lazo
    if useCustomLambda
        L1_vec = max(la1_fun(a, popI), 1e-10);
        L2_vec = max(la2_fun(a, popI), 1e-10);
        Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
                     spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];
    end

    % ---- Búsqueda de r o r fijo
    if ~fix_r
        if S >  crit_S
            rmax = r; r = 0.5*(r + rmin);
        elseif S < -crit_S
            rmin = r; r = 0.5*(r + rmax);
        else
            break
        end

        % actualiza valor inicial para próxima iteración de r
        rr = r*ones(I,1); rr(a<0) = r + theta;
        v0(:,1) = ((((1-0)*zz(:,1) + rr.*a + phi*z1)/(1+tau_c)).^(1-sI)/(1-sI) + Gov)/rho;
        v0(:,2) = ((((1-tau_l)*zz(:,2) + r .*a          )/(1+tau_c)).^(1-sF)/(1-sF) + Gov)/rho;
    else
        break
    end
end

% ---- Resultados fiscales finales
rB = r * Btarget;
PB = Tl + Tc - (Gx + Tr);
BB = PB - rB;
Y  = z1*popI + z2*popF;

% ---- Estadísticas de riqueza y consumo
wI_mean = wmean(a, g(:,1), da);
wF_mean = wmean(a, g(:,2), da);
wT_mean = wmean(a, g(:,1)+g(:,2), da);
wI_med  = wmedian(a, g(:,1), da);
wF_med  = wmedian(a, g(:,2), da);
wT_med  = wmedian(a, g(:,1)+g(:,2), da);
gini_I  = giniW(a, g(:,1), da);
gini_F  = giniW(a, g(:,2), da);
gini_T  = giniW(a, g(:,1)+g(:,2), da);

cI_mean = wmean(c(:,1), g(:,1), da);
cF_mean = wmean(c(:,2), g(:,2), da);
cT_mean = wmean(c(:,1)+c(:,2), g(:,1)+g(:,2), da);
cI_med  = wmedian(c(:,1), g(:,1), da);
cF_med  = wmedian(c(:,2), g(:,2), da);
cT_med  = wmedian([c(:,1); c(:,2)], [g(:,1); g(:,2)], da);

% ---- Prestatarios / prestamistas
idxB = (a<0); idxL = (a>0);
massI = popI; massF = popF;
fracBorrow_I = sum(g(idxB,1))*da / max(massI,1e-12);
fracBorrow_F = sum(g(idxB,2))*da / max(massF,1e-12);
fracLend_I   = sum(g(idxL,1))*da / max(massI,1e-12);
fracLend_F   = sum(g(idxL,2))*da / max(massF,1e-12);

volBorrow_I  = sum(g(idxB,1).*a(idxB))*da;   % < 0
volBorrow_F  = sum(g(idxB,2).*a(idxB))*da;   % < 0
volLend_I    = sum(g(idxL,1).*a(idxL))*da;   % > 0
volLend_F    = sum(g(idxL,2).*a(idxL))*da;   % > 0

% ---- Empaquetar salida
out.r   = r;
out.a   = a;
out.g   = g;
out.c   = c;
out.s   = s;
out.popI = popI;
out.popF = popF;
out.Ctot = Ctot;
out.Y    = Y;

out.fiscal = struct('Tl',Tl,'Tc',Tc,'Tr',Tr,'G',Gx,'B',Btarget,'rB',rB,'PB',PB,'BB',BB);

out.stats  = struct( ...
    'wealth_mean',[wI_mean wF_mean wT_mean], ...
    'wealth_median',[wI_med wF_med wT_med], ...
    'gini',[gini_I gini_F gini_T], ...
    'cons_mean',[cI_mean cF_mean cT_mean], ...
    'cons_median',[cI_med cF_med cT_med], ...
    'p11', p11_rep ...
);

out.borrowers = struct( ...
    'fracBorrow',[fracBorrow_I fracBorrow_F], ...
    'fracLend',  [fracLend_I   fracLend_F], ...
    'volBorrow', [volBorrow_I  volBorrow_F], ...
    'volLend',   [volLend_I    volLend_F] ...
);

end

% ======================== Helpers (ponderados) ===========================
function m = wmean(x,w,da)
x=x(:); w=w(:); if nargin<3, da=1; end
W = sum(w)*da; if W<=0, m=NaN; else, m = sum(x.*w)*da / W; end
end

function med = wmedian(x,w,da)
x=x(:); w=w(:); if nargin<3, da=1; end
W = sum(w)*da; if W<=0, med=NaN; return; end
[xx,ix]=sort(x); ww=w(ix)*da; cw=cumsum(ww);
k = find(cw>=0.5*W,1,'first'); med = xx(k);
end

function gini = giniW(x,w,da)
x=x(:); w=w(:); if nargin<3, da=1; end
W = sum(w)*da; if W<=0, gini=NaN; return; end
xmin = min(0,min(x)); xs = x - xmin + 1e-12;    % desplaza a>0
[xx,ix]=sort(xs); ww=w(ix)*da;
cumw = cumsum(ww)/W; cumxw = cumsum(xx.*ww); if cumxw(end)<=0, gini=NaN; return; end
L = cumxw/cumxw(end); gini = 1 - 2*trapz(cumw,L);
end
