function out = hugguet_no_fiscal_policy(params)
% -------------------------------------------------------------------------
% Modelo Huggett con 2 tipos (Informal/Formal) SIN política fiscal:
%  - tau_l = 0, tau_c = 0, phi = 0, Gov = 0, B(r) = 0
%  - r se determina cerrando el mercado de activos PRIVADO: ∫a g(a) da = 0
%  - Informal paga prima θ cuando a<0
%  - Transiciones ocupacionales: "legacy" con p22 dado y eta_target
%
% ENTRADA (params): RRA_I,RRA_F,rho,theta,z1,z2,I,amin,amax,
%                   r_guess,rmin,rmax,p22_bar,eta_target (mismos nombres que tu base)
% SALIDA (out): estructura compatible con tu main "base", pero con fiscal=0
% -------------------------------------------------------------------------

%% ----- Lectura de parámetros y defaults
arg = @(f,def) (isfield(params,f) && ~isempty(params.(f))) * params.(f) + ...
               ~(isfield(params,f) && ~isempty(params.(f))) * def;

sI    = arg('RRA_I', 3.40);
sF    = arg('RRA_F', 3.40);
rho   = arg('rho',   0.05);
theta = arg('theta', 0.02);

% sin política fiscal
tau_l = 0; tau_c = 0; phi = 0; Gov = 0;

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

p22_bar    = arg('p22_bar', 0.8155);
eta_target = arg('eta_target', 0.64);

useCustomLambda = (isfield(params,'lambdaSpec') && isstruct(params.lambdaSpec) && ...
                   isfield(params.lambdaSpec,'mode') && strcmpi(params.lambdaSpec.mode,'custom'));
if useCustomLambda
    la1_fun = params.lambdaSpec.la1_fun;  % la1_fun(a, eta_guess)
    la2_fun = params.lambdaSpec.la2_fun;  % la2_fun(a, eta_guess)
end

%% ----- Iteración HJB + búsqueda de r
maxit = 100; crit = 1e-6; Delta = 1000;
Ir    = 1000; crit_S = 1e-5;

% Lambda inicial (legacy o custom)
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

% Valor inicial simple
rr = r*ones(I,1); rr(a<0) = r + theta;
v0(:,1) = ((((1-0)*zz(:,1) + rr.*a )).^(1-sI)/(1-sI) + Gov)/rho;
v0(:,2) = ((((1-0)*zz(:,2) + r .*a )).^(1-sF)/(1-sF) + Gov)/rho;

for ir = 1:Ir
    v = v0;
    for n=1:maxit
        V = v;

        % derivadas forward/backward
        dVf = zeros(I,2); dVb = zeros(I,2);
        dVf(1:I-1,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVf(I,1) = ((1-0)*z1 + r*amax)^(-sI);
        dVf(I,2) = ((1-0)*z2 + r*amax)^(-sF);

        dVb(2:I,:) = (V(2:I,:) - V(1:I-1,:))/da;
        dVb(1,1)   = ((1-0)*z1 + (r+theta)*amin)^(-sI);
        dVb(1,2)   = ((1-0)*z2 + r*amin)^(-sF);

        % tasas efectivas (informal paga prima si a<0)
        rr = r*ones(I,1); rr(a<0) = r + theta;

        % “recursos” (sin impuestos ni transferencias)
        res_inf = (1-0)*zz(:,1) + rr.*aa(:,1);
        res_for = (1-0)*zz(:,2) + r .*aa(:,2);

        % consumo (upwind, sin IVA)
        cf = [max(dVf(:,1),1e-10).^(-1/sI), max(dVf(:,2),1e-10).^(-1/sF)];
        cb = [max(dVb(:,1),1e-10).^(-1/sI), max(dVb(:,2),1e-10).^(-1/sF)];
        ssf = [res_inf, res_for] - cf;
        ssb = [res_inf, res_for] - cb;
        c0  = [res_inf, res_for];

        If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
        c  = cf.*If + cb.*Ib + c0.*I0;   % si drift≈0, c = res

        u = [c(:,1).^(1-sI)/(1-sI), c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);

        % generador (suma filas=0)
        X = max(-ssb,0)/da;  Z = max(ssf,0)/da;
        X(1,:) = 0; Z(I,:) = 0;  Y = -(X+Z);
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
        A  = A - spdiags(sum(A,2),0,2*I,2*I);

        % resolver HJB estacionario
        Bmat = (1/Delta + rho)*speye(2*I) - A;
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
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

    % políticas de ahorro s(a)
    s = [res_inf, res_for] - c;

    % ---- Mercado de activos (solo privado)
    popI = sum(g(:,1))*da; popF = sum(g(:,2))*da;
    Ctot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;

    % B(public) = 0  =>   S = ∫ a g(a) da
    S = g(:,1)'*a*da + g(:,2)'*a*da;

    % si λ depende de composición, actualiza
    if useCustomLambda
        L1_vec = max(la1_fun(a, popI), 1e-10);
        L2_vec = max(la2_fun(a, popI), 1e-10);
        Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
                     spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];
    end

    if S >  crit_S
        rmax = r; r = 0.5*(r + rmin);
    elseif S < -crit_S
        rmin = r; r = 0.5*(r + rmax);
    else
        break
    end

    % actualizar v0 para nueva r
    rr = r*ones(I,1); rr(a<0) = r + theta;
    v0(:,1) = ((((1-0)*zz(:,1) + rr.*a )).^(1-sI)/(1-sI) + Gov)/rho;
    v0(:,2) = ((((1-0)*zz(:,2) + r .*a )).^(1-sF)/(1-sF) + Gov)/rho;
end

%% ---- Resultados "fiscales" nulos y agregados
Tl=0; Tc=0; Tr=0; Gx=0; Btarget=0; rB=0; PB=0; BB=0;
Y  = z1*popI + z2*popF;

% Estadísticas
wI_mean = wmean(a, g(:,1), da);  wF_mean = wmean(a, g(:,2), da);
wT_mean = wmean(a, g(:,1)+g(:,2), da);
wI_med  = wmedian(a, g(:,1), da); wF_med  = wmedian(a, g(:,2), da);
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

idxB = (a<0); idxL=(a>0); massI=popI; massF=popF;
fracBorrow_I = sum(g(idxB,1))*da / max(massI,1e-12);
fracBorrow_F = sum(g(idxB,2))*da / max(massF,1e-12);
fracLend_I   = sum(g(idxL,1))*da / max(massI,1e-12);
fracLend_F   = sum(g(idxL,2))*da / max(massF,1e-12);
volBorrow_I  = sum(g(idxB,1).*a(idxB))*da;   % <0
volBorrow_F  = sum(g(idxB,2).*a(idxB))*da;   % <0
volLend_I    = sum(g(idxL,1).*a(idxL))*da;   % >0
volLend_F    = sum(g(idxL,2).*a(idxL))*da;   % >0

% ---- Empaquetar salida
out.r   = r;     out.a   = a;     out.g = g;  out.c = c;  out.s = s;
out.popI = popI; out.popF = popF; out.Ctot = Ctot; out.Y = Y;

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
    'fracLend',  [fracLend_I fracLend_F], ...
    'volBorrow', [volBorrow_I volBorrow_F], ...
    'volLend',   [volLend_I volLend_F] ...
);

end

% ======================== Helpers ===========================
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
xmin = min(0,min(x)); xs = x - xmin + 1e-12;
[xx,ix]=sort(xs); ww=w(ix)*da;
cumw = cumsum(ww)/W; cumxw = cumsum(xx.*ww);
if cumxw(end)<=0, gini=NaN; return; end
L = cumxw/cumxw(end); gini = 1 - 2*trapz(cumw,L);
end
