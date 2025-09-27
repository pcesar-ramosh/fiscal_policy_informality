%% ===================== SAVE NEXT BLOCK AS: solve_two_type_huggett_fiscal.m =====================
% (The solver function starts on the next line; save it as its own file.)

function out = solve_two_type_huggett_fiscal(params)
% -------------------------------------------------------------------------
% SOLVER for a two-type Huggett model with fiscal policy and informality.
% Taxes: tau_l on formal labor income; tau_c as VAT on consumption.
% Transfers: phi * z1 to each informal agent. Public good Gov adds to utility.
% Borrowing premium: theta so that effective rate is r+theta when a<0.
% Government budget closes with endogenous B(r): B = (Tl + Tc - G - Tr)/r.
% Occupational switching uses Poisson hazards λ1 (inf->for), λ2 (for->inf).
% Default: λ2 from p22_bar; λ1 chosen to match eta_target.
% Optional: params.lambdaSpec.mode='custom' with la1_fun(a,eta), la2_fun(a,eta).
%
% Outputs include policies (c,s), distribution g, equilibrium r, fiscal accounts,
% diagnostics (market clearing residual S_residual, HJB residual), and inequality.
% -------------------------------------------------------------------------

%% --- Read params with defaults
arg = @(f,def) (isfield(params,f) && ~isempty(params.(f))) * params.(f) + ...
               ~(isfield(params,f) && ~isempty(params.(f))) * def;

sI    = arg('RRA_I', 3.40);
sF    = arg('RRA_F', 3.40);
rho   = arg('rho',   0.05);
theta = arg('theta', 0.02);

tau_l = arg('tau_l', 0.15);
tau_c = arg('tau_c', 0.18);
phi   = arg('phi',   0.09);
Gov   = arg('Gov',   0.05);

z1    = arg('z1',    0.33);
z2    = arg('z2',    1.00);
z     = [z1, z2];

I     = arg('I',     700);
amin  = arg('amin', -0.30*z1);
amax  = arg('amax',  5.0);
a     = linspace(amin, amax, I)'; da = (amax-amin)/(I-1);
aa    = [a, a]; zz = ones(I,1) * z;

r     = arg('r_guess', 0.03);
rmin  = arg('rmin',     0.005);
rmax  = arg('rmax',     0.08);
fix_r = arg('fix_r',    0);

p22_bar    = arg('p22_bar', 0.8155);
%p22_bar    = arg('p22_bar', 0.8155);
eta_target = arg('eta_target', 0.64);

maxit_V = arg('maxit_V', 100);
crit_V  = arg('crit_V',  1e-6);
Delta   = arg('Delta',   1000);
maxit_r = arg('maxit_r', 1000);
crit_S  = arg('crit_S',  1e-5);

useCustomLambda = (isfield(params,'lambdaSpec') && isstruct(params.lambdaSpec) && ...
                   isfield(params.lambdaSpec,'mode') && strcmpi(params.lambdaSpec.mode,'custom'));
if useCustomLambda
    la1_fun = params.lambdaSpec.la1_fun;  % la1_fun(a, eta_guess)
    la2_fun = params.lambdaSpec.la2_fun;  % la2_fun(a, eta_guess)
end

%% --- Initialize hazards (match eta_target unless custom provided)
if useCustomLambda
    L1_vec = max(la1_fun(a, eta_target), 1e-10);
    L2_vec = max(la2_fun(a, eta_target), 1e-10);
    p11_rep = exp(-mean(L1_vec));
else
    la2 = -log(p22_bar);
    la1 = (1-eta_target)*la2 / max(eta_target,1e-12);
    L1_vec = la1 * ones(I,1);   % inf -> for
    L2_vec = la2 * ones(I,1);   % for -> inf
    p11_rep = exp(-la1);
end
Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
             spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];

%% --- Stable initial value (consuming current resources forever)
rr = r*ones(I,1); rr(a<0) = r + theta;

v0 = zeros(I,2);
v0(:,1) = (u_CRRA( ((1-0)*zz(:,1) + rr.*a + phi*z1)/(1+tau_c), sI ) + Gov)/rho;
v0(:,2) = (u_CRRA( ((1-tau_l)*zz(:,2) + r .*a         )/(1+tau_c), sF ) + Gov)/rho;

%% --- Outer loop on r (bisection unless fix_r)
Ir = maxit_r; if fix_r, Ir = 1; rmin=r; rmax=r; end
for ir = 1:Ir
    v = v0; V = v0; A = [];
    for n=1:maxit_V
        Vprev = V;
        % derivatives (upwind)
        dVf = zeros(I,2); dVb = zeros(I,2);
        dVf(1:I-1,:) = (Vprev(2:I,:) - Vprev(1:I-1,:))/da;
        dVf(I,1) = ((1-0)*z1 + r*amax)^(-sI);              % MU at upper boundary
        dVf(I,2) = ((1-tau_l)*z2 + r*amax)^(-sF);

        dVb(2:I,:) = (Vprev(2:I,:) - Vprev(1:I-1,:))/da;
        dVb(1,1)   = ((1-0)*z1 + (r+theta)*amin)^(-sI);    % MU at lower boundary
        dVb(1,2)   = ((1-tau_l)*z2 + r*amin)^(-sF);

        rr = r*ones(I,1); rr(a<0) = r + theta;             % effective rate

        % resources (before c) and policies via FOC: (1+tau_c)u'(c)=V_a
        res_inf = (1-0)*zz(:,1) + rr.*aa(:,1) + phi*z1;
        res_for = (1-tau_l)*zz(:,2) + r .*aa(:,2);
        cf = [max((1+tau_c)*dVf(:,1),1e-12).^(-1/sI), max((1+tau_c)*dVf(:,2),1e-12).^(-1/sF)];
        cb = [max((1+tau_c)*dVb(:,1),1e-12).^(-1/sI), max((1+tau_c)*dVb(:,2),1e-12).^(-1/sF)];
        ssf = [res_inf, res_for] - (1+tau_c)*cf;           % drift with forward c
        ssb = [res_inf, res_for] - (1+tau_c)*cb;           % drift with backward c
        c0  = [res_inf, res_for];                          % zero drift => c = res/(1+tau_c)
        If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
        c  = max( cf.*If + cb.*Ib + c0.*(I0/(1+tau_c)), 1e-12 );

        u = [u_CRRA(c(:,1), sI), u_CRRA(c(:,2), sF)] + Gov;

        % Generator A (rows sum to zero)
        X = max(-ssb,0)/da;  Z = max(ssf,0)/da; X(1,:) = 0; Z(I,:) = 0;  Y = -(X+Z);
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
        A  = A - spdiags(sum(A,2),0,2*I,2*I);

        % Implicit scheme for stationary HJB
        Bmat = (1/Delta + rho)*speye(2*I) - A;
        u_stacked = [u(:,1);u(:,2)]; V_stacked = [Vprev(:,1);Vprev(:,2)];
        V_stacked = Bmat \ (u_stacked + V_stacked/Delta);
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];

        if max(max(abs(V - Vprev))) < crit_V, break; end
    end

    % Stationary Fokker–Planck for g
    AT = A'; b0 = zeros(2*I,1); i_fix = 1; b0(i_fix)=.1; row = zeros(1,2*I); row(i_fix)=1; AT(i_fix,:) = row;
    gg = AT \ b0; g = [gg(1:I), gg(I+1:2*I)]; g = g / (sum(g(:))*da);

    % Savings drift s(a)
    s_pol = [res_inf, res_for] - (1+tau_c)*c;

    % Government accounts & market clearing residual S
    popI = sum(g(:,1))*da; popF = sum(g(:,2))*da;
    Ctot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;
    Tl   = tau_l * z2 * popF;     % labor tax
    Tc   = tau_c * Ctot;          % VAT revenue
    Tr   = phi   * z1 * popI;     % transfers
    Gx   = Gov * (popI + popF);   % public good provision

    Btarget = (Tl + Tc - (Gx + Tr)) / max(r,1e-9);
    A_priv  = g(:,1)'*a*da + g(:,2)'*a*da;
    S = A_priv - Btarget;         % excess private assets vs public debt

    if ~fix_r
        if S >  crit_S, rmax = r; r = 0.5*(r + rmin);
        elseif S < -crit_S, rmin = r; r = 0.5*(r + rmax);
        else, break; end
        % refresh initial guess values at new r
        rr = r*ones(I,1); rr(a<0) = r + theta;
        v0(:,1) = (u_CRRA( ((1-0)*zz(:,1) + rr.*a + phi*z1)/(1+tau_c), sI ) + Gov)/rho;
        v0(:,2) = (u_CRRA( ((1-tau_l)*zz(:,2) + r .*a         )/(1+tau_c), sF ) + Gov)/rho;
    else
        break
    end
end

% Fiscal results
rB = r * Btarget; PB = Tl + Tc - (Gx + Tr); BB = PB - rB; Y  = z1*popI + z2*popF;

% Inequality & stats
[wI_mean, wF_mean, wT_mean] = deal(wmean(a, g(:,1), da), wmean(a, g(:,2), da), wmean(a, g(:,1)+g(:,2), da));
[wI_med,  wF_med,  wT_med]  = deal(wmedian(a, g(:,1), da), wmedian(a, g(:,2), da), wmedian(a, g(:,1)+g(:,2), da));
[giniW_I, giniW_F, giniW_T] = deal(giniW(a, g(:,1), da), giniW(a, g(:,2), da), giniW(a, g(:,1)+g(:,2), da));
[cI_mean, cF_mean, cT_mean] = deal(wmean(c(:,1), g(:,1), da), wmean(c(:,2), g(:,2), da), wmean(c(:,1)+c(:,2), g(:,1)+g(:,2), da));
[cI_med,  cF_med,  cT_med]  = deal(wmedian(c(:,1), g(:,1), da), wmedian(c(:,2), g(:,2), da), wmedian([c(:,1);c(:,2)], [g(:,1);g(:,2)], da));
[giniC_I, giniC_F, giniC_T] = deal(giniX(c(:,1), g(:,1), da), giniX(c(:,2), g(:,2), da), giniX([c(:,1);c(:,2)], [g(:,1);g(:,2)], da));

% HJB residual diagnostic (∞-norm)
R = rho*[V(:,1);V(:,2)] - [u(:,1);u(:,2)] - A*[V(:,1);V(:,2)];
hjb_res = max(abs(R));

% Borrowers / lenders
idxB = (a<0); idxL = (a>0); massI = popI; massF = popF;
fracBorrow_I = sum(g(idxB,1))*da / max(massI,1e-12);
fracBorrow_F = sum(g(idxB,2))*da / max(massF,1e-12);
fracLend_I   = sum(g(idxL,1))*da / max(massI,1e-12);
fracLend_F   = sum(g(idxL,2))*da / max(massF,1e-12);
volBorrow_I  = sum(g(idxB,1).*a(idxB))*da;   % < 0
volBorrow_F  = sum(g(idxB,2).*a(idxB))*da;   % < 0
volLend_I    = sum(g(idxL,1).*a(idxL))*da;   % > 0
volLend_F    = sum(g(idxL,2).*a(idxL))*da;   % > 0

% Pack output
out.r   = r; out.a = a; out.g = g; out.c = c; out.s = s_pol; out.popI = popI; out.popF = popF; out.Ctot = Ctot; out.Y = Y;
out.fiscal = struct('Tl',Tl,'Tc',Tc,'Tr',Tr,'G',Gx,'B',Btarget,'rB',rB,'PB',PB,'BB',BB);
out.stats  = struct('wealth_mean',[wI_mean wF_mean wT_mean], 'wealth_median',[wI_med wF_med wT_med], ...
                    'giniW',[giniW_I giniW_F giniW_T], 'cons_mean',[cI_mean cF_mean cT_mean], ...
                    'cons_median',[cI_med cF_med cT_med], 'giniC',[giniC_I giniC_F giniC_T], 'p11', p11_rep);
out.borrowers = struct('fracBorrow',[fracBorrow_I fracBorrow_F], 'fracLend',[fracLend_I fracLend_F], ...
                       'volBorrow',[volBorrow_I volBorrow_F], 'volLend',[volLend_I volLend_F]);
out.S_residual = S; out.hjb_residual = hjb_res;

end

%% ========================== Helper functions ============================
function u = u_CRRA(c, sigma)
    c = max(c,1e-12);
    if abs(sigma-1) < 1e-10
        u = log(c);
    else
        u = c.^(1-sigma) / (1-sigma);
    end
end

function m = wmean(x,w,da)
    x=x(:); w=w(:); if nargin<3, da=1; end
    W = sum(w)*da; if W<=0, m=NaN; else, m = sum(x.*w)*da / W; end
end

function med = wmedian(x,w,da)
    x=x(:); w=w(:); if nargin<3, da=1; end
    W = sum(w)*da; if W<=0, med=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cw=cumsum(ww); k = find(cw>=0.5*W,1,'first'); med = xx(k);
end

function gini = giniW(x,w,da)
    x=x(:); w=w(:); if nargin<3, da=1; end
    W = sum(w)*da; if W<=0, gini=NaN; return; end
    xmin = min(0,min(x)); xs = x - xmin + 1e-12;
    [xx,ix]=sort(xs); ww=w(ix)*da; cumw = cumsum(ww)/W; cumxw = cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end
    L = cumxw/cumxw(end); gini = 1 - 2*trapz(cumw,L);
end

function gini = giniX(x,w,da)
    % Gini of a variable x across the population with density weights w
    x=x(:); w=w(:); if nargin<3, da=1; end
    W = sum(w)*da; if W<=0, gini=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cumw = cumsum(ww)/W; cumxw = cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end
    L = cumxw/cumxw(end); gini = 1 - 2*trapz(cumw,L);
end
