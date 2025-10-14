function out = solve_two_type_huggett_fiscal_Bfixed(cfg)
% -------------------------------------------------------------------------
% Dos tipos (informal/formal) estilo Huggett con:
%   - Cierre financiero por bisección en r:  A_priv(r) = B
%   - Deuda pública:  B = Bbar * Y  (B_mode='ratio_to_Y') o nivel ('level')
%   - Cierre fiscal EXACTO por G:     G = Tl + Tc - Tr - rB  => BB = 0
%   - Bien público multiplicativo en utilidad y MU (vía Gpc)
%   - Punto fijo externo en Gpc para coherencia FOC-cuentas
%   - Bracketing robusto de S(r) y fallback al borde si no hay cambio de signo
% -------------------------------------------------------------------------

% ======= Lectura de parámetros =======
arg   = @(f,def) get_arg(cfg,f,def);
sI    = arg('RRA_I',3.40);   sF    = arg('RRA_F',3.40);   rho   = arg('rho',0.08);
tau_l = arg('tau_l',0.15);   tau_c = arg('tau_c',0.18);   phi   = arg('phi',0.09);
z1    = arg('z1',0.33);      z2    = arg('z2',1.00);
theta_I = arg('theta_I',0.06); theta_F = arg('theta_F',0.01);

I     = arg('I',700);  amin = arg('amin',-2.0*z1);  amax = arg('amax',3.0);
a     = linspace(amin,amax,I)';  da = (amax-amin)/(I-1);

r_guess = arg('r_guess',0.03);  rmin = arg('rmin',0.005);  rmax = arg('rmax',0.10);
fix_r   = arg('fix_r',0);       maxit_r = arg('maxit_r',120);  crit_S = arg('crit_S',1e-5);

p22_bar = arg('p22_bar',0.8155); eta_target = arg('eta_target',0.654);
la2 = -log(p22_bar);  la1 = (1-eta_target)*la2/max(eta_target,1e-12);
L1_vec = la1*ones(I,1); L2_vec = la2*ones(I,1); p11_rep = exp(-la1);
Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
             spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];

sigma_a = arg('sigma_a',0.010);  nu = 0.5*(sigma_a^2)/(da^2);

psi_G   = arg('psi_G',0.08);
omegaG  = arg('omegaG',0.50);
report_G = arg('report_G_effects',0);

maxit_V = arg('maxit_V',160);  crit_V = arg('crit_V',1e-6);  Delta = arg('Delta',1400);

B_mode  = arg('B_mode','ratio_to_Y');   % 'ratio_to_Y' | 'level'
Bbar    = arg('Bbar',0.35);
gov_mode= arg('gov_mode','close_by_G'); % 'close_by_G' (recomendado) | 'fixed_G'
G_fixed = arg('G_fixed',0.0);           % usado solo si gov_mode='fixed_G'

% ======= Resolución en r (con bracketing robusto) =======
if fix_r
    out = solve_given_r(r_guess); 
    return;
else
    % evaluar S(rmin), S(rmax); expandir si no hay cambio de signo
    [Smin, ~] = eval_S(rmin);
    [Smax, ~] = eval_S(rmax);
    it_expand = 0;
    while Smin*Smax > 0 && it_expand < 6
        rmin = max(1e-4, 0.7*rmin);
        rmax = 1.4*rmax;
        [Smin, ~] = eval_S(rmin);
        [Smax, ~] = eval_S(rmax);
        it_expand = it_expand + 1;
    end
    if Smin*Smax > 0
        % no hubo cambio de signo; devolvemos el borde "menos malo"
        if abs(Smin) < abs(Smax), out = solve_given_r(rmin); else, out = solve_given_r(rmax); end
        warning('Bisection: no sign change; returning boundary solution.');
        return;
    end

    r_low = rmin; S_low = Smin;
    r_high= rmax; S_high= Smax;
    best = []; r_mid = r_guess;

    for it=1:maxit_r
        r_mid = 0.5*(r_low + r_high);
        [S_mid, sol_mid] = eval_S(r_mid);
        if isempty(best) || abs(S_mid) < abs(best.S_residual), best = sol_mid; end
        if abs(S_mid) < crit_S, break; end
        if S_low * S_mid <= 0
            r_high = r_mid; S_high = S_mid;
        else
            r_low = r_mid;  S_low  = S_mid;
        end
    end
    out = best; 
    return;
end

% =================== subfunciones (nivel 1) ===================
    function [Sval, solr] = eval_S(rr)
        solr = solve_given_r(rr);
        Sval = solr.S_residual;
    end

    function out_r = solve_given_r(r)
        % ----- punto fijo en Gpc -----
        Gpc = 0.0;
        for itG = 1:12
            [V, c, A, rr1, rr2] = solve_HJB_given_Gpc(r, a, da, ...
                sI,sF,rho, z1,z2, tau_l,tau_c, phi, theta_I,theta_F, ...
                psi_G,omegaG, Gpc, sigma_a,nu, Aswitch, maxit_V,crit_V,Delta);

            % Fokker-Planck implícito (pseudo-tiempo)
            g = solve_FP(A, I, da);

            % Agregados
            popI = sum(g(:,1))*da;  popF = sum(g(:,2))*da;  massT = popI + popF;
            Ctot = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;
            Y    = z1*popI + z2*popF;

            % Bono público
            if strcmpi(B_mode,'ratio_to_Y'), B = Bbar*Y; else, B = Bbar; end
            rB = r*B;

            % Fiscal: cierre por G o G fijo
            Tl = tau_l*z2*popF;
            Tc = tau_c*Ctot;
            Tr = phi*z1*popI;
            if strcmpi(gov_mode,'close_by_G')
                Gx  = Tl + Tc - Tr - rB;  % PB = rB, BB = 0
                PB  = rB;  BB = 0.0;
            else
                Gx  = G_fixed;
                PB  = Tl + Tc - Tr - Gx;
                BB  = PB - rB;
            end
            Gpc_new = Gx / max(massT,1e-12);

            if abs(Gpc_new - Gpc) < 1e-9, break; end
            Gpc = 0.6*Gpc + 0.4*Gpc_new;
        end

        % Cierre financiero
        A_priv = g(:,1)'*a*da + g(:,2)'*a*da;
        S = A_priv - B;

        % Políticas de ahorro
        s_pol = [ z1 + rr1.*a + phi*z1 - (1+tau_c)*c(:,1), ...
                  (1-tau_l)*z2 + rr2.*a - (1+tau_c)*c(:,2) ];

        % Stats
        [wI_mean,wF_mean,wT_mean] = deal(wmean(a,g(:,1),da), wmean(a,g(:,2),da), wmean(a,g(:,1)+g(:,2),da));
        [wI_med,wF_med,wT_med]   = deal(wmedian(a,g(:,1),da), wmedian(a,g(:,2),da), wmedian(a,g(:,1)+g(:,2),da));
        [giniW_I,giniW_F,giniW_T]= deal(giniW(a,g(:,1),da),giniW(a,g(:,2),da),giniW(a,g(:,1)+g(:,2),da));
        [cI_mean,cF_mean,cT_mean]= deal(wmean(c(:,1),g(:,1),da),wmean(c(:,2),g(:,2),da),wmean(c(:,1)+c(:,2),g(:,1)+g(:,2),da));
        [cI_med,cF_med,cT_med]   = deal(wmedian(c(:,1),g(:,1),da),wmedian(c(:,2),g(:,2),da),wmedian([c(:,1);c(:,2)],[g(:,1);g(:,2)],da));
        [giniC_I,giniC_F,giniC_T]= deal(giniX(c(:,1),g(:,1),da),giniX(c(:,2),g(:,2),da),giniX([c(:,1);c(:,2)],[g(:,1);g(:,2)],da));

        % Residuo HJB
        u_check = [u_mult(c(:,1),Gpc,sI,psi_G,omegaG), u_mult(c(:,2),Gpc,sF,psi_G,omegaG)];
        R = rho*[V(:,1);V(:,2)] - [u_check(:,1);u_check(:,2)] - A*[V(:,1);V(:,2)];
        hjb_res = max(abs(R));

        % Deudores / Ahorradores
        idxB=(a<0); idxL=(a>0);
        fracBorrow_I=sum(g(idxB,1))*da / max(sum(g(:,1))*da,1e-12);
        fracBorrow_F=sum(g(idxB,2))*da / max(sum(g(:,2))*da,1e-12);
        fracLend_I  =sum(g(idxL,1))*da / max(sum(g(:,1))*da,1e-12);
        fracLend_F  =sum(g(idxL,2))*da / max(sum(g(:,2))*da,1e-12);
        volBorrow_I =sum(g(idxB,1).*a(idxB))*da;
        volBorrow_F =sum(g(idxB,2).*a(idxB))*da;
        volLend_I   =sum(g(idxL,1).*a(idxL))*da;
        volLend_F   =sum(g(idxL,2).*a(idxL))*da;

        % Output
        out_r.r = r; out_r.a = a; out_r.g = g; out_r.c = c; out_r.s = s_pol;
        out_r.popI = popI; out_r.popF = popF; out_r.Ctot = Ctot; out_r.Y = Y;
        out_r.S_residual = S; out_r.hjb_residual = hjb_res; out_r.Gpc = Gpc;

        B_over_Y = B / max(Y,1e-12);
        out_r.fiscal = struct('Tl',Tl,'Tc',Tc,'Tr',Tr,'G',Gx,'B',B,'rB',rB,'PB',PB,'BB',BB, ...
                              'B_over_Y', B_over_Y);

        out_r.stats = struct('wealth_mean',[wI_mean wF_mean wT_mean], ...
                             'wealth_median',[wI_med wF_med wT_med], ...
                             'giniW',[giniW_I giniW_F giniW_T], ...
                             'cons_mean',[cI_mean cF_mean cT_mean], ...
                             'cons_median',[cI_med cF_med cT_med], ...
                             'giniC',[giniC_I giniC_F giniC_T], 'p11',p11_rep);

        out_r.borrowers = struct('fracBorrow',[fracBorrow_I fracBorrow_F], ...
                                 'fracLend',[fracLend_I fracLend_F], ...
                                 'volBorrow',[volBorrow_I volBorrow_F], ...
                                 'volLend',[volLend_I volLend_F]);

        if report_G
            mu_mult_I = (1+psi_G*Gpc)^(omegaG*(1-sI));
            mu_mult_F = (1+psi_G*Gpc)^(omegaG*(1-sF));
            u_mult_fac= (1+psi_G*Gpc)^(omegaG);
            out_r.G_effects = struct('mu_mult_I',mu_mult_I,'mu_mult_F',mu_mult_F,'u_mult',u_mult_fac,'Gpc',Gpc);
        else
            out_r.G_effects = struct();
        end

        out_r.params.public_good = struct('psi_G',psi_G,'omegaG',omegaG);
        out_r.r_path = r; out_r.S_path = S; % trazas simples
    end

% =================== subfunciones (nivel 2) ===================
    function [V, c, A, rr1, rr2] = solve_HJB_given_Gpc(r, a, da, ...
            sI,sF,rho, z1,z2, tau_l,tau_c, phi, theta_I,theta_F, ...
            psi_G,omegaG, Gpc, sigma_a,nu, Aswitch, maxit_V,crit_V,Delta)

        I = numel(a);
        rr1 = r*ones(I,1); rr1(a<0)=r+theta_I;
        rr2 = r*ones(I,1); rr2(a<0)=r+theta_F;

        V = zeros(I,2);
        V(:,1) = u_CRRA(max((z1 + rr1.*a + phi*z1)/(1+tau_c),1e-12), sI)/rho;
        V(:,2) = u_CRRA(max(((1-tau_l)*z2 + rr2.*a)/(1+tau_c),1e-12), sF)/rho;

        facI = (1+psi_G*Gpc)^(omegaG*(1-sI));
        facF = (1+psi_G*Gpc)^(omegaG*(1-sF));

        for it=1:maxit_V
            Vprev = V;

            dVf=zeros(I,2); dVb=zeros(I,2);
            dVf(1:I-1,:) = (Vprev(2:I,:) - Vprev(1:I-1,:))/da;
            dVb(2:I,:)   = (Vprev(2:I,:) - Vprev(1:I-1,:))/da;

            rr1 = r*ones(I,1); rr1(a<0)=r+theta_I;
            rr2 = r*ones(I,1); rr2(a<0)=r+theta_F;

            % Fronteras
            c_inf_max = max((z1 + rr1(end)*a(end) + phi*z1)/(1+tau_c),1e-12);
            c_for_max = max(((1-tau_l)*z2 + rr2(end)*a(end))/(1+tau_c),1e-12);
            c_inf_min = max((z1 + rr1(1)*a(1) + phi*z1)/(1+tau_c),1e-12);
            c_for_min = max(((1-tau_l)*z2 + rr2(1)*a(1))/(1+tau_c),1e-12);

            dVf(I,1) = (1+tau_c)*facI*u_CRRA_prime(c_inf_max,sI);
            dVf(I,2) = (1+tau_c)*facF*u_CRRA_prime(c_for_max,sF);
            dVb(1,1) = (1+tau_c)*facI*u_CRRA_prime(c_inf_min,sI);
            dVb(1,2) = (1+tau_c)*facF*u_CRRA_prime(c_for_min,sF);

            % Recursos previos a c
            res_inf = z1 + rr1.*a + phi*z1;
            res_for = (1-tau_l)*z2 + rr2.*a;

            % Euler (corrigiendo por IVA y multiplicador MU)
            cf = [max(dVf(:,1)./((1+tau_c)*facI),1e-12).^(-1/sI), ...
                  max(dVf(:,2)./((1+tau_c)*facF),1e-12).^(-1/sF)];
            cb = [max(dVb(:,1)./((1+tau_c)*facI),1e-12).^(-1/sI), ...
                  max(dVb(:,2)./((1+tau_c)*facF),1e-12).^(-1/sF)];

            ssf = [res_inf, res_for] - (1+tau_c)*cf;
            ssb = [res_inf, res_for] - (1+tau_c)*cb;
            c0  = [res_inf, res_for];
            If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
            c  = max( cf.*If + cb.*Ib + c0.*(I0/(1+tau_c)), 1e-12 );

            % Utilidad flujo
            u = [u_mult(c(:,1),Gpc,sI,psi_G,omegaG), u_mult(c(:,2),Gpc,sF,psi_G,omegaG)];

            % Generador (drift + difusión + switching)
            X = max(-ssb,0)/da;  Z = max(ssf,0)/da;  X(1,:)=0; Z(I,:)=0; Y=-(X+Z);
            A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);

            if sigma_a>0
                e=nu*ones(I,1); D2=spdiags([e -2*e e],-1:1,I,I);
                D2(1,1)=-nu; D2(1,2)=nu; D2(I,I)=-nu; D2(I,I-1)=nu;
                A1 = A1 + D2;  A2 = A2 + D2;
            end
            A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;
            A = A - spdiags(sum(A,2),0,2*I,2*I);

            Bmat = (1/Delta + rho)*speye(2*I) - A;
            Vst  = Bmat \ ([u(:,1);u(:,2)] + [Vprev(:,1);Vprev(:,2)]/Delta);
            V    = [Vst(1:I), Vst(I+1:2*I)];

            if max(max(abs(V - Vprev))) < crit_V, break; end
        end
    end

    function g = solve_FP(A, I, da)
        AT = A'; g = ones(2*I,1); g = g/(sum(g)*da);
        tauFP = 80*da;
        for k=1:1000
            g_new = (speye(2*I) - tauFP*AT) \ g;
            g_new = max(g_new,0); g_new = g_new/(sum(g_new)*da);
            if norm(g_new - g, inf) < 1e-12, break; end
            g = g_new;
        end
        g = [g(1:I), g(I+1:2*I)];
    end
end

% ====================== Helpers ======================
function v = get_arg(cfg,f,def)
    if isfield(cfg,f) && ~isempty(cfg.(f)), v = cfg.(f); else, v = def; end
end
function u = u_CRRA(c,sigma)
    c=max(c,1e-12); if abs(sigma-1)<1e-10, u=log(c); else, u=c.^(1-sigma)/(1-sigma); end
end
function up = u_CRRA_prime(c,sigma), c=max(c,1e-12); up=c.^(-sigma); end
function u = u_mult(c,Gpc,sigma,psi,omega)
    c_eff = max(c,1e-12).*(1+psi*Gpc).^omega; u = u_CRRA(c_eff,sigma);
end
function m = wmean(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; if W<=0, m=NaN; else, m=sum(x.*w)*da/W; end
end
function med = wmedian(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; if W<=0, med=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cw=cumsum(ww); k=find(cw>=0.5*W,1,'first'); med=xx(k);
end
function gini = giniW(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; if W<=0, gini=NaN; return; end
    xmin=min(0,min(x)); xs=x-xmin+1e-12; [xx,ix]=sort(xs); ww=w(ix)*da;
    cumw=cumsum(ww)/W; cumxw=cumsum(xx.*ww); if cumxw(end)<=0, gini=NaN; return; end
    L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end
function gini = giniX(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; if W<=0, gini=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cumw=cumsum(ww)/W; cumxw=cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end, L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end
