function out = solve_two_type_huggett_fiscal_Bfixed(cfg)
% Huggett de dos tipos con cierres fiscales:
%  (1) 'G_resid_B_fixed'              -> G es residuo, B fijo (nivel o ratio)
%  (2) 'G_fixed_B_fixed_adjust_tax'   -> G fijo, B fijo, ajusta τ_c o τ_l para BB=0
% Sin funciones anidadas (todo son subfunciones) y TODAS con 'end'.

% ------------------ Empaqueta parámetros y grid ------------------
K = build_K(cfg);

switch lower(K.fiscal_mode)
    case 'g_resid_b_fixed'
        % Cierre estándar: G residual con impuestos dados
        out = solve_GE_given_taxes(K, K.tau_l0, K.tau_c0, 'residG');

    case 'g_fixed_b_fixed_adjust_tax'
        % Ajusta 'vat' (τ_c) o 'labor' (τ_l) para BB=0 manteniendo G fijo y B fijo
        adjust = K.adjust_tax;     % 'vat' | 'labor'
        tax_tol = K.tax_tol;

        switch adjust
            case 'vat',   tlo = 0.00; thi = 0.50;   % IVA 0–50%
            case 'labor', tlo = 0.00; thi = 0.60;   % τ_l 0–60%
            otherwise, error('adjust_tax debe ser "vat" o "labor".');
        end

        % Evalúa BB en extremos y fuerza cambio de signo si es posible
        [BBlo, solLo] = residual_BB(K, tlo);
        [BBhi, solHi] = residual_BB(K, thi);
        it=0;
        while BBlo*BBhi>0 && it<8
            tlo = max(0.00, 0.5*tlo);
            thi = min(0.95, 1.2*thi);
            [BBlo, solLo] = residual_BB(K, tlo);
            [BBhi, solHi] = residual_BB(K, thi);
            it=it+1;
        end

        if BBlo*BBhi>0
            warning('No hay cambio de signo en BB(t). Devuelvo el extremo con menor |BB|.');
            if abs(BBlo) < abs(BBhi), out = solLo; else, out = solHi; end
            return;
        end

        % Bisección en el impuesto elegido
        for k=1:50
            tm = 0.5*(tlo+thi);
            [BBm, solM] = residual_BB(K, tm);
            if abs(BBm) < tax_tol
                out = solM; return;
            end
            if BBlo*BBm <= 0
                thi = tm;  BBhi = BBm;
            else
                tlo = tm;  BBlo = BBm;
            end
        end
        out = solM;

    otherwise
        error('fiscal_mode no reconocido.');
end

end % ======= solve_two_type_huggett_fiscal_Bfixed =======


% =====================================================================
% ============================= SUBFUNCIONES ===========================
% =====================================================================

function K = build_K(cfg)
    % Lee argumentos con default
    K.tau_l0 = get_arg(cfg,'tau_l',0.15);
    K.tau_c0 = get_arg(cfg,'tau_c',0.18);
    K.phi    = get_arg(cfg,'phi',0.09);

    K.RRA_I  = get_arg(cfg,'RRA_I',3.40); 
    K.RRA_F  = get_arg(cfg,'RRA_F',3.40); 
    K.rho    = get_arg(cfg,'rho',0.08);

    K.z1     = get_arg(cfg,'z1',0.33);
    K.z2     = get_arg(cfg,'z2',1.00);

    K.theta_I= get_arg(cfg,'theta_I',0.06);
    K.theta_F= get_arg(cfg,'theta_F',0.01);

    K.I      = get_arg(cfg,'I',700);
    amin     = get_arg(cfg,'amin',-2.0*K.z1);
    amax     = get_arg(cfg,'amax',3.0);
    K.a      = linspace(amin,amax,K.I)'; 
    K.da     = (amax-amin)/(K.I-1);

    K.r_guess= get_arg(cfg,'r_guess',0.03);
    K.rmin   = get_arg(cfg,'rmin',0.005);
    K.rmax   = get_arg(cfg,'rmax',0.10);
    K.fix_r  = get_arg(cfg,'fix_r',0);
    K.maxit_r= get_arg(cfg,'maxit_r',80);
    K.crit_S = get_arg(cfg,'crit_S',1e-5);

    K.p22_bar= get_arg(cfg,'p22_bar',0.8155);
    K.eta_target = get_arg(cfg,'eta_target',0.654);
    la2 = -log(K.p22_bar); 
    la1 = (1-K.eta_target)*la2/max(K.eta_target,1e-12);
    K.Aswitch = [ -spdiags(la1*ones(K.I,1),0,K.I,K.I),  spdiags(la1*ones(K.I,1),0,K.I,K.I);
                   spdiags(la2*ones(K.I,1),0,K.I,K.I), -spdiags(la2*ones(K.I,1),0,K.I,K.I) ];
    K.p11_rep = exp(-la1);

    K.sigma_a= get_arg(cfg,'sigma_a',0.007);
    K.nu     = 0.5*(K.sigma_a^2)/(K.da^2);

    K.psi_G  = get_arg(cfg,'psi_G',0.08);
    K.omegaG = get_arg(cfg,'omegaG',0.50);
    K.report_G_effects = get_arg(cfg,'report_G_effects',0);

    K.maxit_V= get_arg(cfg,'maxit_V',160);
    K.crit_V = get_arg(cfg,'crit_V',1e-6);
    K.Delta  = get_arg(cfg,'Delta',1400);

    % --- bloque fiscal ---
    K.fiscal_mode = get_arg(cfg,'fiscal_mode','G_resid_B_fixed');
    K.B_mode  = get_arg(cfg,'B_mode','level');   % 'level' | 'ratio_to_Y'
    K.Bbar    = get_arg(cfg,'Bbar',0.25);
    K.alphaG  = get_arg(cfg,'alphaG',0.50);
    K.clampG0 = get_arg(cfg,'clamp_G_to_zero',true);
    K.G_cap_ratio = get_arg(cfg,'G_cap_ratio',Inf);
    if isfield(cfg,'G_target_ratio'), K.G_target_ratio = cfg.G_target_ratio; else, K.G_target_ratio = []; end
    if isfield(cfg,'G_target_level'), K.G_target_level = cfg.G_target_level; else, K.G_target_level = []; end

    if isfield(cfg,'adjust_tax'), K.adjust_tax = lower(cfg.adjust_tax); else, K.adjust_tax = 'vat'; end
    K.tax_tol = get_arg(cfg,'tax_tol',1e-5);
end

function [BBres, sol] = residual_BB(K, taxval)
    % Evalúa BB(tax) resolviendo el general equilibrium con G fijo y B fijo
    switch K.adjust_tax
        case 'vat',   tau_c = taxval; tau_l = K.tau_l0;
        case 'labor', tau_l = taxval; tau_c = K.tau_c0;
        otherwise, error('adjust_tax debe ser "vat" o "labor".');
    end
    sol   = solve_GE_given_taxes(K, tau_l, tau_c, 'fixedG');
    BBres = sol.fiscal.BB;
end

function out = solve_GE_given_taxes(K, tau_l, tau_c, Gmode)
    % Resuelve r* mediante bisección en S(r)=A_priv - B

    if K.fix_r
        out = solve_given_r(K, K.r_guess, tau_l, tau_c, Gmode);
        out.params.tau_l = tau_l; out.params.tau_c = tau_c;
        return;
    end

    [SL,~] = eval_S_at(K, K.rmin, tau_l, tau_c, Gmode);
    [SH,~] = eval_S_at(K, K.rmax, tau_l, tau_c, Gmode);
    it = 0;
    while SL*SH>0 && it<6
        K.rmin = max(0.003, 0.7*K.rmin);
        K.rmax = 1.4*K.rmax;
        [SL,~] = eval_S_at(K, K.rmin, tau_l, tau_c, Gmode);
        [SH,~] = eval_S_at(K, K.rmax, tau_l, tau_c, Gmode);
        it=it+1;
    end

    if SL*SH>0
        warning('Bisection (r): no sign change; boundary solution.');
        if abs(SL) < abs(SH)
            out = solve_given_r(K, K.rmin, tau_l, tau_c, Gmode);
        else
            out = solve_given_r(K, K.rmax, tau_l, tau_c, Gmode);
        end
        out.params.tau_l = tau_l; out.params.tau_c = tau_c; 
        return;
    end

    rL = K.rmin; rH = K.rmax; best = [];
    for k=1:K.maxit_r
        rm = 0.5*(rL+rH);
        [Sm, solm] = eval_S_at(K, rm, tau_l, tau_c, Gmode);
        if isempty(best) || abs(Sm) < abs(best.S_residual), best = solm; end
        if abs(Sm) < K.crit_S, break; end
        if SL*Sm <= 0
            rH = rm; SH = Sm;
        else
            rL = rm; SL = Sm;
        end
    end
    out = best; out.params.tau_l = tau_l; out.params.tau_c = tau_c;
end

function [Sval, solm] = eval_S_at(K, r, tau_l, tau_c, Gmode)
    solm = solve_given_r(K, r, tau_l, tau_c, Gmode);
    Sval = solm.S_residual;
end

function out = solve_given_r(K, r, tau_l, tau_c, Gmode)
    % ---------- Inicialización coherente ----------
    rr1 = r*ones(K.I,1); rr1(K.a<0) = r + K.theta_I;
    rr2 = r*ones(K.I,1); rr2(K.a<0) = r + K.theta_F;
    V = zeros(K.I,2);
    V(:,1) = u_CRRA(max((K.z1 + rr1.*K.a + K.phi*K.z1)/(1+tau_c),1e-12), K.RRA_I)/K.rho;
    V(:,2) = u_CRRA(max(((1-tau_l)*K.z2 + rr2.*K.a)/(1+tau_c),1e-12), K.RRA_F)/K.rho;

    if K.sigma_a>0
        e=K.nu*ones(K.I,1); D2=spdiags([e -2*e e],-1:1,K.I,K.I);
        D2(1,1)=-K.nu; D2(1,2)=K.nu; D2(K.I,K.I)=-K.nu; D2(K.I,K.I-1)=K.nu;
    else
        D2 = sparse(K.I,K.I);
    end

    % ---------- Iteración HJB ----------
    Gpc = 0;
    for itv=1:K.maxit_V
        Vprev = V;
        dVf=zeros(K.I,2); dVb=zeros(K.I,2);
        dVf(1:K.I-1,:)=(Vprev(2:K.I,:)-Vprev(1:K.I-1,:))/K.da;
        dVb(2:K.I,:)  =(Vprev(2:K.I,:)-Vprev(1:K.I-1,:))/K.da;

        rr1 = r*ones(K.I,1); rr1(K.a<0)=r+K.theta_I;
        rr2 = r*ones(K.I,1); rr2(K.a<0)=r+K.theta_F;

        c_inf_max = max((K.z1 + rr1(end)*K.a(end) + K.phi*K.z1)/(1+tau_c),1e-12);
        c_for_max = max(((1-tau_l)*K.z2 + rr2(end)*K.a(end))/(1+tau_c),1e-12);
        c_inf_min = max((K.z1 + rr1(1)*K.a(1) + K.phi*K.z1)/(1+tau_c),1e-12);
        c_for_min = max(((1-tau_l)*K.z2 + rr2(1)*K.a(1))/(1+tau_c),1e-12);

        facI=(1+K.psi_G*Gpc)^(K.omegaG*(1-K.RRA_I));
        facF=(1+K.psi_G*Gpc)^(K.omegaG*(1-K.RRA_F));
        dVf(K.I,1)=(1+tau_c)*facI*u_CRRA_prime(c_inf_max,K.RRA_I);
        dVf(K.I,2)=(1+tau_c)*facF*u_CRRA_prime(c_for_max,K.RRA_F);
        dVb(1,1)  =(1+tau_c)*facI*u_CRRA_prime(c_inf_min,K.RRA_I);
        dVb(1,2)  =(1+tau_c)*facF*u_CRRA_prime(c_for_min,K.RRA_F);

        res_inf = K.z1 + rr1.*K.a + K.phi*K.z1;
        res_for = (1-tau_l)*K.z2 + rr2.*K.a;

        cf = [max((1+tau_c)*facI*dVf(:,1),1e-12).^(-1/K.RRA_I), ...
              max((1+tau_c)*facF*dVf(:,2),1e-12).^(-1/K.RRA_F)];
        cb = [max((1+tau_c)*facI*dVb(:,1),1e-12).^(-1/K.RRA_I), ...
              max((1+tau_c)*facF*dVb(:,2),1e-12).^(-1/K.RRA_F)];
        ssf = [res_inf,res_for] - (1+tau_c)*cf;
        ssb = [res_inf,res_for] - (1+tau_c)*cb;
        c0  = [res_inf,res_for];
        If = ssf>0; Ib = ssb<0; I0=(1-If-Ib);
        c  = max( cf.*If + cb.*Ib + c0.*(I0/(1+tau_c)), 1e-12 );

        u = [u_mult(c(:,1),Gpc,K.RRA_I,K.psi_G,K.omegaG), u_mult(c(:,2),Gpc,K.RRA_F,K.psi_G,K.omegaG)];

        X = max(-ssb,0)/K.da; Z = max(ssf,0)/K.da; X(1,:)=0; Z(K.I,:)=0; Y=-(X+Z);
        A1=spdiags(Y(:,1),0,K.I,K.I)+spdiags(X(2:K.I,1),-1,K.I,K.I)+spdiags([0;Z(1:K.I-1,1)],1,K.I,K.I);
        A2=spdiags(Y(:,2),0,K.I,K.I)+spdiags(X(2:K.I,2),-1,K.I,K.I)+spdiags([0;Z(1:K.I-1,2)],1,K.I,K.I);
        A1=A1+D2; A2=A2+D2;
        A = [A1,sparse(K.I,K.I); sparse(K.I,K.I),A2] + K.Aswitch;
        A = A - spdiags(sum(A,2),0,2*K.I,2*K.I);

        Bmat = (1/K.Delta + K.rho)*speye(2*K.I) - A;
        Vst  = Bmat \ ([u(:,1);u(:,2)] + [Vprev(:,1);Vprev(:,2)]/K.Delta);
        V    = [Vst(1:K.I) Vst(K.I+1:2*K.I)];

        if max(max(abs(V-Vprev))) < K.crit_V, break; end
    end

    % ---------- Fokker–Planck estable ----------
    AT=A'; g = ones(2*K.I,1); g = g/(sum(g)*K.da); tauFP = 80*K.da;
    for kk=1:800
        g_new = (speye(2*K.I) - tauFP*AT) \ g;
        g_new = max(g_new,0); g_new = g_new/(sum(g_new)*K.da);
        if norm(g_new-g,inf) < 1e-12, break; end
        g = g_new;
    end
    g = [g(1:K.I) g(K.I+1:2*K.I)];

    % ---------- Agregados ----------
    popI = sum(g(:,1))*K.da; popF = sum(g(:,2))*K.da;
    Ctot = sum(c(:,1).*g(:,1))*K.da + sum(c(:,2).*g(:,2))*K.da;
    Y    = K.z1*popI + K.z2*popF;

    % ---------- Fiscal ----------
    if strcmpi(K.B_mode,'ratio_to_y'), B = K.Bbar*Y; else, B = K.Bbar; end
    rB = r*B; Tl = tau_l*K.z2*popF; Tc = tau_c*Ctot; Tr = K.phi*K.z1*popI;

    if strcmpi(Gmode,'residG')
        Gx = Tl + Tc - Tr - rB;
        if K.clampG0, Gx = max(Gx,0); end
        if isfinite(K.G_cap_ratio), Gx = min(Gx, K.G_cap_ratio*Y); end
        PB = Tl + Tc - Tr - Gx; BB = PB - rB;
    else
        if ~isempty(K.G_target_ratio), Gx = K.G_target_ratio*Y;
        elseif ~isempty(K.G_target_level), Gx = K.G_target_level;
        else, Gx = 0.08*Y;
        end
        PB = Tl + Tc - Tr - Gx; BB = PB - rB;
    end

    % ---------- Ahorro/Exceso ----------
    A_priv = g(:,1)'*K.a*K.da + g(:,2)'*K.a*K.da;
    S = A_priv - B;

    % ---------- Políticas de ahorro ----------
    rr1=r*ones(K.I,1); rr1(K.a<0)=r+K.theta_I;
    rr2=r*ones(K.I,1); rr2(K.a<0)=r+K.theta_F;
    s_pol=[ K.z1 + rr1.*K.a + K.phi*K.z1 - (1+tau_c)*c(:,1), ...
            (1-tau_l)*K.z2 + rr2.*K.a - (1+tau_c)*c(:,2) ];

    % ---------- Stats ----------
    [wI_mean,wF_mean,wT_mean]=deal(wmean(K.a,g(:,1),K.da),wmean(K.a,g(:,2),K.da),wmean(K.a,g(:,1)+g(:,2),K.da));
    [wI_med,wF_med,wT_med]=deal(wmedian(K.a,g(:,1),K.da),wmedian(K.a,g(:,2),K.da),wmedian(K.a,g(:,1)+g(:,2),K.da));
    [giniW_I,giniW_F,giniW_T]=deal(giniW(K.a,g(:,1),K.da),giniW(K.a,g(:,2),K.da),giniW(K.a,g(:,1)+g(:,2),K.da));
    [cI_mean,cF_mean,cT_mean]=deal(wmean(c(:,1),g(:,1),K.da),wmean(c(:,2),g(:,2),K.da),wmean(c(:,1)+c(:,2),g(:,1)+g(:,2),K.da));
    [cI_med,cF_med,cT_med]=deal(wmedian(c(:,1),g(:,1),K.da),wmedian(c(:,2),g(:,2),K.da),wmedian([c(:,1);c(:,2)],[g(:,1);g(:,2)],K.da));
    [giniC_I,giniC_F,giniC_T]=deal(giniX(c(:,1),g(:,1),K.da),giniX(c(:,2),g(:,2),K.da),giniX([c(:,1);c(:,2)],[g(:,1);g(:,2)],K.da));

    % ---------- Residuo HJB ----------
    u_check=[u_mult(c(:,1),Gpc,K.RRA_I,K.psi_G,K.omegaG), u_mult(c(:,2),Gpc,K.RRA_F,K.psi_G,K.omegaG)];
    R = K.rho*[V(:,1);V(:,2)] - [u_check(:,1);u_check(:,2)] - A*[V(:,1);V(:,2)];
    hjb_res = max(abs(R));

    % ---------- Borrowers/Lenders ----------
    idxB=(K.a<0); idxL=(K.a>0);
    fracBorrow_I=sum(g(idxB,1))*K.da / max(sum(g(:,1))*K.da,1e-12);
    fracBorrow_F=sum(g(idxB,2))*K.da / max(sum(g(:,2))*K.da,1e-12);
    fracLend_I  =sum(g(idxL,1))*K.da / max(sum(g(:,1))*K.da,1e-12);
    fracLend_F  =sum(g(idxL,2))*K.da / max(sum(g(:,2))*K.da,1e-12);
    volBorrow_I =sum(g(idxB,1).*K.a(idxB))*K.da;
    volBorrow_F =sum(g(idxB,2).*K.a(idxB))*K.da;
    volLend_I   =sum(g(idxL,1).*K.a(idxL))*K.da;
    volLend_F   =sum(g(idxL,2).*K.a(idxL))*K.da;

    % ---------- Empaqueta ----------
    out.r=r; out.a=K.a; out.g=g; out.c=c; out.s=s_pol;
    out.popI=popI; out.popF=popF; out.Ctot=Ctot; out.Y=Y;
    out.S_residual=S; out.hjb_residual=hjb_res;
    out.Gpc = Gx / max(popI+popF,1e-12);
    out.fiscal=struct('Tl',Tl,'Tc',Tc,'Tr',Tr,'G',Gx,'B',B,'rB',rB,'PB',PB,'BB',BB);
    out.stats=struct('wealth_mean',[wI_mean wF_mean wT_mean],'wealth_median',[wI_med wF_med wT_med], ...
                     'giniW',[giniW_I giniW_F giniW_T],'cons_mean',[cI_mean cF_mean cT_mean], ...
                     'cons_median',[cI_med cF_med cT_med],'giniC',[giniC_I giniC_F giniC_T],'p11',K.p11_rep);
    out.borrowers=struct('fracBorrow',[fracBorrow_I fracBorrow_F],'fracLend',[fracLend_I fracLend_F], ...
                         'volBorrow',[volBorrow_I volBorrow_F],'volLend',[volLend_I volLend_F]);
    out.params = struct('tau_l',tau_l,'tau_c',tau_c);
end

% ------------------- Helpers puros (todas con end) -------------------
function v = get_arg(cfg,f,def)
    if isfield(cfg,f) && ~isempty(cfg.(f)), v=cfg.(f); else, v=def; end
end

function u = u_CRRA(c,sigma)
    c=max(c,1e-12); 
    if abs(sigma-1)<1e-10, u=log(c); else, u=c.^(1-sigma)/(1-sigma); end
end

function up= u_CRRA_prime(c,sigma)
    c=max(c,1e-12); up=c.^(-sigma);
end

function u = u_mult(c,Gpc,sigma,psi,omega)
    c_eff = max(c,1e-12).*(1+psi*Gpc).^omega; 
    u = u_CRRA(c_eff,sigma);
end

function m = wmean(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; 
    if W<=0, m=NaN; else, m=sum(x.*w)*da/W; end
end

function med = wmedian(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; 
    if W<=0, med=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cw=cumsum(ww); 
    k=find(cw>=0.5*W,1,'first'); med=xx(k);
end

function gini = giniW(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; 
    if W<=0, gini=NaN; return; end
    xmin=min(0,min(x)); xs=x-xmin+1e-12; 
    [xx,ix]=sort(xs); ww=w(ix)*da; cumw=cumsum(ww)/W; cumxw=cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end
    L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end

function gini = giniX(x,w,da)
    x=x(:); w=w(:); W=sum(w)*da; 
    if W<=0, gini=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cumw=cumsum(ww)/W; cumxw=cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end
    L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end
