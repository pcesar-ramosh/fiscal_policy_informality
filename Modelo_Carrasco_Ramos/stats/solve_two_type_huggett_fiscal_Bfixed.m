function out = solve_two_type_huggett_fiscal_Bfixed(cfg)
% -------------------------------------------------------------------------
% Huggett 2 tipos con B exógeno (en 'level' o 'ratio_to_Y').
% Cierre financiero: bisección en r para A_priv(r)=B.
% Cierre fiscal exacto en SS: G = Tl + Tc - Tr - rB  => BB=0, PB=rB.
% Fokker–Planck: esquema implícito en pseudo-tiempo (estable, sin clamps).
% -------------------------------------------------------------------------

arg = @(f,def) get_arg(cfg,f,def);

% Preferencias e impuestos
sI = arg('RRA_I',3.40); sF = arg('RRA_F',3.40); rho = arg('rho',0.08);
tau_l = arg('tau_l',0.15); tau_c = arg('tau_c',0.18); phi = arg('phi',0.09);

% Ingresos
z1 = arg('z1',0.33); z2 = arg('z2',1.00);

% Spreads
theta_I = arg('theta_I',0.06); theta_F = arg('theta_F',0.01);

% Grid
I    = arg('I',700);
amin = arg('amin',-2.0*z1); amax = arg('amax',3.0);
a    = linspace(amin,amax,I)'; da=(amax-amin)/(I-1);

% Interés / búsqueda
r_guess = arg('r_guess',0.03);
rmin    = arg('rmin',0.003);
rmax    = arg('rmax',0.12);
fix_r   = arg('fix_r',0);
maxit_r = arg('maxit_r',100);
crit_S  = arg('crit_S',1e-5);

% Switching (Markov)
p22_bar = arg('p22_bar',0.8155);
eta_t   = arg('eta_target',0.654);
la2 = -log(p22_bar);
la1 = (1-eta_t)*la2/max(eta_t,1e-12);
L1_vec = la1*ones(I,1); L2_vec = la2*ones(I,1); p11_rep = exp(-la1);
Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I); ...
             spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];

% Difusión
sigma_a = arg('sigma_a',0.007);  nu = 0.5*(sigma_a^2)/(da^2);

% Bien público en utilidad
psi_G  = arg('psi_G',0.08);
omegaG = arg('omegaG',0.50);
report_G = arg('report_G_effects',0);

% HJB numérico
maxit_V = arg('maxit_V',160); crit_V = arg('crit_V',1e-6); Delta = arg('Delta',1400);

% Gobierno / deuda
B_mode  = arg('B_mode','level');   % 'level' | 'ratio_to_Y'
Bbar    = arg('Bbar',0.25);

% ------------------------------ Solver en r -------------------------------
if fix_r
    out = solve_given_r(r_guess);
    return;
else
    [Smin,~] = eval_S(rmin);
    [Smax,~] = eval_S(rmax);
    it_expand=0;
    while Smin*Smax>0 && it_expand<6
        rmin = max(1e-4,0.7*rmin);
        rmax = 1.4*rmax;
        [Smin,~] = eval_S(rmin);
        [Smax,~] = eval_S(rmax);
        it_expand=it_expand+1;
    end
    if Smin*Smax>0
        if abs(Smin)<abs(Smax), out=solve_given_r(rmin); else, out=solve_given_r(rmax); end
        warning('Bisection: no sign change; returning boundary solution.');
        return;
    end
    r_low=rmin; S_low=Smin; r_high=rmax; S_high=Smax;
    best=[]; r_mid=r_guess;
    for it=1:maxit_r
        r_mid = 0.5*(r_low+r_high);
        [S_mid, sol_mid] = eval_S(r_mid);
        if isempty(best) || abs(S_mid) < abs(best.S_residual), best = sol_mid; end
        if abs(S_mid) < crit_S, break; end
        if S_low*S_mid<=0, r_high=r_mid; S_high=S_mid;
        else, r_low=r_mid; S_low=S_mid; end
    end
    out = best; return;
end

% --------------------------- subfunciones locales -------------------------
    function [Sval, solr] = eval_S(rr)
        solr = solve_given_r(rr);
        Sval = solr.S_residual;
    end

    function out_r = solve_given_r(r)
        % V inicial (Gpc=0)
        rr1 = r*ones(I,1); rr1(a<0)=r+theta_I;
        rr2 = r*ones(I,1); rr2(a<0)=r+theta_F;
        V = zeros(I,2);
        V(:,1) = u_CRRA( max((z1+rr1.*a+phi*z1)/(1+tau_c),1e-12), sI )/rho;
        V(:,2) = u_CRRA( max(((1-tau_l)*z2+rr2.*a)/(1+tau_c),1e-12), sF )/rho;

        % Laplaciano
        if sigma_a>0
            e=nu*ones(I,1);
            D2=spdiags([e -2*e e],-1:1,I,I); D2(1,1)=-nu; D2(1,2)=nu; D2(I,I)=-nu; D2(I,I-1)=nu;
        else
            D2=sparse(I,I);
        end

        % HJB iterativo implícito
        Gpc=0;
        for it=1:maxit_V
            Vprev=V;

            dVf=zeros(I,2); dVb=zeros(I,2);
            dVf(1:I-1,:)=(Vprev(2:I,:)-Vprev(1:I-1,:))/da;
            dVb(2:I,:)  =(Vprev(2:I,:)-Vprev(1:I-1,:))/da;

            rr1 = r*ones(I,1); rr1(a<0)=r+theta_I;
            rr2 = r*ones(I,1); rr2(a<0)=r+theta_F;

            c_inf_max = max((z1 + rr1(end)*a(end) + phi*z1)/(1+tau_c),1e-12);
            c_for_max = max(((1-tau_l)*z2 + rr2(end)*a(end))/(1+tau_c),1e-12);
            c_inf_min = max((z1 + rr1(1)*a(1) + phi*z1)/(1+tau_c),1e-12);
            c_for_min = max(((1-tau_l)*z2 + rr2(1)*a(1))/(1+tau_c),1e-12);

            facI = (1+psi_G*Gpc)^(omegaG*(1-sI));
            facF = (1+psi_G*Gpc)^(omegaG*(1-sF));

            dVf(I,1) = (1+tau_c)*facI*u_CRRA_prime(c_inf_max,sI);
            dVf(I,2) = (1+tau_c)*facF*u_CRRA_prime(c_for_max,sF);
            dVb(1,1) = (1+tau_c)*facI*u_CRRA_prime(c_inf_min,sI);
            dVb(1,2) = (1+tau_c)*facF*u_CRRA_prime(c_for_min,sF);

            res_inf = z1 + rr1.*a + phi*z1;
            res_for = (1-tau_l)*z2 + rr2.*a;

            cf = [max((1+tau_c)*facI*dVf(:,1),1e-12).^(-1/sI), ...
                  max((1+tau_c)*facF*dVf(:,2),1e-12).^(-1/sF)];
            cb = [max((1+tau_c)*facI*dVb(:,1),1e-12).^(-1/sI), ...
                  max((1+tau_c)*facF*dVb(:,2),1e-12).^(-1/sF)];
            ssf = [res_inf, res_for] - (1+tau_c)*cf;
            ssb = [res_inf, res_for] - (1+tau_c)*cb;
            c  = max( cf.*(ssf>0) + cb.*(ssb<0) + [res_inf,res_for].*((1-(ssf>0)-(ssb<0))/(1+tau_c)), 1e-12 );

            u = [u_mult(c(:,1),Gpc,sI,psi_G,omegaG), u_mult(c(:,2),Gpc,sF,psi_G,omegaG)];

            X = max(-ssb,0)/da; Z = max(ssf,0)/da; X(1,:)=0; Z(I,:)=0; Y=-(X+Z);
            A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
            A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
            if sigma_a>0, A1=A1+D2; A2=A2+D2; end
            A = [A1,sparse(I,I); sparse(I,I),A2] + Aswitch;
            A = A - spdiags(sum(A,2),0,2*I,2*I);

            Bmat = (1/Delta+rho)*speye(2*I) - A;
            Vst  = Bmat \ ([u(:,1);u(:,2)] + [Vprev(:,1);Vprev(:,2)]/Delta);
            V    = [Vst(1:I), Vst(I+1:2*I)];

            if max(max(abs(V-Vprev)))<crit_V, break; end
        end

        % Fokker–Planck implícito en pseudo-tiempo
        AT = A'; g = ones(2*I,1); g = g/(sum(g)*da);
        tauFP = 80*da;
        for k=1:1000
            g_new = (speye(2*I) - tauFP*AT) \ g;
            g_new = max(g_new,0); g_new = g_new/(sum(g_new)*da);
            if norm(g_new-g,inf)<1e-12, break; end
            g = g_new;
        end
        gg=[g(1:I), g(I+1:2*I)];

        % Agregados
        popI = sum(gg(:,1))*da; popF=sum(gg(:,2))*da; massT=popI+popF;
        Ctot = sum(c(:,1).*gg(:,1))*da + sum(c(:,2).*gg(:,2))*da;
        Y    = z1*popI + z2*popF;

        % Deuda pública
        if strcmpi(B_mode,'ratio_to_Y'), B = Bbar*Y; else, B=Bbar; end
        rB = r*B;

        % Fiscal (cierre exacto)
        Tl = tau_l*z2*popF;
        Tc = tau_c*Ctot;
        Tr = phi*z1*popI;
        Gx = Tl + Tc - Tr - rB;       % PB = rB; BB=0
        Gpc = Gx / max(massT,1e-12);
        PB = rB; BB = 0.0;

        % Cierre financiero
        A_priv = gg(:,1)'*a*da + gg(:,2)'*a*da;
        S = A_priv - B;

        % Políticas de ahorro
        rr1=r*ones(I,1); rr1(a<0)=r+theta_I;
        rr2=r*ones(I,1); rr2(a<0)=r+theta_F;
        s_pol = [ z1 + rr1.*a + phi*z1 - (1+tau_c)*c(:,1), ...
                  (1-tau_l)*z2 + rr2.*a - (1+tau_c)*c(:,2) ];

        % Stats de chequeo
        u_check=[u_mult(c(:,1),Gpc,sI,psi_G,omegaG), u_mult(c(:,2),Gpc,sF,psi_G,omegaG)];
        R = rho*[V(:,1);V(:,2)] - [u_check(:,1);u_check(:,2)] - A*[V(:,1);V(:,2)];
        hjb_res = max(abs(R));

        % Deudores/ahorradores
        idxB=(a<0); idxL=(a>0);
        fracBorrow_I=sum(gg(idxB,1))*da / max(sum(gg(:,1))*da,1e-12);
        fracBorrow_F=sum(gg(idxB,2))*da / max(sum(gg(:,2))*da,1e-12);
        fracLend_I  =sum(gg(idxL,1))*da / max(sum(gg(:,1))*da,1e-12);
        fracLend_F  =sum(gg(idxL,2))*da / max(sum(gg(:,2))*da,1e-12);
        volBorrow_I =sum(gg(idxB,1).*a(idxB))*da;
        volBorrow_F =sum(gg(idxB,2).*a(idxB))*da;
        volLend_I   =sum(gg(idxL,1).*a(idxL))*da;
        volLend_F   =sum(gg(idxL,2).*a(idxL))*da;

        % Output
        out_r.r = r; out_r.a=a; out_r.g=gg; out_r.c=c; out_r.s=s_pol;
        out_r.popI=popI; out_r.popF=popF; out_r.Ctot=Ctot; out_r.Y=Y;
        out_r.S_residual=S; out_r.hjb_residual=hjb_res; out_r.Gpc=Gpc;
        out_r.fiscal=struct('Tl',Tl,'Tc',Tc,'Tr',Tr,'G',Gx,'B',B,'rB',rB,'PB',PB,'BB',BB);
        out_r.stats=struct('p11',p11_rep);  % (dejamos básicos; el main calcula todo)
        out_r.borrowers=struct('fracBorrow',[fracBorrow_I fracBorrow_F], ...
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
        out = out_r;
    end
end

% ============================== Helpers ===================================
function v = get_arg(cfg,f,def)
    if isfield(cfg,f) && ~isempty(cfg.(f)), v = cfg.(f); else, v = def; end
end

function u = u_CRRA(c,sigma)
    c=max(c,1e-12);
    if abs(sigma-1)<1e-10, u = log(c); else, u = c.^(1-sigma)/(1-sigma); end
end

function up = u_CRRA_prime(c,sigma)
    c=max(c,1e-12); up=c.^(-sigma);
end

function u = u_mult(c,Gpc,sigma,psi,omega)
    c_eff = max(c,1e-12).*(1+psi*Gpc).^omega;
    u = u_CRRA(c_eff,sigma);
end
