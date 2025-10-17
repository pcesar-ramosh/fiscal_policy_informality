function out = solve_two_type_huggett_fiscal_VATdifferential(cfg, tau_c)
% -------------------------------------------------------------------------
% Dos tipos Huggett con B exógeno y IVA NO NEUTRAL:
% τc_I = chi_I * τc   (informal, cumplimiento parcial) ;  τc_F = τc (formal).
% Incluye opciones de reciclaje y swap τ_l (apagadas por defecto).
% Cierre financiero: A_priv(r) = B (B = Bbar*Y si ratio_to_Y).
% -------------------------------------------------------------------------

% ======= Helper seguro para defaults
arg = @(f,def) get_arg(cfg,f,def);

% ---- Preferencias e impuestos base
sI  = arg('RRA_I',3.40); sF = arg('RRA_F',3.40); rho = arg('rho',0.06);
tau_l0 = arg('tau_l',0.15); phi = arg('phi',0.09);

% IVA no neutral
chi_I = arg('chi_I',0.40);
if nargin<2 || isempty(tau_c), tau_c = arg('tau_c',0.18); end
tau_c_I = chi_I * tau_c;
tau_c_F = tau_c;

% Opcionales (por defecto apagados)
rebate_rule  = arg('rebate_rule','none');   % {'none','phi_from_VAT','lump_sum_pc'}
rebate_share = arg('rebate_share',0.0);     % ∈[0,1]
taul_swap    = arg('taul_swap',0.0);        % κ: τ_l = τ_l0 - κ*(τ_c - τc_base)
tau_c_base   = arg('tau_c_base',0.18);      % base para el swap
tau_l = max(tau_l0 - taul_swap * (tau_c - tau_c_base), 0.0);

% Ingresos
z1 = arg('z1',0.33); z2 = arg('z2',1.00);

% Primas
theta_I = arg('theta_I',0.06); theta_F = arg('theta_F',0.01);

% Activos
I    = arg('I',700);
amin = arg('amin',-2.0*z1); amax = arg('amax',3.0);
a    = linspace(amin,amax,I)'; da=(amax-amin)/(I-1);

% Interés / búsqueda
r     = arg('r_guess',0.03);
rmin  = arg('rmin',0.005); rmax = arg('rmax',0.10);
fix_r = arg('fix_r',0);
maxit_r = arg('maxit_r',1200); crit_S = arg('crit_S',1e-5);

% Switching informal<->formal
p22_bar=arg('p22_bar',0.8155); eta_target=arg('eta_target',0.654);
la2 = -log(p22_bar);                              % F->I
la1 = (1-eta_target)*la2/max(eta_target,1e-12);   % I->F
L1_vec = la1*ones(I,1); L2_vec = la2*ones(I,1); p11_rep = exp(-la1);
Aswitch = [ -spdiags(L1_vec,0,I,I),  spdiags(L1_vec,0,I,I);
             spdiags(L2_vec,0,I,I), -spdiags(L2_vec,0,I,I) ];

% Difusión
sigma_a = arg('sigma_a',0.010);  nu = 0.5*(sigma_a^2)/(da^2);

% Bien público multiplicativo
psi_G  = arg('psi_G',0.08); omegaG = arg('omegaG',0.50);

% HJB numérico
maxit_V = arg('maxit_V',140); crit_V = arg('crit_V',1e-6); Delta = arg('Delta',1400);

% Fisco: B exógeno, G residual
B_mode  = arg('B_mode','ratio_to_Y'); Bbar = arg('Bbar',0.80);
alphaG  = arg('alphaG',0.50); clampG0 = arg('clamp_G_to_zero',true);
G_cap_ratio = arg('G_cap_ratio',1.0);

% ---- v0 coherente (Gpc=0)
rr1 = r*ones(I,1); rr1(a<0)=r+theta_I;
rr2 = r*ones(I,1); rr2(a<0)=r+theta_F;
V = zeros(I,2);
V(:,1) = u_CRRA( max((z1+rr1.*a+phi*z1)/(1+tau_c_I),1e-12), sI )/rho;
V(:,2) = u_CRRA( max(((1-tau_l)*z2+rr2.*a)/(1+tau_c_F),1e-12), sF )/rho;
v0 = V;

Ir = maxit_r; if fix_r, Ir=1; rmin=r; rmax=r; end
S = NaN; Gx=0; Gx_old=0; B=NaN; Tl=NaN; Tc=NaN; Tr=NaN; rB=NaN; PB=NaN; BB=NaN; Gpc=0;
r_path=[]; S_path=[];

for ir=1:Ir
    % ----------------------- HJB -----------------------------------------
    V=v0;
    for it=1:maxit_V
        Vprev=V;

        dVf=zeros(I,2); dVb=zeros(I,2);
        dVf(1:I-1,:)= (Vprev(2:I,:)-Vprev(1:I-1,:))/da;
        dVb(2:I,:)  = (Vprev(2:I,:)-Vprev(1:I-1,:))/da;

        rr1 = r*ones(I,1); rr1(a<0)=r+theta_I;
        rr2 = r*ones(I,1); rr2(a<0)=r+theta_F;

        % Fronteras coherentes con IVA diferencial
        c_inf_max = max((z1 + rr1(end)*a(end) + phi*z1)/(1+tau_c_I),1e-12);
        c_for_max = max(((1-tau_l)*z2 + rr2(end)*a(end))/(1+tau_c_F),1e-12);
        c_inf_min = max((z1 + rr1(1)*a(1) + phi*z1)/(1+tau_c_I),1e-12);
        c_for_min = max(((1-tau_l)*z2 + rr2(1)*a(1))/(1+tau_c_F),1e-12);

        % FOC multiplicativo con IVA por tipo
        facI = (1+psi_G*Gpc)^(omegaG*(1-sI));
        facF = (1+psi_G*Gpc)^(omegaG*(1-sF));
        dVf(I,1) = (1+tau_c_I)*facI*u_CRRA_prime(c_inf_max,sI);
        dVf(I,2) = (1+tau_c_F)*facF*u_CRRA_prime(c_for_max,sF);
        dVb(1,1) = (1+tau_c_I)*facI*u_CRRA_prime(c_inf_min,sI);
        dVb(1,2) = (1+tau_c_F)*facF*u_CRRA_prime(c_for_min,sF);

        % Recursos antes de c
        res_inf = z1 + rr1.*a + phi*z1;
        res_for = (1-tau_l)*z2 + rr2.*a;

        % Políticas (Godunov) con IVA por tipo
        cf = [max((1+tau_c_I)*facI*dVf(:,1),1e-12).^(-1/sI), ...
              max((1+tau_c_F)*facF*dVf(:,2),1e-12).^(-1/sF)];
        cb = [max((1+tau_c_I)*facI*dVb(:,1),1e-12).^(-1/sI), ...
              max((1+tau_c_F)*facF*dVb(:,2),1e-12).^(-1/sF)];
        ssf = [res_inf, res_for] - [ (1+tau_c_I)*cf(:,1), (1+tau_c_F)*cf(:,2) ];
        ssb = [res_inf, res_for] - [ (1+tau_c_I)*cb(:,1), (1+tau_c_F)*cb(:,2) ];
        c0  = [res_inf/(1+tau_c_I), res_for/(1+tau_c_F)];
        If = ssf>0; Ib = ssb<0; I0 = (1-If-Ib);
        c  = max( cf.*If + cb.*Ib + c0.*(I0), 1e-12 );

        % Utilidad de flujo
        u = [u_mult(c(:,1),Gpc,sI,psi_G,omegaG), u_mult(c(:,2),Gpc,sF,psi_G,omegaG)];

        % Generador drift+difusión+switching
        X = max(-ssb,0)/da; Z = max(ssf,0)/da; X(1,:)=0; Z(I,:)=0; Y=-(X+Z);
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

        if sigma_a>0
            e=nu*ones(I,1); D2=spdiags([e -2*e e],-1:1,I,I);
            D2(1,1)=-nu; D2(1,2)=nu; D2(I,I)=-nu; D2(I,I-1)=nu;
            A1=A1+D2; A2=A2+D2;
        end
        A = [A1,sparse(I,I);sparse(I,I),A2]+Aswitch;
        A = A - spdiags(sum(A,2),0,2*I,2*I);

        % Paso implícito
        Bmat = (1/Delta+rho)*speye(2*I) - A;
        Vst  = Bmat \ ([u(:,1);u(:,2)] + [Vprev(:,1);Vprev(:,2)]/Delta);
        V    = [Vst(1:I), Vst(I+1:2*I)];

        if max(max(abs(V-Vprev)))<crit_V, break; end
    end

    % ----------------------- Fokker-Planck --------------------------------
    AT=A'; b0=zeros(2*I,1); i_fix=1; b0(i_fix)=.1; row=zeros(1,2*I); row(i_fix)=1; AT(i_fix,:)=row;
    gg=AT\b0; g=[gg(1:I), gg(I+1:2*I)]; g=g/(sum(g(:))*da);

    % Agregados
    popI=sum(g(:,1))*da; popF=sum(g(:,2))*da; massT=popI+popF;
    C_I=sum(c(:,1).*g(:,1))*da; C_F=sum(c(:,2).*g(:,2))*da; Ctot = C_I + C_F;
    Y = z1*popI + z2*popF;

    % ----------------------- Fiscal ---------------------------------------
    if strcmpi(B_mode,'ratio_to_Y'), B = Bbar*Y; else, B=Bbar; end
    rB = r*B; Tl = tau_l*z2*popF; Tr = phi*z1*popI;

    % IVA por tipo
    TcI = tau_c_I * C_I;
    TcF = tau_c_F * C_F;
    Tc  = TcI + TcF;

    % Reciclaje (opcional)
    rebate = 0;
    if rebate_share>0 && ~strcmpi(rebate_rule,'none')
        % base: IVA que habría con τc_base y mismo chi_I
        TcI_base = (chi_I * tau_c_base) * C_I;
        TcF_base = (tau_c_base)       * C_F;
        Tc_extra = max(Tc - (TcI_base + TcF_base), 0);
        rebate = rebate_share * Tc_extra;
        switch lower(rebate_rule)
            case 'phi_from_vat'
                Tr = Tr + rebate;        % focaliza a informales
            case 'lump_sum_pc'
                % equivalente: restarle a G y sumar ingreso per cápita
                % aquí lo implementamos vía G (ver más abajo)
            otherwise
        end
    end

    % G residual
    G_resid = Tl + Tc - Tr - rB;
    if clampG0, G_resid=max(0,G_resid); end
    if ~isempty(G_cap_ratio) && G_cap_ratio>0, G_resid=min(G_resid, G_cap_ratio*Y); end

    % si reciclaje per cápita, lo sacamos de G (manteniendo PB)
    if rebate>0 && strcmpi(lower(rebate_rule),'lump_sum_pc')
        G_resid = max(G_resid - rebate, 0);
        % Si quisieras sumar a ingreso individual, se haría en res_inf/res_for (no lo activamos por defecto)
    end

    Gx = (1-alphaG)*Gx_old + alphaG*G_resid; Gx_old=Gx;
    Gpc = Gx / max(massT,1e-12);

    PB = Tl + Tc - Tr - Gx;  BB = PB - rB;

    % ------------------- Cierre financiero --------------------------------
    A_priv = g(:,1)'*a*da + g(:,2)'*a*da;
    S = A_priv - B;
    r_path=[r_path; r]; S_path=[S_path; S];

    if ~fix_r
        if S> +crit_S, rmax=r; r=0.5*(r+rmin);
        elseif S<-crit_S, rmin=r; r=0.5*(r+rmax);
        else, break;
        end
        % reset v0 coherente con nuevo r (Gpc=0)
        rr1=r*ones(I,1); rr1(a<0)=r+theta_I;
        rr2=r*ones(I,1); rr2(a<0)=r+theta_F;
        V(:,1)=u_mult(max((z1+rr1.*a+phi*z1)/(1+tau_c_I),1e-12),0,sI,psi_G,omegaG)/rho;
        V(:,2)=u_mult(max(((1-tau_l)*z2+rr2.*a)/(1+tau_c_F),1e-12),0,sF,psi_G,omegaG)/rho;
        v0=V;
    else
        break;
    end
end

% Políticas de ahorro (con IVA por tipo)
rr1=r*ones(I,1); rr1(a<0)=r+theta_I;
rr2=r*ones(I,1); rr2(a<0)=r+theta_F;
s_pol=[ z1 + rr1.*a + phi*z1 - (1+tau_c_I)*c(:,1), ...
        (1-tau_l)*z2 + rr2.*a - (1+tau_c_F)*c(:,2) ];

% ---------- Estadísticos
[wI_mean,wF_mean,wT_mean]=deal(wmean(a,g(:,1),da),wmean(a,g(:,2),da),wmean(a,g(:,1)+g(:,2),da));
[giniW_I,giniW_F,giniW_T]=deal(giniW(a,g(:,1),da),giniW(a,g(:,2),da),giniW(a,g(:,1)+g(:,2),da));
[cI_mean,cF_mean,cT_mean]=deal(wmean(c(:,1),g(:,1),da),wmean(c(:,2),g(:,2),da),wmean(c(:,1)+c(:,2),g(:,1)+g(:,2),da));
[giniC_I,giniC_F,giniC_T]=deal(giniX(c(:,1),g(:,1),da),giniX(c(:,2),g(:,2),da),giniX([c(:,1);c(:,2)],[g(:,1);g(:,2)],da));

% Residuo HJB
u_check=[u_mult(c(:,1),Gpc,sI,psi_G,omegaG), u_mult(c(:,2),Gpc,sF,psi_G,omegaG)];
R = rho*[V(:,1);V(:,2)] - [u_check(:,1);u_check(:,2)] - A*[V(:,1);V(:,2)];
hjb_res = max(abs(R));

% Deudores/ahorradores
idxB=(a<0); idxL=(a>0);
fracBorrow_I=sum(g(idxB,1))*da / max(sum(g(:,1))*da,1e-12);
fracBorrow_F=sum(g(idxB,2))*da / max(sum(g(:,2))*da,1e-12);
fracLend_I  =sum(g(idxL,1))*da / max(sum(g(:,1))*da,1e-12);
fracLend_F  =sum(g(idxL,2))*da / max(sum(g(:,2))*da,1e-12);
volBorrow_I =sum(g(idxB,1).*a(idxB))*da;
volBorrow_F =sum(g(idxB,2).*a(idxB))*da;
volLend_I   =sum(g(idxL,1).*a(idxL))*da;
volLend_F   =sum(g(idxL,2).*a(idxL))*da;

% ---------- Salida
out.r=r; out.a=a; out.g=g; out.c=c; out.s=s_pol;
out.popI=popI; out.popF=popF; out.Ctot=Ctot; out.Y=Y;
out.S_residual=S; out.hjb_residual=hjb_res; out.Gpc=Gpc;
out.fiscal=struct('Tl',Tl,'Tc',Tc,'TcI',TcI,'TcF',TcF,'Tr',Tr,'G',Gx,'B',B,'rB',rB,'PB',PB,'BB',BB);
out.stats=struct('wealth_mean',[wI_mean wF_mean wT_mean], ...
                 'giniW',[giniW_I giniW_F giniW_T], ...
                 'cons_mean',[cI_mean cF_mean cT_mean], ...
                 'giniC',[giniC_I giniC_F giniC_T], 'p11',p11_rep);
out.borrowers=struct('fracBorrow',[fracBorrow_I fracBorrow_F],'fracLend',[fracLend_I fracLend_F], ...
                     'volBorrow',[volBorrow_I volBorrow_F],'volLend',[volLend_I volLend_F]);
out.params.public_good = struct('psi_G',psi_G,'omegaG',omegaG);
out.params.vat = struct('tau_c',tau_c,'tau_c_I',tau_c_I,'tau_c_F',tau_c_F,'chi_I',chi_I, ...
                        'tau_l',tau_l,'tau_l0',tau_l0,'taul_swap',taul_swap);
out.r_path = r_path; out.S_path = S_path;

end

% ===== Helpers =====
function v = get_arg(cfg,f,def)
    if isfield(cfg,f) && ~isempty(cfg.(f)), v = cfg.(f); else, v = def; end
end
function u = u_CRRA(c,sigma), c=max(c,1e-12); if abs(sigma-1)<1e-10, u=log(c); else, u=c.^(1-sigma)/(1-sigma); end, end
function up= u_CRRA_prime(c,sigma), c=max(c,1e-12); up=c.^(-sigma); end
function u = u_mult(c,Gpc,sigma,psi,omega), c_eff = max(c,1e-12).*(1+psi*Gpc).^omega; u=u_CRRA(c_eff,sigma); end
function m = wmean(x,w,da), x=x(:); w=w(:); W=sum(w)*da; if W<=0, m=NaN; else, m=sum(x.*w)*da/W; end, end
function gini=giniW(x,w,da), x=x(:); w=w(:); W=sum(w)*da; if W<=0, gini=NaN; return; end
    xmin=min(0,min(x)); xs=x-xmin+1e-12; [xx,ix]=sort(xs); ww=w(ix)*da; cumw=cumsum(ww)/W; cumxw=cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end, L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end
function gini=giniX(x,w,da), x=x(:); w=w(:); W=sum(w)*da; if W<=0, gini=NaN; return; end
    [xx,ix]=sort(x); ww=w(ix)*da; cumw=cumsum(ww)/W; cumxw=cumsum(xx.*ww);
    if cumxw(end)<=0, gini=NaN; return; end, L=cumxw/cumxw(end); gini=1-2*trapz(cumw,L);
end
