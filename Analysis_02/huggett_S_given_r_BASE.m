function [S_by_agent, S_total, g_out, c_out, a_grid, G_val, Ec_total] = ...
    huggett_S_given_r_BASE(share_I, sI_vec, sF_vec, trans, cfg, r_fixed)
% HJB+FP con r = r_fixed y política fiscal BASE:
%  tau_l solo FORMAL; tau_c ambos; Transfer INFORMAL = phi*z1.
%  Bien público G se calcula por presupuesto público y entra a utilidad
%  como +xiG*log(G+eps), pero NO afecta políticas (aditivo).

% --- Unpack ---
I = cfg.I; amax = cfg.amax; amin = cfg.amin_abs;
rho = cfg.rho; z1 = cfg.z1; z2 = cfg.z2;
tau_l = cfg.tau_l; tau_c = cfg.tau_c; theta = cfg.theta;
xiG = cfg.xiG; epsG = cfg.epsG;

% Ocupación intensities consistent with share_I
la1 = trans.la1; la2 = trans.la2;   % from main (consistent with share_I)
Aswitch = @(I) [-speye(I)*la1,  speye(I)*la1; ...
                 speye(I)*la2, -speye(I)*la2];

% Grid
a_grid = linspace(amin, amax, I)'; da = (amax-amin)/(I-1);

J = numel(sI_vec);
S_by_agent = zeros(J,1);
g_out = cell(J,1); c_out = cell(J,1);

% Nota: Como G no afecta políticas, resolvemos HJB sin G y lo añadimos a u ex-post
for jj=1:J
    sI = sI_vec(jj); sF = sF_vec(jj);
    A_sw = Aswitch(I);

    % Inicialización de V
    r = r_fixed; rr = r + theta*(a_grid<0);
    cI0 = max((1-0)*z1 + rr.*a_grid + cfg.phi*z1, 1e-10);
    cF0 = max((1-tau_l)*z2 + r*a_grid,            1e-10);
    V = zeros(I,2);
    V(:,1) = (cI0.^(1-sI))/(1-sI)/rho;
    V(:,2) = (cF0.^(1-sF))/(1-sF)/rho;

    % Iteración HJB
    Delta = 50; maxit = 200; crit = 1e-6;
    for n=1:maxit
        dVf = zeros(I,2); dVb=zeros(I,2);

        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        cI_sup = max((1-0)*z1 + r*amax + cfg.phi*z1, 1e-10);
        cF_sup = max((1-tau_l)*z2 + r*amax,          1e-10);
        dVf(I,1) = cI_sup^(-sI)/(1+tau_c);
        dVf(I,2) = cF_sup^(-sF)/(1+tau_c);

        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        cI_inf = max((1-0)*z1 + (r+theta)*amin + cfg.phi*z1, 1e-10);
        cF_inf = max((1-tau_l)*z2 + r*amin,                  1e-10);
        dVb(1,1) = cI_inf^(-sI)/(1+tau_c);
        dVb(1,2) = cF_inf^(-sF)/(1+tau_c);

        rr = r + theta*(a_grid<0);
        resI = (1-0)*z1 + rr.*a_grid + cfg.phi*z1;
        resF = (1-tau_l)*z2 + r*a_grid;

        cf = [ ((1+tau_c)*max(dVf(:,1),1e-10)).^(-1/sI) , ...
               ((1+tau_c)*max(dVf(:,2),1e-10)).^(-1/sF) ];
        cb = [ ((1+tau_c)*max(dVb(:,1),1e-10)).^(-1/sI) , ...
               ((1+tau_c)*max(dVb(:,2),1e-10)).^(-1/sF) ];

        ssf = [resI,resF] - cf;
        ssb = [resI,resF] - cb;
        If = ssf>0; Ib = ssb<0; I0 = ~(If|Ib);
        c  = cf.*If + cb.*Ib + [resI,resF].*I0;

        U1 = (max(c(:,1),1e-12).^(1-sI))/(1-sI);
        U2 = (max(c(:,2),1e-12).^(1-sF))/(1-sF);
        u  = [U1,U2];  % (sin G; G se añade ex-post a bienestar)

        X = -min(ssb,0)/da;  Y = -max(ssf,0)/da + min(ssb,0)/da;  Z = max(ssf,0)/da;
        A1 = spdiags(Y(:,1),0,I,I) + spdiags(X(2:I,1),-1,I,I) + spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I) + spdiags(X(2:I,2),-1,I,I) + spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1, sparse(I,I); sparse(I,I), A2] + A_sw;

        B = (1/Delta + rho)*speye(2*I) - A;
        Vst = B \ ([u(:,1);u(:,2)] + [V(:,1);V(:,2)]/Delta);
        Vn  = [Vst(1:I), Vst(I+1:2*I)];

        if max(max(abs(Vn-V)))<crit, V=Vn; break; end
        V=Vn;
    end

    % FP con normalización directa
    AT = A'; bb = zeros(2*I,1);
    AT(1,:) = ones(1, 2*I)*da; bb(1)=1;
    gg = AT \ bb; gg = max(gg,0); gg = gg/(sum(gg)*da);
    g  = [gg(1:I), gg(I+1:2*I)];

    % Guardar
    g_out{jj}=g; c_out{jj}=c;

    % Exceso por agente
    S_by_agent(jj) = (g(:,1)'*a_grid + g(:,2)'*a_grid)*da;
end

S_total = sum(S_by_agent);

% --- G y E[c] agregados (para reporte y utilidad ex-post) ---
Ec_total = 0;
for j=1:J
    g  = g_out{j}; c = c_out{j};
    Ec_total = Ec_total + sum(g(:,1).*c(:,1) + g(:,2).*c(:,2))*da;
end
pi_I = la2/(la1+la2); pi_F = la1/(la1+la2);
revenues = tau_l*z2*pi_F + tau_c*Ec_total;
transfs  = cfg.phi*z1*pi_I;
G_val    = max(0.0, revenues - transfs); %#ok<NASGU>  % entra en utilidad, no en políticas
% (Si quieres medir bienestar, añade +xiG*log(G_val+epsG) a u al final.)
end
