% ========================= main_vat_nonneutral.m =========================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------------------- BASE configuration -------------------------------
cfg0 = struct();

% Preferencias e impuestos base
cfg0.RRA_I = 3.40; cfg0.RRA_F = 3.40;
cfg0.rho   = 0.08;
cfg0.tau_l = 0.15;          % laboral (puede moverse con swap si activas abajo)
cfg0.tau_c = 0.18;          % sólo informativo; el solver usará el del escenario
cfg0.phi   = 0.09;

% Ingresos
cfg0.z1 = 0.33; cfg0.z2 = 1.00;

% Primas de endeudamiento
cfg0.theta_I = 0.06; cfg0.theta_F = 0.01;

% Mercado laboral informal (objetivo de masa)
cfg0.eta_target = 0.654;   cfg0.p22_bar = 0.8155;

% Grid de activos y numérico
cfg0.I = 700; cfg0.amax = 3.0; cfg0.amin = -2.0*cfg0.z1;
cfg0.r_guess=0.03; cfg0.rmin=0.005; cfg0.rmax=0.10;
cfg0.maxit_V=140; cfg0.crit_V=1e-6; cfg0.Delta=1400;
cfg0.maxit_r=1500; cfg0.crit_S=1e-5; cfg0.fix_r=0;

% Bien público (multiplicativo)
cfg0.psi_G = 0.08; cfg0.omegaG = 0.50; cfg0.report_G_effects = 0;

% Difusión (Fokker-Planck)
cfg0.sigma_a = 0.010;

% Gobierno / Deuda
cfg0.B_mode  = 'ratio_to_Y';  % o 'level'
cfg0.Bbar    = 0.80;          % B/Y
cfg0.alphaG  = 0.50;
cfg0.clamp_G_to_zero = true;
cfg0.G_cap_ratio = 1;         % ~ sin tope efectivo

% -------------------- NO NEUTRALIDAD DEL IVA (CLAVE) ---------------------
cfg0.chi_I = 0.40;      % Cumplimiento informal: τc_I = chi_I * τc   (chi_I<1)
% Opcionales (apagados por defecto):
cfg0.rebate_rule  = 'none';   % {'none','phi_from_VAT','lump_sum_pc'}
cfg0.rebate_share = 0.0;      % fracción del IVA extra reciclada
cfg0.taul_swap    = 0.0;      % κ para bajar τ_l cuando sube τ_c (swap recaudación-neutral aprox.)

% --------------------------- Escenarios τ_c -------------------------------
tau_vec = [0.18 0.20 0.22 0.26];  % puedes ajustar
labels  = arrayfun(@(x)sprintf('\\tau_c=%.2f',x), tau_vec, 'uni', 0);

S = struct('name',{},'tau_c',{},'cfg',{},'sol',{});
for k=1:numel(tau_vec)
    S(k).name  = sprintf('VATnn_%02d', round(100*tau_vec(k)));
    S(k).tau_c = tau_vec(k);
    S(k).cfg   = cfg0;
end

% ----------------------------- Resolver ----------------------------------
for k=1:numel(S)
    tic;
    S(k).sol = solve_two_type_huggett_fiscal_VATdifferential(S(k).cfg, S(k).tau_c);
    s = S(k).sol;
    fprintf('[%s] τc=%.2f | r*=%.4f | η=%.3f | B=%.3f | TcI=%.3f TcF=%.3f | GiniW=%.3f GiniC=%.3f (%.2fs)\n', ...
        S(k).name, S(k).tau_c, s.r, s.popI/(s.popI+s.popF), ...
        s.fiscal.B, s.fiscal.TcI, s.fiscal.TcF, ...
        s.stats.giniW(3), s.stats.giniC(3), toc);

    % ---- CSVs por escenario ----
    a=s.a; g=s.g; da=a(2)-a(1);
    A_I=sum(g(:,1).*a)*da; A_F=sum(g(:,2).*a)*da; A_priv=A_I+A_F;
    scen=string(S(k).name);

    T_fiscal = table(scen, s.fiscal.Tl, s.fiscal.Tc, s.fiscal.TcI, s.fiscal.TcF, ...
        s.fiscal.Tl+s.fiscal.Tc, s.fiscal.Tr, s.fiscal.G, s.fiscal.rB, ...
        s.fiscal.PB, s.fiscal.B, s.fiscal.BB, ...
        'VariableNames', {'scenario','Tl','Tc','TcI','TcF','Rev_total','Tr','G','rB','PB','B','BB'});
    writetable(T_fiscal, fullfile(outdir_tabs, sprintf('fiscal_%s.csv', S(k).name)));

    T_assets = table(scen, A_priv, s.fiscal.B, ...
        'VariableNames', {'scenario','A_private','B_public'});
    writetable(T_assets, fullfile(outdir_tabs, sprintf('assets_%s.csv', S(k).name)));

    st=s.stats;
    T_stats = table(scen, s.r, s.popI, s.popF, s.Y, s.Ctot, ...
        st.wealth_mean(1),st.wealth_mean(2),st.wealth_mean(3), ...
        st.giniW(1),st.giniW(2),st.giniW(3), ...
        st.cons_mean(1),st.cons_mean(2),st.cons_mean(3), ...
        st.giniC(1),st.giniC(2),st.giniC(3), st.p11, ...
        'VariableNames', {'scenario','r','popI','popF','Y','C', ...
        'wmean_I','wmean_F','wmean_T','giniW_I','giniW_F','giniW_T', ...
        'cmean_I','cmean_F','cmean_T','giniC_I','giniC_F','giniC_T','p11_rep'});
    writetable(T_stats, fullfile(outdir_tabs, sprintf('hh_stats_%s.csv', S(k).name)));
end

% ------------------- Resumen y figuras clave -----------------------------
RES = cell(numel(S), 20);
for k=1:numel(S)
    s=S(k).sol; da=s.a(2)-s.a(1);
    A_priv = sum((s.g(:,1)+s.g(:,2)).*s.a)*da;
    RES(k,:) = { labels{k}, S(k).tau_c, s.r, A_priv, s.fiscal.B, s.Y, s.Ctot, ...
                 s.fiscal.TcI, s.fiscal.TcF, s.fiscal.Tc, s.fiscal.Tl, s.fiscal.Tr, ...
                 s.fiscal.G, s.fiscal.PB, s.fiscal.BB, ...
                 s.stats.giniW(1), s.stats.giniW(2), s.stats.giniW(3), ...
                 s.stats.giniC(3), s.popI/(s.popI+s.popF) };
end
Tsum = cell2table(RES, 'VariableNames', ...
 {'scenario','tau_c','r','A_priv','B','Y','C','TcI','TcF','Tc','Tl','Tr','G','PB','BB', ...
  'giniW_I','giniW_F','giniW_T','giniC_T','eta_real'});
writetable(Tsum, fullfile(outdir_tabs,'vat_nonneutral_summary.csv'));
disp(Tsum);

% -------- Lorenz total (riqueza) y consumo
cols = lines(numel(S));
fig=figure('Name','Lorenz wealth (total) - nonneutral VAT'); hold on;
for k=1:numel(S)
    [~,Lw,cumPop]=lorenz_from_density(S(k).sol.a, S(k).sol.g);
    plot(cumPop,Lw,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share');
legend(Tsum.scenario,'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'nn_lorenz_wealth_total'));

fig=figure('Name','Lorenz consumption (total) - nonneutral VAT'); hold on;
for k=1:numel(S)
    [~,Lc,cumPopC]=lorenz_from_values(S(k).sol.a, S(k).sol.g, S(k).sol.c(:,1)+S(k).sol.c(:,2));
    plot(cumPopC,Lc,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share');
legend(Tsum.scenario,'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'nn_lorenz_consumption_total'));

% -------- Políticas s(a) (para ver separación I vs F)
fig=figure('Name','Policies: savings by τ_c (nonneutral)'); hold on;
for k=1:numel(S)
    a=S(k).sol.a; s=S(k).sol.s;
    plot(a,s(:,1),'LineWidth',2,'Color',cols(k,:),'LineStyle','-');   % I
    plot(a,s(:,2),'LineWidth',2,'Color',cols(k,:),'LineStyle','--');  % F
end
grid on; xlabel('Assets a'); ylabel('Savings s(a)');
title('I sólido / F punteado'); legend(Tsum.scenario,'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'nn_policy_savings'));

% =============================== Helpers =================================
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
end
function [gT,L,cumPop]=lorenz_from_density(a,g)
    da=a(2)-a(1); gT=g(:,1)+g(:,2); W=sum(gT)*da; [as,ix]=sort(a); gs=gT(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as-min(0,min(as))+1e-12; cumWealth=cumsum(wealth_pos.*gs); L=cumWealth/max(cumWealth(end),1e-12);
end
function [vals,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w); vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix);
    cumPop=cumsum(w_s)/W; cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end
