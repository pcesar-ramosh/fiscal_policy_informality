% =================== main_VAT_compare_B_endog.m ===================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------- Parámetros “limpios” ----------
cfg = struct();
cfg.RRA_I = 3.40; cfg.RRA_F = 3.40;
cfg.rho   = 0.08;
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;         % IVA base
cfg.phi   = 0.09;
cfg.z1 = 0.33; cfg.z2 = 1.00;

cfg.theta_I = 0.06; cfg.theta_F = 0.01;

cfg.eta_target = 0.654; cfg.p22_bar = 0.8155;

cfg.I = 700; cfg.amin = -2.0*cfg.z1; cfg.amax = 3.0;

cfg.r_guess = 0.03; cfg.rmin = 0.005; cfg.rmax = 0.10;
cfg.maxit_r = 80; cfg.crit_S = 1e-5;

cfg.maxit_V = 160; cfg.crit_V = 1e-6; cfg.Delta = 1400;
cfg.sigma_a = 0.007;

cfg.psi_G   = 0.08; cfg.omegaG = 0.50; cfg.report_G_effects = 0;

% ---- Nuevo: cierre fiscal ----
%  Modo 1 (referencia): primero resolvemos con G residual para medir gY
cfg_ref = cfg; cfg_ref.fiscal_mode = 'G_resid_B_fixed';
cfg_ref.B_mode = 'level'; cfg_ref.Bbar = 0.25;      % usado SOLO en el pase de referencia
ref = solve_two_type_huggett_fiscal_Bfixed(cfg_ref);
gY = ref.fiscal.G / max(ref.Y,1e-12);               % ratio de G/PIB estable a igualar

%  Modo 2 (comparativo): G fijo a gY y B endógeno
cfg.fiscal_mode = 'G_fixed_B_endog';
cfg.G_target_ratio = gY;                             % MISMO gY en ambos escenarios

paper_style();

% --------------------- Escenario BASE ---------------------
base = solve_two_type_huggett_fiscal_Bfixed(cfg);

% --------------------- Escenario VAT-UP -------------------
d_tau_c = 0.05;                 % +5 pp de IVA (ajusta aquí)
cfgVAT = cfg; cfgVAT.tau_c = cfg.tau_c + d_tau_c;
vat = solve_two_type_huggett_fiscal_Bfixed(cfgVAT);

% ------------------- Reporte en consola -------------------
disp('== BASE (G fixed, B endog) ==');   print_summary(base);
disp('== VAT↑  (G fixed, B endog) ==');  print_summary(vat);

% --------------------------- TABLAS ------------------------
% Fiscal
T_fiscal = table( ...
    ["BASE"; "VAT_UP"], ...
    [base.fiscal.Tl; vat.fiscal.Tl], [base.fiscal.Tc; vat.fiscal.Tc], ...
    [base.fiscal.Tl+base.fiscal.Tc; vat.fiscal.Tl+vat.fiscal.Tc], ...
    [base.fiscal.Tr; vat.fiscal.Tr], [base.fiscal.G; vat.fiscal.G], ...
    [base.fiscal.rB; vat.fiscal.rB], [base.fiscal.PB; vat.fiscal.PB], ...
    [base.fiscal.B; vat.fiscal.B],   [base.fiscal.BB; vat.fiscal.BB], ...
'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers','public_good','debt_serv','primary_bal','debt_stock','global_bal'});
writetable(T_fiscal, fullfile(outdir_tabs,'vat_Bendog_fiscal.csv'));

% Mercado de activos (al equilibrio)
daB = base.a(2)-base.a(1);
A_I_base = sum(base.g(:,1).*base.a)*daB;  A_F_base = sum(base.g(:,2).*base.a)*daB;
daV = vat.a(2)-vat.a(1);
A_I_vat  = sum(vat.g(:,1).*vat.a)*daV;    A_F_vat  = sum(vat.g(:,2).*vat.a)*daV;

T_assets = table( ...
    ["BASE"; "VAT_UP"], ...
    [A_I_base; A_I_vat], [A_F_base; A_F_vat], [A_I_base + A_F_base; A_I_vat + A_F_vat], ...
    [base.fiscal.B; vat.fiscal.B], ...
'VariableNames', {'scenario','A_I','A_F','A_private','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'vat_Bendog_asset_market.csv'));

% Estadísticas
Sb = base.stats; Sv = vat.stats;
T_stats = table( ...
    ["BASE"; "VAT_UP"], ...
    [base.r; vat.r], [base.popI; vat.popI], [base.popF; vat.popF], [base.Y; vat.Y], [base.Ctot; vat.Ctot], ...
    [Sb.wealth_mean(1); Sv.wealth_mean(1)], [Sb.wealth_mean(2); Sv.wealth_mean(2)], [Sb.wealth_mean(3); Sv.wealth_mean(3)], ...
    [Sb.wealth_median(1); Sv.wealth_median(1)], [Sb.wealth_median(2); Sv.wealth_median(2)], [Sb.wealth_median(3); Sv.wealth_median(3)], ...
    [Sb.giniW(1); Sv.giniW(1)], [Sb.giniW(2); Sv.giniW(2)], [Sb.giniW(3); Sv.giniW(3)], ...
    [Sb.cons_mean(1); Sv.cons_mean(1)], [Sb.cons_mean(2); Sv.cons_mean(2)], [Sb.cons_mean(3); Sv.cons_mean(3)], ...
    [Sb.cons_median(1); Sv.cons_median(1)], [Sb.cons_median(2); Sv.cons_median(2)], [Sb.cons_median(3); Sv.cons_median(3)], ...
    [Sb.giniC(1); Sv.giniC(1)], [Sb.giniC(2); Sv.giniC(2)], [Sb.giniC(3); Sv.giniC(3)], ...
'VariableNames', {'scenario','r','popI','popF','Y','Ctot', ...
    'wealth_mean_I','wealth_mean_F','wealth_mean_T', ...
    'wealth_med_I','wealth_med_F','wealth_med_T', ...
    'giniW_I','giniW_F','giniW_T', ...
    'cons_mean_I','cons_mean_F','cons_mean_T', ...
    'cons_med_I','cons_med_F','cons_med_T', ...
    'giniC_I','giniC_F','giniC_T'});
writetable(T_stats, fullfile(outdir_tabs,'vat_Bendog_stats.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

% --------------------------- GRÁFICOS ----------------------
paper_style();

% Consumo / Ahorro
fig=figure('Name','c(a) por tipo: BASE vs VAT↑');
plot(base.a,base.c(:,1),'LineWidth',2); hold on; plot(base.a,base.c(:,2),'LineWidth',2);
plot(vat.a, vat.c(:,1),'--','LineWidth',2); plot(vat.a, vat.c(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('c(a)'); legend({'Inf BASE','For BASE','Inf VAT↑','For VAT↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'Bendog_policy_c'));

fig=figure('Name','s(a) por tipo: BASE vs VAT↑');
plot(base.a,base.s(:,1),'LineWidth',2); hold on; plot(base.a,base.s(:,2),'LineWidth',2);
plot(vat.a, vat.s(:,1),'--','LineWidth',2); plot(vat.a, vat.s(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('s(a)'); legend({'Inf BASE','For BASE','Inf VAT↑','For VAT↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'Bendog_policy_s'));

% Densidades g(a)
fig=figure('Name','g(a) - Informal'); plot(base.a,base.g(:,1),'LineWidth',2); hold on; plot(vat.a,vat.g(:,1),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_1(a)'); legend({'BASE','VAT↑'}); export_fig(fig, fullfile(outdir_figs,'Bendog_g_informal'));
fig=figure('Name','g(a) - Formal');   plot(base.a,base.g(:,2),'LineWidth',2); hold on; plot(vat.a,vat.g(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_2(a)'); legend({'BASE','VAT↑'}); export_fig(fig, fullfile(outdir_figs,'Bendog_g_formal'));

% Cierre financiero (barras)
cats = categorical({'A_{priv}','B (public)'}); cats=reordercats(cats,{'A_{priv}','B (public)'});
Yb = [ (A_I_base + A_F_base) base.fiscal.B ; (A_I_vat + A_F_vat) vat.fiscal.B ];
fig=figure('Name','Cierre financiero: BASE vs VAT↑');
bar(cats, Yb.','grouped'); grid on; ylabel('level'); legend({'BASE','VAT↑'},'Location','best');
title(sprintf('r* BASE=%.4f | r* VAT↑=%.4f', base.r, vat.r));
export_fig(fig, fullfile(outdir_figs,'Bendog_asset_closure'));

% Curvas financieras A_I(r), A_F(r), A_priv(r), B(r)
[rB,AIb,AFb,APb,Bb,Sb] = financial_curves(cfg,   base.r);
[rV,AIv,AFv,APv,Bv,Sv] = financial_curves(cfgVAT, vat.r);
fig=figure('Name','Demanda privada vs oferta pública');
plot(rB,AIb,'LineWidth',2); hold on; plot(rB,AFb,'LineWidth',2); plot(rB,APb,'LineWidth',2); plot(rB,Bb,'k--','LineWidth',1.6);
plot(rV,AIv,'--','LineWidth',2);    plot(rV,AFv,'--','LineWidth',2);    plot(rV,APv,'--','LineWidth',2);    plot(rV,Bv,'k:','LineWidth',1.6);
xline(base.r,'b:'); xline(vat.r,'r:'); grid on; xlabel('r'); ylabel('level');
legend({'A_I BASE','A_F BASE','A_{priv} BASE','B BASE','A_I VAT↑','A_F VAT↑','A_{priv} VAT↑','B VAT↑','r* BASE','r* VAT↑'},'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'Bendog_financial_curves'));

fig=figure('Name','Exceso S(r)');
plot(rB,Sb,'LineWidth',2); hold on; plot(rV,Sv,'--','LineWidth',2);
yline(0,'k--'); xline(base.r,'b:'); xline(vat.r,'r:'); grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
legend({'BASE','VAT↑','S=0','r* BASE','r* VAT↑'},'Location','best'); export_fig(fig, fullfile(outdir_figs,'Bendog_excess_supply'));

% Fiscal agrupado
catsR = categorical({'Labor','VAT','Total'}); catsR=reordercats(catsR,{'Labor','VAT','Total'});
YR = [ base.fiscal.Tl base.fiscal.Tc base.fiscal.Tl+base.fiscal.Tc; ...
       vat.fiscal.Tl  vat.fiscal.Tc  vat.fiscal.Tl + vat.fiscal.Tc ];
fig=figure('Name','Revenues'); bar(catsR,YR,'grouped'); grid on; ylabel('amount'); legend({'BASE','VAT↑'}); export_fig(fig, fullfile(outdir_figs,'Bendog_fiscal_rev'));

catsE = categorical({'Debt serv','G','Transfers'}); catsE=reordercats(catsE,{'Debt serv','G','Transfers'});
YE = [ base.fiscal.rB base.fiscal.G base.fiscal.Tr; vat.fiscal.rB  vat.fiscal.G  vat.fiscal.Tr ];
fig=figure('Name','Expenditures'); bar(catsE,YE,'grouped'); grid on; ylabel('amount'); legend({'BASE','VAT↑'}); export_fig(fig, fullfile(outdir_figs,'Bendog_fiscal_exp'));

catsB = categorical({'Primary','Overall'}); catsB=reordercats(catsB,{'Primary','Overall'});
YB = [ base.fiscal.PB base.fiscal.BB; vat.fiscal.PB vat.fiscal.BB ];
fig=figure('Name','Balances'); bar(catsB,YB,'grouped'); grid on; ylabel('amount'); legend({'BASE','VAT↑'}); yline(0,'k--'); export_fig(fig, fullfile(outdir_figs,'Bendog_fiscal_bal'));

% Lorenz total y por tipo (riqueza y consumo) — igual que antes
% Total riqueza
fig=figure('Name','Lorenz wealth (total)'); [~,LwB,cPB]=lorenz_from_density(base.a,base.g); [~,LwV,cPV]=lorenz_from_density(vat.a,vat.g);
plot(cPB,LwB,'LineWidth',2); hold on; plot(cPV,LwV,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Gini_W^T BASE=%.3f | VAT=%.3f', base.stats.giniW(3), vat.stats.giniW(3)));
legend({'BASE','VAT↑','45°'},'Location','southeast'); export_fig(fig, fullfile(outdir_figs,'Bendog_lorenzW_total'));

% Por tipo riqueza
fig=figure('Name','Lorenz wealth by type');
[LwI_B,cI_B]=lorenz_from_assets_single(base.a, base.g(:,1)); [LwF_B,cF_B]=lorenz_from_assets_single(base.a, base.g(:,2));
[LwI_V,cI_V]=lorenz_from_assets_single(vat.a,  vat.g(:,1));  [LwF_V,cF_V]=lorenz_from_assets_single(vat.a,  vat.g(:,2));
subplot(1,2,1); plot(cI_B,LwI_B,'LineWidth',2); hold on; plot(cI_V,LwI_V,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share'); title(sprintf('Informal: Gini %.3f | %.3f', base.stats.giniW(1), vat.stats.giniW(1)));
subplot(1,2,2); plot(cF_B,LwF_B,'LineWidth',2); hold on; plot(cF_V,LwF_V,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share'); title(sprintf('Formal: Gini %.3f | %.3f', base.stats.giniW(2), vat.stats.giniW(2)));
export_fig(fig, fullfile(outdir_figs,'Bendog_lorenzW_bytype'));

% Consumo
fig=figure('Name','Lorenz consumption (total)'); [~,LcB,cCB]=lorenz_from_values(base.a,base.g,base.c(:,1)+base.c(:,2));
[~,LcV,cCV]=lorenz_from_values(vat.a,vat.g,vat.c(:,1)+vat.c(:,2));
plot(cCB,LcB,'LineWidth',2); hold on; plot(cCV,LcV,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share'); title(sprintf('Gini_C^T %.3f | %.3f', base.stats.giniC(3), vat.stats.giniC(3)));
legend({'BASE','VAT↑','45°'},'Location','southeast'); export_fig(fig, fullfile(outdir_figs,'Bendog_lorenzC_total'));

fig=figure('Name','Lorenz consumption by type');
[LcI_B, cIB] = lorenz_from_values_single(base.c(:,1), base.g(:,1), base.a);
[LcF_B, cFB] = lorenz_from_values_single(base.c(:,2), base.g(:,2), base.a);
[LcI_V, cIV] = lorenz_from_values_single(vat.c(:,1),  vat.g(:,1),  vat.a);
[LcF_V, cFV] = lorenz_from_values_single(vat.c(:,2),  vat.g(:,2),  vat.a);
subplot(1,2,1); plot(cIB,LcI_B,'LineWidth',2); hold on; plot(cIV,LcI_V,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share'); title(sprintf('Informal: Gini %.3f | %.3f', base.stats.giniC(1), vat.stats.giniC(1)));
subplot(1,2,2); plot(cFB,LcF_B,'LineWidth',2); hold on; plot(cFV,LcF_V,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share'); title(sprintf('Formal: Gini %.3f | %.3f', base.stats.giniC(2), vat.stats.giniC(2)));
export_fig(fig, fullfile(outdir_figs,'Bendog_lorenzC_bytype'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% --------------------- Helpers gráficos/tablas ---------------------
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
end
function print_summary(sol)
    fb=sol.fiscal;
    fprintf('r*=%.5f | Y=%.6f | C=%.6f | A_priv=%.6f | B=%.6f | S=%.2e\n', sol.r, sol.Y, sol.Ctot, ...
        sum(sol.g(:,1).*sol.a)*(sol.a(2)-sol.a(1))+sum(sol.g(:,2).*sol.a)*(sol.a(2)-sol.a(1)), fb.B, sol.S_residual);
    fprintf('Tl=%.6f Tc=%.6f Tr=%.6f G=%.6f | PB=%.6f rB=%.6f BB=%.1e\n', fb.Tl,fb.Tc,fb.Tr,fb.G,fb.PB,fb.rB,fb.BB);
end
function [r_span,Ai,Af,Apriv,Bg,Sg] = financial_curves(cfg_in, r_star)
    r_lo = max(cfg_in.rmin, max(0.5*r_star, 0.003));
    r_hi = min(cfg_in.rmax, max(1.8*r_star, r_star + 0.02));
    if r_hi <= r_lo, r_hi = r_lo*(1+0.05); end
    r_span = linspace(r_lo,r_hi,35);
    Ai=nan(size(r_span)); Af=Ai; Apriv=Ai; Bg=Ai; Sg=Ai;
    for k=1:numel(r_span)
        cfr = cfg_in; cfr.fix_r=1; cfr.r_guess=r_span(k); cfr.rmin=r_span(k); cfr.rmax=r_span(k);
        solr = solve_two_type_huggett_fiscal_Bfixed(cfr);
        a=solr.a; g=solr.g; da=a(2)-a(1);
        Ai(k)=sum(g(:,1).*a)*da; Af(k)=sum(g(:,2).*a)*da; Apriv(k)=Ai(k)+Af(k);
        Bg(k)=solr.fiscal.B; Sg(k)=Apriv(k)-Bg(k);
    end
end
function [gT,L,cumPop]=lorenz_from_density(a,g)
    da=a(2)-a(1); gT=g(:,1)+g(:,2); W=sum(gT)*da; [as,ix]=sort(a); gs=gT(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as-min(0,min(as))+1e-12; cumWealth=cumsum(wealth_pos.*gs); L=cumWealth/max(cumWealth(end),1e-12);
end
function [vals,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w); vals=x(:); [xs,ix]=sort(vals); ws=w(ix);
    cumPop=cumsum(ws)/W; cumx=cumsum(xs.*ws); L=cumx/max(cumx(end),1e-12);
end
function [L,cumPop]=lorenz_from_assets_single(a, gk)
    da=a(2)-a(1); W=sum(gk)*da; [as,ix]=sort(a); gs=gk(ix)*da; cumPop=cumsum(gs)/max(W,1e-12);
    wealth_pos = as - min(0,min(as)) + 1e-12;  cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
function [L,cumPop]=lorenz_from_values_single(xk, gk, a)
    da=a(2)-a(1); w = gk(:)*da; W=sum(w); [xs,ix]=sort(xk(:)); ws=w(ix);
    cumPop = cumsum(ws)/max(W,1e-12);  cumx = cumsum(xs.*ws);  L = cumx / max(cumx(end),1e-12);
end
