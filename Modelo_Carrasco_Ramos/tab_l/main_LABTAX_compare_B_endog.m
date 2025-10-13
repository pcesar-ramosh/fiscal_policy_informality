% ================= main_LABTAX_compare_B_endog.m =================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------- Parámetros base limpios ----------
cfg = struct();
cfg.RRA_I = 3.40; cfg.RRA_F = 3.40;
cfg.rho   = 0.08;
cfg.tau_l = 0.15;                 % impuesto laboral (solo formales)
cfg.tau_c = 0.18;                 % IVA
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

% ===== Paso 0: referenciar G/PIB con el cierre antiguo para fijar gY =====
cfg_ref = cfg; cfg_ref.fiscal_mode = 'G_resid_B_fixed';
cfg_ref.B_mode = 'level'; cfg_ref.Bbar = 0.25;       % usado SOLO aquí
ref = solve_two_type_huggett_fiscal_Bfixed(cfg_ref);
gY  = ref.fiscal.G / max(ref.Y,1e-12);               % ratio objetivo de G

% ===== Cierre de trabajo: G fijo (gY) y B endógeno =====
cfg.fiscal_mode    = 'G_fixed_B_endog';
cfg.G_target_ratio = gY;

paper_style();

% --------------------- Escenario BASE ---------------------
base = solve_two_type_huggett_fiscal_Bfixed(cfg);

% --------------------- Escenario LABTAX-UP ----------------
d_tau_l = 0.05;                 % +5 pp al impuesto laboral
cfgLAB = cfg; cfgLAB.tau_l = cfg.tau_l + d_tau_l;
lab = solve_two_type_huggett_fiscal_Bfixed(cfgLAB);

% ------------------- Reporte en consola -------------------
disp('== BASE (G fixed, B endog) ==');   print_summary(base);
disp('== LABTAX↑ (G fixed, B endog) ==');print_summary(lab);

% --------------------------- TABLAS ------------------------
% Fiscal
T_fiscal = table( ...
    ["BASE"; "LABTAX_UP"], ...
    [base.fiscal.Tl; lab.fiscal.Tl], [base.fiscal.Tc; lab.fiscal.Tc], ...
    [base.fiscal.Tl+base.fiscal.Tc; lab.fiscal.Tl+lab.fiscal.Tc], ...
    [base.fiscal.Tr; lab.fiscal.Tr], [base.fiscal.G; lab.fiscal.G], ...
    [base.fiscal.rB; lab.fiscal.rB], [base.fiscal.PB; lab.fiscal.PB], ...
    [base.fiscal.B; lab.fiscal.B],   [base.fiscal.BB; lab.fiscal.BB], ...
'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers','public_good','debt_serv','primary_bal','debt_stock','global_bal'});
writetable(T_fiscal, fullfile(outdir_tabs,'labtax_Bendog_fiscal.csv'));

% Mercado de activos (equilibrio)
daB = base.a(2)-base.a(1);
A_I_base = sum(base.g(:,1).*base.a)*daB;  A_F_base = sum(base.g(:,2).*base.a)*daB;
daL = lab.a(2)-lab.a(1);
A_I_lab  = sum(lab.g(:,1).*lab.a)*daL;    A_F_lab  = sum(lab.g(:,2).*lab.a)*daL;

T_assets = table( ...
    ["BASE"; "LABTAX_UP"], ...
    [A_I_base; A_I_lab], [A_F_base; A_F_lab], [A_I_base + A_F_base; A_I_lab + A_F_lab], ...
    [base.fiscal.B; lab.fiscal.B], ...
'VariableNames', {'scenario','A_I','A_F','A_private','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'labtax_Bendog_asset_market.csv'));

% Estadísticas
Sb = base.stats; Sl = lab.stats;
T_stats = table( ...
    ["BASE"; "LABTAX_UP"], ...
    [base.r; lab.r], [base.popI; lab.popI], [base.popF; lab.popF], [base.Y; lab.Y], [base.Ctot; lab.Ctot], ...
    [Sb.wealth_mean(1); Sl.wealth_mean(1)], [Sb.wealth_mean(2); Sl.wealth_mean(2)], [Sb.wealth_mean(3); Sl.wealth_mean(3)], ...
    [Sb.wealth_median(1); Sl.wealth_median(1)], [Sb.wealth_median(2); Sl.wealth_median(2)], [Sb.wealth_median(3); Sl.wealth_median(3)], ...
    [Sb.giniW(1); Sl.giniW(1)], [Sb.giniW(2); Sl.giniW(2)], [Sb.giniW(3); Sl.giniW(3)], ...
    [Sb.cons_mean(1); Sl.cons_mean(1)], [Sb.cons_mean(2); Sl.cons_mean(2)], [Sb.cons_mean(3); Sl.cons_mean(3)], ...
    [Sb.cons_median(1); Sl.cons_median(1)], [Sb.cons_median(2); Sl.cons_median(2)], [Sb.cons_median(3); Sl.cons_median(3)], ...
    [Sb.giniC(1); Sl.giniC(1)], [Sb.giniC(2); Sl.giniC(2)], [Sb.giniC(3); Sl.giniC(3)], ...
'VariableNames', {'scenario','r','popI','popF','Y','Ctot', ...
    'wealth_mean_I','wealth_mean_F','wealth_mean_T', ...
    'wealth_med_I','wealth_med_F','wealth_med_T', ...
    'giniW_I','giniW_F','giniW_T', ...
    'cons_mean_I','cons_mean_F','cons_mean_T', ...
    'cons_med_I','cons_med_F','cons_med_T', ...
    'giniC_I','giniC_F','giniC_T'});
writetable(T_stats, fullfile(outdir_tabs,'labtax_Bendog_stats.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

% --------------------------- GRÁFICOS ----------------------
paper_style();

% c(a) y s(a)
fig=figure('Name','c(a): BASE vs LABTAX↑');
plot(base.a,base.c(:,1),'LineWidth',2); hold on; plot(base.a,base.c(:,2),'LineWidth',2);
plot(lab.a, lab.c(:,1),'--','LineWidth',2); plot(lab.a, lab.c(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('c(a)'); legend({'Inf BASE','For BASE','Inf LAB↑','For LAB↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'LABTAX_policy_c'));

fig=figure('Name','s(a): BASE vs LABTAX↑');
plot(base.a,base.s(:,1),'LineWidth',2); hold on; plot(base.a,base.s(:,2),'LineWidth',2);
plot(lab.a, lab.s(:,1),'--','LineWidth',2); plot(lab.a, lab.s(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('s(a)'); legend({'Inf BASE','For BASE','Inf LAB↑','For LAB↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'LABTAX_policy_s'));

% g(a) por tipo
fig=figure('Name','g(a) - Informal'); plot(base.a,base.g(:,1),'LineWidth',2); hold on; plot(lab.a,lab.g(:,1),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_1(a)'); legend({'BASE','LABTAX↑'}); export_fig(fig, fullfile(outdir_figs,'LABTAX_g_informal'));
fig=figure('Name','g(a) - Formal');   plot(base.a,base.g(:,2),'LineWidth',2); hold on; plot(lab.a,lab.g(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_2(a)'); legend({'BASE','LABTAX↑'}); export_fig(fig, fullfile(outdir_figs,'LABTAX_g_formal'));

% Cierre financiero (barras)
cats = categorical({'A_{priv}','B (public)'}); cats=reordercats(cats,{'A_{priv}','B (public)'});
Yb = [ (A_I_base + A_F_base) base.fiscal.B ; (A_I_lab + A_F_lab) lab.fiscal.B ];
fig=figure('Name','Cierre financiero: BASE vs LABTAX↑');
bar(cats, Yb.','grouped'); grid on; ylabel('level'); legend({'BASE','LABTAX↑'},'Location','best');
title(sprintf('r* BASE=%.4f | r* LAB↑=%.4f', base.r, lab.r));
export_fig(fig, fullfile(outdir_figs,'LABTAX_asset_closure'));

% Curvas financieras A_I(r), A_F(r), A_priv(r), B(r)
[rB,AIb,AFb,APb,Bb,Sb] = financial_curves(cfg,    base.r);
[rL,AIl,AFl,APl,Bl,Sl] = financial_curves(cfgLAB,  lab.r);
fig=figure('Name','Demanda privada vs oferta pública');
plot(rB,AIb,'LineWidth',2); hold on; plot(rB,AFb,'LineWidth',2); plot(rB,APb,'LineWidth',2); plot(rB,Bb,'k--','LineWidth',1.6);
plot(rL,AIl,'--','LineWidth',2);    plot(rL,AFl,'--','LineWidth',2);    plot(rL,APl,'--','LineWidth',2);    plot(rL,Bl,'k:','LineWidth',1.6);
xline(base.r,'b:'); xline(lab.r,'r:'); grid on; xlabel('r'); ylabel('level');
legend({'A_I BASE','A_F BASE','A_{priv} BASE','B BASE','A_I LAB↑','A_F LAB↑','A_{priv} LAB↑','B LAB↑','r* BASE','r* LAB↑'},'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'LABTAX_financial_curves'));

fig=figure('Name','Exceso S(r)');
plot(rB,Sb,'LineWidth',2); hold on; plot(rL,Sl,'--','LineWidth',2);
yline(0,'k--'); xline(base.r,'b:'); xline(lab.r,'r:'); grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
legend({'BASE','LABTAX↑','S=0','r* BASE','r* LAB↑'},'Location','best'); export_fig(fig, fullfile(outdir_figs,'LABTAX_excess_supply'));

% Fiscal agrupado
catsR = categorical({'Labor','VAT','Total'}); catsR=reordercats(catsR,{'Labor','VAT','Total'});
YR = [ base.fiscal.Tl base.fiscal.Tc base.fiscal.Tl+base.fiscal.Tc; ...
       lab.fiscal.Tl  lab.fiscal.Tc  lab.fiscal.Tl + lab.fiscal.Tc ];
fig=figure('Name','Revenues'); bar(catsR,YR,'grouped'); grid on; ylabel('amount'); legend({'BASE','LABTAX↑'}); export_fig(fig, fullfile(outdir_figs,'LABTAX_fiscal_rev'));

catsE = categorical({'Debt serv','G','Transfers'}); catsE=reordercats(catsE,{'Debt serv','G','Transfers'});
YE = [ base.fiscal.rB base.fiscal.G base.fiscal.Tr; lab.fiscal.rB  lab.fiscal.G  lab.fiscal.Tr ];
fig=figure('Name','Expenditures'); bar(catsE,YE,'grouped'); grid on; ylabel('amount'); legend({'BASE','LABTAX↑'}); export_fig(fig, fullfile(outdir_figs,'LABTAX_fiscal_exp'));

catsB = categorical({'Primary','Overall'}); catsB=reordercats(catsB,{'Primary','Overall'});
YB = [ base.fiscal.PB base.fiscal.BB; lab.fiscal.PB lab.fiscal.BB ];
fig=figure('Name','Balances'); bar(catsB,YB,'grouped'); grid on; ylabel('amount'); legend({'BASE','LABTAX↑'}); yline(0,'k--'); export_fig(fig, fullfile(outdir_figs,'LABTAX_fiscal_bal'));

% Lorenz (riqueza total)
fig=figure('Name','Lorenz wealth (total)');
[~,LwB,cPB]=lorenz_from_density(base.a,base.g); [~,LwL,cPL]=lorenz_from_density(lab.a,lab.g);
plot(cPB,LwB,'LineWidth',2); hold on; plot(cPL,LwL,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Gini_W^T BASE=%.3f | LAB↑=%.3f', base.stats.giniW(3), lab.stats.giniW(3)));
legend({'BASE','LABTAX↑','45°'},'Location','southeast'); export_fig(fig, fullfile(outdir_figs,'LABTAX_lorenzW_total'));

% Lorenz por tipo (riqueza)
fig=figure('Name','Lorenz wealth by type');
[LwI_B,cI_B]=lorenz_from_assets_single(base.a, base.g(:,1)); [LwF_B,cF_B]=lorenz_from_assets_single(base.a, base.g(:,2));
[LwI_L,cI_L]=lorenz_from_assets_single(lab.a,  lab.g(:,1));  [LwF_L,cF_L]=lorenz_from_assets_single(lab.a,  lab.g(:,2));
subplot(1,2,1); plot(cI_B,LwI_B,'LineWidth',2); hold on; plot(cI_L,LwI_L,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share'); title(sprintf('Informal: Gini %.3f | %.3f', base.stats.giniW(1), lab.stats.giniW(1)));
subplot(1,2,2); plot(cF_B,LwF_B,'LineWidth',2); hold on; plot(cF_L,LwF_L,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Wealth share'); title(sprintf('Formal: Gini %.3f | %.3f', base.stats.giniW(2), lab.stats.giniW(2)));
export_fig(fig, fullfile(outdir_figs,'LABTAX_lorenzW_bytype'));

% Lorenz (consumo total y por tipo)
fig=figure('Name','Lorenz consumption (total)');
[~,LcB,cCB]=lorenz_from_values(base.a,base.g,base.c(:,1)+base.c(:,2));
[~,LcL,cCL]=lorenz_from_values(lab.a,lab.g,lab.c(:,1)+lab.c(:,2));
plot(cCB,LcB,'LineWidth',2); hold on; plot(cCL,LcL,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share'); title(sprintf('Gini_C^T %.3f | %.3f', base.stats.giniC(3), lab.stats.giniC(3)));
legend({'BASE','LABTAX↑','45°'},'Location','southeast'); export_fig(fig, fullfile(outdir_figs,'LABTAX_lorenzC_total'));

fig=figure('Name','Lorenz consumption by type');
[LcI_B, cIB] = lorenz_from_values_single(base.c(:,1), base.g(:,1), base.a);
[LcF_B, cFB] = lorenz_from_values_single(base.c(:,2), base.g(:,2), base.a);
[LcI_L, cIL] = lorenz_from_values_single(lab.c(:,1),  lab.g(:,1),  lab.a);
[LcF_L, cFL] = lorenz_from_values_single(lab.c(:,2),  lab.g(:,2),  lab.a);
subplot(1,2,1); plot(cIB,LcI_B,'LineWidth',2); hold on; plot(cIL,LcI_L,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share'); title(sprintf('Informal: Gini %.3f | %.3f', base.stats.giniC(1), lab.stats.giniC(1)));
subplot(1,2,2); plot(cFB,LcF_B,'LineWidth',2); hold on; plot(cFL,LcF_L,'--','LineWidth',2); plot([0,1],[0,1],'k--'); axis square; grid on;
xlabel('Population share'); ylabel('Consumption share'); title(sprintf('Formal: Gini %.3f | %.3f', base.stats.giniC(2), lab.stats.giniC(2)));
export_fig(fig, fullfile(outdir_figs,'LABTAX_lorenzC_bytype'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ------------------------ Helpers ------------------------
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
    fb=sol.fiscal; da=sol.a(2)-sol.a(1);
    Apriv = sum(sol.g(:,1).*sol.a)*da + sum(sol.g(:,2).*sol.a)*da;
    fprintf('r*=%.5f | Y=%.6f | C=%.6f | A_priv=%.6f | B=%.6f | S=%.2e\n', sol.r, sol.Y, sol.Ctot, Apriv, fb.B, sol.S_residual);
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
