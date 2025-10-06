% ===== main_compare_income_shock_FAST.m =====
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ----------------------- FAST switches -----------------------------------
FAST              = true;     % si false: modo completo
I_grid            = FAST * 500 + (~FAST) * 700;   % puntos en malla de activos
MAXITV_BASE       = FAST * 110 + (~FAST) * 140;   % HJB base (escenarios)
MAXITV_DIAG       = FAST * 50  + (~FAST) * 100;   % HJB en diagnósticos r-span
N_DIAG            = FAST * 11  + (~FAST) * 41;    % puntos en r-span
EXPORT_PDF        = false;     % true si quieres también PDF (más lento)
DIAG_WINDOW_MULT  = 0.9;       % ventana alrededor de r* (menor = más rápido)

% ----------------------- Baseline config (CURRENT policy) ----------------
cfg0 = struct();
cfg0.RRA_I = 3.40; cfg0.RRA_F = 3.40;
cfg0.rho   = 0.08;
cfg0.tau_l = 0.15; cfg0.tau_c = 0.18; cfg0.phi = 0.09;
cfg0.z1 = 0.33; cfg0.z2 = 1.00;
cfg0.theta_I = 0.06;  cfg0.theta_F = 0.01;

cfg0.eta_target = 0.654; cfg0.p22_bar = 0.8155;
cfg0.I = I_grid; cfg0.amax = 3.0; cfg0.amin = -2.0*cfg0.z1;

cfg0.r_guess = 0.03; cfg0.rmin = 0.005; cfg0.rmax = 0.10;
cfg0.maxit_V = MAXITV_BASE; cfg0.crit_V = 1e-6; cfg0.Delta = 1400;
cfg0.maxit_r = 1500; cfg0.crit_S = 1e-5; cfg0.fix_r = 0;

% Bien público (suave, últimos refinamientos)
cfg0.psi_G = 0.08; cfg0.omegaG = 0.50; cfg0.report_G_effects = 0;
cfg0.sigma_a = 0.010;

% Gobierno / Deuda
cfg0.B_mode  = 'ratio_to_Y'; cfg0.Bbar = 0.35;
cfg0.alphaG  = 0.50; cfg0.clamp_G_to_zero = true; cfg0.G_cap_ratio = 0.08;

% -------------------- Income drop scenarios ------------------------------
% (b) Caída de ingresos: -35% informal (z1), -25% formal (z2)
cfgINC = cfg0;
cfgINC.z1 = 0.65 * cfg0.z1;
cfgINC.z2 = 0.75 * cfg0.z2;

% (c) Caída de ingresos + incremento en transferencias (+15%)
cfgINC_TR = cfgINC;
cfgINC_TR.phi = 1.15 * cfg0.phi;

% --------------------------- Escenarios ----------------------------------
sc = struct('name',{},'cfg',{},'sol',{});
sc(1).name = 'BASE';        sc(1).cfg = cfg0;
sc(2).name = 'INCOME_DROP'; sc(2).cfg = cfgINC;
sc(3).name = 'INCOME_DROP_TR'; sc(3).cfg = cfgINC_TR;

% Nombres bonitos
pretty = @(nm) strrep(strrep(strrep(nm,'BASE','Current policy'), ...
                              'INCOME_DROP_TR','Income drop + transfers +15%'), ...
                              'INCOME_DROP','Income drop');

% ---------------------------- Resolver -----------------------------------
for k=1:numel(sc)
    tic;
    sc(k).sol = solve_two_type_huggett_fiscal_Bfixed(sc(k).cfg);
    eta_real = sc(k).sol.popI/(sc(k).sol.popI+sc(k).sol.popF);
    fprintf('[%s] r*=%.4f | eta_real=%.3f | B=%.3f | G=%.3f | Y=%.3f | C=%.3f  (%.2fs)\n', ...
        sc(k).name, sc(k).sol.r, eta_real, sc(k).sol.fiscal.B, sc(k).sol.fiscal.G, ...
        sc(k).sol.Y, sc(k).sol.Ctot, toc);

    % ---- CSVs por escenario ----
    a  = sc(k).sol.a; g = sc(k).sol.g; da = a(2)-a(1);
    A_I = sum(g(:,1).*a)*da; A_F = sum(g(:,2).*a)*da; A_priv = A_I + A_F;

    scen = string(pretty(sc(k).name));
    T_fiscal = table(scen, sc(k).sol.fiscal.Tl, sc(k).sol.fiscal.Tc, ...
        sc(k).sol.fiscal.Tl+sc(k).sol.fiscal.Tc, sc(k).sol.fiscal.Tr, ...
        sc(k).sol.fiscal.G, sc(k).sol.fiscal.rB, sc(k).sol.fiscal.PB, ...
        sc(k).sol.fiscal.B, sc(k).sol.fiscal.BB, ...
        'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers', ...
                          'public_good','debt_serv','primary_bal','debt_stock','global_bal'});
    writetable(T_fiscal, fullfile(outdir_tabs, sprintf('fiscal_%s.csv', sc(k).name)));

    T_assets = table(scen, A_priv, sc(k).sol.fiscal.B, ...
        'VariableNames', {'scenario','A_private','B_public'});
    writetable(T_assets, fullfile(outdir_tabs, sprintf('assets_%s.csv', sc(k).name)));

    S = sc(k).sol.stats;
    T_stats = table(scen, sc(k).sol.r, sc(k).sol.popI, sc(k).sol.popF, sc(k).sol.Y, sc(k).sol.Ctot, ...
        S.wealth_mean(1),S.wealth_mean(2),S.wealth_mean(3), ...
        S.wealth_median(1),S.wealth_median(2),S.wealth_median(3), ...
        S.giniW(1),S.giniW(2),S.giniW(3), ...
        S.cons_mean(1),S.cons_mean(2),S.cons_mean(3), ...
        S.cons_median(1),S.cons_median(2),S.cons_median(3), ...
        S.giniC(1),S.giniC(2),S.giniC(3), S.p11, ...
        'VariableNames', {'scenario','r','popI','popF','Y','Ctot', ...
        'wealth_mean_I','wealth_mean_F','wealth_mean_T', ...
        'wealth_med_I','wealth_med_F','wealth_med_T', ...
        'giniW_I','giniW_F','giniW_T', ...
        'cons_mean_I','cons_mean_F','cons_mean_T', ...
        'cons_med_I','cons_med_F','cons_med_T', ...
        'giniC_I','giniC_F','giniC_T','p11_rep'});
    writetable(T_stats, fullfile(outdir_tabs, sprintf('hh_stats_%s.csv', sc(k).name)));
end

% -------- Tabla comparativa resumen -------------------------------------
RES = cell(numel(sc), 15);
for k=1:numel(sc)
    s = sc(k).sol; da = s.a(2)-s.a(1);
    A_priv = sum((s.g(:,1)+s.g(:,2)).*s.a)*da;
    RES(k,:) = { pretty(sc(k).name), sc(k).cfg.eta_target, s.r, s.fiscal.B, A_priv, ...
                 s.popI/(s.popI+s.popF), s.Y, s.Ctot, s.fiscal.G, ...
                 s.fiscal.PB, s.fiscal.BB, s.stats.giniW(1), s.stats.giniW(2), ...
                 s.stats.giniW(3), s.stats.giniC(3) };
end
Tsum = cell2table(RES, 'VariableNames', ...
 {'scenario','eta_target','r','B','A_priv','eta_real','Y','C','G','PB','BB','giniW_I','giniW_F','giniW_T','giniC_T'});
writetable(Tsum, fullfile(outdir_tabs,'compare_income_shock_summary.csv'));
disp(Tsum);

% --------------------------- Figuras (mismo set) -------------------------
paper_style();
cols = lines(numel(sc));
names_char = arrayfun(@(j) char(pretty(sc(j).name)), 1:numel(sc), 'uni', 0);
mk = @(suffix) cellfun(@(s)[s suffix], names_char, 'uni', 0);

% Policies: c(a)
fig=figure('Name','Policies: consumption by scenario'); hold on;
hI=gobjects(numel(sc),1); hF=hI;
for k=1:numel(sc)
    a=sc(k).sol.a; c=sc(k).sol.c;
    hI(k)=plot(a,c(:,1),'Color',cols(k,:),'LineWidth',2,'LineStyle','-');
    hF(k)=plot(a,c(:,2),'Color',cols(k,:),'LineWidth',2,'LineStyle','--');
end
grid on; xlabel('Assets a'); ylabel('Consumption c(a)');
legend([hI;hF],[mk(' - I'), mk(' - F')],'Location','bestoutside');
title('Consumption policies (I solid, F dashed)');
export_fig(fig, fullfile(outdir_figs,'inc_policy_consumption'), EXPORT_PDF);

% Policies: s(a)
fig=figure('Name','Policies: savings by scenario'); hold on;
hI=gobjects(numel(sc),1); hF=hI;
for k=1:numel(sc)
    a=sc(k).sol.a; s=sc(k).sol.s;
    hI(k)=plot(a,s(:,1),'Color',cols(k,:),'LineWidth',2,'LineStyle','-');
    hF(k)=plot(a,s(:,2),'Color',cols(k,:),'LineWidth',2,'LineStyle','--');
end
grid on; xlabel('Assets a'); ylabel('Savings s(a)');
legend([hI;hF],[mk(' - I'), mk(' - F')],'Location','bestoutside');
title('Savings policies (I solid, F dashed)');
export_fig(fig, fullfile(outdir_figs,'inc_policy_savings'), EXPORT_PDF);

% Wealth densities overlay
fig=figure('Name','Wealth densities by type (overlay)');
subplot(2,1,1); hold on;
for k=1:numel(sc), plot(sc(k).sol.a, sc(k).sol.g(:,1), 'LineWidth',2,'Color',cols(k,:)); end
grid on; xlabel('Assets a'); ylabel('Density g_I(a)'); title('Informal density'); legend(names_char,'Location','best');
subplot(2,1,2); hold on;
for k=1:numel(sc), plot(sc(k).sol.a, sc(k).sol.g(:,2), 'LineWidth',2,'Color',cols(k,:)); end
grid on; xlabel('Assets a'); ylabel('Density g_F(a)'); title('Formal density'); legend(names_char,'Location','best');
export_fig(fig, fullfile(outdir_figs,'inc_wealth_distributions'), EXPORT_PDF);

% Asset market closure (stacked)
fig=figure('Name','Asset market closure by scenario (stacked)');
tiledlayout(1,numel(sc),'Padding','compact','TileSpacing','compact');
for k=1:numel(sc)
    nexttile;
    a=sc(k).sol.a; g=sc(k).sol.g; da=a(2)-a(1);
    A_I=sum(g(:,1).*a)*da; A_F=sum(g(:,2).*a)*da; Bval=sc(k).sol.fiscal.B;
    Ybars = [A_I A_F 0; 0 0 Bval];
    bar(categorical({'Private demand','Public debt B'}), Ybars, 'stacked');
    grid on; ylabel('Assets level'); title(sprintf('%s (r^*=%.4f)', pretty(sc(k).name), sc(k).sol.r));
    legend({'A_I','A_F','B'},'Location','northeast');
    xt=[1 2]; yt=[A_I+A_F, Bval];
    for j=1:2, text(xt(j), yt(j), sprintf(' %.3f',yt(j)), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10); end
end
export_fig(fig, fullfile(outdir_figs,'inc_asset_market_closure'), EXPORT_PDF);

% A_priv(r) and B(r) — FAST alrededor de r*
fig=figure('Name','A_{priv}(r) & B(r) by scenario'); hold on;
hAp=gobjects(numel(sc),1); hB=hAp; hr=hAp;
for k=1:numel(sc)
    rstar = sc(k).sol.r;
    rmin = max(0.005, rstar*(1-DIAG_WINDOW_MULT));
    rmax = min(0.10,  rstar*(1+DIAG_WINDOW_MULT));
    r_span = linspace(rmin,rmax,N_DIAG);

    alt=sc(k).cfg; alt.fix_r=1; alt.maxit_r=1; alt.maxit_V = MAXITV_DIAG;
    Apriv=nan(size(r_span)); Bgrid=Apriv;

    for i=1:numel(r_span)
        alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal_Bfixed(alt);
        da_i=tmp.a(2)-tmp.a(1);
        Apriv(i)=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
        Bgrid(i)=tmp.fiscal.B;
    end
    hAp(k)=plot(r_span,Apriv,'LineWidth',2,'Color',cols(k,:));
    hB(k)=plot(r_span,Bgrid,'--','LineWidth',1.6,'Color',cols(k,:));
    hr(k)=xline(sc(k).sol.r,':','Color',cols(k,:),'LineWidth',1.6);
end
grid on; xlabel('Interest rate r'); ylabel('Assets / Bonds level');
labels = [ cellfun(@(s)[s ' - A_{priv}(r)'], names_char,'uni',0), ...
           cellfun(@(s)[s ' - B(r)'],        names_char,'uni',0), ...
           cellfun(@(s)[s ' - r^*'],         names_char,'uni',0) ];
legend([hAp;hB;hr], labels, 'Location','bestoutside');
title('Private demand A_{priv}(r) and public supply B(r)');
export_fig(fig, fullfile(outdir_figs,'inc_asset_demand_vs_B'), EXPORT_PDF);

% Excess S(r) — FAST
fig=figure('Name','Excess S(r) by scenario'); hold on;
hS=gobjects(numel(sc),1);
for k=1:numel(sc)
    rstar = sc(k).sol.r;
    rmin = max(0.005, rstar*(1-DIAG_WINDOW_MULT));
    rmax = min(0.10,  rstar*(1+DIAG_WINDOW_MULT));
    r_span = linspace(rmin,rmax,N_DIAG);

    alt=sc(k).cfg; alt.fix_r=1; alt.maxit_r=1; alt.maxit_V = MAXITV_DIAG;
    Sgrid=nan(size(r_span));
    for i=1:numel(r_span)
        alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal_Bfixed(alt);
        da_i=tmp.a(2)-tmp.a(1);
        Apriv_i = sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
        Sgrid(i) = Apriv_i - tmp.fiscal.B;
    end
    hS(k)=plot(r_span,Sgrid,'LineWidth',2,'Color',cols(k,:));
    xline(sc(k).sol.r,':','Color',cols(k,:),'LineWidth',1.3);
end
yline(0,'k--'); grid on; xlabel('Interest rate r'); ylabel('S(r)=A_{priv}-B');
legend(hS, names_char, 'Location','best'); title('Excess-asset curve by scenario');
export_fig(fig, fullfile(outdir_figs,'inc_excess_supply_curve'), EXPORT_PDF);

% Fiscal accounts (revenues/expenditures/balances)
fig=figure('Name','Fiscal accounts and balances (compare)');
subplot(1,3,1);
vals=zeros(2,numel(sc));
for k=1:numel(sc), vals(:,k)=[sc(k).sol.fiscal.Tl; sc(k).sol.fiscal.Tc]; end
h=bar(categorical({'Labor tax','VAT'}), vals, 'grouped'); grid on; ylabel('Revenues'); title('Revenues'); legend(names_char,'Location','best');
add_grouped_labels(h,'%.3f');
subplot(1,3,2);
vals=zeros(3,numel(sc));
for k=1:numel(sc), vals(:,k)=[sc(k).sol.fiscal.rB; sc(k).sol.fiscal.G; sc(k).sol.fiscal.Tr]; end
h=bar(categorical({'Debt service','Public good','Transfers'}), vals, 'grouped'); grid on; ylabel('Expenditures'); title('Expenditures'); legend(names_char,'Location','best');
add_grouped_labels(h,'%.3f');
subplot(1,3,3);
vals=zeros(2,numel(sc));
for k=1:numel(sc), vals(:,k)=[sc(k).sol.fiscal.PB; sc(k).sol.fiscal.BB]; end
h=bar(categorical({'Primary balance','Overall (global)'}), vals, 'grouped'); yline(0,'k--');
grid on; ylabel('Balance'); title('Balances'); legend(names_char,'Location','best');
add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'inc_fiscal_accounts'), EXPORT_PDF);

% Borrowers / Lenders
fig=figure('Name','Borrowers-Lenders compare');
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
subplot(1,2,1);
Borrow=zeros(2,numel(sc)); Lend=zeros(2,numel(sc));
for k=1:numel(sc)
    Borrow(:,k)=sc(k).sol.borrowers.fracBorrow(:);
    Lend(:,k)  =sc(k).sol.borrowers.fracLend(:);
end
h=bar(catX, Borrow, 'grouped'); grid on; ylabel('Share'); title('Borrowers (a<0)'); legend(names_char,'Location','best');
add_grouped_labels(h,'%.3f');
subplot(1,2,2);
VolB=zeros(2,numel(sc)); VolL=zeros(2,numel(sc));
for k=1:numel(sc)
    VolB(:,k)=abs(sc(k).sol.borrowers.volBorrow(:));
    VolL(:,k)=sc(k).sol.borrowers.volLend(:);
end
h=bar(catX, VolB, 'grouped'); grid on; ylabel('Volume'); title('|Debt| volumes'); legend(names_char,'Location','best');
add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'inc_borrowers_lenders'), EXPORT_PDF);

% Lorenz (Wealth total)
fig=figure('Name','Lorenz wealth (total) compare'); hold on;
hW=gobjects(numel(sc),1);
for k=1:numel(sc)
    [~,Lw,cumPop]=lorenz_from_density(sc(k).sol.a, sc(k).sol.g);
    hW(k)=plot(cumPop,Lw,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
labels = arrayfun(@(k) sprintf('%s (Gini=%.3f)', pretty(sc(k).name), sc(k).sol.stats.giniW(3)), 1:numel(sc), 'uni',0);
legend(hW, labels, 'Location','southeast'); title('Lorenz curve (Wealth, total)');
export_fig(fig, fullfile(outdir_figs,'inc_lorenz_wealth_total'), EXPORT_PDF);

% Lorenz (Wealth por tipo)
fig=figure('Name','Lorenz wealth by type compare'); hold on;
hI=gobjects(numel(sc),1); hF=hI;
for k=1:numel(sc)
    [LwI,cI]=lorenz_from_assets_single(sc(k).sol.a, sc(k).sol.g(:,1));
    hI(k)=plot(cI,LwI,'LineWidth',2,'Color',cols(k,:),'LineStyle','-');
end
for k=1:numel(sc)
    [LwF,cF]=lorenz_from_assets_single(sc(k).sol.a, sc(k).sol.g(:,2));
    hF(k)=plot(cF,LwF,'LineWidth',2,'Color',cols(k,:),'LineStyle','--');
end
plot([0,1],[0,1],'k:'); grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
labelsI=arrayfun(@(k) sprintf('%s - I (Gini=%.3f)', pretty(sc(k).name), sc(k).sol.stats.giniW(1)), 1:numel(sc),'uni',0);
labelsF=arrayfun(@(k) sprintf('%s - F (Gini=%.3f)', pretty(sc(k).name), sc(k).sol.stats.giniW(2)), 1:numel(sc),'uni',0);
legend([hI;hF],[labelsI, labelsF],'Location','southeast');
title('Lorenz curve (Wealth) by type (I solid, F dashed)');
export_fig(fig, fullfile(outdir_figs,'inc_lorenz_wealth_bytype'), EXPORT_PDF);

% Lorenz (Consumption total)
fig=figure('Name','Lorenz consumption (total) compare'); hold on;
hC=gobjects(numel(sc),1);
for k=1:numel(sc)
    [~,Lc,cumPopC]=lorenz_from_values(sc(k).sol.a, sc(k).sol.g, sc(k).sol.c(:,1)+sc(k).sol.c(:,2));
    hC(k)=plot(cumPopC,Lc,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
labels=arrayfun(@(k) sprintf('%s (Gini=%.3f)', pretty(sc(k).name), sc(k).sol.stats.giniC(3)), 1:numel(sc), 'uni',0);
legend(hC, labels, 'Location','southeast'); title('Lorenz curve (Consumption, total)');
export_fig(fig, fullfile(outdir_figs,'inc_lorenz_consumption_total'), EXPORT_PDF);

% Lorenz (Consumption por tipo)
fig=figure('Name','Lorenz consumption by type compare'); hold on;
hI=gobjects(numel(sc),1); hF=hI;
for k=1:numel(sc)
    [LcI,cI]=lorenz_from_values_single(sc(k).sol.c(:,1), sc(k).sol.g(:,1), sc(k).sol.a);
    hI(k)=plot(cI,LcI,'LineWidth',2,'Color',cols(k,:), 'LineStyle','-');
end
for k=1:numel(sc)
    [LcF,cF]=lorenz_from_values_single(sc(k).sol.c(:,2), sc(k).sol.g(:,2), sc(k).sol.a);
    hF(k)=plot(cF,LcF,'LineWidth',2,'Color',cols(k,:), 'LineStyle','--');
end
plot([0,1],[0,1],'k:'); grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
labelsI=arrayfun(@(k) sprintf('%s - I (Gini=%.3f)', pretty(sc(k).name), sc(k).sol.stats.giniC(1)), 1:numel(sc),'uni',0);
labelsF=arrayfun(@(k) sprintf('%s - F (Gini=%.3f)', pretty(sc(k).name), sc(k).sol.stats.giniC(2)), 1:numel(sc),'uni',0);
legend([hI;hF],[labelsI, labelsF],'Location','southeast');
title('Lorenz curve (Consumption) by type (I solid, F dashed)');
export_fig(fig, fullfile(outdir_figs,'inc_lorenz_consumption_bytype'), EXPORT_PDF);

% MPC(a) — rápida
fig=figure('Name','MPC(a) compare'); hold on;
hI=gobjects(numel(sc),1); hF=hI;
for k=1:numel(sc)
    cfg=sc(k).cfg; r=sc(k).sol.r; a=sc(k).sol.a;
    alt=cfg; alt.fix_r=1; alt.r_guess=r; alt.rmin=r; alt.rmax=r;
    alt.maxit_V = MAXITV_DIAG;
    eps_z=0.01; alt.z1=cfg.z1*(1+eps_z); alt.z2=cfg.z2*(1+eps_z);
    solP=solve_two_type_huggett_fiscal_Bfixed(alt);
    dres1=(cfg.z1*(1+eps_z)-cfg.z1) * (1 + cfg.phi);
    dres2=(cfg.z2*(1+eps_z)-cfg.z2) * (1 - cfg.tau_l);
    MPC1=(solP.c(:,1)-sc(k).sol.c(:,1))/max(dres1,1e-12);
    MPC2=(solP.c(:,2)-sc(k).sol.c(:,2))/max(dres2,1e-12);
    hI(k)=plot(a,MPC1,'Color',cols(k,:),'LineWidth',2,'LineStyle','-');
    hF(k)=plot(a,MPC2,'Color',cols(k,:),'LineWidth',2,'LineStyle','--');
end
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)');
legend([hI;hF],[mk(' - I'), mk(' - F')],'Location','bestoutside');
title('MPC along asset grid (I solid, F dashed)');
export_fig(fig, fullfile(outdir_figs,'inc_mpc_by_assets'), EXPORT_PDF);

fprintf('Figures saved to %s\n', outdir_figs);

% =============================== Helpers =================================
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath, export_pdf)
    if nargin<1||isempty(fig), fig=gcf; end
    if nargin<3, export_pdf=true; end
    print(fig,[basepath '.png'],'-dpng','-r220');
    if export_pdf, print(fig,[basepath '.pdf'],'-dpdf'); end
end
function [gT,L,cumPop]=lorenz_from_density(a,g)
    da=a(2)-a(1); gT=g(:,1)+g(:,2); W=sum(gT)*da; [as,ix]=sort(a); gs=gT(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as-min(0,min(as))+1e-12; cumWealth=cumsum(wealth_pos.*gs); L=cumWealth/max(cumWealth(end),1e-12);
end
function [L,cumPop]=lorenz_from_assets_single(a, gk)
    da=a(2)-a(1); W=sum(gk)*da; [as,ix]=sort(a); gs=gk(ix)*da; cumPop=cumsum(gs)/max(W,1e-12);
    wealth_pos = as - min(0,min(as)) + 1e-12;  cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
function [vals,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w); vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix);
    cumPop=cumsum(w_s)/W; cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end
function [L,cumPop]=lorenz_from_values_single(xk, gk, a)
    da=a(2)-a(1); w = gk(:)*da; W=sum(w); [xs,ix]=sort(xk(:)); ws=w(ix);
    cumPop = cumsum(ws)/max(W,1e-12); cumx = cumsum(xs.*ws); L = cumx / max(cumx(end),1e-12);
end
function add_grouped_labels(h, fmt)
    if nargin<2, fmt = '%.3f'; end
    for i=1:numel(h)
        x = h(i).XEndPoints; y = h(i).YEndPoints;
        for j=1:numel(x)
            text(x(j), y(j), sprintf([' ',fmt], y(j)), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
        end
    end
end
