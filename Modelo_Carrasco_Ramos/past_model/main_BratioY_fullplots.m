% ====================== main_BratioY_fullplots.m ======================
clear; clc; close all;

outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------- Parámetros base ----------
cfg = struct();

% Preferencias e ingresos
cfg.RRA_I = 3.40;  cfg.RRA_F = 3.40;
cfg.rho   = 0.08;
cfg.z1    = 0.33;  cfg.z2    = 1.00;

% Impuestos y transferencias
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;
cfg.phi   = 0.09;   % transf a informales como fracción de z1

% Primas por tipo (a<0)
cfg.theta_I = 0.06;
cfg.theta_F = 0.01;

% Ocupaciones / objetivo de informalidad
cfg.eta_target = 0.654;
cfg.p22_bar    = 0.8155;

% Grid de activos y numérico
cfg.I      = 700;
cfg.amin   = -2.0*cfg.z1;
cfg.amax   =  3.0;
cfg.r_guess= 0.03;  cfg.rmin=0.005; cfg.rmax=0.10;
cfg.maxit_V= 160;   cfg.crit_V=1e-6; cfg.Delta=1400;
cfg.maxit_r= 120;   cfg.crit_S=1e-5; cfg.fix_r=0;

% Bien público (en utilidad y MU)
cfg.psi_G  = 0.08;
cfg.omegaG = 0.50;
cfg.report_G_effects = 1;

% Difusión en activos (suaviza FP)
cfg.sigma_a = 0.010;

% ---- Gobierno / deuda: B = Bbar * Y  ----
cfg.B_mode = 'ratio_to_Y';   % deuda como % del PIB
cfg.Bbar   = 0.35;           % objetivo B/Y

% ------------------ Resolver ------------------
sol = solve_two_type_huggett_fiscal_Bfixed(cfg);

% ------------------ Reporte -------------------
fb = sol.fiscal;
fprintf('\n== EQUILIBRIO (B/Y=%.3f) ==\n', fb.B_over_Y);
fprintf('r* = %.4f   |   S(r*) = %.2e\n', sol.r, sol.S_residual);
fprintf('popI=%.4f  popF=%.4f  (eta=%.3f)\n', sol.popI, sol.popF, sol.popI/(sol.popI+sol.popF));
fprintf('Y=%.6f  Ctot=%.6f  Gpc=%.6f\n', sol.Y, sol.Ctot, sol.Gpc);
fprintf('Tl=%.6f  Tc=%.6f  Tr=%.6f  G=%.6f  rB=%.6f\n', fb.Tl, fb.Tc, fb.Tr, fb.G, fb.rB);
fprintf('PB=%.6e  BB=%.6e  (BB debe ser 0)\n', fb.PB, fb.BB);

% ------------------ CSVs ----------------------
da = sol.a(2)-sol.a(1);
A_I = sum(sol.g(:,1).*sol.a)*da;  A_F = sum(sol.g(:,2).*sol.a)*da;  A_priv = A_I + A_F;

T_fiscal = table("BASE", fb.B_over_Y, fb.B, fb.rB, fb.Tl, fb.Tc, fb.Tr, fb.G, fb.PB, fb.BB, ...
    'VariableNames', {'scenario','B_over_Y','B_level','rB','Tl','Tc','Tr','G','PB','BB'});
writetable(T_fiscal, fullfile(outdir_tabs,'fiscal_breakdown.csv'));

T_assets = table("BASE", A_I, A_F, A_priv, fb.B, ...
    'VariableNames', {'scenario','A_I','A_F','A_priv','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'asset_market.csv'));

S = sol.stats;
T_stats = table("BASE", sol.r, sol.popI, sol.popF, sol.Y, sol.Ctot, ...
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
writetable(T_stats, fullfile(outdir_tabs,'household_stats.csv'));

Borr = sol.borrowers;
T_borr = table("BASE", Borr.fracBorrow(1),Borr.fracBorrow(2),Borr.fracLend(1),Borr.fracLend(2), ...
               Borr.volBorrow(1),Borr.volBorrow(2),Borr.volLend(1),Borr.volLend(2), ...
               'VariableNames', {'scenario','fracBorrow_I','fracBorrow_F','fracLend_I','fracLend_F', ...
                                 'volBorrow_I','volBorrow_F','volLend_I','volLend_F'});
writetable(T_borr, fullfile(outdir_tabs,'borrowers_lenders.csv'));
fprintf('CSV exportados en %s\n', outdir_tabs);

% ------------------ GRÁFICOS ------------------
paper_style();

% 1) Políticas: Consumo c(a)
fig=figure('Name','Policy: Consumption');
plot(sol.a,sol.c(:,1),'LineWidth',2); hold on; plot(sol.a,sol.c(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('Consumption c(a)'); legend({'Informal','Formal'},'Location','best');
title('Consumption Policies');
export_fig(fig, fullfile(outdir_figs,'policy_consumption'));

% 2) Políticas: Ahorro s(a)
fig=figure('Name','Policy: Savings');
plot(sol.a,sol.s(:,1),'LineWidth',2); hold on; plot(sol.a,sol.s(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('Savings s(a)'); legend({'Informal','Formal'},'Location','best');
title('Savings Policies');
export_fig(fig, fullfile(outdir_figs,'policy_savings'));

% 3) Densidades de riqueza por tipo
fig=figure('Name','Wealth distributions');
subplot(2,1,1); bar(sol.a,sol.g(:,1),'FaceAlpha',0.9,'EdgeColor','none'); grid on;
xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density');
subplot(2,1,2); bar(sol.a,sol.g(:,2),'FaceAlpha',0.9,'EdgeColor','none'); grid on;
xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density');
export_fig(fig, fullfile(outdir_figs,'wealth_distributions'));

% 4) Cierre financiero (barras)
fig=figure('Name','Asset market closure (aggregate)');
bar(categorical({'Private demand','Public debt B'}), [A_priv, fb.B]);
ylabel('Assets level'); grid on; title(sprintf('Asset market closes at r* = %.4f', sol.r));
export_fig(fig, fullfile(outdir_figs,'asset_market_closure'));

% 5) Cierre por tipo (apilado)
fig=figure('Name','Asset market closure by type (stacked)');
X = categorical({'Private demand','Public debt B'}); X=reordercats(X,{'Private demand','Public debt B'});
Ybars = [A_I A_F 0; 0 0 fb.B];
bar(X, Ybars, 'stacked'); grid on; ylabel('Assets level');
legend({'A_I (informal)','A_F (formal)','B (public)'},'Location','best');
title(sprintf('Composition at r* = %.4f', sol.r));
export_fig(fig, fullfile(outdir_figs,'asset_market_closure_by_type'));

% 6) Curvas S(r) y A(r), B(r)
r_span = linspace(max(0.005,cfg.r_guess*0.4), cfg.r_guess*2.2, 41);
Sgrid  = nan(size(r_span)); Apriv_grid = Sgrid; Apriv_I_grid = Sgrid; Apriv_F_grid = Sgrid;
B_grid = Sgrid; Y_grid = Sgrid;
alt = cfg; alt.fix_r = 1; alt.maxit_r = 1;
for i=1:numel(r_span)
    alt.r_guess = r_span(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
    tmp = solve_two_type_huggett_fiscal_Bfixed(alt);
    da_i = tmp.a(2)-tmp.a(1);
    Ai = sum(tmp.g(:,1).*tmp.a)*da_i; Af = sum(tmp.g(:,2).*tmp.a)*da_i;
    Apriv_I_grid(i) = Ai;  Apriv_F_grid(i) = Af;  Apriv_grid(i) = Ai + Af;
    Y_grid(i) = tmp.Y;     B_grid(i) = tmp.fiscal.B;
    Sgrid(i)  = Apriv_grid(i) - B_grid(i);
end
r_star = sol.r; ix = find(diff(sign(Sgrid))~=0,1,'first');
if ~isempty(ix), r1=r_span(ix); r2=r_span(ix+1); s1=Sgrid(ix); s2=Sgrid(ix+1);
    r_star = r1 - s1*(r2-r1)/max(s2-s1,1e-12);
end

fig=figure('Name','Excess asset supply S(r)');
plot(r_span,Sgrid,'LineWidth',2); hold on; yline(0,'k--'); xline(r_star,'r:','LineWidth',1.5);
grid on; xlabel('Interest rate r'); ylabel('S(r)=A_{priv}(r)-B(r)'); title('Asset market diagnostic');
export_fig(fig, fullfile(outdir_figs,'excess_supply_curve'));

fig=figure('Name','Private demand vs Public supply');
plot(r_span,Apriv_I_grid,'LineWidth',2); hold on; plot(r_span,Apriv_F_grid,'LineWidth',2);
plot(r_span,Apriv_grid,'LineWidth',2);   plot(r_span,B_grid,'k--','LineWidth',1.8);
xline(r_star,'r:','LineWidth',1.5);
grid on; xlabel('Interest rate r'); ylabel('Assets / Bonds level');
legend({'A_I(r)','A_F(r)','A_{priv}(r)','B(r)','r^*'},'Location','northwest');
title('Private demand for bonds vs Public bond supply');
export_fig(fig, fullfile(outdir_figs,'asset_demand_vs_B_function_of_r'));

% 7) Fiscal: ingresos/egresos y balances
fig=figure('Name','Fiscal accounts and balances');
subplot(2,2,1); bar(categorical({'Labor tax','VAT'}), [fb.Tl, fb.Tc]); ylabel('Revenues'); title('Revenues'); grid on;
subplot(2,2,2); bar(categorical({'Debt service','Public good','Transfers'}), [fb.rB, fb.G, fb.Tr]); ylabel('Expenditures'); title('Expenditures'); grid on;
subplot(2,2,3); bar(categorical({'Primary balance','Overall (global)'}), [fb.PB, fb.BB]); yline(0,'k--'); ylabel('Balance'); title('Balances'); grid on;
subplot(2,2,4); bar(categorical({'Revenues total','Expenditures total'}), [fb.Tl+fb.Tc, fb.rB+fb.G+fb.Tr]); title('Totals'); grid on;
export_fig(fig, fullfile(outdir_figs,'fiscal_accounts_balances'));

% 8) Borrowers/Lenders
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
fig=figure('Name','Borrowers-Lenders: shares & volumes');
subplot(1,2,1); bar(catX, [sol.borrowers.fracBorrow(:) sol.borrowers.fracLend(:)]); grid on;
ylabel('Share'); legend({'Borrowers (a<0)','Lenders (a>0)'},'Location','best'); title('Shares by type');
subplot(1,2,2); bar(catX, [abs(sol.borrowers.volBorrow(:)) sol.borrowers.volLend(:)]); grid on;
ylabel('Volume'); legend({'|Debt|','Savings'},'Location','best'); title('Volumes by type');
export_fig(fig, fullfile(outdir_figs,'borrowers_lenders'));

% 9) Lorenz de riqueza (total y por tipo)
fig=figure('Name','Lorenz wealth (total)');
[~,Lw,cumPop]=lorenz_from_density(sol.a,sol.g);
plot(cumPop,Lw,'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Wealth share');
text(0.60,0.12,sprintf('Gini_W^T = %.3f', sol.stats.giniW(3)),'FontSize',12);
title('Lorenz curve (Wealth, total)');
export_fig(fig, fullfile(outdir_figs,'lorenz_wealth_total'));

fig=figure('Name','Lorenz wealth by type');
[LwI,cumI] = lorenz_from_assets_single(sol.a, sol.g(:,1));
[LwF,cumF] = lorenz_from_assets_single(sol.a, sol.g(:,2));
plot(cumI,LwI,'LineWidth',2); hold on; plot(cumF,LwF,'LineWidth',2); plot([0,1],[0,1],'k--');
grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
legend({sprintf('Informal (Gini=%.3f)', sol.stats.giniW(1)), ...
        sprintf('Formal (Gini=%.3f)',   sol.stats.giniW(2)), '45°'}, 'Location','southeast');
title('Lorenz curve (Wealth) by type');
export_fig(fig, fullfile(outdir_figs,'lorenz_wealth_bytype'));

% 10) Lorenz de consumo (total y por tipo)
fig=figure('Name','Lorenz consumption (total)');
[~,Lc,cumPopC]=lorenz_from_values(sol.a,sol.g,sol.c(:,1)+sol.c(:,2));
plot(cumPopC,Lc,'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Consumption share');
text(0.60,0.12,sprintf('Gini_C^T = %.3f', sol.stats.giniC(3)),'FontSize',12);
title('Lorenz curve (Consumption, total)');
export_fig(fig, fullfile(outdir_figs,'lorenz_consumption_total'));

fig=figure('Name','Lorenz consumption by type');
[LcI,cI] = lorenz_from_values_single(sol.c(:,1), sol.g(:,1), sol.a);
[LcF,cF] = lorenz_from_values_single(sol.c(:,2), sol.g(:,2), sol.a);
plot(cI, LcI, 'LineWidth',2); hold on; plot(cF, LcF, 'LineWidth',2); plot([0,1],[0,1],'k--');
grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
legend({sprintf('Informal (Gini=%.3f)', sol.stats.giniC(1)), ...
        sprintf('Formal (Gini=%.3f)',   sol.stats.giniC(2)), '45°'}, 'Location','southeast');
title('Lorenz curve (Consumption) by type');
export_fig(fig, fullfile(outdir_figs,'lorenz_consumption_bytype'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ================= Helpers locales =================
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
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
function [L,cumPop]=lorenz_from_assets_single(a, gk)
    da=a(2)-a(1); W=sum(gk)*da; [as,ix]=sort(a); gs=gk(ix)*da; cumPop=cumsum(gs)/max(W,1e-12);
    wealth_pos = as - min(0,min(as)) + 1e-12;  cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
function [L,cumPop]=lorenz_from_values_single(xk, gk, a)
    da=a(2)-a(1); w = gk(:)*da; W=sum(w); [xs,ix]=sort(xk(:)); ws=w(ix);
    cumPop = cumsum(ws)/max(W,1e-12);  cumx = cumsum(xs.*ws);  L = cumx / max(cumx(end),1e-12);
end
