% ================== main_Bfixed_income_shocks_scenarios.m =================
clear; clc; close all;

% --- carpetas de salida ---
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% -------- Parámetros BASE --------
cfg = struct();
cfg.RRA_I = 3.40; cfg.RRA_F = 3.40;
cfg.rho   = 0.08;
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;
cfg.phi   = 0.09;
cfg.z1 = 0.33; cfg.z2 = 1.00;

cfg.theta_I = 0.06;
cfg.theta_F = 0.01;

cfg.eta_target = 0.654; cfg.p22_bar = 0.8155;

cfg.I = 700; cfg.amax = 3.0;
cfg.amin = -2.0*cfg.z1;

cfg.r_guess = 0.03; cfg.rmin = 0.005; cfg.rmax = 0.10;
cfg.maxit_V = 140; cfg.crit_V = 1e-6; cfg.Delta = 1400;
cfg.maxit_r = 1500; cfg.crit_S = 1e-5; cfg.fix_r = 0;

% ===== Bien público =====
cfg.psi_G   = 0.08;
cfg.omegaG  = 0.50;
cfg.report_G_effects = 1;

% Difusión explícita
cfg.sigma_a = 0.010;

% ---- Gobierno / deuda ----
cfg.B_mode  = 'ratio_to_Y';
cfg.Bbar    = 0.35;
cfg.alphaG  = 0.50;
cfg.clamp_G_to_zero = true;
cfg.G_cap_ratio = 0.08;

% -------- Control: ¿mantener r fijo ajustando B? --------
hold_r_fixed = false;     % << pon true para comparar sólo composición (cierra S≈0 con B/Y)
r_target     = 0.025;

% --------- Estilo figuras ---------
paper_style();

% =============== Escenarios: BASE / caída / caída+Tr↑ =====================
phi_up_rel = 0.15;   % +15% relativo en transferencias
z1_drop    = 0.85;   % -15% informales
z2_drop    = 0.75;   % -25% formales

base_amin  = cfg.amin;   % mismo grid de activos en todos
base_amax  = cfg.amax;

SC = struct('label',{},'cfg',{},'sol',{},'color',{},'ls',{});

% (a) BASE
SC(1).label = 'BASE';
SC(1).cfg   = cfg;
SC(1).color = [0 0 0];
SC(1).ls    = '-';

% (b) INCOME_DROP
cfg_b = cfg;
cfg_b.z1 = cfg.z1 * z1_drop;
cfg_b.z2 = cfg.z2 * z2_drop;
cfg_b.amin = base_amin;  cfg_b.amax = base_amax;
SC(2).label = 'INCOME_DROP';
SC(2).cfg   = cfg_b;
SC(2).color = [0.10 0.45 0.90];
SC(2).ls    = '--';

% (c) INCOME_DROP_TRUP  (income drop + phi ↑ 15%)
cfg_c = cfg_b;
cfg_c.phi = cfg.phi * (1 + phi_up_rel);
SC(3).label = 'INCOME_DROP_TRUP';
SC(3).cfg   = cfg_c;
SC(3).color = [0.85 0.33 0.10];
SC(3).ls    = '-.';

% ========================= Resolver escenarios ============================
for k=1:numel(SC)
    if hold_r_fixed
        [SC(k).cfg, SC(k).sol] = tune_to_r_target_Bbar(SC(k).cfg, r_target, 30, 1e-4);
    else
        SC(k).sol = solve_two_type_huggett_fiscal_Bfixed(SC(k).cfg);
    end
    
    s  = SC(k).sol; fb = s.fiscal;
    fprintf('\n== %s ==\n', SC(k).label);
    fprintf('r* = %.4f  | S_excess = %.3e\n', s.r, s.S_residual);
    fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', s.popI, s.popF, s.popI/(s.popI+s.popF));
    fprintf('Y=%.6f  Ctot=%.6f  Gpc=%.6f  Bbar=%.3f  (mode: %s)\n', s.Y, s.Ctot, s.Gpc, SC(k).cfg.Bbar, SC(k).cfg.B_mode);
    fprintf('Tl=%.6f  Tc=%.6f  Tr=%.6f  G=%.6f  rB=%.6f  -> PB=%.6f  BB=%.6e\n', ...
        fb.Tl, fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.BB);
    fprintf('GiniW: I=%.3f F=%.3f T=%.3f | GiniC: I=%.3f F=%.3f T=%.3f\n', ...
        s.stats.giniW(1),s.stats.giniW(2),s.stats.giniW(3), ...
        s.stats.giniC(1),s.stats.giniC(2),s.stats.giniC(3));
    fprintf('||HJB residual||_∞ ≈ %.3e | clampG=%d | capG=%d\n', s.hjb_residual, fb.clamp_active, fb.cap_active);
end

% ============================ CSV combinados ==============================
TF = []; TA = []; TH = []; TB = [];
for k=1:numel(SC)
    s = SC(k).sol; fb = s.fiscal;
    T_fiscal = table(string(SC(k).label), fb.Tl, fb.Tc, fb.Tl+fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.B, fb.BB, ...
        'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers','public_good','debt_serv','primary_bal','debt_stock','global_bal'});
    TF = [TF; T_fiscal];

    da = s.a(2)-s.a(1);
    A_I = sum(s.g(:,1).*s.a)*da;  A_F = sum(s.g(:,2).*s.a)*da;  A_priv = A_I + A_F;
    T_assets = table(string(SC(k).label), A_I, A_F, A_priv, fb.B, 'VariableNames', ...
        {'scenario','A_I','A_F','A_private','B_public'});
    TA = [TA; T_assets];

    St = s.stats;
    T_stats = table(string(SC(k).label), s.r, s.popI, s.popF, s.Y, s.Ctot, ...
        St.wealth_mean(1),St.wealth_mean(2),St.wealth_mean(3), ...
        St.wealth_median(1),St.wealth_median(2),St.wealth_median(3), ...
        St.giniW(1),St.giniW(2),St.giniW(3), ...
        St.cons_mean(1),St.cons_mean(2),St.cons_mean(3), ...
        St.cons_median(1),St.cons_median(2),St.cons_median(3), ...
        St.giniC(1),St.giniC(2),St.giniC(3), St.p11, ...
        'VariableNames', {'scenario','r','popI','popF','Y','Ctot', ...
        'wealth_mean_I','wealth_mean_F','wealth_mean_T', ...
        'wealth_med_I','wealth_med_F','wealth_med_T', ...
        'giniW_I','giniW_F','giniW_T', ...
        'cons_mean_I','cons_mean_F','cons_mean_T', ...
        'cons_med_I','cons_med_F','cons_med_T', ...
        'giniC_I','giniC_F','giniC_T','p11_rep'});
    TH = [TH; T_stats];

    Borr = s.borrowers;
    T_borr = table(string(SC(k).label), Borr.fracBorrow(1),Borr.fracBorrow(2),Borr.fracLend(1),Borr.fracLend(2), ...
                   Borr.volBorrow(1),Borr.volBorrow(2),Borr.volLend(1),Borr.volLend(2), ...
                   'VariableNames', {'scenario','fracBorrow_I','fracBorrow_F','fracLend_I','fracLend_F', ...
                                     'volBorrow_I','volBorrow_F','volLend_I','volLend_F'});
    TB = [TB; T_borr];
end
writetable(TF, fullfile(outdir_tabs,'scen_income_fiscal_breakdown.csv'));
writetable(TA, fullfile(outdir_tabs,'scen_income_asset_market.csv'));
writetable(TH, fullfile(outdir_tabs,'scen_income_household_stats.csv'));
writetable(TB, fullfile(outdir_tabs,'scen_income_borrowers_lenders.csv'));
fprintf('CSV combinados exportados en %s\n', outdir_tabs);

% =============================== Figuras ==================================
labels = {SC.label};
cats   = categorical(labels);
cats   = reordercats(cats, labels);

% --- Policy c(a) y s(a), overlay por escenario y tipo ---
a = SC(1).sol.a;
fig = figure('Name','Policy: Consumption (compare)');
subplot(1,2,1); hold on; grid on; title('c_I(a) Informal'); xlabel('Assets a'); ylabel('c(a)');
for k=1:numel(SC), plot(a,SC(k).sol.c(:,1),'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2); end
legend({SC.label},'Location','best');
subplot(1,2,2); hold on; grid on; title('c_F(a) Formal'); xlabel('Assets a'); ylabel('c(a)');
for k=1:numel(SC), plot(a,SC(k).sol.c(:,2),'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2); end
legend({SC.label},'Location','best');
save_fig(fig, fullfile(outdir_figs,'policy_consumption_income_compare'));

fig = figure('Name','Policy: Savings (compare)');
subplot(1,2,1); hold on; grid on; title('s_I(a) Informal'); xlabel('Assets a'); ylabel('s(a)');
for k=1:numel(SC), plot(a,SC(k).sol.s(:,1),'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2); end
legend({SC.label},'Location','best');
subplot(1,2,2); hold on; grid on; title('s_F(a) Formal'); xlabel('Assets a'); ylabel('s(a)');
for k=1:numel(SC), plot(a,SC(k).sol.s(:,2),'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2); end
legend({SC.label},'Location','best');
save_fig(fig, fullfile(outdir_figs,'policy_savings_income_compare'));

% --- Densidades (líneas) ---
fig = figure('Name','Wealth density by type (compare)');
subplot(2,1,1); hold on; grid on; title('Informal density g_1(a)'); xlabel('Assets a'); ylabel('density');
for k=1:numel(SC), plot(a,SC(k).sol.g(:,1),'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2); end
legend({SC.label},'Location','best');
subplot(2,1,2); hold on; grid on; title('Formal density g_2(a)'); xlabel('Assets a'); ylabel('density');
for k=1:numel(SC), plot(a,SC(k).sol.g(:,2),'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2); end
legend({SC.label},'Location','best');
save_fig(fig, fullfile(outdir_figs,'wealth_distributions_income_compare'));

% --- Cierre financiero (agregado) ---
fig = figure('Name','Asset market closure (aggregate) compare');
A_priv = zeros(numel(SC),1); B_pub = zeros(numel(SC),1);
for k=1:numel(SC)
    da = SC(k).sol.a(2)-SC(k).sol.a(1);
    A_priv(k) = sum(SC(k).sol.g(:,1).*SC(k).sol.a)*da + sum(SC(k).sol.g(:,2).*SC(k).sol.a)*da;
    B_pub(k)  = SC(k).sol.fiscal.B;
end
bar(cats, [A_priv B_pub]); grid on; ylabel('Assets level');
legend({'A_{priv}','B (public)'},'Location','best');
title(sprintf('Asset market closure: r^* BASE=%.4f | DROP=%.4f | DROP+Tr=%.4f', SC(1).sol.r, SC(2).sol.r, SC(3).sol.r));
save_fig(fig, fullfile(outdir_figs,'asset_market_closure_income_compare'));

% --- Cierre por tipo (stacked) ---
fig = figure('Name','Asset market closure by type (stacked) compare'); hold on;
A_I = zeros(numel(SC),1); A_F = zeros(numel(SC),1);
for k=1:numel(SC)
    da = SC(k).sol.a(2)-SC(k).sol.a(1);
    A_I(k) = sum(SC(k).sol.g(:,1).*SC(k).sol.a)*da;
    A_F(k) = sum(SC(k).sol.g(:,2).*SC(k).sol.a)*da;
end
bar(cats, [A_I A_F B_pub],'stacked'); grid on; ylabel('Assets level');
legend({'A_I','A_F','B'},'Location','best');
title('Composition of closure by scenario');
save_fig(fig, fullfile(outdir_figs,'asset_market_closure_by_type_income_compare'));

% --- Curvas S(r) y A_priv/B vs r (sin memoria de G) ---
r_span = linspace(0.005, 0.10, 31);
figS = figure('Name','Excess asset supply S(r) compare'); hold on; grid on;
figAB= figure('Name','Apriv vs B across r (compare)'); hold on; grid on;

for k=1:numel(SC)
    alt = SC(k).cfg;
    alt.fix_r = 1; alt.maxit_r = 1;
    alt.reset_G_each_call = true;
    alt.G_init = 0;  alt.alphaG = 1.0;

    Sgrid  = nan(size(r_span)); Apriv_grid = Sgrid; B_grid = Sgrid;
    for i=1:numel(r_span)
        alt.r_guess = r_span(i); alt.rmin = r_span(i); alt.rmax = r_span(i);
        tmp = solve_two_type_huggett_fiscal_Bfixed(alt);
        da_i = tmp.a(2)-tmp.a(1);
        Ai = sum(tmp.g(:,1).*tmp.a)*da_i; Af = sum(tmp.g(:,2).*tmp.a)*da_i;
        Apriv_grid(i) = Ai+Af;  B_grid(i) = tmp.fiscal.B;
        Sgrid(i)      = Apriv_grid(i) - B_grid(i);
    end
    figure(figS); plot(r_span,Sgrid,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    figure(figAB); plot(r_span,Apriv_grid,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    figure(figAB); plot(r_span,B_grid,  'LineStyle',':','Color',SC(k).color,'LineWidth',2);
end
figure(figS); yline(0,'k--'); xlabel('Interest rate r'); ylabel('S(r)=A_{priv}-B'); 
legend({SC.label},'Location','best'); title('Asset market diagnostic: S(r) by scenario');
save_fig(figS, fullfile(outdir_figs,'excess_supply_curve_income_compare'));
figure(figAB); xlabel('Interest rate r'); ylabel('Assets / Bonds level');
legend_entries = [strcat("A_{priv} ",string({SC.label})) strcat("B ",string({SC.label}))];
legend(legend_entries,'Location','best'); title('A_{priv}(r) vs B(r) by scenario');
save_fig(figAB, fullfile(outdir_figs,'asset_demand_vs_B_income_compare'));

% --- Fiscal (revenues / exp / balances) --- (¡AJUSTADO!)
fig = figure('Name','Fiscal accounts and balances (compare)');

subplot(1,3,1); hold on; grid on; title('Revenues'); ylabel('amount');
bar(cats, [TF{:,'labor_tax'} TF{:,'vat_tax'} TF{:,'rev_total'}]);
legend({'Labor','VAT','Total'},'Location','best');

subplot(1,3,2); hold on; grid on; title('Expenditures'); ylabel('amount');
bar(cats, [TF{:,'debt_serv'} TF{:,'public_good'} TF{:,'transfers'}]);
legend({'Debt serv','G','Tr'},'Location','best');

subplot(1,3,3); hold on; grid on; title('Balances'); ylabel('amount');
bar(cats, [TF{:,'primary_bal'} TF{:,'global_bal'}]);   % primero bar(...)
yline(0,'k--');                                        % luego yline
legend({'Primary','Overall'},'Location','best');

save_fig(fig, fullfile(outdir_figs,'fiscal_accounts_and_balances_income_compare'));

% --- Borrowers/Lenders (agrupado) ---
shares = zeros(numel(SC),4);  % [B_I, L_I, B_F, L_F]
vols   = zeros(numel(SC),4);  % [|Debt|_I, Save_I, |Debt|_F, Save_F]
for k=1:numel(SC)
    bk = SC(k).sol.borrowers;
    shares(k,:) = [bk.fracBorrow(1), bk.fracLend(1), bk.fracBorrow(2), bk.fracLend(2)];
    vols(k,:)   = [abs(bk.volBorrow(1)), bk.volLend(1), abs(bk.volBorrow(2)), bk.volLend(2)];
end
fig = figure('Name','Borrowers-Lenders shares & volumes (compare)');
subplot(1,2,1); bar(cats, shares, 'grouped'); grid on;
ylabel('Share'); title('Borrowers/Lenders by scenario');
legend({'B\_I','L\_I','B\_F','L\_F'},'Location','best');
subplot(1,2,2); bar(cats, vols, 'grouped'); grid on;
ylabel('Volume'); title('Borrowers/Lenders Volumes');
legend({'|Debt|\_I','Save\_I','|Debt|\_F','Save\_F'},'Location','best');
save_fig(fig, fullfile(outdir_figs,'borrowers_lenders_income_compare'));

% --- Lorenz (riqueza) total ---
fig = figure('Name','Lorenz wealth (total) compare'); hold on; grid on; axis square;
legW = cell(1,numel(SC));
for k=1:numel(SC)
    [~,Lw,cumPop]=lorenz_from_density(SC(k).sol.a,SC(k).sol.g);
    plot(cumPop,Lw,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    legW{k} = sprintf('%s (Gini=%.3f)', SC(k).label, SC(k).sol.stats.giniW(3));
end
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Wealth share');
title('Lorenz curve (Wealth, total)'); legend(legW,'Location','southeast');
save_fig(fig, fullfile(outdir_figs,'lorenz_wealth_total_income_compare'));

% --- Lorenz (riqueza) por tipo: dos subplots (I y F) con su Gini ---
fig = figure('Name','Lorenz wealth by type compare'); 
subplot(1,2,1); hold on; grid on; axis square; title('Informal (Wealth)');
legI = cell(1,numel(SC));
for k=1:numel(SC)
    [LwI,cumI] = lorenz_from_assets_single(SC(k).sol.a, SC(k).sol.g(:,1));
    plot(cumI,LwI,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    legI{k} = sprintf('%s (Gini=%.3f)', SC(k).label, SC(k).sol.stats.giniW(1));
end
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Wealth share'); legend(legI,'Location','southeast');

subplot(1,2,2); hold on; grid on; axis square; title('Formal (Wealth)');
legF = cell(1,numel(SC));
for k=1:numel(SC)
    [LwF,cumF] = lorenz_from_assets_single(SC(k).sol.a, SC(k).sol.g(:,2));
    plot(cumF,LwF,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    legF{k} = sprintf('%s (Gini=%.3f)', SC(k).label, SC(k).sol.stats.giniW(2));
end
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Wealth share'); legend(legF,'Location','southeast');
save_fig(fig, fullfile(outdir_figs,'lorenz_wealth_bytype_income_compare'));

% --- Lorenz (consumo) total ---
fig = figure('Name','Lorenz consumption (total) compare'); hold on; grid on; axis square;
legC = cell(1,numel(SC));
for k=1:numel(SC)
    [~,Lc,cumPopC]=lorenz_from_values(SC(k).sol.a,SC(k).sol.g, SC(k).sol.c(:,1)+SC(k).sol.c(:,2));
    plot(cumPopC,Lc,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    legC{k} = sprintf('%s (Gini=%.3f)', SC(k).label, SC(k).sol.stats.giniC(3));
end
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Consumption share');
title('Lorenz curve (Consumption, total)'); legend(legC,'Location','southeast');
save_fig(fig, fullfile(outdir_figs,'lorenz_consumption_total_income_compare'));

% --- Lorenz (consumo) por tipo (I y F) con su Gini ---
fig = figure('Name','Lorenz consumption by type compare');
subplot(1,2,1); hold on; grid on; axis square; title('Informal (Consumption)');
legCI = cell(1,numel(SC));
for k=1:numel(SC)
    [LcI,cI] = lorenz_from_values_single(SC(k).sol.c(:,1), SC(k).sol.g(:,1), SC(k).sol.a);
    plot(cI,LcI,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    legCI{k} = sprintf('%s (Gini=%.3f)', SC(k).label, SC(k).sol.stats.giniC(1));
end
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Consumption share'); legend(legCI,'Location','southeast');

subplot(1,2,2); hold on; grid on; axis square; title('Formal (Consumption)');
legCF = cell(1,numel(SC));
for k=1:numel(SC)
    [LcF,cF] = lorenz_from_values_single(SC(k).sol.c(:,2), SC(k).sol.g(:,2), SC(k).sol.a);
    plot(cF,LcF,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2);
    legCF{k} = sprintf('%s (Gini=%.3f)', SC(k).label, SC(k).sol.stats.giniC(2));
end
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Consumption share'); legend(legCF,'Location','southeast');
save_fig(fig, fullfile(outdir_figs,'lorenz_consumption_bytype_income_compare'));

% --- Barras de Gini por escenario (Wealth & Consumption, I/F/T) ---
GW = [TH{:,'giniW_I'}, TH{:,'giniW_F'}, TH{:,'giniW_T'}];
GC = [TH{:,'giniC_I'}, TH{:,'giniC_F'}, TH{:,'giniC_T'}];

fig = figure('Name','Gini indices (Wealth & Consumption)'); 
subplot(1,2,1); bar(cats, GW, 'grouped'); grid on; title('Wealth Gini'); ylabel('Gini');
legend({'I','F','Total'},'Location','best');
subplot(1,2,2); bar(cats, GC, 'grouped'); grid on; title('Consumption Gini'); ylabel('Gini');
legend({'I','F','Total'},'Location','best');
save_fig(fig, fullfile(outdir_figs,'gini_indices_income_compare'));

% --- MPC(a): neto de IVA y con r fijo por escenario ---
eps_z=0.01;
fig = figure('Name','MPC by assets (compare)');
subplot(1,2,1); hold on; grid on; title('MPC_I(a)'); xlabel('Assets a'); ylabel('MPC');
for k=1:numel(SC)
    alt = SC(k).cfg; s0 = SC(k).sol;
    alt.fix_r=1; alt.r_guess=s0.r; alt.rmin=s0.r; alt.rmax=s0.r;
    alt.z1 = SC(k).cfg.z1*(1+eps_z); alt.z2 = SC(k).cfg.z2*(1+eps_z);
    solP = solve_two_type_huggett_fiscal_Bfixed(alt);
    dres1 = ((SC(k).cfg.z1*(1+eps_z)-SC(k).cfg.z1)*(1 + SC(k).cfg.phi)) / (1 + SC(k).cfg.tau_c);
    MPC1 = (solP.c(:,1)-s0.c(:,1))/max(dres1,1e-12);
    plot(a,MPC1,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2,'DisplayName',SC(k).label);
end
legend('Location','best');

subplot(1,2,2); hold on; grid on; title('MPC_F(a)'); xlabel('Assets a'); ylabel('MPC');
for k=1:numel(SC)
    alt = SC(k).cfg; s0 = SC(k).sol;
    alt.fix_r=1; alt.r_guess=s0.r; alt.rmin=s0.r; alt.rmax=s0.r;
    alt.z1 = SC(k).cfg.z1*(1+eps_z); alt.z2 = SC(k).cfg.z2*(1+eps_z);
    solP = solve_two_type_huggett_fiscal_Bfixed(alt);
    dres2 = ((SC(k).cfg.z2*(1+eps_z)-SC(k).cfg.z2)*(1 - SC(k).cfg.tau_l)) / (1 + SC(k).cfg.tau_c);
    MPC2 = (solP.c(:,2)-s0.c(:,2))/max(dres2,1e-12);
    plot(a,MPC2,'LineStyle',SC(k).ls,'Color',SC(k).color,'LineWidth',2,'DisplayName',SC(k).label);
end
legend('Location','best');
save_fig(fig, fullfile(outdir_figs,'mpc_by_assets_income_compare'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ================================ HELPERS ================================
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function save_fig(fig,basepath)
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

% ========= Auto-tuning: Bbar para alcanzar r_target =========
function [cfg_out, sol_best] = tune_to_r_target_Bbar(cfg_in, r_target, maxit, tol)
    if nargin<3, maxit=25; end, if nargin<4, tol=1e-4; end
    cfg_out = cfg_in;
    % Bracketing inicial
    BLo = max(0.05, 0.25);  BHi = 1.50;   % deuda/PIB
    [rLo,solLo] = r_given_Bbar(cfg_out,BLo);
    [rHi,solHi] = r_given_Bbar(cfg_out,BHi);
    % Amplía el bracket si es necesario
    it_widen=0;
    while ~((rLo<=r_target && r_target<=rHi) || (rHi<=r_target && r_target<=rLo))
        it_widen=it_widen+1; if it_widen>8, break; end
        BHi = BHi*1.5; [rHi,solHi] = r_given_Bbar(cfg_out,BHi);
        if r_target < min(rLo,rHi), BLo = max(0.02, BLo*0.6); [rLo,solLo] = r_given_Bbar(cfg_out,BLo); end
    end
    % Bisección
    for it=1:maxit
        Bmid = 0.5*(BLo+BHi);
        [rmid,solmid] = r_given_Bbar(cfg_out,Bmid);
        if abs(rmid - r_target) < tol
            cfg_out.Bbar = Bmid; sol_best = solmid; return;
        end
        if (rLo <= rHi && rmid < r_target) || (rHi < rLo && rmid > r_target)
            BLo = Bmid; rLo = rmid; solLo = solmid;
        else
            BHi = Bmid; rHi = rmid; solHi = solmid;
        end
    end
    % Salida
    if abs(rLo - r_target) < abs(rHi - r_target)
        cfg_out.Bbar = BLo; sol_best = solLo;
    else
        cfg_out.Bbar = BHi; sol_best = solHi;
    end
end
function [r_star, sol] = r_given_Bbar(cfg_in,Bbar)
    cfg_tmp = cfg_in; cfg_tmp.Bbar = Bbar;
    sol = solve_two_type_huggett_fiscal_Bfixed(cfg_tmp);
    r_star = sol.r;
end
