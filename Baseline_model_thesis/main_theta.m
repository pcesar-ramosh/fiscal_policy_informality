%% ============================ FILE: scenario_theta_compare.m ============================
% Compares steady states for LOWER, BASELINE, and HIGHER borrowing spread (theta)
% - theta affects the effective rate on a<0: r_eff = r + theta (informals)
% Reuses solve_two_type_huggett_fiscal.m. Exports paper-ready figures & a summary CSV.
%
% OUTPUTS
%   ./figures/compare_theta_*   (PNG+PDF)
%   ./tables/compare_theta_summary.csv

clear; clc; close all;

% ---------- 1) Base settings shared across scenarios ----------
base = struct();
base.RRA_I = 3.40; base.RRA_F = 3.40;
base.rho   = 0.05;
base.theta = 0.02;                 % BASELINE borrowing spread
base.tau_l = 0.15;
base.tau_c = 0.18;
base.Gov   = 0.05;
base.phi   = 0.09;
base.z1    = 0.33; base.z2 = 1.00;
base.p22_bar = 0.8155;
base.I = 700; base.amax = 5.0; base.amin = -0.30*base.z1;
base.r_guess = 0.03; base.rmin = 0.005; base.rmax = 0.08;
base.maxit_V = 100; base.crit_V = 1e-6; base.Delta = 1000;
base.maxit_r = 1000; base.crit_S = 1e-5; base.fix_r = 0;
base.eta_target = 0.654;           % as in your BASE paper

% Three theta scenarios (adjust if you prefer other values)
theta_low  = 0.01;
theta_mid  = 0.02;   % baseline
theta_high = 0.05;

cfgL = base; cfgL.theta = theta_low;
cfgM = base; cfgM.theta = theta_mid;
cfgH = base; cfgH.theta = theta_high;

outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end
paper_style();

% ---------- 2) Solve steady states ----------
solL = solve_two_type_huggett_fiscal(cfgL);
solM = solve_two_type_huggett_fiscal(cfgM);
solH = solve_two_type_huggett_fiscal(cfgH);

% Aliases
L = solL; M = solM; H = solH; a = L.a; da = a(2)-a(1);
labL = sprintf('Lower spread (theta=%.2f)', theta_low);
labM = sprintf('Baseline (theta=%.2f)',    theta_mid);
labH = sprintf('Higher spread (theta=%.2f)', theta_high);

% ---------- 3) Summary CSV ----------
S = table;
S.scenario = ["theta_low"; "theta_mid"; "theta_high"];
S.theta    = [theta_low;  theta_mid;  theta_high];
S.r        = [L.r; M.r; H.r];
S.popI     = [L.popI; M.popI; H.popI];
S.popF     = [L.popF; M.popF; H.popF];
S.Y        = [L.Y; M.Y; H.Y];
S.Ctot     = [L.Ctot; M.Ctot; H.Ctot];
S.Tl       = [L.fiscal.Tl; M.fiscal.Tl; H.fiscal.Tl];
S.Tc       = [L.fiscal.Tc; M.fiscal.Tc; H.fiscal.Tc];
S.Tr       = [L.fiscal.Tr; M.fiscal.Tr; H.fiscal.Tr];
S.G        = [L.fiscal.G;  M.fiscal.G;  H.fiscal.G];
S.B        = [L.fiscal.B;  M.fiscal.B;  H.fiscal.B];
S.PB       = [L.fiscal.PB; M.fiscal.PB; H.fiscal.PB];
S.rB       = [L.fiscal.rB; M.fiscal.rB; H.fiscal.rB];
S.BB       = [L.fiscal.BB; M.fiscal.BB; H.fiscal.BB];
S.GiniW_T  = [L.stats.giniW(3); M.stats.giniW(3); H.stats.giniW(3)];
S.GiniC_T  = [L.stats.giniC(3); M.stats.giniC(3); H.stats.giniC(3)];
S.p11_rep  = [L.stats.p11; M.stats.p11; H.stats.p11];
writetable(S, fullfile('./tables','compare_theta_summary.csv'));

% ---------- 4) Policy functions (consumption & savings) ----------
fig = figure('Name','Policy: Consumption vs theta');
subplot(1,2,1);
plot(a, L.c(:,1),'LineWidth',2); hold on;
plot(a, M.c(:,1),'--','LineWidth',2);
plot(a, H.c(:,1),'LineWidth',2);
xlabel('Assets a'); ylabel('Consumption c(a)'); title('Informal'); grid on;
legend({labL,labM,labH},'Location','best');
subplot(1,2,2);
plot(a, L.c(:,2),'LineWidth',2); hold on;
plot(a, M.c(:,2),'--','LineWidth',2);
plot(a, H.c(:,2),'LineWidth',2);
xlabel('Assets a'); ylabel('Consumption c(a)'); title('Formal'); grid on;
legend({labL,labM,labH},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_theta_policy_consumption'));

fig = figure('Name','Policy: Savings vs theta');
subplot(1,2,1);
plot(a, L.s(:,1),'LineWidth',2); hold on;
plot(a, M.s(:,1),'--','LineWidth',2);
plot(a, H.s(:,1),'LineWidth',2);
xlabel('Assets a'); ylabel('Savings s(a)'); title('Informal'); grid on;
legend({labL,labM,labH},'Location','best');
subplot(1,2,2);
plot(a, L.s(:,2),'LineWidth',2); hold on;
plot(a, M.s(:,2),'--','LineWidth',2);
plot(a, H.s(:,2),'LineWidth',2);
xlabel('Assets a'); ylabel('Savings s(a)'); title('Formal'); grid on;
legend({labL,labM,labH},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_theta_policy_savings'));

% ---------- 5) Wealth distributions g(a) ----------
fig = figure('Name','Wealth density vs theta');
subplot(2,1,1);
plot(a, L.g(:,1),'LineWidth',2); hold on;
plot(a, M.g(:,1),'--','LineWidth',2);
plot(a, H.g(:,1),'LineWidth',2);
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density'); grid on;
legend({labL,labM,labH},'Location','best');
subplot(2,1,2);
plot(a, L.g(:,2),'LineWidth',2); hold on;
plot(a, M.g(:,2),'--','LineWidth',2);
plot(a, H.g(:,2),'LineWidth',2);
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density'); grid on;
legend({labL,labM,labH},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_theta_wealth_distributions'));

% ---------- 6) Asset market closure ----------
A_priv_L = sum( (L.g(:,1)+L.g(:,2)).*a ) * da;
A_priv_M = sum( (M.g(:,1)+M.g(:,2)).*a ) * da;
A_priv_H = sum( (H.g(:,1)+H.g(:,2)).*a ) * da;
Mmat = [A_priv_L, A_priv_M, A_priv_H;
        L.fiscal.B, M.fiscal.B, H.fiscal.B];
fig = figure('Name','Asset market closure vs theta');
bar(Mmat); grid on; ylabel('Assets level');
set(gca,'XTickLabel',{'Private demand','Public supply B(r)'});
legend({labL,labM,labH},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_theta_asset_market'));

% ---------- 7) Fiscal accounts ----------
fig = figure('Name','Fiscal accounts vs theta');
subplot(1,2,1);
bar([L.fiscal.Tl, M.fiscal.Tl, H.fiscal.Tl;
     L.fiscal.Tc, M.fiscal.Tc, H.fiscal.Tc]);
set(gca,'XTickLabel',{'Labor tax','VAT'}); ylabel('Revenues'); grid on;
legend({labL,labM,labH},'Location','best'); title('Revenues');
subplot(1,2,2);
bar([L.fiscal.Tr, M.fiscal.Tr, H.fiscal.Tr;
     L.fiscal.G,  M.fiscal.G,  H.fiscal.G;
     L.fiscal.rB, M.fiscal.rB, H.fiscal.rB]);
set(gca,'XTickLabel',{'Transfers','Public good','Debt service'});
ylabel('Expenditures'); grid on; legend({labL,labM,labH},'Location','best');
title('Expenditures');
export_fig(fig, fullfile(outdir_figs,'compare_theta_fiscal_accounts'));

% ---------- 8) Borrowers/Lenders ----------
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});
fig = figure('Name','Borrowers-Lenders shares vs theta');
subplot(1,2,1);
bar(catX,[L.borrowers.fracBorrow(:), M.borrowers.fracBorrow(:), H.borrowers.fracBorrow(:)]);
ylabel('Fraction borrowers (a<0)'); title('Borrowers'); grid on; legend({labL,labM,labH},'Location','best');
subplot(1,2,2);
bar(catX,[L.borrowers.fracLend(:),   M.borrowers.fracLend(:),   H.borrowers.fracLend(:)]);
ylabel('Fraction lenders (a>0)');   title('Lenders');   grid on; legend({labL,labM,labH},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_theta_borrowers_lenders_shares'));

fig = figure('Name','Borrowers-Lenders volumes vs theta');
subplot(1,2,1);
bar(catX,[abs(L.borrowers.volBorrow(:)), abs(M.borrowers.volBorrow(:)), abs(H.borrowers.volBorrow(:))]);
ylabel('|Aggregate debt|'); title('Debt volume (a<0)'); grid on; legend({labL,labM,labH},'Location','best');
subplot(1,2,2);
bar(catX,[L.borrowers.volLend(:), M.borrowers.volLend(:), H.borrowers.volLend(:)]);
ylabel('Aggregate savings'); title('Savings volume (a>0)'); grid on; legend({labL,labM,labH},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_theta_borrowers_lenders_volumes'));

% ---------- 9) Lorenz curves ----------
[~, LwL, cumPopL] = lorenz_from_density(a, L.g);
[~, LwM, cumPopM] = lorenz_from_density(a, M.g);
[~, LwH, cumPopH] = lorenz_from_density(a, H.g);
fig = figure('Name','Lorenz wealth vs theta');
plot(cumPopL,LwL,'LineWidth',2); hold on;
plot(cumPopM,LwM,'--','LineWidth',2);
plot(cumPopH,LwH,'LineWidth',2);
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Wealth share'); grid on; axis square;
legend({sprintf('%s (Gini=%.3f)',labL,L.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labM,M.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labH,H.stats.giniW(3))}, 'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'compare_theta_lorenz_wealth'));

[~, LcL, cumPopCL] = lorenz_from_values(a, L.g, L.c(:,1)+L.c(:,2));
[~, LcM, cumPopCM] = lorenz_from_values(a, M.g, M.c(:,1)+M.c(:,2));
[~, LcH, cumPopCH] = lorenz_from_values(a, H.g, H.c(:,1)+H.c(:,2));
fig = figure('Name','Lorenz consumption vs theta');
plot(cumPopCL,LcL,'LineWidth',2); hold on;
plot(cumPopCM,LcM,'--','LineWidth',2);
plot(cumPopCH,LcH,'LineWidth',2);
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Consumption share'); grid on; axis square;
legend({sprintf('%s (Gini=%.3f)',labL,L.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labM,M.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labH,H.stats.giniC(3))}, 'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'compare_theta_lorenz_consumption'));

% ---------- 10) MPC(a) overlay by type ----------
% Same tiny income bump as in previous analyses (keeps comparability).
DO_MPC = true;
if DO_MPC
    eps_z = 0.01;
    % low-theta equilibrium
    altL = cfgL; altL.fix_r = 1; altL.r_guess = L.r; altL.rmin=L.r; altL.rmax=L.r;
    altL.z1 = base.z1*(1+eps_z); altL.z2 = base.z2*(1+eps_z);
    solPL = solve_two_type_huggett_fiscal(altL);
    dres1L = (base.z1*eps_z) * (1 + base.phi);
    dres2L = (base.z2*eps_z) * (1 - base.tau_l);
    MPC1L = (solPL.c(:,1) - L.c(:,1)) / max(dres1L,1e-12);
    MPC2L = (solPL.c(:,2) - L.c(:,2)) / max(dres2L,1e-12);

    % baseline theta
    altM = cfgM; altM.fix_r = 1; altM.r_guess = M.r; altM.rmin=M.r; altM.rmax=M.r;
    altM.z1 = base.z1*(1+eps_z); altM.z2 = base.z2*(1+eps_z);
    solPM = solve_two_type_huggett_fiscal(altM);
    MPC1M = (solPM.c(:,1) - M.c(:,1)) / max(dres1L,1e-12);
    MPC2M = (solPM.c(:,2) - M.c(:,2)) / max(dres2L,1e-12);

    % high-theta equilibrium
    altH = cfgH; altH.fix_r = 1; altH.r_guess = H.r; altH.rmin=H.r; altH.rmax=H.r;
    altH.z1 = base.z1*(1+eps_z); altH.z2 = base.z2*(1+eps_z);
    solPH = solve_two_type_huggett_fiscal(altH);
    MPC1H = (solPH.c(:,1) - H.c(:,1)) / max(dres1L,1e-12);
    MPC2H = (solPH.c(:,2) - H.c(:,2)) / max(dres2L,1e-12);

    fig = figure('Name','MPC by assets vs theta');
    subplot(1,2,1);
    plot(a, MPC1L,'LineWidth',2); hold on;
    plot(a, MPC1M,'--','LineWidth',2);
    plot(a, MPC1H,'LineWidth',2);
    xlabel('Assets a'); ylabel('MPC(a)'); title('Informal'); grid on;
    legend({labL,labM,labH},'Location','best');

    subplot(1,2,2);
    plot(a, MPC2L,'LineWidth',2); hold on;
    plot(a, MPC2M,'--','LineWidth',2);
    plot(a, MPC2H,'LineWidth',2);
    xlabel('Assets a'); ylabel('MPC(a)'); title('Formal'); grid on;
    legend({labL,labM,labH},'Location','best');

    export_fig(fig, fullfile(outdir_figs,'compare_theta_mpc_by_assets'));
end

% ---------- 11) Excess supply S(r) curves ----------
DO_S_CURVE = true;
if DO_S_CURVE
    rgridL = linspace(max(1e-3, L.r*0.5), L.r*1.8, 15);
    rgridM = linspace(max(1e-3, M.r*0.5), M.r*1.8, 15);
    rgridH = linspace(max(1e-3, H.r*0.5), H.r*1.8, 15);
    SgridL = nan(size(rgridL)); SgridM = nan(size(rgridM)); SgridH = nan(size(rgridH));
    alt = base; alt.maxit_r = 1; alt.fix_r = 1;
    for i=1:numel(rgridL)
        alt.theta = theta_low;  alt.eta_target = base.eta_target;
        alt.r_guess = rgridL(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp = solve_two_type_huggett_fiscal(alt); SgridL(i) = tmp.S_residual;
    end
    for i=1:numel(rgridM)
        alt.theta = theta_mid;  alt.eta_target = base.eta_target;
        alt.r_guess = rgridM(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp = solve_two_type_huggett_fiscal(alt); SgridM(i) = tmp.S_residual;
    end
    for i=1:numel(rgridH)
        alt.theta = theta_high; alt.eta_target = base.eta_target;
        alt.r_guess = rgridH(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp = solve_two_type_huggett_fiscal(alt); SgridH(i) = tmp.S_residual;
    end
    fig = figure('Name','Excess supply S(r) vs theta');
    plot(rgridL, SgridL,'LineWidth',2); hold on;
    plot(rgridM, SgridM,'--','LineWidth',2);
    plot(rgridH, SgridH,'LineWidth',2);
    yline(0,'k--'); grid on;
    xlabel('Interest rate r'); ylabel('Excess assets S(r)');
    legend({labL,labM,labH},'Location','best');
    export_fig(fig, fullfile(outdir_figs,'compare_theta_excess_supply_curve'));
end

% ---------- 12) Bond demand vs supply (by r) ----------
rgridL = linspace(max(1e-3, L.r*0.5), L.r*1.8, 41);
rgridM = linspace(max(1e-3, M.r*0.5), M.r*1.8, 41);
rgridH = linspace(max(1e-3, H.r*0.5), H.r*1.8, 41);
[~, A_L, B_L] = bond_curves_by_r(cfgL, rgridL(1), rgridL(end), numel(rgridL));
[~, A_M, B_M] = bond_curves_by_r(cfgM, rgridM(1), rgridM(end), numel(rgridM));
[~, A_H, B_H] = bond_curves_by_r(cfgH, rgridH(1), rgridH(end), numel(rgridH));
fig = figure('Name','Bond demand vs supply (by interest rate, theta)');
tiledlayout(1,2);

nexttile;
plot(rgridL, A_L,'LineWidth',2); hold on;
plot(rgridL, B_L,'LineWidth',2);
plot(rgridM, A_M,'--','LineWidth',2);
plot(rgridM, B_M,'--','LineWidth',2);
xline(L.r,'k--'); yline(L.fiscal.B,'k:'); scatter(L.r,L.fiscal.B,40,'k','filled');
grid on; xlabel('Interest rate r'); ylabel('Assets / Debt level');
title(labL);
legend({'A_{priv}(low \theta)','B(low \theta)','A_{priv}(base)','B(base)','r^*_{low}','B^*_{low}','Equilibrium'},'Location','best');

nexttile;
plot(rgridH, A_H,'LineWidth',2); hold on;
plot(rgridH, B_H,'LineWidth',2);
plot(rgridM, A_M,'--','LineWidth',2);
plot(rgridM, B_M,'--','LineWidth',2);
xline(H.r,'k--'); yline(H.fiscal.B,'k:'); scatter(H.r,H.fiscal.B,40,'k','filled');
grid on; xlabel('Interest rate r'); ylabel('Assets / Debt level');
title(labH);
legend({'A_{priv}(high \theta)','B(high \theta)','A_{priv}(base)','B(base)','r^*_{high}','B^*_{high}','Equilibrium'},'Location','best');

export_fig(fig, fullfile(outdir_figs,'compare_theta_bond_demand_supply'));

% ---------- 12b) Bond demand by type vs public supply ----------
rgridL = linspace(max(1e-3, L.r*0.6), L.r*1.6, 41);
rgridM = linspace(max(1e-3, M.r*0.6), M.r*1.6, 41);
rgridH = linspace(max(1e-3, H.r*0.6), H.r*1.6, 41);
[~, A1_L, A2_L, B_L] = bond_curves_by_r_types(cfgL, rgridL(1), rgridL(end), numel(rgridL));
[~, A1_M, A2_M, B_M] = bond_curves_by_r_types(cfgM, rgridM(1), rgridM(end), numel(rgridM));
[~, A1_H, A2_H, B_H] = bond_curves_by_r_types(cfgH, rgridH(1), rgridH(end), numel(rgridH));

fig = figure('Name','Bond demand by type vs public supply (by r, theta)');
tiledlayout(1,2);

nexttile;
plot(rgridL,A1_L,'LineWidth',2); hold on;
plot(rgridL,A2_L,'LineWidth',2);
plot(rgridL,B_L, 'LineWidth',2);
plot(rgridM,A1_M,'--','LineWidth',2);
plot(rgridM,A2_M,'--','LineWidth',2);
plot(rgridM,B_M, '--','LineWidth',2);
xline(L.r,'k--'); yline(L.fiscal.B,'k:'); scatter(L.r,L.fiscal.B,40,'k','filled');
grid on; xlabel('Interest rate r'); ylabel('Assets / Debt level');
title(labL);
legend({'A_I(low \theta)','A_F(low \theta)','B(low \theta)', ...
        'A_I(base)','A_F(base)','B(base)','r^*_{low}','B^*_{low}','Equilibrium'},'Location','best');

nexttile;
plot(rgridH,A1_H,'LineWidth',2); hold on;
plot(rgridH,A2_H,'LineWidth',2);
plot(rgridH,B_H, 'LineWidth',2);
plot(rgridM,A1_M,'--','LineWidth',2);
plot(rgridM,A2_M,'--','LineWidth',2);
plot(rgridM,B_M, '--','LineWidth',2);
xline(H.r,'k--'); yline(H.fiscal.B,'k:'); scatter(H.r,H.fiscal.B,40,'k','filled');
grid on; xlabel('Interest rate r'); ylabel('Assets / Debt level');
title(labH);
legend({'A_I(high \theta)','A_F(high \theta)','B(high \theta)', ...
        'A_I(base)','A_F(base)','B(base)','r^*_{high}','B^*_{high}','Equilibrium'},'Location','best');

export_fig(fig, fullfile(outdir_figs,'compare_theta_bond_demand_by_type'));

fprintf('All comparative figures saved to %s\n', outdir_figs);

% ============================ Helpers (keep at the very end) =========================
function paper_style()
    set(groot,'defaulttextinterpreter','tex');
    set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex');
    set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.8);
    set(groot,'DefaultFigureColor','w');
end

function export_fig(fig, basepath)
    if nargin<1 || isempty(fig), fig=gcf; end
    print(fig, [basepath '.png'], '-dpng', '-r300');
    print(fig, [basepath '.pdf'], '-dpdf');
end

function [gT, L, cumPop] = lorenz_from_density(a, g)
    da = a(2)-a(1); gT = g(:,1)+g(:,2);
    W = sum(gT)*da; [as, ix] = sort(a); gs = gT(ix)*da; cumPop = cumsum(gs)/W;
    wealth_pos = as - min(0,min(as)) + 1e-12; cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end

function [vals, L, cumPop] = lorenz_from_values(a, g, x)
    da = a(2)-a(1); w = (g(:,1)+g(:,2))*da; W = sum(w);
    vals = x(:); [vals_s, ix] = sort(vals); w_s = w(ix); cumPop = cumsum(w_s)/W;
    cumx = cumsum(vals_s .* w_s); L = cumx / max(cumx(end),1e-12);
end

function [rgrid, Agrid, Bgrid] = bond_curves_by_r(cfg, r_lo, r_hi, n)
    rgrid = linspace(r_lo, r_hi, n);
    Agrid = nan(size(rgrid)); Bgrid = nan(size(rgrid));
    alt = cfg; alt.maxit_r = 1; alt.fix_r = 1;  % solve once per r (no bisection)
    for i = 1:numel(rgrid)
        alt.r_guess = rgrid(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
        sol = solve_two_type_huggett_fiscal(alt);
        a = sol.a; da = a(2)-a(1);
        Agrid(i) = sum((sol.g(:,1)+sol.g(:,2)).*a)*da;
        PB = sol.fiscal.PB; Bgrid(i) = PB / max(rgrid(i),1e-12);
    end
end

function [rgrid, A1grid, A2grid, Bgrid] = bond_curves_by_r_types(cfg, r_lo, r_hi, n)
    rgrid  = linspace(r_lo, r_hi, n);
    A1grid = nan(size(rgrid)); A2grid = nan(size(rgrid)); Bgrid = nan(size(rgrid));
    alt = cfg; alt.maxit_r = 1; alt.fix_r = 1;
    for i = 1:numel(rgrid)
        alt.r_guess = rgrid(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
        sol = solve_two_type_huggett_fiscal(alt);
        a  = sol.a; da = a(2)-a(1);
        A1grid(i) = sum(sol.g(:,1).*a)*da;   % informal
        A2grid(i) = sum(sol.g(:,2).*a)*da;   % formal
        PB = sol.fiscal.PB;
        Bgrid(i) = PB / max(rgrid(i), 1e-12);
    end
end
