%% ============================ FILE: main_base_paper.m ============================
% Paper-ready BASE run (no shocks) for a two-type Huggett model with fiscal policy.
% SAVE THIS BLOCK AS: main_base_paper.m
% - Exports tidy CSVs (fiscal, asset market, household stats, borrowers/lenders)
% - Produces publication-grade figures (English labels) and saves to ./figures
% - Adds diagnostics (market clearing residual, HJB residual) and optional MPC(a)
%
% HOW TO RUN
%   1) Save this block as main_base_paper.m
%   2) Save the next block (solver) as solve_two_type_huggett_fiscal.m
%   3) Ensure subfolders ./tables and ./figures exist (code will create them).
%   4) Run main_base_paper.m
%
% NOTES
%   * Continuous-time HJB solved by implicit upwind FD; stationary Fokker–Planck for g(a).
%   * Government budget closes with endogenous B(r). Asset market clears A_private = B(r).
%   * All figure labels/titles are in ENGLISH, suitable for journal submission.

clear; clc; close all;

%% 1) Base parameters (replicate your baseline, minor renames for clarity)
cfg = struct();
cfg.RRA_I = 3.40;                  % CRRA (informal)
cfg.RRA_F = 3.40;                  % CRRA (formal)
cfg.rho   = 0.05;                  % subjective discount rate
cfg.theta = 0.02;                  % borrowing premium (only if a<0)
cfg.tau_l = 0.15;                  % labor income tax (formal only)
cfg.tau_c = 0.18;                  % VAT (both types)
cfg.Gov   = 0.05;                  % public good in utility (additive flow)
cfg.phi   = 0.09;                  % transfers to informal (phi * z1 per informal)
cfg.z1    = 0.33;                  % informal income
cfg.z2    = 1.00;                  % formal income

cfg.eta_target = 0.654;            % target informality share
cfg.p22_bar    = 0.8155;           % persistence formal (implies λ2)

cfg.I    = 700;                    % grid points
cfg.amax = 5.0;
cfg.amin = -0.30*cfg.z1;           % borrowing limit (ad hoc, tighter than natural)

% interest rate search (bisection) and initial guess
cfg.r_guess = 0.03; cfg.rmin = 0.005; cfg.rmax = 0.08;

% Solver controls
cfg.maxit_V = 100;                 % HJB value iterations
cfg.crit_V  = 1e-6;
cfg.Delta   = 1000;                % implicit scheme weight
cfg.maxit_r = 1000;                % bisection iterations on r
cfg.crit_S  = 1e-5;                % market clearing tolerance
cfg.fix_r   = 0;                   % =1 to hold r fixed at r_guess (used for MPC & curves)

% I/O
outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

%% 2) Solve BASE
base = solve_two_type_huggett_fiscal(cfg);

% Aliases
a = base.a; g = base.g; c = base.c; s = base.s;
popI = base.popI; popF = base.popF; r = base.r; Y = base.Y; Ctot = base.Ctot;
fb = base.fiscal; S_res = base.S_residual; hjb_res = base.hjb_residual;

fprintf('\n== BASE (paper) ==\n');
fprintf('r = %.6f  | S_excess = %.3e (market clears ~ 0)\n', r, S_res);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y = %.6f,   Ctot = %.6f\n', Y, Ctot);
fprintf('PB = %.6f,  rB = %.6f,  BB = %.6f (≈0 in steady state)\n', fb.PB, fb.rB, fb.BB);
fprintf('||HJB residual||_∞ ≈ %.3e\n', hjb_res);

%% 3) Export tidy CSVs (English variable names)
% 3.1 Fiscal breakdown
T_fiscal = table;
T_fiscal.scenario   = "BASE";
T_fiscal.labor_tax  = fb.Tl;
T_fiscal.vat_tax    = fb.Tc;
T_fiscal.rev_total  = fb.Tl + fb.Tc;
T_fiscal.transfers  = fb.Tr;
T_fiscal.public_good= fb.G;
T_fiscal.debt_serv  = fb.rB;
T_fiscal.exp_total  = fb.Tr + fb.G + fb.rB;
T_fiscal.primary_bal= fb.PB;
T_fiscal.debt_stock = fb.B;
T_fiscal.global_bal = fb.BB;
writetable(T_fiscal, fullfile(outdir_tabs,'base_fiscal_breakdown.csv'));

% 3.2 Asset market (private demand vs public supply)
A_priv = sum( (g(:,1)+g(:,2)).*a ) * (a(2)-a(1));
T_assets = table("BASE", A_priv, fb.B, 'VariableNames', {'scenario','A_private','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'base_asset_market.csv'));

% 3.3 Household statistics
S = base.stats;
T_stats = table;
T_stats.scenario = "BASE";
T_stats.r = r; T_stats.popI = popI; T_stats.popF = popF; T_stats.Y = Y; T_stats.Ctot = Ctot;
T_stats.wealth_mean_I = S.wealth_mean(1);
T_stats.wealth_mean_F = S.wealth_mean(2);
T_stats.wealth_mean_T = S.wealth_mean(3);
T_stats.wealth_med_I  = S.wealth_median(1);
T_stats.wealth_med_F  = S.wealth_median(2);
T_stats.wealth_med_T  = S.wealth_median(3);
T_stats.giniW_I       = S.giniW(1);
T_stats.giniW_F       = S.giniW(2);
T_stats.giniW_T       = S.giniW(3);
T_stats.cons_mean_I   = S.cons_mean(1);
T_stats.cons_mean_F   = S.cons_mean(2);
T_stats.cons_mean_T   = S.cons_mean(3);
T_stats.cons_med_I    = S.cons_median(1);
T_stats.cons_med_F    = S.cons_median(2);
T_stats.cons_med_T    = S.cons_median(3);
T_stats.giniC_I       = S.giniC(1);
T_stats.giniC_F       = S.giniC(2);
T_stats.giniC_T       = S.giniC(3);
T_stats.p11_rep       = S.p11;  % informal persistence proxy
writetable(T_stats, fullfile(outdir_tabs,'base_household_stats.csv'));

% 3.4 Borrowers / lenders
Borr = base.borrowers;
T_borr = table;
T_borr.scenario   = "BASE";
T_borr.fracBorrow_I = Borr.fracBorrow(1);
T_borr.fracBorrow_F = Borr.fracBorrow(2);
T_borr.fracLend_I   = Borr.fracLend(1);
T_borr.fracLend_F   = Borr.fracLend(2);
T_borr.volBorrow_I  = Borr.volBorrow(1);  % < 0
T_borr.volBorrow_F  = Borr.volBorrow(2);
T_borr.volLend_I    = Borr.volLend(1);    % > 0
T_borr.volLend_F    = Borr.volLend(2);
writetable(T_borr, fullfile(outdir_tabs,'base_borrowers_lenders.csv'));

fprintf('Exported CSVs to %s\n', outdir_tabs);

%% 4) Figures (publication style)
paper_style(); da = a(2)-a(1);

% 4.1 Consumption policy c(a) by type
fig = figure('Name','Policy: Consumption');
plot(a, c(:,1), 'LineWidth',2); hold on;
plot(a, c(:,2), 'LineWidth',2);
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('Consumption c(a)');
legend({'Informal','Formal'},'Location','best'); title('Consumption Policies');
export_fig(fig, fullfile(outdir_figs,'policy_consumption'));

% 4.2 Savings policy s(a) by type
fig = figure('Name','Policy: Savings');
plot(a, s(:,1), 'LineWidth',2); hold on;
plot(a, s(:,2), 'LineWidth',2);
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('Savings s(a)');
legend({'Informal','Formal'},'Location','best'); title('Savings Policies');
export_fig(fig, fullfile(outdir_figs,'policy_savings'));

% 4.3 Wealth distributions g(a)
fig = figure('Name','Wealth distributions');
subplot(2,1,1);
bar(a, g(:,1), 'FaceAlpha',0.7, 'EdgeColor','none'); grid on;
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density');
subplot(2,1,2);
bar(a, g(:,2), 'FaceAlpha',0.7, 'EdgeColor','none'); grid on;
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density');
export_fig(fig, fullfile(outdir_figs,'wealth_distributions'));

% 4.4 Asset market: Private demand vs Public supply at equilibrium r
fig = figure('Name','Asset market closure');
bar(categorical({'Private demand','Public supply B(r)'}), [A_priv, fb.B]);
ylabel('Assets level'); grid on; title(sprintf('Asset market closes at r = %.4f', r));
export_fig(fig, fullfile(outdir_figs,'asset_market_closure'));

% 4.5 Fiscal accounts
fig = figure('Name','Fiscal accounts');
subplot(1,2,1);
bar(categorical({'Labor tax','VAT'}), [fb.Tl, fb.Tc]); ylabel('Revenues'); title('Revenues'); grid on;
subplot(1,2,2);
bar(categorical({'Transfers','Public good','Debt service'}), [fb.Tr, fb.G, fb.rB]); ylabel('Expenditures'); title('Expenditures'); grid on;
export_fig(fig, fullfile(outdir_figs,'fiscal_accounts'));

% 4.6 Borrowers/lenders: fractions and volumes
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});
fig = figure('Name','Borrowers-Lenders shares');
subplot(1,2,1); bar(catX, [Borr.fracBorrow(:)]); ylabel('Fraction borrowers (a<0)'); title('Borrowers'); grid on;
subplot(1,2,2); bar(catX, [Borr.fracLend(:)]);   ylabel('Fraction lenders (a>0)');   title('Lenders'); grid on;
export_fig(fig, fullfile(outdir_figs,'borrowers_lenders_shares'));

fig = figure('Name','Borrowers-Lenders volumes');
subplot(1,2,1); bar(catX, abs([Borr.volBorrow(:)])); ylabel('|Aggregate debt|');  title('Debt volume (a<0)'); grid on;
subplot(1,2,2); bar(catX, [Borr.volLend(:)]);     ylabel('Aggregate savings'); title('Savings volume (a>0)'); grid on;
export_fig(fig, fullfile(outdir_figs,'borrowers_lenders_volumes'));

% 4.7 Lorenz curve (WEALTH)
fig = figure('Name','Lorenz wealth');
[gT, Lw, cumPop] = lorenz_from_density(a, g);
plot(cumPop, Lw, 'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Wealth share');
text(0.6,0.1,sprintf('Gini_W = %.3f', base.stats.giniW(3)),'FontSize',12);
export_fig(fig, fullfile(outdir_figs,'lorenz_wealth'));

% 4.7b Lorenz curve (WEALTH) by type
fig = figure('Name','Lorenz wealth by type');
[LwI, cumPopI] = lorenz_from_density_single(a, g(:,1));
[LwF, cumPopF] = lorenz_from_density_single(a, g(:,2));
plot(cumPopI, LwI, 'LineWidth',2); hold on;
plot(cumPopF, LwF, 'LineWidth',2);
plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Wealth share');
legend({'Informal','Formal','45°'},'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'lorenz_wealth_by_type'));

% 4.8 Lorenz curve (CONSUMPTION)
fig = figure('Name','Lorenz consumption');
[~, Lc, cumPopC] = lorenz_from_values(a, g, c(:,1)+c(:,2));
plot(cumPopC, Lc, 'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Consumption share');
text(0.6,0.1,sprintf('Gini_C = %.3f', base.stats.giniC(3)),'FontSize',12);
export_fig(fig, fullfile(outdir_figs,'lorenz_consumption'));

% 4.9 Excess asset supply curve S(r) around equilibrium (diagnostic)
r_span = linspace(max(1e-3,cfg.r_guess*0.4), cfg.r_guess*2.0, 21);
Sgrid  = nan(size(r_span));
alt = cfg; alt.fix_r = 1; alt.maxit_r = 1;  % disable bisection
for i=1:numel(r_span)
    alt.r_guess = r_span(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
    tmp = solve_two_type_huggett_fiscal(alt);
    Sgrid(i) = tmp.S_residual;  % A_private - B(r)
end
fig = figure('Name','Excess asset supply S(r)');
plot(r_span, Sgrid, 'LineWidth',2); yline(0,'k--'); grid on;
xlabel('Interest rate r'); ylabel('Excess assets S(r)'); title('Asset market diagnostic');
export_fig(fig, fullfile(outdir_figs,'excess_supply_curve'));

% 4.10 Optional: MPC(a) by type via small income shock (finite diff, r fixed)
DO_MPC = true;
if DO_MPC
    eps_z = 0.01; alt = cfg; alt.fix_r = 1; alt.r_guess = r; alt.rmin=r; alt.rmax=r;  % hold r
    alt.z1 = cfg.z1 * (1+eps_z); alt.z2 = cfg.z2 * (1+eps_z);
    solP = solve_two_type_huggett_fiscal(alt);
    % resources change by type
    dres1 = (cfg.z1*(1+eps_z) - cfg.z1) * (1 + cfg.phi);         % dz1*(1+phi)
    dres2 = (cfg.z2*(1+eps_z) - cfg.z2) * (1 - cfg.tau_l);       % dz2*(1-tau_l)
    MPC1 = (solP.c(:,1) - c(:,1)) / max(dres1,1e-12);
    MPC2 = (solP.c(:,2) - c(:,2)) / max(dres2,1e-12);
    fig = figure('Name','MPC by assets');
    plot(a, MPC1, 'LineWidth',2); hold on; plot(a, MPC2, 'LineWidth',2);
    yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)');
    legend({'Informal','Formal'},'Location','best'); title('MPC along asset grid (finite-diff)');
    export_fig(fig, fullfile(outdir_figs,'mpc_by_assets'));
end

fprintf('Saved figures to %s\n', outdir_figs);

%% ======= Local helpers (style, export, Lorenz) =======
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

function [L, cumPop] = lorenz_from_density_single(a, gcol)
    % Lorenz curve of wealth for a single-type density gcol(a)
    da = a(2)-a(1); W = sum(gcol)*da; W = max(W, 1e-12);
    [as, ix] = sort(a); gs = gcol(ix)*da; cumPop = cumsum(gs)/W;
    wealth_pos = as - min(0,min(as)) + 1e-12; cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end

function [vals, L, cumPop] = lorenz_from_values(a, g, x)
    da = a(2)-a(1); w = (g(:,1)+g(:,2))*da; W = sum(w);
    vals = x(:); [vals_s, ix] = sort(vals); w_s = w(ix); cumPop = cumsum(w_s)/W;
    cumx = cumsum(vals_s .* w_s); L = cumx / max(cumx(end),1e-12);
end


