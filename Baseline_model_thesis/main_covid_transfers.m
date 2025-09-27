%% ===================== FILE: scenario_covid_transfers_compare.m =====================
% COVID shock: common fall in incomes (z1, z2) + three transfer policies:
%  (1) Fixed transfers (baseline phi)
%  (2) Higher transfers phi = 0.13
%  (3) Higher transfers phi = 0.15
%
% Reuses solve_two_type_huggett_fiscal.m
% Outputs: ./figures/compare_covid_phi_*  and ./tables/compare_covid_phi_summary.csv

clear; clc; close all;

% ---------- 1) Baseline settings ----------
base = struct();
base.RRA_I = 3.40; base.RRA_F = 3.40;
base.rho   = 0.05;
base.theta = 0.02;                 % extra spread for a<0 (informal)
base.tau_l = 0.15;                 % labor tax
base.tau_c = 0.18;                 % VAT
base.Gov   = 0.05;                 % public good (additive in utility)
base.phi   = 0.09;                 % baseline transfers to informals (phi*z1)
base.z1    = 0.33; base.z2 = 1.00; % pre-COVID productivities/incomes
base.p22_bar = 0.8155;             % formal persistence (legacy)
base.I = 700; base.amax = 5.0; base.amin = -0.30*base.z1;
base.r_guess = 0.03; base.rmin = 0.005; base.rmax = 0.08;
base.maxit_V = 100; base.crit_V = 1e-6; base.Delta = 1000;
base.maxit_r = 1000; base.crit_S = 1e-5; base.fix_r = 0;
base.eta_target = 0.654;           % informality share target

% ---------- COVID income shock ----------
% Fraction of pre-COVID income retained (same for z1 and z2).
% Adjust this to change the severity of the shock.
y_shock = 0.85;     % => 15% drop in incomes (informal & formal)
z1_covid = base.z1 * y_shock;
z2_covid = base.z2 * y_shock;

% ---------- Transfer policies under COVID ----------
phi_fixed = base.phi;   % fixed baseline transfers despite income fall
phi_013   = 0.13;       % higher transfers 13%
phi_015   = 0.15;       % higher transfers 15%

% Build the three configs (all with COVID incomes)
cfgF = base; cfgF.z1=z1_covid; cfgF.z2=z2_covid; cfgF.phi=phi_fixed;
cfgH = base; cfgH.z1=z1_covid; cfgH.z2=z2_covid; cfgH.phi=phi_013;
cfgX = base; cfgX.z1=z1_covid; cfgX.z2=z2_covid; cfgX.phi=phi_015;

outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end
paper_style();

% ---------- 2) Solve steady states ----------
solF = solve_two_type_huggett_fiscal(cfgF);   % fixed phi
solH = solve_two_type_huggett_fiscal(cfgH);   % phi=0.13
solX = solve_two_type_huggett_fiscal(cfgX);   % phi=0.15

% Aliases
F = solF; H = solH; X = solX; a = F.a; da = a(2)-a(1);
labF = sprintf('Fixed transfers (phi=%.2f)', phi_fixed);
labH = sprintf('Higher transfers (phi=%.2f)', phi_013);
labX = sprintf('Higher transfers (phi=%.2f)', phi_015);

% ---------- 3) Summary CSV ----------
S = table;
S.scenario = ["phi_fixed"; "phi_013"; "phi_015"];
S.y_shock  = [y_shock; y_shock; y_shock];
S.phi      = [phi_fixed; phi_013; phi_015];
S.r        = [F.r; H.r; X.r];
S.popI     = [F.popI; H.popI; X.popI];
S.popF     = [F.popF; H.popF; X.popF];
S.Y        = [F.Y; H.Y; X.Y];
S.Ctot     = [F.Ctot; H.Ctot; X.Ctot];
S.Tl       = [F.fiscal.Tl; H.fiscal.Tl; X.fiscal.Tl];
S.Tc       = [F.fiscal.Tc; H.fiscal.Tc; X.fiscal.Tc];
S.Tr       = [F.fiscal.Tr; H.fiscal.Tr; X.fiscal.Tr];
S.G        = [F.fiscal.G;  H.fiscal.G;  X.fiscal.G];
S.B        = [F.fiscal.B;  H.fiscal.B;  X.fiscal.B];
S.PB       = [F.fiscal.PB; H.fiscal.PB; X.fiscal.PB];
S.rB       = [F.fiscal.rB; H.fiscal.rB; X.fiscal.rB];
S.BB       = [F.fiscal.BB; H.fiscal.BB; X.fiscal.BB];
S.GiniW_T  = [F.stats.giniW(3); H.stats.giniW(3); X.stats.giniW(3)];
S.GiniC_T  = [F.stats.giniC(3); H.stats.giniC(3); X.stats.giniC(3)];
S.p11_rep  = [F.stats.p11; H.stats.p11; X.stats.p11];
S.Tr_over_Y = [F.fiscal.Tr/F.Y; H.fiscal.Tr/H.Y; X.fiscal.Tr/X.Y];
S.B_over_Y  = [F.fiscal.B/F.Y;  H.fiscal.B/H.Y;  X.fiscal.B/X.Y];
writetable(S, fullfile('./tables','compare_covid_phi_summary.csv'));

% ---------- 4) Policy functions ----------
fig = figure('Name','Policy: Consumption under COVID (by transfers)');
subplot(1,2,1);
plot(a,F.c(:,1),'LineWidth',2); hold on; plot(a,H.c(:,1),'--','LineWidth',2); plot(a,X.c(:,1),'LineWidth',2);
xlabel('Assets a'); ylabel('Consumption c(a)'); title('Informal'); grid on; legend({labF,labH,labX},'Location','best');
subplot(1,2,2);
plot(a,F.c(:,2),'LineWidth',2); hold on; plot(a,H.c(:,2),'--','LineWidth',2); plot(a,X.c(:,2),'LineWidth',2);
xlabel('Assets a'); ylabel('Consumption c(a)'); title('Formal'); grid on; legend({labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_policy_consumption'));

fig = figure('Name','Policy: Savings under COVID (by transfers)');
subplot(1,2,1);
plot(a,F.s(:,1),'LineWidth',2); hold on; plot(a,H.s(:,1),'--','LineWidth',2); plot(a,X.s(:,1),'LineWidth',2);
xlabel('Assets a'); ylabel('Savings s(a)'); title('Informal'); grid on; legend({labF,labH,labX},'Location','best');
subplot(1,2,2);
plot(a,F.s(:,2),'LineWidth',2); hold on; plot(a,H.s(:,2),'--','LineWidth',2); plot(a,X.s(:,2),'LineWidth',2);
xlabel('Assets a'); ylabel('Savings s(a)'); title('Formal'); grid on; legend({labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_policy_savings'));

% ---------- 5) Wealth densities ----------
fig = figure('Name','Wealth density under COVID (by transfers)');
subplot(2,1,1);
plot(a,F.g(:,1),'LineWidth',2); hold on; plot(a,H.g(:,1),'--','LineWidth',2); plot(a,X.g(:,1),'LineWidth',2);
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density'); grid on; legend({labF,labH,labX},'Location','best');
subplot(2,1,2);
plot(a,F.g(:,2),'LineWidth',2); hold on; plot(a,H.g(:,2),'--','LineWidth',2); plot(a,X.g(:,2),'LineWidth',2);
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density'); grid on; legend({labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_wealth_distributions'));

% ---------- 6) Asset market closure ----------
A_priv_F = sum( (F.g(:,1)+F.g(:,2)).*a ) * da;
A_priv_H = sum( (H.g(:,1)+H.g(:,2)).*a ) * da;
A_priv_X = sum( (X.g(:,1)+X.g(:,2)).*a ) * da;
Mmat = [A_priv_F, A_priv_H, A_priv_X;
        F.fiscal.B, H.fiscal.B, X.fiscal.B];
fig = figure('Name','Asset market closure under COVID (by transfers)');
bar(Mmat); grid on; ylabel('Assets level');
set(gca,'XTickLabel',{'Private demand','Public supply B(r)'});
legend({labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_asset_market'));

% ---------- 7) Fiscal accounts ----------
fig = figure('Name','Fiscal accounts under COVID (by transfers)');
subplot(1,2,1);
bar([F.fiscal.Tl, H.fiscal.Tl, X.fiscal.Tl;
     F.fiscal.Tc, H.fiscal.Tc, X.fiscal.Tc]);
set(gca,'XTickLabel',{'Labor tax','VAT'}); ylabel('Revenues'); grid on;
legend({labF,labH,labX},'Location','best'); title('Revenues');
subplot(1,2,2);
bar([F.fiscal.Tr, H.fiscal.Tr, X.fiscal.Tr;
     F.fiscal.G,  H.fiscal.G,  X.fiscal.G;
     F.fiscal.rB, H.fiscal.rB, X.fiscal.rB]);
set(gca,'XTickLabel',{'Transfers','Public good','Debt service'});
ylabel('Expenditures'); grid on; legend({labF,labH,labX},'Location','best');
title('Expenditures');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_fiscal_accounts'));

% ---------- 8) Borrowers / Lenders ----------
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});
fig = figure('Name','Borrowers-Lenders shares under COVID');
subplot(1,2,1);
bar(catX,[F.borrowers.fracBorrow(:), H.borrowers.fracBorrow(:), X.borrowers.fracBorrow(:)]);
ylabel('Fraction borrowers (a<0)'); title('Borrowers'); grid on; legend({labF,labH,labX},'Location','best');
subplot(1,2,2);
bar(catX,[F.borrowers.fracLend(:),   H.borrowers.fracLend(:),   X.borrowers.fracLend(:)]);
ylabel('Fraction lenders (a>0)');   title('Lenders');   grid on; legend({labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_borrowers_lenders_shares'));

fig = figure('Name','Borrowers-Lenders volumes under COVID');
subplot(1,2,1);
bar(catX,[abs(F.borrowers.volBorrow(:)), abs(H.borrowers.volBorrow(:)), abs(X.borrowers.volBorrow(:))]);
ylabel('|Aggregate debt|'); title('Debt volume (a<0)'); grid on; legend({labF,labH,labX},'Location','best');
subplot(1,2,2);
bar(catX,[F.borrowers.volLend(:), H.borrowers.volLend(:), X.borrowers.volLend(:)]);
ylabel('Aggregate savings'); title('Savings volume (a>0)'); grid on; legend({labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_borrowers_lenders_volumes'));

% ---------- 9) Lorenz curves ----------
[~, LwF, cumPopF] = lorenz_from_density(a, F.g);
[~, LwH, cumPopH] = lorenz_from_density(a, H.g);
[~, LwX, cumPopX] = lorenz_from_density(a, X.g);
fig = figure('Name','Lorenz wealth under COVID');
plot(cumPopF,LwF,'LineWidth',2); hold on;
plot(cumPopH,LwH,'--','LineWidth',2);
plot(cumPopX,LwX,'LineWidth',2);
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Wealth share'); grid on; axis square;
legend({sprintf('%s (Gini=%.3f)',labF,F.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labH,H.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labX,X.stats.giniW(3))}, 'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_lorenz_wealth'));

[~, LcF, cumPopCF] = lorenz_from_values(a, F.g, F.c(:,1)+F.c(:,2));
[~, LcH, cumPopCH] = lorenz_from_values(a, H.g, H.c(:,1)+H.c(:,2));
[~, LcX, cumPopCX] = lorenz_from_values(a, X.g, X.c(:,1)+X.c(:,2));
fig = figure('Name','Lorenz consumption under COVID');
plot(cumPopCF,LcF,'LineWidth',2); hold on;
plot(cumPopCH,LcH,'--','LineWidth',2);
plot(cumPopCX,LcX,'LineWidth',2);
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Consumption share'); grid on; axis square;
legend({sprintf('%s (Gini=%.3f)',labF,F.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labH,H.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labX,X.stats.giniC(3))}, 'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_lorenz_consumption'));

% ---------- 10) MPC(a) overlay (tiny income bump around each equilibrium) ----------
DO_MPC = true;
if DO_MPC
    eps_z = 0.01;
    % F
    altF = cfgF; altF.fix_r = 1; altF.r_guess=F.r; altF.rmin=F.r; altF.rmax=F.r;
    altF.z1 = cfgF.z1*(1+eps_z); altF.z2 = cfgF.z2*(1+eps_z);
    solPF = solve_two_type_huggett_fiscal(altF);
    dres1 = (cfgF.z1*eps_z) * (1 + cfgF.phi);
    dres2 = (cfgF.z2*eps_z) * (1 - cfgF.tau_l);
    MPC1F = (solPF.c(:,1)-F.c(:,1))/max(dres1,1e-12);
    MPC2F = (solPF.c(:,2)-F.c(:,2))/max(dres2,1e-12);
    % H
    altH = cfgH; altH.fix_r = 1; altH.r_guess=H.r; altH.rmin=H.r; altH.rmax=H.r;
    altH.z1 = cfgH.z1*(1+eps_z); altH.z2 = cfgH.z2*(1+eps_z);
    solPH = solve_two_type_huggett_fiscal(altH);
    MPC1H = (solPH.c(:,1)-H.c(:,1))/max(dres1,1e-12);
    MPC2H = (solPH.c(:,2)-H.c(:,2))/max(dres2,1e-12);
    % X
    altX = cfgX; altX.fix_r = 1; altX.r_guess=X.r; altX.rmin=X.r; altX.rmax=X.r;
    altX.z1 = cfgX.z1*(1+eps_z); altX.z2 = cfgX.z2*(1+eps_z);
    solPX = solve_two_type_huggett_fiscal(altX);
    MPC1X = (solPX.c(:,1)-X.c(:,1))/max(dres1,1e-12);
    MPC2X = (solPX.c(:,2)-X.c(:,2))/max(dres2,1e-12);

    fig = figure('Name','MPC by assets under COVID (by transfers)');
    subplot(1,2,1);
    plot(a,MPC1F,'LineWidth',2); hold on;
    plot(a,MPC1H,'--','LineWidth',2);
    plot(a,MPC1X,'LineWidth',2);
    xlabel('Assets a'); ylabel('MPC(a)'); title('Informal'); grid on; legend({labF,labH,labX},'Location','best');

    subplot(1,2,2);
    plot(a,MPC2F,'LineWidth',2); hold on;
    plot(a,MPC2H,'--','LineWidth',2);
    plot(a,MPC2X,'LineWidth',2);
    xlabel('Assets a'); ylabel('MPC(a)'); title('Formal'); grid on; legend({labF,labH,labX},'Location','best');

    export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_mpc_by_assets'));
end

% ---------- 11) Excess supply S(r) ----------
DO_S_CURVE = true;
if DO_S_CURVE
    rgridF = linspace(max(1e-3,F.r*0.5), F.r*1.8, 15);
    rgridH = linspace(max(1e-3,H.r*0.5), H.r*1.8, 15);
    rgridX = linspace(max(1e-3,X.r*0.5), X.r*1.8, 15);
    SgridF = nan(size(rgridF)); SgridH = nan(size(rgridH)); SgridX = nan(size(rgridX));
    alt = base; alt.maxit_r = 1; alt.fix_r = 1;
    for i=1:numel(rgridF)
        alt.z1=z1_covid; alt.z2=z2_covid; alt.phi=phi_fixed;
        alt.r_guess=rgridF(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal(alt); SgridF(i)=tmp.S_residual;
    end
    for i=1:numel(rgridH)
        alt.z1=z1_covid; alt.z2=z2_covid; alt.phi=phi_013;
        alt.r_guess=rgridH(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal(alt); SgridH(i)=tmp.S_residual;
    end
    for i=1:numel(rgridX)
        alt.z1=z1_covid; alt.z2=z2_covid; alt.phi=phi_015;
        alt.r_guess=rgridX(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal(alt); SgridX(i)=tmp.S_residual;
    end
    fig = figure('Name','Excess supply S(r) under COVID');
    plot(rgridF,SgridF,'LineWidth',2); hold on;
    plot(rgridH,SgridH,'--','LineWidth',2);
    plot(rgridX,SgridX,'LineWidth',2);
    yline(0,'k--'); grid on; xlabel('Interest rate r'); ylabel('Excess assets S(r)');
    legend({labF,labH,labX},'Location','best');
    export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_excess_supply_curve'));
end

% ---------- 12) Bond demand vs supply (by r) ----------
rgridF = linspace(max(1e-3,F.r*0.5), F.r*1.8, 41);
rgridH = linspace(max(1e-3,H.r*0.5), H.r*1.8, 41);
rgridX = linspace(max(1e-3,X.r*0.5), X.r*1.8, 41);
[~, A_F, B_F] = bond_curves_by_r(cfgF, rgridF(1), rgridF(end), numel(rgridF));
[~, A_H, B_H] = bond_curves_by_r(cfgH, rgridH(1), rgridH(end), numel(rgridH));
[~, A_X, B_X] = bond_curves_by_r(cfgX, rgridX(1), rgridX(end), numel(rgridX));

fig = figure('Name','Bond demand vs supply (COVID, by transfers)');
tiledlayout(1,2);

% Left: fixed vs phi=0.13
nexttile;
plot(rgridF,A_F,'LineWidth',2); hold on; plot(rgridF,B_F,'LineWidth',2);
plot(rgridH,A_H,'--','LineWidth',2); plot(rgridH,B_H,'--','LineWidth',2);
xline(F.r,'k--'); yline(F.fiscal.B,'k:'); scatter(F.r,F.fiscal.B,40,'k','filled');
grid on; xlabel('Interest rate r'); ylabel('Assets / Debt level');
title('Fixed vs 13%'); legend({'A_{priv}(fixed)','B(fixed)','A_{priv}(13%)','B(13%)','r^*','B^*','Eq.'},'Location','best');

% Right: fixed vs phi=0.15
nexttile;
plot(rgridF,A_F,'LineWidth',2); hold on; plot(rgridF,B_F,'LineWidth',2);
plot(rgridX,A_X,'--','LineWidth',2); plot(rgridX,B_X,'--','LineWidth',2);
xline(X.r,'k--'); yline(X.fiscal.B,'k:'); scatter(X.r,X.fiscal.B,40,'k','filled');
grid on; xlabel('Interest rate r'); ylabel('Assets / Debt level');
title('Fixed vs 15%'); legend({'A_{priv}(fixed)','B(fixed)','A_{priv}(15%)','B(15%)','r^*','B^*','Eq.'},'Location','best');

export_fig(fig, fullfile(outdir_figs,'compare_covid_phi_bond_demand_supply'));

% ============================ Helpers (end of file) ============================
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
    alt = cfg; alt.maxit_r = 1; alt.fix_r = 1;
    for i = 1:numel(rgrid)
        alt.r_guess = rgrid(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
        sol = solve_two_type_huggett_fiscal(alt);
        a = sol.a; da = a(2)-a(1);
        Agrid(i) = sum((sol.g(:,1)+sol.g(:,2)).*a)*da;
        PB = sol.fiscal.PB; Bgrid(i) = PB / max(rgrid(i),1e-12);
    end
end
