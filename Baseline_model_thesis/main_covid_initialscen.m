%% ===================== FILE: main_covid_initialscen.m =====================
% Compara PRE-COVID vs COVID (caída asimétrica de ingresos) bajo 3 políticas de transferencias:
%   - Pre-COVID: baseline (z1,z2; phi=baseline)
%   - COVID-Fijo: z1,z2 caen (z1*0.65, z2*0.75); phi = baseline (fijo)
%   - COVID-0.13: idem, pero phi = 0.13
%   - COVID-0.15: idem, pero phi = 0.15
%
% Reproduce las figuras paper-ready + un dashboard de %Δ vs Pre-COVID.
% Requiere: solve_two_type_huggett_fiscal.m en el path.

clear; clc; close all;

% ---------- 1) Parámetros base ----------
base = struct();
base.RRA_I = 3.40; base.RRA_F = 3.40;
base.rho   = 0.05;
base.theta = 0.02;                 % prima a<0 (informales)
base.tau_l = 0.15; base.tau_c = 0.18;
base.Gov   = 0.05; base.phi = 0.09;
base.z1    = 0.33; base.z2 = 1.00;
base.p22_bar = 0.8155;
base.I = 700; base.amax = 5.0; base.amin = -0.30*base.z1;
base.r_guess = 0.03; base.rmin = 0.005; base.rmax = 0.08;
base.maxit_V = 100; base.crit_V = 1e-6; base.Delta = 1000;
base.maxit_r = 1000; base.crit_S = 1e-5; base.fix_r = 0;
base.eta_target = 0.654;

% ---------- 2) Shock COVID en ingresos (ASIMÉTRICO) ----------
y_shock_I = 0.65;   % Informales caen 35%
y_shock_F = 0.75;   % Formales caen 25%
z1_covid  = base.z1 * y_shock_I;
z2_covid  = base.z2 * y_shock_F;

% ---------- 3) Políticas de transferencias ----------
phi_pre = base.phi;      % baseline (pre y covid fijo)
phi_013 = 0.13;
phi_015 = 0.15;

% ---------- 4) Configs ----------
cfgPRE = base;                           % Pre-COVID (sin shock)
cfgF   = base; cfgF.z1=z1_covid; cfgF.z2=z2_covid; cfgF.phi=phi_pre;   % COVID fijo
cfgH   = base; cfgH.z1=z1_covid; cfgH.z2=z2_covid; cfgH.phi=phi_013;   % COVID 13%
cfgX   = base; cfgX.z1=z1_covid; cfgX.z2=z2_covid; cfgX.phi=phi_015;   % COVID 15%

outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end
paper_style();

% ---------- 5) Resolver SS ----------
solPRE = solve_two_type_huggett_fiscal(cfgPRE);
solF   = solve_two_type_huggett_fiscal(cfgF);
solH   = solve_two_type_huggett_fiscal(cfgH);
solX   = solve_two_type_huggett_fiscal(cfgX);

P = solPRE; F = solF; H = solH; X = solX;
a = P.a; da = a(2)-a(1);

labP = sprintf('Pre-COVID (phi=%.2f)', phi_pre);
labF = sprintf('COVID fixed (phi=%.2f)', phi_pre);
labH = sprintf('COVID (phi=%.2f)', phi_013);
labX = sprintf('COVID (phi=%.2f)', phi_015);

% ---------- 6) Summary CSV: niveles + %Δ vs Pre ----------
S = table;
S.scenario = ["pre"; "covid_fixed"; "covid_phi013"; "covid_phi015"];
S.y_shock_I  = [1; y_shock_I; y_shock_I; y_shock_I];
S.y_shock_F  = [1; y_shock_F; y_shock_F; y_shock_F];
S.phi        = [phi_pre; phi_pre; phi_013; phi_015];
S.r          = [P.r; F.r; H.r; X.r];
S.Y          = [P.Y; F.Y; H.Y; X.Y];
S.Ctot       = [P.Ctot; F.Ctot; H.Ctot; X.Ctot];
S.B          = [P.fiscal.B; F.fiscal.B; H.fiscal.B; X.fiscal.B];
S.Tl         = [P.fiscal.Tl; F.fiscal.Tl; H.fiscal.Tl; X.fiscal.Tl];
S.Tc         = [P.fiscal.Tc; F.fiscal.Tc; H.fiscal.Tc; X.fiscal.Tc];
S.Tr         = [P.fiscal.Tr; F.fiscal.Tr; H.fiscal.Tr; X.fiscal.Tr];
S.rB         = [P.fiscal.rB; F.fiscal.rB; H.fiscal.rB; X.fiscal.rB];
S.GiniW      = [P.stats.giniW(3); F.stats.giniW(3); H.stats.giniW(3); X.stats.giniW(3)];
S.GiniC      = [P.stats.giniC(3); F.stats.giniC(3); H.stats.giniC(3); X.stats.giniC(3)];
S.fracBorrow_I = [fracB(P,1,da); fracB(F,1,da); fracB(H,1,da); fracB(X,1,da)];
S.fracBorrow_F = [fracB(P,2,da); fracB(F,2,da); fracB(H,2,da); fracB(X,2,da)];
S.debt_I       = [volB(P,1,da);  volB(F,1,da);  volB(H,1,da);  volB(X,1,da)];
S.debt_F       = [volB(P,2,da);  volB(F,2,da);  volB(H,2,da);  volB(X,2,da)];

% Cambios % vs pre (usa la función local pct(x))
S.dY_pct    = pct(S.Y);
S.dC_pct    = pct(S.Ctot);
S.dB_pct    = pct(S.B);
S.drB_pct   = pct(S.rB);

% Graba CSV
writetable(S, fullfile(outdir_tabs,'compare_covid_with_pre_summary.csv'));

% ---------- 7) Políticas: Consumo ----------
fig = figure('Name','Policy: Consumption (Pre vs COVID transfers)');
subplot(1,2,1);
plot(a,P.c(:,1),'k-','LineWidth',2); hold on;
plot(a,F.c(:,1),'LineWidth',2); plot(a,H.c(:,1),'--','LineWidth',2); plot(a,X.c(:,1),'LineWidth',2);
xlabel('Assets a'); ylabel('Consumption c(a)'); title('Informal'); grid on; legend({labP,labF,labH,labX},'Location','best');
subplot(1,2,2);
plot(a,P.c(:,2),'k-','LineWidth',2); hold on;
plot(a,F.c(:,2),'LineWidth',2); plot(a,H.c(:,2),'--','LineWidth',2); plot(a,X.c(:,2),'LineWidth',2);
xlabel('Assets a'); ylabel('Consumption c(a)'); title('Formal'); grid on; legend({labP,labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'covid_pre_policy_consumption'));

% ---------- 8) Políticas: Ahorro ----------
fig = figure('Name','Policy: Savings (Pre vs COVID transfers)');
subplot(1,2,1);
plot(a,P.s(:,1),'k-','LineWidth',2); hold on;
plot(a,F.s(:,1),'LineWidth',2); plot(a,H.s(:,1),'--','LineWidth',2); plot(a,X.s(:,1),'LineWidth',2);
xlabel('Assets a'); ylabel('Savings s(a)'); title('Informal'); grid on; legend({labP,labF,labH,labX},'Location','best');
subplot(1,2,2);
plot(a,P.s(:,2),'k-','LineWidth',2); hold on;
plot(a,F.s(:,2),'LineWidth',2); plot(a,H.s(:,2),'--','LineWidth',2); plot(a,X.s(:,2),'LineWidth',2);
xlabel('Assets a'); ylabel('Savings s(a)'); title('Formal'); grid on; legend({labP,labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'covid_pre_policy_savings'));

% ---------- 9) Densidades de riqueza ----------
fig = figure('Name','Wealth density (Pre vs COVID transfers)');
subplot(2,1,1);
plot(a,P.g(:,1),'k-','LineWidth',2); hold on;
plot(a,F.g(:,1),'LineWidth',2); plot(a,H.g(:,1),'--','LineWidth',2); plot(a,X.g(:,1),'LineWidth',2);
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density'); grid on; legend({labP,labF,labH,labX},'Location','best');
subplot(2,1,2);
plot(a,P.g(:,2),'k-','LineWidth',2); hold on;
plot(a,F.g(:,2),'LineWidth',2); plot(a,H.g(:,2),'--','LineWidth',2); plot(a,X.g(:,2),'LineWidth',2);
xlim([min(a) min(1.0,max(a))]); xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density'); grid on; legend({labP,labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'covid_pre_wealth_distributions'));

% ---------- 10) Cierre del mercado de activos ----------
A_priv = @(S_) sum( (S_.g(:,1)+S_.g(:,2)).*a ) * da;
Mmat = [A_priv(P), A_priv(F), A_priv(H), A_priv(X);
        P.fiscal.B, F.fiscal.B, H.fiscal.B, X.fiscal.B];
fig = figure('Name','Asset market closure (Pre vs COVID transfers)');
bar(Mmat); grid on; ylabel('Assets level');
set(gca,'XTickLabel',{'Private demand','Public supply B(r)'}); 
legend({labP,labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'covid_pre_asset_market'));

% ---------- 11) Cuentas fiscales ----------
fig = figure('Name','Fiscal accounts (Pre vs COVID transfers)');
subplot(1,2,1);
bar([P.fiscal.Tl, F.fiscal.Tl, H.fiscal.Tl, X.fiscal.Tl;
     P.fiscal.Tc, F.fiscal.Tc, H.fiscal.Tc, X.fiscal.Tc]);
set(gca,'XTickLabel',{'Labor tax','VAT'}); ylabel('Revenues'); grid on;
legend({labP,labF,labH,labX},'Location','best'); title('Revenues');
subplot(1,2,2);
bar([P.fiscal.Tr, F.fiscal.Tr, H.fiscal.Tr, X.fiscal.Tr;
     P.fiscal.G,  F.fiscal.G,  H.fiscal.G,  X.fiscal.G;
     P.fiscal.rB, F.fiscal.rB, H.fiscal.rB, X.fiscal.rB]);
set(gca,'XTickLabel',{'Transfers','Public good','Debt service'});
ylabel('Expenditures'); grid on; legend({labP,labF,labH,labX},'Location','best');
title('Expenditures');
export_fig(fig, fullfile(outdir_figs,'covid_pre_fiscal_accounts'));

% ---------- 12) Prestatarios/Prestamistas ----------
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});
fig = figure('Name','Borrowers-Lenders shares (Pre vs COVID transfers)');
subplot(1,2,1);
bar(catX,[P.borrowers.fracBorrow(:), F.borrowers.fracBorrow(:), H.borrowers.fracBorrow(:), X.borrowers.fracBorrow(:)]);
ylabel('Fraction borrowers (a<0)'); title('Borrowers'); grid on; legend({labP,labF,labH,labX},'Location','best');
subplot(1,2,2);
bar(catX,[P.borrowers.fracLend(:),   F.borrowers.fracLend(:),   H.borrowers.fracLend(:),   X.borrowers.fracLend(:)]);
ylabel('Fraction lenders (a>0)');   title('Lenders');   grid on; legend({labP,labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'covid_pre_borrowers_lenders_shares'));

fig = figure('Name','Borrowers-Lenders volumes (Pre vs COVID transfers)');
subplot(1,2,1);
bar(catX,[abs(P.borrowers.volBorrow(:)), abs(F.borrowers.volBorrow(:)), abs(H.borrowers.volBorrow(:)), abs(X.borrowers.volBorrow(:))]);
ylabel('|Aggregate debt|'); title('Debt volume (a<0)'); grid on; legend({labP,labF,labH,labX},'Location','best');
subplot(1,2,2);
bar(catX,[P.borrowers.volLend(:), F.borrowers.volLend(:), H.borrowers.volLend(:), X.borrowers.volLend(:)]);
ylabel('Aggregate savings'); title('Savings volume (a>0)'); grid on; legend({labP,labF,labH,labX},'Location','best');
export_fig(fig, fullfile(outdir_figs,'covid_pre_borrowers_lenders_volumes'));

% ---------- 13) Curvas de Lorenz ----------
[~, LwP, cumPopP] = lorenz_from_density(a, P.g);
[~, LwF, cumPopF] = lorenz_from_density(a, F.g);
[~, LwH, cumPopH] = lorenz_from_density(a, H.g);
[~, LwX, cumPopX] = lorenz_from_density(a, X.g);
fig = figure('Name','Lorenz wealth (Pre vs COVID transfers)');
plot(cumPopP,LwP,'k-','LineWidth',2); hold on;
plot(cumPopF,LwF,'LineWidth',2); plot(cumPopH,LwH,'--','LineWidth',2); plot(cumPopX,LwX,'LineWidth',2);
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Wealth share'); grid on; axis square;
legend({sprintf('%s (Gini=%.3f)',labP,P.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labF,F.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labH,H.stats.giniW(3)), ...
        sprintf('%s (Gini=%.3f)',labX,X.stats.giniW(3))}, 'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'covid_pre_lorenz_wealth'));

[~, LcP, cumPopCP] = lorenz_from_values(a, P.g, P.c(:,1)+P.c(:,2));
[~, LcF, cumPopCF] = lorenz_from_values(a, F.g, F.c(:,1)+F.c(:,2));
[~, LcH, cumPopCH] = lorenz_from_values(a, H.g, H.c(:,1)+H.c(:,2));
[~, LcX, cumPopCX] = lorenz_from_values(a, X.g, X.c(:,1)+X.c(:,2));
fig = figure('Name','Lorenz consumption (Pre vs COVID transfers)');
plot(cumPopCP,LcP,'k-','LineWidth',2); hold on;
plot(cumPopCF,LcF,'LineWidth',2); plot(cumPopCH,LcH,'--','LineWidth',2); plot(cumPopCX,LcX,'LineWidth',2);
plot([0,1],[0,1],'k--'); xlabel('Population share'); ylabel('Consumption share'); grid on; axis square;
legend({sprintf('%s (Gini=%.3f)',labP,P.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labF,F.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labH,H.stats.giniC(3)), ...
        sprintf('%s (Gini=%.3f)',labX,X.stats.giniC(3))}, 'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'covid_pre_lorenz_consumption'));

% ---------- 14) MPC(a) (Pre vs COVID políticas) ----------
DO_MPC = true;
if DO_MPC
    eps_z = 0.01;
    % PRE
    altP = cfgPRE; altP.fix_r=1; altP.r_guess=P.r; altP.rmin=P.r; altP.rmax=P.r;
    altP.z1 = altP.z1*(1+eps_z); altP.z2 = altP.z2*(1+eps_z);
    solPP = solve_two_type_huggett_fiscal(altP);
    dres1P = (cfgPRE.z1*eps_z) * (1 + cfgPRE.phi);
    dres2P = (cfgPRE.z2*eps_z) * (1 - cfgPRE.tau_l);
    MPC1P = (solPP.c(:,1)-P.c(:,1))/max(dres1P,1e-12);
    MPC2P = (solPP.c(:,2)-P.c(:,2))/max(dres2P,1e-12);

    % COVID fijo
    altF = cfgF; altF.fix_r=1; altF.r_guess=F.r; altF.rmin=F.r; altF.rmax=F.r;
    altF.z1 = altF.z1*(1+eps_z); altF.z2 = altF.z2*(1+eps_z);
    solPF = solve_two_type_huggett_fiscal(altF);
    MPC1F = (solPF.c(:,1)-F.c(:,1))/max(dres1P,1e-12);
    MPC2F = (solPF.c(:,2)-F.c(:,2))/max(dres2P,1e-12);

    % COVID 0.13
    altH = cfgH; altH.fix_r=1; altH.r_guess=H.r; altH.rmin=H.r; altH.rmax=H.r;
    altH.z1 = altH.z1*(1+eps_z); altH.z2 = altH.z2*(1+eps_z);
    solPH = solve_two_type_huggett_fiscal(altH);
    MPC1H = (solPH.c(:,1)-H.c(:,1))/max(dres1P,1e-12);
    MPC2H = (solPH.c(:,2)-H.c(:,2))/max(dres2P,1e-12);

    % COVID 0.15
    altX = cfgX; altX.fix_r=1; altX.r_guess=X.r; altX.rmin=X.r; altX.rmax=X.r;
    altX.z1 = altX.z1*(1+eps_z); altX.z2 = altX.z2*(1+eps_z);
    solPX = solve_two_type_huggett_fiscal(altX);
    MPC1X = (solPX.c(:,1)-X.c(:,1))/max(dres1P,1e-12);
    MPC2X = (solPX.c(:,2)-X.c(:,2))/max(dres2P,1e-12);

    fig = figure('Name','MPC by assets (Pre vs COVID transfers)');
    subplot(1,2,1);
    plot(a,MPC1P,'k-','LineWidth',2); hold on;
    plot(a,MPC1F,'LineWidth',2); plot(a,MPC1H,'--','LineWidth',2); plot(a,MPC1X,'LineWidth',2);
    xlabel('Assets a'); ylabel('MPC(a)'); title('Informal'); grid on; legend({labP,labF,labH,labX},'Location','best');
    subplot(1,2,2);
    plot(a,MPC2P,'k-','LineWidth',2); hold on;
    plot(a,MPC2F,'LineWidth',2); plot(a,MPC2H,'--','LineWidth',2); plot(a,MPC2X,'LineWidth',2);
    xlabel('Assets a'); ylabel('MPC(a)'); title('Formal'); grid on; legend({labP,labF,labH,labX},'Location','best');
    export_fig(fig, fullfile(outdir_figs,'covid_pre_mpc_by_assets'));
end

% ---------- 15) Dashboard: %Δ vs Pre-COVID ----------
dR_pct  = pct([S.r(1);  S.r(2);  S.r(3);  S.r(4)]);
dTl_pct = pct([S.Tl(1); S.Tl(2); S.Tl(3); S.Tl(4)]);
dTc_pct = pct([S.Tc(1); S.Tc(2); S.Tc(3); S.Tc(4)]);
dTr_pct = pct([S.Tr(1); S.Tr(2); S.Tr(3); S.Tr(4)]);

valsA = [ S.dY_pct(2:4),  S.dC_pct(2:4),  dR_pct(2:4),  S.dB_pct(2:4),  dTl_pct(2:4) ];
namesA = {'Y','Ctot','r','B','Tl'};

valsB = [ dTc_pct(2:4),  dTr_pct(2:4),  S.drB_pct(2:4), ...
          (S.GiniW(2:4)-S.GiniW(1))*100, (S.GiniC(2:4)-S.GiniC(1))*100 ];
namesB = {'Tc','Tr','rB','GiniW(pp*100)','GiniC(pp*100)'};

fig = figure('Name','Δ vs Pre-COVID (percent changes)');
tiledlayout(2,1);

nexttile;
bar(valsA,'grouped'); grid on; title('Macro & prices');
set(gca,'XTickLabel',namesA);
legend({'COVID fixed','\phi=0.13','\phi=0.15'},'Location','best');

nexttile;
bar(valsB,'grouped'); grid on; title('Fiscal & inequality');
set(gca,'XTickLabel',namesB);
legend({'COVID fixed','\phi=0.13','\phi=0.15'},'Location','best');

export_fig(fig, fullfile(outdir_figs,'covid_pre_percent_changes_dashboard'));

fprintf('All figures written to %s\n', outdir_figs);

% ============================ Helpers (DEBEN IR AL FINAL) ============================
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

function f = fracB(S, j, da)
    a = S.a;
    f = sum(S.g(a<0,j))*da / max(sum(S.g(:,j))*da,1e-12);
end

function v = volB(S, j, da)
    a = S.a;
    v = sum(S.g(a<0,j).*a(a<0))*da;
end

function y = pct(x)
    x = x(:);
    y = 100*(x - x(1)) / max(1e-12, x(1));
end
