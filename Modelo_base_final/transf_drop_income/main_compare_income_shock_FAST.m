% ================= main_compare_income_shock_FAST.m =================
% Escenarios:
% (A) BASE
% (B) Caída de ingresos: z1 -> 0.65 z1 ; z2 -> 0.75 z2
% (C) (B) + +15% en transferencias: phi -> 1.15*phi
%
% Requiere: solve_two_type_huggett_fiscal_Bfixed.m en el path.

clear; clc; close all;

outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

paper_style();

% ---------- PARÁMETROS BASE ----------
base = struct();
% Preferencias
base.RRA_I = 3.40; base.RRA_F = 3.40; base.rho = 0.08;
% Impuestos/transferencias
base.tau_l = 0.15; base.tau_c = 0.18; base.phi = 0.09;
% Ingresos (endowment)
base.z1 = 0.33; base.z2 = 1.00;
% Primas por endeudamiento
base.theta_I = 0.06; base.theta_F = 0.01;
% Objetivo de informalidad y persistencia
base.eta_target = 0.654; base.p22_bar = 0.8155;
% Malla de activos
base.I = 500; base.amax = 3.0; base.amin = -2.0*base.z1;
% Interés y búsqueda
base.r_guess = 0.03; base.rmin = 0.005; base.rmax = 0.08;
base.maxit_V = 120; base.crit_V = 1e-6; base.Delta = 1200;
base.maxit_r = 1500; base.crit_S = 1e-5; base.fix_r = 0;
% Difusión y bien público
base.sigma_a = 0.010; base.psi_G = 0.30; base.omegaG = 1.0;
% Deuda pública y regla de G
base.B_mode = 'ratio_to_Y'; % 'ratio_to_Y' o 'level'
base.Bbar   = 0.35;         % si ratio_to_Y: B = 0.35*Y ; si level: B=0.35
base.alphaG  = 0.50;        % suavizado de G
base.clamp_G_to_zero = true;
base.G_cap_ratio = 0.08;    % CAP opcional G <= 0.08*Y

% ---------- ESCENARIOS ----------
sc = struct([]);
% (A) BASE
sc(1).name = "BASE";
sc(1).cfg  = base;

% (B) Caída de ingresos
sc(2).name = "LOW INCOME";
cfgB = base; cfgB.z1 = 0.65*base.z1; cfgB.z2 = 0.75*base.z2;
sc(2).cfg  = cfgB;

% (C) Caída de ingresos + +15% transferencias
sc(3).name = "LOW INC + 15% TR";
cfgC = cfgB; cfgC.phi = 1.15*base.phi;
sc(3).cfg  = cfgC;

% ---------- RESOLVER ----------
for k=1:numel(sc)
    t0 = tic;
    sc(k).sol = solve_two_type_huggett_fiscal_Bfixed(sc(k).cfg);
    t1 = toc(t0);
    BY = sc(k).sol.fiscal.B / max(sc(k).sol.Y,1e-12);
    fprintf('[%s] r*=%.4f | B=%.3f | Y=%.3f | B/Y=%.3f | G=%.3f | C=%.3f  (%.2fs)\n', ...
        sc(k).name, sc(k).sol.r, sc(k).sol.fiscal.B, sc(k).sol.Y, BY, ...
        sc(k).sol.fiscal.G, sc(k).sol.Ctot, t1);
end

% ---------- TABLA RESUMEN ----------
names      = string(arrayfun(@(q) q.name, sc, 'UniformOutput', false)).';
eta_target = arrayfun(@(q) q.cfg.eta_target, sc).';
r_star     = arrayfun(@(q) q.sol.r, sc).';
B_level    = arrayfun(@(q) q.sol.fiscal.B, sc).';
A_priv     = arrayfun(@(q) sum((q.sol.g(:,1)+q.sol.g(:,2)).*q.sol.a)*(q.sol.a(2)-q.sol.a(1)), sc).';
eta_real   = arrayfun(@(q) q.sol.stats.p11, sc).';   % p11 reportado
Y_out      = arrayfun(@(q) q.sol.Y, sc).';
BY_ratio   = B_level ./ max(Y_out,1e-12);
C_tot      = arrayfun(@(q) q.sol.Ctot, sc).';
G_spend    = arrayfun(@(q) q.sol.fiscal.G, sc).';
PB_bal     = arrayfun(@(q) q.sol.fiscal.PB, sc).';
BB_bal     = arrayfun(@(q) q.sol.fiscal.BB, sc).';
giniW_I    = arrayfun(@(q) q.sol.stats.giniW(1), sc).';
giniW_F    = arrayfun(@(q) q.sol.stats.giniW(2), sc).';
giniW_T    = arrayfun(@(q) q.sol.stats.giniW(3), sc).';
giniC_T    = arrayfun(@(q) q.sol.stats.giniC(3), sc).';

T = table(names, eta_target, r_star, B_level, BY_ratio, A_priv, eta_real, ...
          Y_out, C_tot, G_spend, PB_bal, BB_bal, ...
          giniW_I, giniW_F, giniW_T, giniC_T, ...
          'VariableNames', {'scenario','eta_target','r','B','B_over_Y','A_priv','eta_real','Y','C','G','PB','BB','giniW_I','giniW_F','giniW_T','giniC_T'});
disp(T);
writetable(T, fullfile(outdir_tabs,'compare_income_shock_fast.csv'));

% ---------- GRÁFICOS ----------
% 1) Lorenz (riqueza) con Gini
fig=figure('Name','Lorenz wealth compare');
hold on; grid on; axis square;
for k=1:numel(sc)
    [~,Lw,cPop]=lorenz_from_density(sc(k).sol.a,sc(k).sol.g);
    plot(cPop,Lw,'LineWidth',2,'DisplayName',sprintf('%s (Gini=%.3f)',sc(k).name,sc(k).sol.stats.giniW(3)));
end
plot([0,1],[0,1],'k--','HandleVisibility','off');
xlabel('Population share'); ylabel('Wealth share');
title('Lorenz curve (Wealth)'); legend('Location','southeast');
export_fig(fig, fullfile(outdir_figs,'lorenz_wealth_compare'));

% 2) Cierre de mercado de activos (por escenario)
fig=figure('Name','Asset market closure by scenario');
for k=1:numel(sc)
    subplot(1,numel(sc),k);
    a=sc(k).sol.a; g=sc(k).sol.g; da=a(2)-a(1);
    A_I = sum(g(:,1).*a)*da; A_F = sum(g(:,2).*a)*da; B = sc(k).sol.fiscal.B;
    bar([1 2],[A_I+A_F, B]); grid on;
    set(gca,'XTickLabel',{'Private demand (A_I+A_F)','Public debt B'});
    title(sprintf('%s (r^*=%.4f)', sc(k).name, sc(k).sol.r));
    ylim([min(0,min([A_I+A_F,B])) max([A_I+A_F,B])*1.25]);
end
export_fig(fig, fullfile(outdir_figs,'asset_closure_compare'));

% 3) Borrowers/Lenders: armo datos primero y grafico una sola vez
cats = categorical(string({sc.name}));
% shares totales de prestatarios por escenario (I+F)
borrow_shares = arrayfun(@(q) sum(q.sol.borrowers.fracBorrow(:)), sc).';
% volúmenes de deuda (valor absoluto)
borrow_vols   = arrayfun(@(q) abs(sum(q.sol.borrowers.volBorrow(:))), sc).';

fig=figure('Name','Borrowers-Lenders compare');
subplot(1,2,1); grid on; bar(cats, borrow_shares, 'FaceAlpha',0.9);
ylabel('Share'); title('Borrowers (a<0)');

subplot(1,2,2); grid on; bar(cats, borrow_vols, 'FaceAlpha',0.9);
ylabel('Volume'); title('|Debt| volumes');
export_fig(fig, fullfile(outdir_figs,'borrowers_lenders_compare'));

fprintf('Tablas en %s | Figuras en %s\n', outdir_tabs, outdir_figs);

% ---- Helpers locales ----
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.8); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
end
function [gT,L,cumPop]=lorenz_from_density(a,g)
    da=a(2)-a(1); gT=g(:,1)+g(:,2); W=sum(gT)*da; [as,ix]=sort(a); gs=gT(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as-min(0,min(as))+1e-12; cumWealth=cumsum(wealth_pos.*gs); L=cumWealth/max(cumWealth(end),1e-12);
end
