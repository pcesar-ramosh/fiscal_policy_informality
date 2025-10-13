% ================= main_INCOMEshock_TRANSF_compare.m =================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

%% -------------------- Parámetros base (pre-shock) ------------------------
cfg = struct();

% Preferencias e impuestos (tasas constantes)
cfg.RRA_I = 2.40; 
cfg.RRA_F = 3.40; 
cfg.rho   = 0.08;
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;
cfg.phi   = 0.09;   % nivel base de transferencias (por unidad de z1)

% Ingresos BASE (antes de shock, sólo para referencia)
cfg.z1 = 0.33;      % informal
cfg.z2 = 1.00;      % formal

% Spreads (a<0)
cfg.theta_I = 0.06;
cfg.theta_F = 0.01;

% Markov tipos
cfg.eta_target = 0.654; 
cfg.p22_bar    = 0.8155;

% Grid de activos
cfg.I    = 700;
cfg.amin = -2.0*cfg.z1; 
cfg.amax = 3.0;

% Interés
cfg.r_guess = 0.03; 
cfg.rmin    = 0.005; 
cfg.rmax    = 0.10;
cfg.maxit_r = 80;     
cfg.crit_S  = 1e-5;

% HJB
cfg.maxit_V = 160; 
cfg.crit_V  = 1e-6; 
cfg.Delta   = 1400;

% Difusión estabilizadora
cfg.sigma_a = 0.007;

% Bien público en utilidad
cfg.psi_G  = 0.08;
cfg.omegaG = 0.50;
cfg.report_G_effects = 1;

% Gobierno y bonos (nivel)
cfg.B_mode = 'level';
cfg.Bbar   = 0.25;       % será recalibrado vía pre-pase

paper_style();

%% ------------------ Definición del SHOCK de ingresos ---------------------
drop_inf = 0.20;          % -20% z1
drop_for = 0.10;          % -10% z2
z1_shock = (1-drop_inf)*cfg.z1;
z2_shock = (1-drop_for)*cfg.z2;

% Escala de transferencias en el escenario de apoyo
transf_factor = 1.25;     % phi_up = 1.25*phi_base (ajusta si deseas)

%% ---------- Pre-pase para calibrar B con el shock y phi=0 ----------------
cfgNT = cfg;                   % "No Transfers" bajo el shock
cfgNT.z1   = z1_shock;
cfgNT.z2   = z2_shock;
cfgNT.phi  = 0.0;              % SIN transferencias
cfgNT.Bbar = 0.0;              % temporal para pre-pase

cfg_lo = cfgNT; cfg_lo.fix_r=1; cfg_lo.r_guess=cfgNT.rmin; cfg_lo.rmin=cfgNT.rmin; cfg_lo.rmax=cfgNT.rmin;
cfg_hi = cfgNT; cfg_hi.fix_r=1; cfg_hi.r_guess=cfgNT.rmax; cfg_hi.rmin=cfgNT.rmax; cfg_hi.rmax=cfgNT.rmax;

sol_lo = solve_two_type_huggett_fiscal_Bfixed(cfg_lo);
sol_hi = solve_two_type_huggett_fiscal_Bfixed(cfg_hi);

da = sol_lo.a(2)-sol_lo.a(1);
Apriv_lo = sum(sol_lo.g(:,1).*sol_lo.a)*da + sum(sol_lo.g(:,2).*sol_lo.a)*da;
Apriv_hi = sum(sol_hi.g(:,1).*sol_hi.a)*da + sum(sol_hi.g(:,2).*sol_hi.a)*da;

B_pre = max(0.05, min(1.50, 0.5*(Apriv_lo + Apriv_hi)));
cfgNT.Bbar = B_pre;

fprintf('Pre-pase (shock, phi=0): Apriv[rmin]=%.4f, Apriv[rmax]=%.4f => B=%.4f\n', Apriv_lo, Apriv_hi, B_pre);

%% ---------------- Equilibrio 1: SHOCK + SIN transferencias ---------------
NT = solve_two_type_huggett_fiscal_Bfixed(cfgNT);

%% -------------- Equilibrio 2: SHOCK + transferencias ↑ -------------------
cfgTR = cfgNT;
cfgTR.phi  = transf_factor * cfg.phi;  % incremento relativo al phi-base
cfgTR.Bbar = B_pre;                    % mismo B
TR = solve_two_type_huggett_fiscal_Bfixed(cfgTR);

%% --------------------- CSVs (foco fiscal y agregados) --------------------
cats = ["SHOCK_NoTrans"; "SHOCK_TransUp"];
T_fiscal = table(cats, ...
   [NT.fiscal.Tl; TR.fiscal.Tl], [NT.fiscal.Tc; TR.fiscal.Tc], ...
   [NT.fiscal.Tr; TR.fiscal.Tr], [NT.fiscal.G; TR.fiscal.G], ...
   [NT.fiscal.rB; TR.fiscal.rB], [NT.fiscal.PB; TR.fiscal.PB], ...
   [NT.fiscal.B ; TR.fiscal.B ], [NT.fiscal.BB; TR.fiscal.BB], ...
   'VariableNames', {'scenario','Tl','Tc','Tr','G','rB','PB','B','BB'});
writetable(T_fiscal, fullfile(outdir_tabs,'shock_transfers_fiscal.csv'));

T_macro = table(cats, ...
   [NT.r; TR.r], [NT.popI; TR.popI], [NT.popF; TR.popF], ...
   [NT.Y; TR.Y], [NT.Ctot; TR.Ctot], ...
   [NT.stats.giniW(1); TR.stats.giniW(1)], [NT.stats.giniW(2); TR.stats.giniW(2)], [NT.stats.giniW(3); TR.stats.giniW(3)], ...
   [NT.stats.giniC(1); TR.stats.giniC(1)], [NT.stats.giniC(2); TR.stats.giniC(2)], [NT.stats.giniC(3); TR.stats.giniC(3)], ...
   'VariableNames', {'scenario','r','popI','popF','Y','Ctot','giniW_I','giniW_F','giniW_T','giniC_I','giniC_F','giniC_T'});
writetable(T_macro, fullfile(outdir_tabs,'shock_transfers_macro.csv'));

% Activos por tipo
A_I0 = sum(NT.g(:,1).*NT.a)*da;  A_F0 = sum(NT.g(:,2).*NT.a)*da;  A_priv0=A_I0+A_F0;
A_I1 = sum(TR.g(:,1).*TR.a)*da;  A_F1 = sum(TR.g(:,2).*TR.a)*da;  A_priv1=A_I1+A_F1;
T_assets = table(cats, [A_I0;A_I1], [A_F0;A_F1], [A_priv0;A_priv1], [NT.fiscal.B; TR.fiscal.B], ...
   'VariableNames', {'scenario','A_I','A_F','A_priv','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'shock_transfers_assets.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

%% ------------------------------ FIGURAS ----------------------------------
a=NT.a; g0=NT.g; g1=TR.g; c0=NT.c; c1=TR.c; s0=NT.s; s1=TR.s;

% (A) Políticas c(a) y s(a)
fig=figure('Name','Policies under income shock: NoTrans vs TransUp');
subplot(1,2,1);
plot(a,c0(:,1),'b-','LineWidth',2); hold on; plot(a,c0(:,2),'r-','LineWidth',2);
plot(a,c1(:,1),'b--','LineWidth',2); plot(a,c1(:,2),'r--','LineWidth',2);
grid on; xlabel('Assets a'); ylabel('c(a)'); title('Consumption');
legend({'Inf NoTr','For NoTr','Inf Trans↑','For Trans↑'},'Location','best');
subplot(1,2,2);
plot(a,s0(:,1),'b-','LineWidth',2); hold on; plot(a,s0(:,2),'r-','LineWidth',2);
plot(a,s1(:,1),'b--','LineWidth',2); plot(a,s1(:,2),'r--','LineWidth',2);
grid on; xlabel('Assets a'); ylabel('s(a)'); title('Savings');
legend({'Inf NoTr','For NoTr','Inf Trans↑','For Trans↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'shock_transfers_policies'));

% (B) Densidades
fig=figure('Name','Wealth densities (shock): NoTrans vs TransUp');
subplot(2,1,1);
plot(a,g0(:,1),'b-','LineWidth',2); hold on; plot(a,g1(:,1),'b--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_1(a)'); title('Informal'); legend({'NoTr','Trans↑'},'Location','best');
subplot(2,1,2);
plot(a,g0(:,2),'r-','LineWidth',2); hold on; plot(a,g1(:,2),'r--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_2(a)'); title('Formal');   legend({'NoTr','Trans↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'shock_transfers_densities'));

% (C) Agregados Y y C (enfasis ingreso/gasto privado)
fig=figure('Name','Macro aggregates: Y and C');
catsY = categorical({'Y','C_{tot}'}); catsY=reordercats(catsY,{'Y','C_{tot}'});
bar(catsY, [NT.Y NT.Ctot; TR.Y TR.Ctot]);
legend({'NoTr','Trans↑'},'Location','best'); grid on; ylabel('level');
title('Output and total consumption (income shock)');
export_fig(fig, fullfile(outdir_figs,'shock_transfers_Y_C'));

% (D) Cierre financiero (barras lado a lado)
fig=figure('Name','Asset market closure (income shock)');
catsC = categorical({'A_{priv}','B (public)'}); catsC=reordercats(catsC,{'A_{priv}','B (public)'});
bar(catsC, [A_priv0 NT.fiscal.B; A_priv1 TR.fiscal.B]);
legend({'NoTr','Trans↑'},'Location','best'); ylabel('level'); grid on;
title(sprintf('r^* NoTr=%.4f | r^* Trans↑=%.4f', NT.r, TR.r));
export_fig(fig, fullfile(outdir_figs,'shock_transfers_asset_closure'));

% (E) Curva S(r) por escenario (misma B)
r_lo = max(cfg.rmin, max(0.5*min(NT.r,TR.r), 0.002));
r_hi = min(cfg.rmax, max(1.8*max(NT.r,TR.r), max(NT.r,TR.r)+0.02));
if r_hi <= r_lo, r_hi = r_lo*(1+0.05); end
r_span = linspace(r_lo, r_hi, 33);
S0 = nan(size(r_span)); S1 = S0;
for i=1:numel(r_span)
    cf = cfgNT; cf.fix_r=1; cf.r_guess=r_span(i); cf.rmin=r_span(i); cf.rmax=r_span(i);
    sol0 = solve_two_type_huggett_fiscal_Bfixed(cf); S0(i)=sol0.S_residual;

    ct = cfgTR; ct.fix_r=1; ct.r_guess=r_span(i); ct.rmin=r_span(i); ct.rmax=r_span(i);
    sol1 = solve_two_type_huggett_fiscal_Bfixed(ct); S1(i)=sol1.S_residual;
end
fig=figure('Name','Excess asset supply S(r): NoTrans vs TransUp');
plot(r_span,S0,'b-','LineWidth',2); hold on; plot(r_span,S1,'r--','LineWidth',2);
yline(0,'k--'); xline(NT.r,':b','LineWidth',1.2); xline(TR.r,':r','LineWidth',1.2);
grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
legend({'NoTr','Trans↑','S=0','r^* NoTr','r^* Trans↑'},'Location','best');
title('Asset market clearing (income shock)');
export_fig(fig, fullfile(outdir_figs,'shock_transfers_Scurve'));

% (F) Cuentas fiscales: ingresos, gasto y balances (énfasis)
fig=figure('Name','Fiscal accounts (income shock)');
subplot(1,3,1); hold on; grid on; title('Revenues'); ylabel('amount');
bar(categorical({'Labor','VAT','Total'}), ...
    [NT.fiscal.Tl NT.fiscal.Tc NT.fiscal.Tl+NT.fiscal.Tc; ...
     TR.fiscal.Tl TR.fiscal.Tc TR.fiscal.Tl+TR.fiscal.Tc]);
legend({'NoTr','Trans↑'},'Location','best');
subplot(1,3,2); hold on; grid on; title('Expenditures'); ylabel('amount');
bar(categorical({'Debt serv','G','Transfers'}), ...
    [NT.fiscal.rB NT.fiscal.G NT.fiscal.Tr; ...
     TR.fiscal.rB TR.fiscal.G TR.fiscal.Tr]);
legend({'NoTr','Trans↑'},'Location','best');
subplot(1,3,3); hold on; grid on; title('Balances'); ylabel('amount');
bar(categorical({'Primary','Overall'}), ...
    [NT.fiscal.PB NT.fiscal.BB; TR.fiscal.PB TR.fiscal.BB]); yline(0,'k--');
legend({'NoTr','Trans↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'shock_transfers_fiscal'));

% (G) Prestatarios/Prestamistas (shares y volúmenes)
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
fig=figure('Name','Borrowers/Lenders (income shock)');
subplot(1,2,1);
bar(catX, [NT.borrowers.fracBorrow(:) NT.borrowers.fracLend(:); ...
           TR.borrowers.fracBorrow(:) TR.borrowers.fracLend(:)]);
legend({'Borrowers (a<0)','Lenders (a>0)'},'Location','best'); grid on; ylabel('Share'); title('Shares by type');
subplot(1,2,2);
bar(catX, [abs(NT.borrowers.volBorrow(:)) NT.borrowers.volLend(:); ...
           abs(TR.borrowers.volBorrow(:)) TR.borrowers.volLend(:)]);
legend({'|Debt|','Savings'},'Location','best'); grid on; ylabel('Volume'); title('Volumes by type');
export_fig(fig, fullfile(outdir_figs,'shock_transfers_borrowers'));

% (H) Lorenz riqueza por tipo
fig=figure('Name','Lorenz wealth by type (income shock)');
[LwI0,cI0]=lorenz_from_assets_single(a, g0(:,1));
[LwF0,cF0]=lorenz_from_assets_single(a, g0(:,2));
[LwI1,cI1]=lorenz_from_assets_single(a, g1(:,1));
[LwF1,cF1]=lorenz_from_assets_single(a, g1(:,2));
subplot(1,2,1); plot(cI0,LwI0,'b-','LineWidth',2); hold on; plot(cI1,LwI1,'b--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Informal: Gini NoTr=%.3f | Trans↑=%.3f', NT.stats.giniW(1), TR.stats.giniW(1)));
subplot(1,2,2); plot(cF0,LwF0,'r-','LineWidth',2); hold on; plot(cF1,LwF1,'r--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Formal: Gini NoTr=%.3f | Trans↑=%.3f', NT.stats.giniW(2), TR.stats.giniW(2)));
export_fig(fig, fullfile(outdir_figs,'shock_transfers_lorenz_wealth'));

% (I) Lorenz consumo por tipo
fig=figure('Name','Lorenz consumption by type (income shock)');
[LcI0,cIp0]=lorenz_from_values_single(c0(:,1), g0(:,1), a);
[LcF0,cFp0]=lorenz_from_values_single(c0(:,2), g0(:,2), a);
[LcI1,cIp1]=lorenz_from_values_single(c1(:,1), g1(:,1), a);
[LcF1,cFp1]=lorenz_from_values_single(c1(:,2), g1(:,2), a);
subplot(1,2,1); plot(cIp0,LcI0,'b-','LineWidth',2); hold on; plot(cIp1,LcI1,'b--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Informal: GiniC NoTr=%.3f | Trans↑=%.3f', NT.stats.giniC(1), TR.stats.giniC(1)));
subplot(1,2,2); plot(cFp0,LcF0,'r-','LineWidth',2); hold on; plot(cFp1,LcF1,'r--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Formal: GiniC NoTr=%.3f | Trans↑=%.3f', NT.stats.giniC(2), TR.stats.giniC(2)));
export_fig(fig, fullfile(outdir_figs,'shock_transfers_lorenz_consumption'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

%% ------------------------------- HELPERS ---------------------------------
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
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
% =======================================================================
