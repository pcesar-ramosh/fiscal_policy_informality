% ====================== main_TRANSF_up_compare_Bfixed.m =====================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

%% --------------------- Parámetros BASE (limpios) --------------------------
cfg = struct();
% Preferencias e impuestos
cfg.RRA_I = 2.40; 
cfg.RRA_F = 3.40; 
cfg.rho   = 0.08;
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;
cfg.phi   = 0.09;     % <- transfer per-cápita para informales (lump-sum, proporcional a z1)

% Ingresos por tipo
cfg.z1 = 0.33;  % informal
cfg.z2 = 1.00;  % formal

% Spreads (a<0)
cfg.theta_I = 0.06;
cfg.theta_F = 0.01;

% Markov de tipos
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

% Difusión
cfg.sigma_a = 0.007;

% Bien público en utilidad
cfg.psi_G  = 0.08;
cfg.omegaG = 0.50;
cfg.report_G_effects = 1;

% Gobierno y bonos (calibraremos B con un pre-pase)
cfg.B_mode = 'level';
cfg.Bbar   = 0.25;

paper_style();

%% ----------- PRE-PASE: fija B para que S(r) cambie de signo ---------------
cfg_lo = cfg; cfg_lo.fix_r=1; cfg_lo.r_guess=cfg.rmin; cfg_lo.rmin=cfg.rmin; cfg_lo.rmax=cfg.rmin; cfg_lo.Bbar=0.0;
cfg_hi = cfg; cfg_hi.fix_r=1; cfg_hi.r_guess=cfg.rmax; cfg_hi.rmin=cfg.rmax; cfg_hi.rmax=cfg.rmax; cfg_hi.Bbar=0.0;

sol_lo = solve_two_type_huggett_fiscal_Bfixed(cfg_lo);
sol_hi = solve_two_type_huggett_fiscal_Bfixed(cfg_hi);

da = sol_lo.a(2)-sol_lo.a(1);
Apriv_lo = sum(sol_lo.g(:,1).*sol_lo.a)*da + sum(sol_lo.g(:,2).*sol_lo.a)*da;
Apriv_hi = sum(sol_hi.g(:,1).*sol_hi.a)*da + sum(sol_hi.g(:,2).*sol_hi.a)*da;

B_pre = max(0.05, min(1.50, 0.5*(Apriv_lo + Apriv_hi)));
cfg.Bbar = B_pre;

fprintf('Pre-pase BASE: Apriv(rmin)=%.4f, Apriv(rmax)=%.4f => Bbar=%.4f\n', Apriv_lo, Apriv_hi, cfg.Bbar);

%% --------------------------- Equilibrio BASE ------------------------------
base = solve_two_type_huggett_fiscal_Bfixed(cfg);

% Guardar B* para mantenerlo fijo en el comparativo
B_fixed = base.fiscal.B;

%% --------------------- Escenario: TRANSFERENCIAS ↑ ------------------------
phi_up = 1.15 * cfg.phi;    % +15% (ajústalo si quieres otro delta)
cfgT = cfg;
cfgT.phi    = phi_up;       % más transferencias a informales
cfgT.Bbar   = B_fixed;      % MISMO stock de deuda en nivel
cfgT.B_mode = 'level';

trup = solve_two_type_huggett_fiscal_Bfixed(cfgT);

%% ---------------------------- Consolidados CSV ---------------------------
catsc = ["BASE";"TRANSF_UP"];
T_fiscal = table(catsc, ...
    [base.fiscal.Tl; trup.fiscal.Tl], ...
    [base.fiscal.Tc; trup.fiscal.Tc], ...
    [base.fiscal.Tl+base.fiscal.Tc; trup.fiscal.Tl+trup.fiscal.Tc], ...
    [base.fiscal.Tr; trup.fiscal.Tr], ...
    [base.fiscal.G;  trup.fiscal.G], ...
    [base.fiscal.rB; trup.fiscal.rB], ...
    [base.fiscal.PB; trup.fiscal.PB], ...
    [base.fiscal.B;  trup.fiscal.B], ...
    [base.fiscal.BB; trup.fiscal.BB], ...
 'VariableNames', {'scenario','Tl','Tc','Rev_total','Transfers','G','rB','Primary','B_stock','Global'});
writetable(T_fiscal, fullfile(outdir_tabs,'compare_transfers_fiscal.csv'));

T_stats = table(catsc, ...
   [base.r; trup.r], ...
   [base.popI; trup.popI], [base.popF; trup.popF], ...
   [base.Y; trup.Y], [base.Ctot; trup.Ctot], ...
   [base.stats.giniW(1); trup.stats.giniW(1)], [base.stats.giniW(2); trup.stats.giniW(2)], [base.stats.giniW(3); trup.stats.giniW(3)], ...
   [base.stats.giniC(1); trup.stats.giniC(1)], [base.stats.giniC(2); trup.stats.giniC(2)], [base.stats.giniC(3); trup.stats.giniC(3)], ...
   [base.S_residual; trup.S_residual], ...
   'VariableNames', {'scenario','r','popI','popF','Y','Ctot','giniW_I','giniW_F','giniW_T','giniC_I','giniC_F','giniC_T','S_excess'});
writetable(T_stats, fullfile(outdir_tabs,'compare_transfers_stats.csv'));

% Ddescomposición A por tipo
da = base.a(2)-base.a(1);
A_I0 = sum(base.g(:,1).*base.a)*da;  A_F0 = sum(base.g(:,2).*base.a)*da;  A_priv0=A_I0+A_F0;
A_I1 = sum(trup.g(:,1).*trup.a)*da;  A_F1 = sum(trup.g(:,2).*trup.a)*da;  A_priv1=A_I1+A_F1;

T_assets = table(catsc, [A_I0;A_I1], [A_F0;A_F1], [A_priv0;A_priv1], [base.fiscal.B; trup.fiscal.B], ...
   'VariableNames', {'scenario','A_I','A_F','A_priv','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'compare_transfers_assets.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

%% ------------------------------ FIGURAS ----------------------------------
% Atajos
a=base.a; gb=base.g; cb=base.c; sb=base.s; rt=trup.r; rb=base.r;
gt=trup.g; ct=trup.c; st=trup.s;

% 1) Políticas c(a) y s(a) (mismas figuras, overlay)
fig=figure('Name','c(a): BASE vs TRANS↑ (tax adj)');
subplot(1,2,1);
plot(a,cb(:,1),'b-','LineWidth',2); hold on; plot(a,cb(:,2),'r-','LineWidth',2);
plot(a,ct(:,1),'b--','LineWidth',2); plot(a,ct(:,2),'r--','LineWidth',2);
grid on; xlabel('a'); ylabel('c(a)');
legend({'Inf BASE','For BASE','Inf TR↑','For TR↑'},'Location','best'); title('Consumption');
subplot(1,2,2);
plot(a,sb(:,1),'b-','LineWidth',2); hold on; plot(a,sb(:,2),'r-','LineWidth',2);
plot(a,st(:,1),'b--','LineWidth',2); plot(a,st(:,2),'r--','LineWidth',2);
grid on; xlabel('a'); ylabel('s(a)');
legend({'Inf BASE','For BASE','Inf TR↑','For TR↑'},'Location','best'); title('Savings');
export_fig(fig, fullfile(outdir_figs,'transf_compare_policies'));

% 2) Densidades (líneas overlay)
fig=figure('Name','Wealth densities: BASE vs TRANS↑');
subplot(2,1,1);
plot(a,gb(:,1),'b-','LineWidth',2); hold on; plot(a,gt(:,1),'b--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_1(a)'); title('Informal density');
legend({'BASE','TRANS↑'},'Location','best');
subplot(2,1,2);
plot(a,gb(:,2),'r-','LineWidth',2); hold on; plot(a,gt(:,2),'r--','LineWidth',2);
grid on; xlabel('a'); ylabel('g_2(a)'); title('Formal density');
legend({'BASE','TRANS↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'transf_compare_densities'));

% 3) Cierre financiero (aggregate, barras lado a lado)
fig=figure('Name','Asset market closure (aggregate, compare)');
cats = categorical({'A_{priv}','B (public)'}); cats=reordercats(cats,{'A_{priv}','B (public)'});
bar(cats, [A_priv0 base.fiscal.B; A_priv1 trup.fiscal.B]); 
legend({'BASE','TRANS↑'},'Location','best'); ylabel('level'); grid on;
title(sprintf('Closes at r^* base=%.4f | r^* transf=%.4f', base.r, trup.r));
export_fig(fig, fullfile(outdir_figs,'transf_compare_asset_closure'));

% 4) Curva S(r) por escenario (overlay)
r_lo = max(cfg.rmin, max(0.5*min(base.r,trup.r), 0.002));
r_hi = min(cfg.rmax, max(1.8*max(base.r,trup.r), max(base.r,trup.r)+0.02));
if r_hi <= r_lo, r_hi = r_lo*(1+0.05); end
r_span = linspace(r_lo, r_hi, 33);

S0 = nan(size(r_span)); S1 = S0;
for i=1:numel(r_span)
    % BASE con B fijo
    cfgr = cfg;  cfgr.fix_r=1; cfgr.r_guess=r_span(i); cfgr.rmin=r_span(i); cfgr.rmax=r_span(i);
    cfgr.Bbar = B_fixed; cfgr.B_mode='level';
    solr0 = solve_two_type_huggett_fiscal_Bfixed(cfgr); S0(i)=solr0.S_residual;

    % TRANS↑ con B fijo y phi_up
    cfrt = cfgr; cfrt.phi = phi_up;
    solr1 = solve_two_type_huggett_fiscal_Bfixed(cfrt); S1(i)=solr1.S_residual;
end
fig=figure('Name','Excess asset supply by scenario');
plot(r_span,S0,'b-','LineWidth',2); hold on;
plot(r_span,S1,'r--','LineWidth',2);
yline(0,'k--'); xline(base.r,':b','LineWidth',1.2); xline(trup.r,':r','LineWidth',1.2);
grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
legend({'BASE','TRANS↑','S=0','r^* BASE','r^* TRANS↑'},'Location','best');
title('Excess asset supply by scenario');
export_fig(fig, fullfile(outdir_figs,'transf_compare_Scurve'));

% 5) Cuentas fiscales (tres paneles, overlay)
fig=figure('Name','Fiscal accounts: BASE vs TRANS↑');
subplot(1,3,1); hold on; grid on; title('Revenues'); ylabel('amount');
bar(categorical({'Labor','VAT','Total'}), [base.fiscal.Tl base.fiscal.Tc base.fiscal.Tl+base.fiscal.Tc; ...
                                           trup.fiscal.Tl trup.fiscal.Tc trup.fiscal.Tl+trup.fiscal.Tc]);
legend({'BASE','TRANS↑'},'Location','best');
subplot(1,3,2); hold on; grid on; title('Expenditures'); ylabel('amount');
bar(categorical({'Debt serv','G','Transfers'}), [base.fiscal.rB base.fiscal.G base.fiscal.Tr; ...
                                                trup.fiscal.rB trup.fiscal.G trup.fiscal.Tr]);
legend({'BASE','TRANS↑'},'Location','best');
subplot(1,3,3); hold on; grid on; title('Balances'); ylabel('amount');
bar(categorical({'Primary','Overall'}), [base.fiscal.PB base.fiscal.BB; ...
                                         trup.fiscal.PB  trup.fiscal.BB]); yline(0,'k--');
legend({'BASE','TRANS↑'},'Location','best');
export_fig(fig, fullfile(outdir_figs,'transf_compare_fiscal'));

% 6) Prestatarios / Prestamistas (shares y volúmenes)
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
fig=figure('Name','Borrowers-Lenders: BASE vs TRANS↑');
subplot(1,2,1); 
bar(catX, [base.borrowers.fracBorrow(:) base.borrowers.fracLend(:); ...
           trup.borrowers.fracBorrow(:) trup.borrowers.fracLend(:)]);
legend({'Borrowers (a<0)','Lenders (a>0)'},'Location','best'); grid on; ylabel('Share');
title('Shares by type'); 
subplot(1,2,2);
bar(catX, [abs(base.borrowers.volBorrow(:)) base.borrowers.volLend(:); ...
           abs(trup.borrowers.volBorrow(:))  trup.borrowers.volLend(:)]);
legend({'|Debt|','Savings'},'Location','best'); grid on; ylabel('Volume');
title('Volumes by type');
export_fig(fig, fullfile(outdir_figs,'transf_compare_borrowers'));

% 7) Lorenz riqueza por tipo (overlay)
fig=figure('Name','Lorenz wealth by type: BASE vs TRANS↑');
[LwI0,cI0]=lorenz_from_assets_single(a, gb(:,1));
[LwF0,cF0]=lorenz_from_assets_single(a, gb(:,2));
[LwI1,cI1]=lorenz_from_assets_single(a, gt(:,1));
[LwF1,cF1]=lorenz_from_assets_single(a, gt(:,2));
subplot(1,2,1);
plot(cI0,LwI0,'b-','LineWidth',2); hold on; plot(cI1,LwI1,'b--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Informal: Gini BASE=%.3f | TR↑=%.3f', base.stats.giniW(1), trup.stats.giniW(1)));
subplot(1,2,2);
plot(cF0,LwF0,'r-','LineWidth',2); hold on; plot(cF1,LwF1,'r--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Formal: Gini BASE=%.3f | TR↑=%.3f', base.stats.giniW(2), trup.stats.giniW(2)));
export_fig(fig, fullfile(outdir_figs,'transf_compare_lorenz_wealth'));

% 8) Lorenz consumo por tipo (overlay)
fig=figure('Name','Lorenz consumption by type: BASE vs TRANS↑');
[LcI0,cIp0]=lorenz_from_values_single(cb(:,1), gb(:,1), a);
[LcF0,cFp0]=lorenz_from_values_single(cb(:,2), gb(:,2), a);
[LcI1,cIp1]=lorenz_from_values_single(ct(:,1), gt(:,1), a);
[LcF1,cFp1]=lorenz_from_values_single(ct(:,2), gt(:,2), a);
subplot(1,2,1);
plot(cIp0,LcI0,'b-','LineWidth',2); hold on; plot(cIp1,LcI1,'b--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Informal: GiniC BASE=%.3f | TR↑=%.3f', base.stats.giniC(1), trup.stats.giniC(1)));
subplot(1,2,2);
plot(cFp0,LcF0,'r-','LineWidth',2); hold on; plot(cFp1,LcF1,'r--','LineWidth',2); plot([0,1],[0,1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Formal: GiniC BASE=%.3f | TR↑=%.3f', base.stats.giniC(2), trup.stats.giniC(2)));
export_fig(fig, fullfile(outdir_figs,'transf_compare_lorenz_cons'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

%% ------------------------------- Helpers ----------------------------------
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
% ==========================================================================
