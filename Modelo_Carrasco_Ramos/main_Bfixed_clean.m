% ============================ main_Bfixed_clean.m =========================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% --------------------- Parámetros base (limpios) --------------------------
cfg = struct();

% Preferencias e impuestos
cfg.RRA_I = 2.40; 
cfg.RRA_F = 3.40; 
cfg.rho   = 0.08;
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;
cfg.phi   = 0.09;

% Ingresos por tipo
cfg.z1 = 0.33;  % informal
cfg.z2 = 1.00;  % formal

% Primas en deuda
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

% Difusión moderada
cfg.sigma_a = 0.007;

% Bien público en utilidad (sin complicar)
cfg.psi_G  = 0.08;
cfg.omegaG = 0.50;
cfg.report_G_effects = 1;

% Gobierno y bonos
cfg.B_mode = 'level';
cfg.Bbar   = 0.25;       % valor provisional; lo ajustamos abajo

paper_style();

% -------------------------------------------------------------------------
% PRE-PASE: elije B dentro del rango de A_priv(rmin) y A_priv(rmax)
% para garantizar cambio de signo en S(r) y evitar r en el borde.
% -------------------------------------------------------------------------
cfg_lo = cfg; cfg_lo.fix_r=1; cfg_lo.r_guess=cfg.rmin; cfg_lo.rmin=cfg.rmin; cfg_lo.rmax=cfg.rmin;
cfg_hi = cfg; cfg_hi.fix_r=1; cfg_hi.r_guess=cfg.rmax; cfg_hi.rmin=cfg.rmax; cfg_hi.rmax=cfg.rmax;

% Usamos Bbar=0 en el pre-pase para no afectar mucho decisiones vía rB/G
cfg_lo.Bbar = 0.0; cfg_hi.Bbar = 0.0;

sol_lo = solve_two_type_huggett_fiscal_Bfixed(cfg_lo);
sol_hi = solve_two_type_huggett_fiscal_Bfixed(cfg_hi);

da = sol_lo.a(2)-sol_lo.a(1);
Apriv_lo = sum(sol_lo.g(:,1).*sol_lo.a)*da + sum(sol_lo.g(:,2).*sol_lo.a)*da;
Apriv_hi = sum(sol_hi.g(:,1).*sol_hi.a)*da + sum(sol_hi.g(:,2).*sol_hi.a)*da;

% Coloca B en el medio del rango (y evita valores degenerados)
B_pre = 0.5*(Apriv_lo + Apriv_hi);
B_pre = max(0.05, min(1.50, B_pre));   % límites razonables
cfg.Bbar = B_pre;

fprintf('Pre-pase: Apriv(rmin)=%.4f, Apriv(rmax)=%.4f => Bbar=%.4f\n', Apriv_lo, Apriv_hi, cfg.Bbar);

% ----------------------------- Resolver equilibrio -----------------------
sol = solve_two_type_huggett_fiscal_Bfixed(cfg);

% ------------------------------- Reportes ---------------------------------
a   = sol.a; g = sol.g; c = sol.c; s = sol.s; r = sol.r;
popI = sol.popI; popF = sol.popF; Y = sol.Y; Ctot = sol.Ctot; fb = sol.fiscal;

fprintf('\n== EQUILIBRIO LIMPIO ==\n');
fprintf('r* = %.6f  | S_excess = %.3e\n', r, sol.S_residual);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y=%.6f  Ctot=%.6f  B=%.6f\n', Y, Ctot, fb.B);
fprintf('Tl=%.6f  Tc=%.6f  Tr=%.6f  G=%.6f  rB=%.6f  -> PB=%.6f  BB=%.6e (debe ser 0)\n', ...
        fb.Tl, fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.BB);
if isfield(sol,'G_effects') && ~isempty(sol.G_effects)
    ge = sol.G_effects;
    fprintf('[G effects] u_mult=(1+psi*Gpc)^omega=%.3f | MU_I=%.3f MU_F=%.3f\n', ...
        ge.u_mult, ge.mu_mult_I, ge.mu_mult_F);
end
fprintf('||HJB residual||_∞ ≈ %.3e\n', sol.hjb_residual);

% ------------------------------- CSVs -------------------------------------
T_fiscal = table("EQUILIBRIO_LIMPIO", fb.Tl, fb.Tc, fb.Tl+fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.B, fb.BB, ...
 'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers','public_good','debt_serv','primary_bal','debt_stock','global_bal'});
writetable(T_fiscal, fullfile(outdir_tabs,'clean_fiscal_breakdown.csv'));

A_I = sum(g(:,1).*a)*da;  A_F = sum(g(:,2).*a)*da;  A_priv = A_I + A_F;
T_assets = table("EQUILIBRIO_LIMPIO", A_I, A_F, A_priv, fb.B, 'VariableNames', ...
    {'scenario','A_I','A_F','A_private','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'clean_asset_market.csv'));

S = sol.stats;
T_stats = table("EQUILIBRIO_LIMPIO", r, popI, popF, Y, Ctot, ...
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
writetable(T_stats, fullfile(outdir_tabs,'clean_household_stats.csv'));

Borr = sol.borrowers;
T_borr = table("EQUILIBRIO_LIMPIO", Borr.fracBorrow(1),Borr.fracBorrow(2),Borr.fracLend(1),Borr.fracLend(2), ...
               Borr.volBorrow(1),Borr.volBorrow(2),Borr.volLend(1),Borr.volLend(2), ...
 'VariableNames', {'scenario','fracBorrow_I','fracBorrow_F','fracLend_I','fracLend_F', ...
                   'volBorrow_I','volBorrow_F','volLend_I','volLend_F'});
writetable(T_borr, fullfile(outdir_tabs,'clean_borrowers_lenders.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

% -------------------------------- Figuras ---------------------------------
% c(a) y s(a)
fig=figure('Name','Policies: Consumption & Savings'); 
subplot(1,2,1);
plot(a,c(:,1),'LineWidth',2); hold on; plot(a,c(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('c(a)'); legend({'Informal','Formal'},'Location','best'); title('Consumption');
subplot(1,2,2);
plot(a,s(:,1),'LineWidth',2); hold on; plot(a,s(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('s(a)'); legend({'Informal','Formal'},'Location','best'); title('Savings');
export_fig(fig, fullfile(outdir_figs,'clean_policies'));

% Densidades (líneas)
fig=figure('Name','Wealth densities');
subplot(2,1,1); plot(a,g(:,1),'LineWidth',2); grid on; xlabel('a'); ylabel('g_1(a)'); title('Informal density');
subplot(2,1,2); plot(a,g(:,2),'LineWidth',2); grid on; xlabel('a'); ylabel('g_2(a)'); title('Formal density');
export_fig(fig, fullfile(outdir_figs,'clean_wealth_densities'));

% Cierre financiero (barras)
fig=figure('Name','Asset market closure (aggregate)');
cats = categorical({'A_{priv}','B (public)'}); cats=reordercats(cats,{'A_{priv}','B (public)'});
bar(cats, [A_priv, fb.B]); ylabel('level'); grid on; 
title(sprintf('Closes at r* = %.4f', r));
export_fig(fig, fullfile(outdir_figs,'clean_asset_market_closure'));

% Curva S(r)
r_span = linspace(cfg.rmin, cfg.rmax, 31);
Sgrid  = nan(size(r_span));
for i=1:numel(r_span)
    cfg_r = cfg; cfg_r.fix_r=1; cfg_r.r_guess=r_span(i); cfg_r.rmin=r_span(i); cfg_r.rmax=r_span(i);
    % usa el Bbar ya calibrado arriba
    solr = solve_two_type_huggett_fiscal_Bfixed(cfg_r);
    Sgrid(i) = solr.S_residual;
end
fig=figure('Name','Excess asset supply S(r)');
plot(r_span,Sgrid,'LineWidth',2); hold on; yline(0,'k--'); xline(r,'r:','LineWidth',1.5);
grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B'); title('Asset market diagnostic');
export_fig(fig, fullfile(outdir_figs,'clean_excess_supply_curve'));

% Fiscal (ingresos/egresos y balances)
fig=figure('Name','Fiscal accounts and balances (clean)');
subplot(1,3,1); hold on; grid on; title('Revenues'); ylabel('amount');
bar(categorical({'Labor','VAT','Total'}), [fb.Tl, fb.Tc, fb.Tl+fb.Tc]);
subplot(1,3,2); hold on; grid on; title('Expenditures'); ylabel('amount');
bar(categorical({'Debt serv','G','Transfers'}), [fb.rB, fb.G, fb.Tr]);
subplot(1,3,3); hold on; grid on; title('Balances'); ylabel('amount');
bar(categorical({'Primary','Overall'}), [fb.PB, fb.BB]); yline(0,'k--');
export_fig(fig, fullfile(outdir_figs,'clean_fiscal_accounts'));

% Lorenz riqueza (total y por tipo)
fig=figure('Name','Lorenz wealth (total)'); 
[~,Lw,cumPop]=lorenz_from_density(a,g);
plot(cumPop,Lw,'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Wealth share');
title(sprintf('Gini_W^T = %.3f', sol.stats.giniW(3)));
export_fig(fig, fullfile(outdir_figs,'clean_lorenz_wealth_total'));

fig=figure('Name','Lorenz wealth by type');
[LwI,cumI] = lorenz_from_assets_single(a, g(:,1));
[LwF,cumF] = lorenz_from_assets_single(a, g(:,2));
plot(cumI,LwI,'LineWidth',2); hold on; plot(cumF,LwF,'LineWidth',2); plot([0,1],[0,1],'k--');
grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
legend({sprintf('Informal (Gini=%.3f)', sol.stats.giniW(1)), ...
        sprintf('Formal (Gini=%.3f)',   sol.stats.giniW(2)), '45°'}, 'Location','southeast');
title('Lorenz curve (Wealth) by type');
export_fig(fig, fullfile(outdir_figs,'clean_lorenz_wealth_bytype'));

% Lorenz consumo (total y por tipo)
fig=figure('Name','Lorenz consumption (total)'); 
[~,Lc,cumPopC]=lorenz_from_values(a,g,c(:,1)+c(:,2));
plot(cumPopC,Lc,'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Gini_C^T = %.3f', sol.stats.giniC(3)));
export_fig(fig, fullfile(outdir_figs,'clean_lorenz_consumption_total'));

fig=figure('Name','Lorenz consumption by type');
[LcI,cI] = lorenz_from_values_single(c(:,1), g(:,1), a);
[LcF,cF] = lorenz_from_values_single(c(:,2), g(:,2), a);
plot(cI, LcI, 'LineWidth',2); hold on; plot(cF, LcF, 'LineWidth',2); plot([0,1],[0,1],'k--');
grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
legend({sprintf('Informal (Gini=%.3f)', sol.stats.giniC(1)), ...
        sprintf('Formal (Gini=%.3f)',   sol.stats.giniC(2)), '45°'}, 'Location','southeast');
title('Lorenz curve (Consumption) by type');
export_fig(fig, fullfile(outdir_figs,'clean_lorenz_consumption_bytype'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ------------------------------- Helpers ----------------------------------
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
