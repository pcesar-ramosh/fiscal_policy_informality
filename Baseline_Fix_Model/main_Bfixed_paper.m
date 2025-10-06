% ============================ main_Bfixed_paper.m ============================
clear; clc; close all;
outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% -------- Parámetros pensados para que r ~ 3% y haya intersección --------
cfg = struct();
cfg.RRA_I = 3.40; cfg.RRA_F = 3.40;
cfg.rho   = 0.08;                     % más impacientes => menos A_priv
cfg.tau_l = 0.15; cfg.tau_c = 0.18; cfg.phi = 0.09;
cfg.z1 = 0.33; cfg.z2 = 1.00;

cfg.theta_I = 0.06;                   % prima alta informal
cfg.theta_F = 0.01;                   % prima baja formal

cfg.eta_target = 0.654; cfg.p22_bar = 0.8155;

cfg.I = 700; cfg.amax = 3.0; cfg.amin = -2.0*cfg.z1; % (reduce A_priv)

cfg.r_guess = 0.03; cfg.rmin = 0.005; cfg.rmax = 0.10;
cfg.maxit_V = 140; cfg.crit_V = 1e-6; cfg.Delta = 1400;
cfg.maxit_r = 1500; cfg.crit_S = 1e-5; cfg.fix_r = 0;

cfg.sigma_a = 0.010;       % suaviza FP
cfg.psi_G   = 0.30;        % G multiplicativo (intensidad)
cfg.omegaG  = 1.0;
cfg.B_mode  = 'ratio_to_Y';
cfg.Bbar    = 0.45;        % Perú-like deuda/PIB
cfg.alphaG  = 0.50; cfg.clamp_G_to_zero = true;
cfg.G_cap_ratio = 0.08;    % CAP opcional para que G no sea descomunal

% ---------------------------- Resolver modelo -----------------------------
base = solve_two_type_huggett_fiscal_Bfixed(cfg);

a   = base.a; g = base.g; c = base.c; s = base.s; r = base.r;
popI = base.popI; popF = base.popF; Y = base.Y; Ctot = base.Ctot; Gpc = base.Gpc;
fb = base.fiscal; S_res = base.S_residual; hjb_res = base.hjb_residual;

fprintf('\n== BASE (B fixed, G multiplicativo) ==\n');
fprintf('r = %.4f  | S_excess = %.3e\n', r, S_res);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y = %.6f,   Ctot = %.6f,  Gpc = %.6f\n', Y, Ctot, Gpc);
fprintf('Tl=%.6f  Tc=%.6f  Tr=%.6f  G=%.6f  rB=%.6f  -> PB=%.6f  BB=%.6e\n', ...
        fb.Tl, fb.Tc, fb.Tl+fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.BB);
fprintf('||HJB residual||_∞ ≈ %.3e\n', hjb_res);

% ------------------------------- CSVs -------------------------------------
T_fiscal = table("BASE_Bfixed", fb.Tl, fb.Tc, fb.Tl+fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.B, fb.BB, ...
    'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers','public_good','debt_serv','primary_bal','debt_stock','global_bal'});
writetable(T_fiscal, fullfile(outdir_tabs,'base_fiscal_breakdown.csv'));

da = a(2)-a(1);
A_I = sum(g(:,1).*a)*da;  % <-- desagregado por tipo
A_F = sum(g(:,2).*a)*da;
A_priv = A_I + A_F;

T_assets = table("BASE_Bfixed", A_priv, fb.B, 'VariableNames', {'scenario','A_private','B_public'});
writetable(T_assets, fullfile(outdir_tabs,'base_asset_market.csv'));

S = base.stats;
T_stats = table("BASE_Bfixed", r, popI, popF, Y, Ctot, ...
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
writetable(T_stats, fullfile(outdir_tabs,'base_household_stats.csv'));

Borr = base.borrowers;
T_borr = table("BASE_Bfixed", Borr.fracBorrow(1),Borr.fracBorrow(2),Borr.fracLend(1),Borr.fracLend(2), ...
               Borr.volBorrow(1),Borr.volBorrow(2),Borr.volLend(1),Borr.volLend(2), ...
               'VariableNames', {'scenario','fracBorrow_I','fracBorrow_F','fracLend_I','fracLend_F', ...
                                 'volBorrow_I','volBorrow_F','volLend_I','volLend_F'});
writetable(T_borr, fullfile(outdir_tabs,'base_borrowers_lenders.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

% ------------------------------- Figuras ----------------------------------
paper_style();

% c(a)
fig=figure('Name','Policy: Consumption'); plot(a,c(:,1),'LineWidth',2); hold on; plot(a,c(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('Consumption c(a)'); legend({'Informal','Formal'},'Location','best'); title('Consumption Policies');
export_fig(fig, fullfile(outdir_figs,'policy_consumption'));

% s(a)
fig=figure('Name','Policy: Savings'); plot(a,s(:,1),'LineWidth',2); hold on; plot(a,s(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('Savings s(a)'); legend({'Informal','Formal'},'Location','best'); title('Savings Policies');
export_fig(fig, fullfile(outdir_figs,'policy_savings'));

% Densidades de riqueza por tipo
fig=figure('Name','Wealth distributions');
subplot(2,1,1); bar(a,g(:,1),'FaceAlpha',0.9,'EdgeColor','none'); grid on; xlim([min(a) min(1.2,max(a))]);
xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density');
subplot(2,1,2); bar(a,g(:,2),'FaceAlpha',0.9,'EdgeColor','none'); grid on; xlim([min(a) min(5.0,max(a))]);
xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density');
export_fig(fig, fullfile(outdir_figs,'wealth_distributions'));

% Cierre financiero (barras simples)
fig=figure('Name','Asset market closure (aggregate)');
bar(categorical({'Private demand','Public debt B'}), [A_priv, fb.B]);
ylabel('Assets level'); grid on; title(sprintf('Asset market closes at r = %.4f', r));
export_fig(fig, fullfile(outdir_figs,'asset_market_closure'));

% === NUEVO: Cierre financiero desagregado (stacked vs B) ===
fig=figure('Name','Asset market closure by type (stacked)');
X = categorical({'Private demand','Public debt B'}); X=reordercats(X,{'Private demand','Public debt B'});
Y = [A_I A_F 0; 0 0 fb.B]; % stacked: (I,F) para privados; B en barra aparte
bar(X, Y, 'stacked'); grid on;
ylabel('Assets level');
legend({'A_I (informal)','A_F (formal)','B (public)'},'Location','best');
title(sprintf('Composition at r^* = %.4f', r));
export_fig(fig, fullfile(outdir_figs,'asset_market_closure_by_type'));

% Curva S(r) y A_priv(r) vs B
r_span = linspace(max(0.005,cfg.r_guess*0.4), cfg.r_guess*2.2, 31);
Sgrid  = nan(size(r_span)); Apriv_grid = nan(size(r_span));
Apriv_I_grid = nan(size(r_span)); Apriv_F_grid = nan(size(r_span));  % <-- NUEVO
alt = cfg; alt.fix_r = 1; alt.maxit_r = 1;
for i=1:numel(r_span)
    alt.r_guess = r_span(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
    tmp = solve_two_type_huggett_fiscal_Bfixed(alt);
    da_i = tmp.a(2)-tmp.a(1);
    Sgrid(i) = tmp.S_residual;
    Apriv_I_grid(i) = sum(tmp.g(:,1).*tmp.a)*da_i;
    Apriv_F_grid(i) = sum(tmp.g(:,2).*tmp.a)*da_i;
    Apriv_grid(i)   = Apriv_I_grid(i) + Apriv_F_grid(i);
end

fig=figure('Name','Excess asset supply S(r)');
plot(r_span,Sgrid,'LineWidth',2); yline(0,'k--');
grid on; xlabel('Interest rate r'); ylabel('Excess assets S(r)'); title('Asset market diagnostic (B fixed)'); xline(r,'r:','LineWidth',1.5);
export_fig(fig, fullfile(outdir_figs,'excess_supply_curve'));

fig=figure('Name','A_priv(r) vs B (by type)');
plot(r_span,Apriv_I_grid,'LineWidth',2); hold on;
plot(r_span,Apriv_F_grid,'LineWidth',2);
plot(r_span,Apriv_grid,'LineWidth',2);
yline(fb.B,'k--','LineWidth',1.5); xline(r,'r:','LineWidth',1.5);
grid on; xlabel('Interest rate r'); ylabel('Assets level');
legend({'A_I(r)','A_F(r)','A_{priv}(r)','B (supply)','r^*'},'Location','best');
title('Private demand for assets by type vs Public supply');
export_fig(fig, fullfile(outdir_figs,'asset_demand_vs_B_by_type'));

%=== Graf Financial Market ======

% --------- Curvas de mercado: A_I(r), A_F(r), A_priv(r) y B(r) ----------
r_span = linspace(max(0.005,cfg.r_guess*0.4), cfg.r_guess*2.2, 41);

Sgrid  = nan(size(r_span));
Apriv_grid = nan(size(r_span));
Apriv_I_grid = nan(size(r_span));
Apriv_F_grid = nan(size(r_span));
B_grid = nan(size(r_span));
Y_grid = nan(size(r_span));

alt = cfg; alt.fix_r = 1; alt.maxit_r = 1;

for i=1:numel(r_span)
    alt.r_guess = r_span(i); alt.rmin = alt.r_guess; alt.rmax = alt.r_guess;
    tmp = solve_two_type_huggett_fiscal_Bfixed(alt);

    da_i = tmp.a(2)-tmp.a(1);
    Ai = sum(tmp.g(:,1).*tmp.a)*da_i;
    Af = sum(tmp.g(:,2).*tmp.a)*da_i;
    Apriv_I_grid(i) = Ai;
    Apriv_F_grid(i) = Af;
    Apriv_grid(i)   = Ai + Af;

    % OJO: en 'ratio_to_Y' el B de cada solución ya depende de Y(r)
    Y_grid(i) = tmp.Y;
    B_grid(i) = tmp.fiscal.B;

    Sgrid(i) = Apriv_grid(i) - B_grid(i);
end

% Hallar (aprox) el r* del grid (cero de S)
r_star = base.r; A_star = A_priv; B_star = fb.B;  % defaults del equilibrio base
ix = find(diff(sign(Sgrid))~=0,1,'first');
if ~isempty(ix)
    % Interpolación lineal simple en el segmento [ix, ix+1]
    r1=r_span(ix); r2=r_span(ix+1); s1=Sgrid(ix); s2=Sgrid(ix+1);
    r_star = r1 - s1*(r2-r1)/max(s2-s1,1e-12);
    % Interpolar niveles de A_priv y B en r_star
    A_star = interp1(r_span, Apriv_grid, r_star, 'linear','extrap');
    B_star = interp1(r_span, B_grid,     r_star, 'linear','extrap');
end

% --- Gráfico 1: Exceso S(r) con cruce en cero ---
fig=figure('Name','Excess asset supply S(r)');
plot(r_span,Sgrid,'LineWidth',2); hold on;
yline(0,'k--'); xline(r_star,'r:','LineWidth',1.5);
grid on; xlabel('Interest rate r'); ylabel('Excess assets S(r)=A_{priv}(r)-B(r)');
title('Asset market diagnostic (B depends on Y when ratio\_to\_Y)');
export_fig(fig, fullfile(outdir_figs,'excess_supply_curve'));

% --- Gráfico 2: Demanda privada vs Oferta pública, con cruce ---
fig=figure('Name','A_priv(r) & B(r): Market equilibrium');
plot(r_span,Apriv_I_grid,'LineWidth',2); hold on;
plot(r_span,Apriv_F_grid,'LineWidth',2);
plot(r_span,Apriv_grid,'LineWidth',2);
plot(r_span,B_grid,'k--','LineWidth',1.8);                     % B(r) endógeno
xline(r_star,'r:','LineWidth',1.5);
plot(r_star,A_star,'ko','MarkerFaceColor','k','MarkerSize',6);  % punto de equilibrio

% (opcional) línea gris con B del equilibrio base como referencia
yline(fb.B,'Color',[.5 .5 .5],'LineStyle',':','LineWidth',1);

grid on; xlabel('Interest rate r'); ylabel('Assets / Bonds level');
legend({'A_I(r)','A_F(r)','A_{priv}(r)','B(r) (public supply)','r^*','Equilibrium point','B(base) ref'}, ...
       'Location','northwest');
title('Private demand for bonds vs Public bond supply');
export_fig(fig, fullfile(outdir_figs,'asset_demand_vs_B_function_of_r'));

% Guardar tabla de las curvas (útil para inspección)
T_curves = table(r_span(:), Apriv_I_grid(:), Apriv_F_grid(:), Apriv_grid(:), B_grid(:), Y_grid(:), Sgrid(:), ...
    'VariableNames', {'r','A_I','A_F','A_priv','B_of_r','Y_of_r','S'});
writetable(T_curves, fullfile(outdir_tabs,'curves_market_equilibrium.csv'));


% === NUEVO: Trayectoria de ajuste de r (bisección) ===
if isfield(base,'r_path') && numel(base.r_path)>=1
    fig=figure('Name','Ajuste de r: path de bisección');
    plot(base.r_path, base.S_path,'-o','LineWidth',1.8); hold on;
    yline(0,'k--'); xline(r,'r:','LineWidth',1.2);
    grid on; xlabel('r (iteraciones)'); ylabel('S(r) = A_{priv}-B');
    title(sprintf('Búsqueda de r: converge a r^*=%.4f (n=%d)', r, numel(base.r_path)));
    export_fig(fig, fullfile(outdir_figs,'r_adjustment_path'));
end

% Fiscal (ingresos/egresos y balances)
fig=figure('Name','Fiscal accounts and balances');
subplot(1,3,1); bar(categorical({'Labor tax','VAT'}), [fb.Tl, fb.Tc]); ylabel('Revenues'); title('Revenues'); grid on;
subplot(1,3,2); bar(categorical({'Debt service','Public good','Transfers'}), [fb.rB, fb.G, fb.Tr]); ylabel('Expenditures'); title('Expenditures'); grid on;
subplot(1,3,3); bar(categorical({'Primary balance','Overall (global)'}), [fb.PB, fb.BB]); yline(0,'k--'); ylabel('Balance'); title('Balances'); grid on;
export_fig(fig, fullfile(outdir_figs,'fiscal_accounts_and_balances'));

% Borrowers/Lenders
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
fig=figure('Name','Borrowers-Lenders shares & volumes');
subplot(1,2,1); bar(catX, [base.borrowers.fracBorrow(:) base.borrowers.fracLend(:)]); grid on;
ylabel('Share'); legend({'Borrowers (a<0)','Lenders (a>0)'},'Location','best'); title('Shares by type');
subplot(1,2,2); bar(catX, [abs(base.borrowers.volBorrow(:)) base.borrowers.volLend(:)]); grid on;
ylabel('Volume'); legend({'|Debt|','Savings'},'Location','best'); title('Volumes by type');
export_fig(fig, fullfile(outdir_figs,'borrowers_lenders'));

% Lorenz (riqueza total, existente)
fig=figure('Name','Lorenz wealth (total)'); [~,Lw,cumPop]=lorenz_from_density(a,g);
plot(cumPop,Lw,'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Wealth share'); text(0.60,0.12,sprintf('Gini_W^T = %.3f', base.stats.giniW(3)),'FontSize',12);
title('Lorenz curve (Wealth, total)'); export_fig(fig, fullfile(outdir_figs,'lorenz_wealth_total'));

% === NUEVO: Lorenz riqueza por tipo ===
fig=figure('Name','Lorenz wealth by type');
[LwI,cumI] = lorenz_from_assets_single(a, g(:,1));
[LwF,cumF] = lorenz_from_assets_single(a, g(:,2));
plot(cumI,LwI,'LineWidth',2); hold on;
plot(cumF,LwF,'LineWidth',2);
plot([0,1],[0,1],'k--');
grid on; axis square;
xlabel('Population share'); ylabel('Wealth share');
legend({sprintf('Informal (Gini=%.3f)', base.stats.giniW(1)), ...
        sprintf('Formal (Gini=%.3f)',   base.stats.giniW(2)), ...
        '45°'}, 'Location','southeast');
title('Lorenz curve (Wealth) by type');
export_fig(fig, fullfile(outdir_figs,'lorenz_wealth_bytype'));

% Lorenz (consumo total, existente)
fig=figure('Name','Lorenz consumption (total)'); [~,Lc,cumPopC]=lorenz_from_values(a,g,c(:,1)+c(:,2));
plot(cumPopC,Lc,'LineWidth',2); hold on; plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Consumption share'); text(0.60,0.12,sprintf('Gini_C^T = %.3f', base.stats.giniC(3)),'FontSize',12);
title('Lorenz curve (Consumption, total)'); export_fig(fig, fullfile(outdir_figs,'lorenz_consumption_total'));

% === NUEVO: Lorenz consumo por tipo ===
fig=figure('Name','Lorenz consumption by type');
[LcI,cI] = lorenz_from_values_single(c(:,1), g(:,1), a);
[LcF,cF] = lorenz_from_values_single(c(:,2), g(:,2), a);
plot(cI, LcI, 'LineWidth',2); hold on;
plot(cF, LcF, 'LineWidth',2);
plot([0,1],[0,1],'k--');
grid on; axis square;
xlabel('Population share'); ylabel('Consumption share');
legend({sprintf('Informal (Gini=%.3f)', base.stats.giniC(1)), ...
        sprintf('Formal (Gini=%.3f)',   base.stats.giniC(2)), ...
        '45°'}, 'Location','southeast');
title('Lorenz curve (Consumption) by type');
export_fig(fig, fullfile(outdir_figs,'lorenz_consumption_bytype'));

% MPC(a)
eps_z=0.01; alt=cfg; alt.fix_r=1; alt.r_guess=r; alt.rmin=r; alt.rmax=r;
alt.z1 = cfg.z1*(1+eps_z); alt.z2 = cfg.z2*(1+eps_z);
solP = solve_two_type_huggett_fiscal_Bfixed(alt);
dres1 = (cfg.z1*(1+eps_z)-cfg.z1) * (1 + cfg.phi);
dres2 = (cfg.z2*(1+eps_z)-cfg.z2) * (1 - cfg.tau_l);
MPC1 = (solP.c(:,1)-c(:,1))/max(dres1,1e-12);
MPC2 = (solP.c(:,2)-c(:,2))/max(dres2,1e-12);
fig=figure('Name','MPC by assets'); plot(a,MPC1,'LineWidth',2); hold on; plot(a,MPC2,'LineWidth',2);
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)'); legend({'Informal','Formal'},'Location','best');
title('MPC along asset grid (finite-diff)'); export_fig(fig, fullfile(outdir_figs,'mpc_by_assets'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ---- Helpers ----
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end, print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
end
function [gT,L,cumPop]=lorenz_from_density(a,g)
    da=a(2)-a(1); gT=g(:,1)+g(:,2); W=sum(gT)*da; [as,ix]=sort(a); gs=gT(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as-min(0,min(as))+1e-12; cumWealth=cumsum(wealth_pos.*gs); L=cumWealth/max(cumWealth(end),1e-12);
end
function [vals,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w); vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix);
    cumPop=cumsum(w_s)/W; cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end

% ===== NUEVOS helpers de Lorenz por tipo =====
function [L,cumPop]=lorenz_from_assets_single(a, gk)
    da=a(2)-a(1); W=sum(gk)*da; [as,ix]=sort(a); gs=gk(ix)*da; cumPop=cumsum(gs)/max(W,1e-12);
    wealth_pos = as - min(0,min(as)) + 1e-12;  % desplaza a>=0
    cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
function [L,cumPop]=lorenz_from_values_single(xk, gk, a)
    da=a(2)-a(1); w = gk(:)*da; W=sum(w); [xs,ix]=sort(xk(:)); ws=w(ix);
    cumPop = cumsum(ws)/max(W,1e-12);
    cumx   = cumsum(xs.*ws);
    L      = cumx / max(cumx(end),1e-12);
end
