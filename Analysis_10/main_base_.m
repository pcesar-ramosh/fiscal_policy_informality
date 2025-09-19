%% main_base_dashboard.m  (BASE sin shocks)
% Genera CSVs y gráficos:
% - c(a), s(a), g(a)
% - Curvas A_priv(r) vs B(r) y equilibrio r*
% - Barras de Gini (I/F/Total)
% - MPC(a) + c(a) por tipo (r fijo)
% Guarda CSV en out/ y figuras en fig/

clear; clc; close all;

if ~exist('out','dir'), mkdir out; end
if ~exist('fig','dir'), mkdir fig; end

%% ===== 1) PARÁMETROS BASE =====
RRA   = 2.30;         % RRA_I = RRA_F
rho   = 0.05;
theta = 0.02;
tau_l = 0.15;
tau_c = 0.15;
Gov   = 0.05;
phi   = 0.10;
z1    = 0.33;
z2    = 1.00;

eta_target = 0.54;
p22_bar    = 0.8155;

Igrid = 700;
amax  = 5.0;
amin  = -0.30*z1;

r_guess = 0.03; rmin = 0.005; rmax = 0.08;

paramsBase = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov,'phi',phi, ...
    'z1',z1,'z2',z2,'I',Igrid,'amax',amax,'amin',amin, ...
    'r_guess',r_guess,'rmin',rmin,'rmax',rmax, ...
    'p22_bar',p22_bar,'eta_target',eta_target, ...
    'fix_r',0);  % r endógeno (clearing)

%% ===== 2) RESOLVER ECONOMÍA BASE =====
base = huggett_base_function(paramsBase);

% Short-hands
a   = base.a;   g = base.g;   c = base.c;   s = base.s;
popI = base.popI; popF = base.popF;
Y    = base.Y; Ctot = base.Ctot; r_star = base.r; fb = base.fiscal;

fprintf('\n== BASE ==\n');
fprintf('r*        = %.6f\n', r_star);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y         = %.6f,   Ctot = %.6f\n', Y, Ctot);
fprintf('PB        = %.6f,   rB   = %.6f,   BB = %.6f (≈0)\n', fb.PB, fb.rB, fb.BB);

%% ===== 3) EXPORTAR CSVs =====
% 3.1 Fiscal breakdown
T_fiscal = table;
T_fiscal.scenario   = "BASE";
T_fiscal.r          = r_star;
T_fiscal.Tl         = fb.Tl;
T_fiscal.Tc         = fb.Tc;
T_fiscal.Ingresos   = fb.Tl + fb.Tc;
T_fiscal.Tr         = fb.Tr;
T_fiscal.G          = fb.G;
T_fiscal.rB         = fb.rB;
T_fiscal.Gastos     = fb.Tr + fb.G + fb.rB;
T_fiscal.PB         = fb.PB;
T_fiscal.B          = fb.B;
T_fiscal.BB         = fb.BB;
writetable(T_fiscal, fullfile('out','base_fiscal_breakdown.csv'));

% 3.2 Mercado de activos
da = a(2)-a(1);
A_priv = sum( (g(:,1)+g(:,2)).*a ) * da;
T_assets = table("BASE", r_star, A_priv, fb.B, 'VariableNames', ...
                {'scenario','r','A_private','B_public'});
writetable(T_assets, fullfile('out','base_asset_market.csv'));

% 3.3 Stats de hogar
S = base.stats;
T_stats = table;
T_stats.scenario = "BASE";
T_stats.r = r_star; T_stats.popI=popI; T_stats.popF=popF; T_stats.Y=Y; T_stats.Ctot=Ctot;
T_stats.wealth_mean_I = S.wealth_mean(1);
T_stats.wealth_mean_F = S.wealth_mean(2);
T_stats.wealth_mean_T = S.wealth_mean(3);
T_stats.wealth_med_I  = S.wealth_median(1);
T_stats.wealth_med_F  = S.wealth_median(2);
T_stats.wealth_med_T  = S.wealth_median(3);
T_stats.gini_I = S.gini(1);
T_stats.gini_F = S.gini(2);
T_stats.gini_T = S.gini(3);
T_stats.cons_mean_I = S.cons_mean(1);
T_stats.cons_mean_F = S.cons_mean(2);
T_stats.cons_mean_T = S.cons_mean(3);
T_stats.cons_med_I  = S.cons_median(1);
T_stats.cons_med_F  = S.cons_median(2);
T_stats.cons_med_T  = S.cons_median(3);
writetable(T_stats, fullfile('out','base_household_stats.csv'));

% 3.4 Prestatarios / prestamistas
Borr = base.borrowers;
T_borr = table;
T_borr.scenario = "BASE";
T_borr.fracBorrow_I = Borr.fracBorrow(1);
T_borr.fracBorrow_F = Borr.fracBorrow(2);
T_borr.fracLend_I   = Borr.fracLend(1);
T_borr.fracLend_F   = Borr.fracLend(2);
T_borr.volBorrow_I  = Borr.volBorrow(1);
T_borr.volBorrow_F  = Borr.volBorrow(2);
T_borr.volLend_I    = Borr.volLend(1);
T_borr.volLend_F    = Borr.volLend(2);
writetable(T_borr, fullfile('out','base_borrowers_lenders.csv'));

%% ===== 4) GRÁFICOS BÁSICOS =====
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% 4.1 Consumo por tipo
f1 = figure('Name','Consumo por tipo (BASE)');
plot(a, c(:,1), 'LineWidth',1.6); hold on;
plot(a, c(:,2), 'LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('c(a)');
legend('Informal','Formal','Location','best'); title('Políticas de consumo'); set(gcf,'Color','w');
exportgraphics(f1, fullfile('fig','consumo_base.png'), 'Resolution', 200);

% 4.2 Ahorro por tipo
f2 = figure('Name','Ahorro por tipo (BASE)');
plot(a, s(:,1), 'LineWidth',1.6); hold on;
plot(a, s(:,2), 'LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('s(a)');
legend('Informal','Formal','Location','best'); title('Políticas de ahorro'); set(gcf,'Color','w');
exportgraphics(f2, fullfile('fig','ahorro_base.png'), 'Resolution', 200);

% 4.3 Distribución g(a)
f3 = figure('Name','Distribucion de riqueza g(a)');
subplot(2,1,1);
bar(a, g(:,1), 'FaceAlpha',0.6, 'EdgeColor','none'); grid on;
xlim([min(a) 1.0]); xlabel('a'); ylabel('g_1(a)'); title('Informales'); set(gcf,'Color','w');
subplot(2,1,2);
bar(a, g(:,2), 'FaceAlpha',0.6, 'EdgeColor','none'); grid on;
xlim([min(a) 1.0]); xlabel('a'); ylabel('g_2(a)'); title('Formales');
exportgraphics(f3, fullfile('fig','distribucion_base.png'), 'Resolution', 200);

%% ===== 5) CURVAS A_priv(r) y B(r) =====
r_grid = linspace(max(0.5*r_star,0.002), min(1.5*r_star,0.12), 9);
A_priv_vec = zeros(size(r_grid));
B_pub_vec  = zeros(size(r_grid));
PB_vec     = zeros(size(r_grid));

for ii = 1:numel(r_grid)
    r_fix = r_grid(ii);
    parR = paramsBase;
    parR.r_guess = r_fix; parR.rmin = r_fix; parR.rmax = r_fix;
    parR.fix_r = 1;
    solR = huggett_base_function(parR);

    daR = solR.a(2) - solR.a(1);
    A_priv_vec(ii) = sum((solR.g(:,1) + solR.g(:,2)) .* solR.a) * daR;
    PB_vec(ii)     = solR.fiscal.PB;
    B_pub_vec(ii)  = solR.fiscal.B;   % = PB/r_fix
end

f4 = figure('Name','Mercado de activos: A_{priv}(r) vs B(r)');
plot(r_grid, A_priv_vec, 'LineWidth',1.8); hold on;
plot(r_grid, B_pub_vec,  'LineWidth',1.8);
yline(0,'k:'); grid on; xlabel('tasa de interés r'); ylabel('activos agregados');
legend('A_{priv}(r)','B(r)','Location','best');
title(sprintf('Cruce en r*=%.4f (BASE)', r_star)); set(gcf,'Color','w');
exportgraphics(f4, fullfile('fig','mercado_activos_curvas.png'), 'Resolution', 200);

% Estilo Walras (r en vertical)
[~, idx_star] = min(abs(r_grid - r_star));
f5 = figure('Name','Equilibrio de activos (estilo Walras)');
plot(A_priv_vec, r_grid, 'LineWidth',1.8); hold on;
plot(B_pub_vec,  r_grid, 'LineWidth',1.8);
plot([A_priv_vec(idx_star) A_priv_vec(idx_star)], [min(r_grid) max(r_grid)], 'k--');
plot([min([A_priv_vec B_pub_vec],[],'all') max([A_priv_vec B_pub_vec],[],'all')], [r_star r_star], 'k--');
grid on; ylabel('r'); xlabel('activos agregados');
legend('A_{priv}(r)','B(r)','r* vertical','r* horizontal','Location','best');
set(gcf,'Color','w');
exportgraphics(f5, fullfile('fig','equilibrio_walras.png'), 'Resolution', 200);

%% ===== 6) Gini en barras =====
gI = S.gini(1); gF = S.gini(2); gT = S.gini(3);
f6 = figure('Name','Gini (BASE)');
bar(categorical({'Informal','Formal','Total'}), [gI, gF, gT]);
ylabel('Gini'); title('Índice de Gini por tipo y total'); grid on; set(gcf,'Color','w');
exportgraphics(f6, fullfile('fig','gini_barras.png'), 'Resolution', 200);

%% ===== 7) MPC(a) y c(a) por tipo (r fijo) =====
eps_z = 0.01; % +1% de z1,z2
parM = paramsBase; parM.fix_r = 1; parM.r_guess = r_star; parM.rmin = r_star; parM.rmax = r_star;
parM.z1 = paramsBase.z1 * (1 + eps_z);
parM.z2 = paramsBase.z2 * (1 + eps_z);
solP = huggett_base_function(parM);

% Cambios en recursos por tipo
dres1 = (paramsBase.z1*(1+eps_z) - paramsBase.z1) * (1 + paramsBase.phi);   % dz1*(1+phi)
dres2 = (paramsBase.z2*(1+eps_z) - paramsBase.z2) * (1 - paramsBase.tau_l); % dz2*(1-tau_l)

MPC1 = (solP.c(:,1) - base.c(:,1)) / max(dres1, 1e-12);
MPC2 = (solP.c(:,2) - base.c(:,2)) / max(dres2, 1e-12);

% Informales
f7 = figure('Name','Consumo y MPC por activos (Informal)');
yyaxis left
plot(a, base.c(:,1), 'LineWidth',1.6); ylabel('c_I(a)');
yyaxis right
plot(a, MPC1, 'LineWidth',1.6);
grid on; xlabel('Activos a'); ylabel('MPC_I(a)');
legend('c_I(a)','MPC_I(a)','Location','best'); title('Informales'); set(gcf,'Color','w');
exportgraphics(f7, fullfile('fig','mpc_informal.png'), 'Resolution', 200);

% Formales
f8 = figure('Name','Consumo y MPC por activos (Formal)');
yyaxis left
plot(a, base.c(:,2), 'LineWidth',1.6); ylabel('c_F(a)');
yyaxis right
plot(a, MPC2, 'LineWidth',1.6);
grid on; xlabel('Activos a'); ylabel('MPC_F(a)');
legend('c_F(a)','MPC_F(a)','Location','best'); title('Formales'); set(gcf,'Color','w');
exportgraphics(f8, fullfile('fig','mpc_formal.png'), 'Resolution', 200);

% MPC promedio ponderado (reporte en consola)
gI = base.g(:,1); gF = base.g(:,2);
WI = sum(gI)*da; WF = sum(gF)*da;
MPC1_bar = sum(MPC1 .* gI) * da / max(WI,1e-12);
MPC2_bar = sum(MPC2 .* gF) * da / max(WF,1e-12);
fprintf('MPC promedio (I/F): %.4f / %.4f\n', MPC1_bar, MPC2_bar);

disp('Listo. CSVs en out/ y figuras en fig/');
