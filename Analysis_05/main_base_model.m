%% =============================================================
%  main_base_model.m
%  BASELINE (formal vs. informal) con:
%   - IVA tau_c (ambos), impuesto laboral tau_l (solo formales)
%   - Transferencias Tr_i = phi * y_i (solo informales)
%   - Bien público Gov en utilidad (exógeno baseline)
%   - r endógeno por clearing con B(r) del presupuesto público
%   - Exporta CSVs y grafica políticas, bonos por tipo, y r vs. informalidad
%
%  Requiere:
%     huggett_Equi_RRA_function_transfer.m  (versión de abajo)
% =============================================================

clear; clc; close all;

%% -------- Parámetros para la grilla de RRA (heterogeneidad en informales)
n_agents  = 19;                 % número de puntos para RRA(informal)
s_min     = 3.15;
s_max     = 5.30;

eta_vector = 0.64 * ones(1, n_agents);        % tamaño de informalidad
sI_vector1 = linspace(s_min, s_max, n_agents); % RRA (informal) heterog.
sF_vector1 = 5.30 * ones(1, n_agents);        % RRA (formal) constante

%% -------- Resolver el GE para cada combinación (vectorizado dentro)
[r_vec, ir, pop1_vector, statsMatrix, statsCMatrix, GDistribution, a, ...
 Distribution, Fiscal, Cpolicies, Spolicies] = ...
    huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1);

% Baseline representativo (último jj)
jj0 = n_agents;
g   = Distribution{jj0};           % [I x 2] (informal, formal) densidad en activos
aa  = a(:);
da  = (aa(end)-aa(1))/(numel(aa)-1);

%% -------- Fracciones de prestatarios (baseline)
isBorrow = (aa < 0);
massI = sum(g(:,1))*da; massF = sum(g(:,2))*da;
fracBorrow_total = sum((g(:,1)+g(:,2)).*isBorrow)*da;
fracBorrow_inf   = sum(g(:,1).*isBorrow)*da / max(massI,eps);
fracBorrow_for   = sum(g(:,2).*isBorrow)*da / max(massF,eps);

%% -------- Oferta/Demanda de bonos por tipo (baseline)
IaB = sum(g(:,1).*isBorrow.*aa)*da;         % informal prestatarios (a<0)
IaL = sum(g(:,1).*(aa>0).*aa)*da;           % informal prestamistas   (a>0)
FaB = sum(g(:,2).*isBorrow.*aa)*da;         % formal prestatarios
FaL = sum(g(:,2).*(aa>0).*aa)*da;

totalSupply = IaB + FaB;    % total borrowing  (<0)
totalDemand = IaL + FaL;    % total lending    (>0)

%% -------- Políticas (consumo y "ahorro/flujo de activos") por tipo
c_baseline   = Cpolicies{jj0};      % [I x 2]  c(a): col1 informal, col2 formal
sflow_baseln = Spolicies{jj0};      % [I x 2]  \dot a(a): col1 informal, col2 formal

%% ===================== GRAFICOS BASE =====================

% (1) Distribución de activos por tipo (masa por bin)
figure;
edges = linspace(min(aa), max(aa), 81);   % 80 bins -> 81 edges
nb = numel(edges) - 1;
idx = discretize(aa, edges);

mass_total = (g(:,1) + g(:,2)) * da;
mass_inf   = (g(:,1)) * da;
mass_for   = (g(:,2)) * da;

valid = ~isnan(idx);
counts_total = accumarray(idx(valid), mass_total(valid), [nb,1], @sum, 0);
counts_inf   = accumarray(idx(valid), mass_inf(valid),   [nb,1], @sum, 0);
counts_for   = accumarray(idx(valid), mass_for(valid),   [nb,1], @sum, 0);

histogram('BinEdges',edges,'BinCounts',counts_total, 'DisplayStyle','bar'); hold on;
histogram('BinEdges',edges,'BinCounts',counts_inf,   'DisplayStyle','bar');
histogram('BinEdges',edges,'BinCounts',counts_for,   'DisplayStyle','bar');
legend('Total','Informal','Formal'); xlabel('Activos a'); ylabel('Masa');
title('Distribución de riqueza (activos)');
set(gcf,'Color','w');

% (2) Fracciones de prestatarios
figure; 
bar([fracBorrow_total, fracBorrow_inf, fracBorrow_for]*100);
set(gca,'XTickLabel',{'Total','Informal','Formal'});
ylabel('Prestatarios (%)'); title('Fracción de prestatarios (a<0)'); grid on;
set(gcf,'Color','w');

% (3) Oferta vs Demanda (niveles) por tipo y total
figure;
barData = [abs(IaB) IaL; abs(FaB) FaL];  % usar abs para oferta (borrowing)
bar(barData,'stacked'); 
set(gca,'XTickLabel',{'Informal','Formal'});
legend('Oferta (|a<0|)','Demanda (a>0)','Location','best');
ylabel('Activos agregados'); title('Mercado de bonos: oferta/demanda por tipo'); grid on;
set(gcf,'Color','w');

% (4) Políticas: Consumo c(a)
figure;
plot(aa, c_baseline(:,1), '-', 'LineWidth', 1.6); hold on;
plot(aa, c_baseline(:,2), '--', 'LineWidth', 1.6);
xlabel('Activos a'); ylabel('Consumo c(a)'); grid on;
legend('Informal','Formal','Location','best');
title('Política de consumo c(a) por tipo');
set(gcf,'Color','w');

% (5) Políticas: Ahorro/flujo de activos \dot a(a)
figure;
plot(aa, sflow_baseln(:,1), '-', 'LineWidth', 1.6); hold on;
plot(aa, sflow_baseln(:,2), '--', 'LineWidth', 1.6);
yline(0,'k:');
xlabel('Activos a'); ylabel('\dot{a}(a)'); grid on;
legend('Informal','Formal','Location','best');
title('Política de ahorro (flujo de activos) por tipo');
set(gcf,'Color','w');

%% ===================== ANALISIS r vs. INFORMALIDAD =====================
% Barrido de eta en [0.40, 0.80] con RRA fija (punto medio)
eta_grid = linspace(0.40, 0.80, 9);
sI_mid   = median(sI_vector1);
sF_mid   = median(sF_vector1);

[r_eta, ~, pop_eta, ~, ~, ~, ~, ~, Fiscal_eta] = ...
    huggett_Equi_RRA_function_transfer(eta_grid, sI_mid*ones(size(eta_grid)), sF_mid*ones(size(eta_grid)));

figure;
plot(pop_eta, r_eta, '-o','LineWidth',1.6);
xlabel('Informalidad \eta'); ylabel('Tasa de interés r'); grid on;
title('Relación r(\eta) en el baseline (con B(r) endógeno)');
set(gcf,'Color','w');

%% ===================== EXPORTAR CSVs =====================

% (A) Estadísticas de riqueza por simulación (incluye Gini por tipo y total)
% statsMatrix cols:
% [gmean_inf gmean_for gmed_inf gmed_for gmean_tot gmed_tot gini_inf gini_for gini_tot p11]
T_wealth = table( ...
    eta_vector(:), sI_vector1(:), sF_vector1(:), r_vec(:), ...
    statsMatrix(:,1), statsMatrix(:,2), statsMatrix(:,3), statsMatrix(:,4), ...
    statsMatrix(:,5), statsMatrix(:,6), statsMatrix(:,7), statsMatrix(:,8), ...
    statsMatrix(:,9), statsMatrix(:,10), ...
    'VariableNames', {'eta','sI','sF','r', ...
      'gmean_inf','gmean_for','gmed_inf','gmed_for', ...
      'gmean_tot','gmed_tot','gini_inf','gini_for','gini_tot','p11'} );
writetable(T_wealth, 'stats_wealth.csv');

% (B) Estadísticas de consumo por simulación
T_cons = table( ...
    eta_vector(:), sI_vector1(:), sF_vector1(:), r_vec(:), ...
    statsCMatrix(:,1), statsCMatrix(:,2), statsCMatrix(:,3), ...
    statsCMatrix(:,4), statsCMatrix(:,5), statsCMatrix(:,6), ...
    'VariableNames', {'eta','sI','sF','r', ...
                      'Cmean_inf','Cmean_for','Cmed_inf','Cmed_for', ...
                      'Cmean_tot','Cmed_tot'} );
writetable(T_cons, 'stats_consumption.csv');

% (C) Resumen fiscal por simulación
% Fiscal{jj} = [Tc, Tl, Tr, G, rB, Gap, popI, popF, Y, Btarget_eq];
F = cell2mat(Fiscal(:));
T_fiscal = table( ...
    eta_vector(:), sI_vector1(:), sF_vector1(:), r_vec(:), ...
    F(:,1), F(:,2), F(:,3), F(:,4), F(:,5), F(:,6), F(:,7), F(:,8), F(:,9), F(:,10), ...
    'VariableNames', {'eta','sI','sF','r','Tc','Tl','Tr','G','rB','Gap','popI','popF','Y','Btarget'} );
writetable(T_fiscal, 'fiscal_summary.csv');

% (D) Políticas baseline (c y \dot a) por tipo
T_policy = table(aa, c_baseline(:,1), c_baseline(:,2), sflow_baseln(:,1), sflow_baseln(:,2), ...
    'VariableNames', {'a','c_informal','c_formal','sflow_informal','sflow_formal'});
writetable(T_policy, 'baseline_policies.csv');

% (E) Descomposición de bonos baseline
T_bonds = table(IaB, IaL, FaB, FaL, totalSupply, totalDemand, ...
    'VariableNames', {'IaB','IaL','FaB','FaL','totalSupply','totalDemand'});
writetable(T_bonds, 'bond_decomposition_baseline.csv');

% (F) r vs. informalidad
Fe = cell2mat(Fiscal_eta(:));
T_r_eta = table(pop_eta(:), r_eta(:), Fe(:,7), Fe(:,8), Fe(:,9), Fe(:,10), ...
    'VariableNames', {'eta','r','popI','popF','Y','Btarget'});
writetable(T_r_eta, 'r_vs_eta.csv');

fprintf('\nCSV exportados:\n - stats_wealth.csv\n - stats_consumption.csv\n - fiscal_summary.csv\n - baseline_policies.csv\n - bond_decomposition_baseline.csv\n - r_vs_eta.csv\n');
