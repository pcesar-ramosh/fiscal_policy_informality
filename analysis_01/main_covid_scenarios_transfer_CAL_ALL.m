%% ==============================================================
%  main_covid_scenarios_transfer_CAL_ALL.m
%  Benchmark vs. COVID shock (incomes down; transfers go UP)
%  Versión completa: r*_j por agente + r*_common agregado
% Author: Hamilton Galindo Gil
% Contribuidores: Ian Carrasco & Cesar Ramos
%  Sept 2025
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Configuración general
%% -----------------------
n_agents  = 20;                 % # de puntos para RRA informal
s_min     = 3.15;               % RRA informal min
s_max     = 5.30;               % RRA informal max

eta_vector = 0.75 * ones(1, n_agents);            % tamaño sector informal
sI_vector1 = linspace(s_min, s_max, n_agents);     % RRA informal heterogénea
sF_vector1 = 5.30 * ones(1, n_agents);             % RRA formal constante

%% -----------------------
%  Escenario 1: Benchmark (calibrado para r* ≈ 2–3%)
%% -----------------------
cfg_base = struct();
cfg_base.scenario              = "baseline";
cfg_base.psi_I                 = 1.00;
cfg_base.psi_F                 = 1.00;
cfg_base.transfer_multiplier   = 1.00;
cfg_base.keep_transfers_level  = true;

% === Claves económicas para elevar r* y vaciar mercado agregado ===
cfg_base.amin_policy           = "absolute";  % 'absolute' | 'relative' | 'baseline' | 'shocked'
cfg_base.amin_abs              = -1.00;       % a_min más laxo (antes ~ -0.1)
cfg_base.theta                 = 0.005;       % menor spread de deuda
cfg_base.phi                   = 0.05;        % transfer base moderada

% Impuestos
cfg_base.taxF                  = 0.10;
cfg_base.taxI                  = 0.00;
cfg_base.tauc                  = 0.00;

% Búsqueda de r* por agente (bracket razonable)
cfg_base.r0                    = 0.02;
cfg_base.rmin                  = -0.01;
cfg_base.rmax                  = 0.06;

[r_b, ir_b, pop_b, a_b, g_b, c_b] = ...
    huggett_Equi_RRA_function_transfer_on_CAL(eta_vector, sI_vector1, sF_vector1, cfg_base);

%% -----------------------
%  Escenario 2: COVID shock (ingresos ↓, transferencias ↑ moderadas)
%% -----------------------
cfg_covid = cfg_base;
cfg_covid.scenario             = "covid_uptransfer";
cfg_covid.psi_I                = 0.75;        % caída 25% ingreso informal
cfg_covid.psi_F                = 0.85;        % caída 15% ingreso formal
cfg_covid.transfer_multiplier  = 1.25;        % aumento 25% (no 50%)
cfg_covid.keep_transfers_level = true;

[r_c, ir_c, pop_c, a_c, g_c, c_c] = ...
    huggett_Equi_RRA_function_transfer_on_CAL(eta_vector, sI_vector1, sF_vector1, cfg_covid);

%% --------------------------------------------------------------
%  r*_j y verificación de clearing por agente
%% --------------------------------------------------------------
fprintf('\n=== Equilibrium interest rate by agent (baseline) ===\n');
disp(r_b);

fprintf('\n[Check] Market clearing by agent (baseline):\n');
da_b = a_b(2) - a_b(1);
for j=1:n_agents
    G = g_b{j};
    S_j = (G(:,1)'*a_b + G(:,2)'*a_b)*da_b;
    fprintf('  j=%02d  r*_j=%.4f   S_j=%.3e\n', j, r_b(j), S_j);
end

fprintf('\n[Check] Market clearing by agent (covid):\n');
da_c = a_c(2) - a_c(1);
for j=1:n_agents
    G = g_c{j};
    S_j = (G(:,1)'*a_c + G(:,2)'*a_c)*da_c;
    fprintf('  j=%02d  r*_j=%.4f   S_j=%.3e\n', j, r_c(j), S_j);
end

fprintf('\n=== Mean r* across agents ===\n');
fprintf('Baseline: %.6f | COVID+T↑: %.6f | Δ: %.6f\n', mean(r_b), mean(r_c), mean(r_c - r_b));

%% --------------------------------------------------------------
%  Fracción en la restricción (condicional por ocupación)
%% --------------------------------------------------------------
frac_const_base = zeros(n_agents, 2);
frac_const_cvd  = zeros(n_agents, 2);

for j = 1:n_agents
    Gb = g_b{j}; Gc = g_c{j};
    mass_inf_b = sum(Gb(:,1))*da_b;  mass_for_b = sum(Gb(:,2))*da_b;
    mass_inf_c = sum(Gc(:,1))*da_c;  mass_for_c = sum(Gc(:,2))*da_c;

    frac_const_base(j,1) = (Gb(1,1)*da_b)/max(mass_inf_b, 1e-14);
    frac_const_base(j,2) = (Gb(1,2)*da_b)/max(mass_for_b, 1e-14);
    frac_const_cvd(j,1)  = (Gc(1,1)*da_c)/max(mass_inf_c, 1e-14);
    frac_const_cvd(j,2)  = (Gc(1,2)*da_c)/max(mass_for_c, 1e-14);
end

fprintf('\n=== Fraction at borrowing limit (mean, conditional by occupation) ===\n');
fprintf('Informal: base=%.4f | covid+T↑=%.4f | Δ=%.4f\n', ...
    mean(frac_const_base(:,1)), mean(frac_const_cvd(:,1)), mean(frac_const_cvd(:,1)-frac_const_base(:,1)));
fprintf('Formal  : base=%.4f | covid+T↑=%.4f | Δ=%.4f\n', ...
    mean(frac_const_base(:,2)), mean(frac_const_cvd(:,2)), mean(frac_const_cvd(:,2)-frac_const_base(:,2)));

%% --------------------------------------------------------------
%  Gráficos: políticas y distribuciones (bandas)
%% --------------------------------------------------------------
% (1) Políticas de consumo
figure('Name','Consumption Policies: Base vs COVID+Transfers↑','Color','w');
subplot(1,2,1); hold on;
for j = 1:n_agents
    plot(a_b, c_b{j}(:,1), 'r-', 'LineWidth', 0.8); % informal
    plot(a_b, c_b{j}(:,2), 'b-', 'LineWidth', 0.8); % formal
end
title('Benchmark'); xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);

subplot(1,2,2); hold on;
for j = 1:n_agents
    plot(a_c, c_c{j}(:,1), 'r-', 'LineWidth', 0.8);
    plot(a_c, c_c{j}(:,2), 'b-', 'LineWidth', 0.8);
end
title('COVID shock (incomes↓, transfers↑)'); xlabel('Assets (a)'); grid on; xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast');

% (2) Distribuciones de riqueza
figure('Name','Wealth Distributions: Base vs COVID+Transfers↑','Color','w');
subplot(1,2,1); hold on;
for j = 1:n_agents
    g1 = g_b{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_b{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_b, g1, 'r-', 'LineWidth', 0.8);
    plot(a_b, g2, 'b-', 'LineWidth', 0.8);
end
title('Benchmark'); xlabel('Assets (a)'); ylabel('Wealth distribution'); grid on; xlim([0.01 0.60]);

subplot(1,2,2); hold on;
for j = 1:n_agents
    g1 = g_c{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_c{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_c, g1, 'r-', 'LineWidth', 0.8);
    plot(a_c, g2, 'b-', 'LineWidth', 0.8);
end
title('COVID shock (incomes↓, transfers↑)'); xlabel('Assets (a)'); grid on; xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast');

%% --------------------------------------------------------------
%  ONE FIGURE: Agente típico — consumo y riqueza
%% --------------------------------------------------------------
[~, j_star] = min(abs(sI_vector1 - median(sI_vector1)));
a_base = a_b;  a_covid = a_c;
c_base = c_b{j_star};  c_covid = c_c{j_star};
g_base = g_b{j_star};  g_covid = g_c{j_star};

g_base_plot  = g_base;   g_base_plot(g_base_plot < 1e-5)   = NaN;
g_covid_plot = g_covid;  g_covid_plot(g_covid_plot < 1e-5) = NaN;

figure('Name','Typical agent: Before vs After','Color','w','Position',[100 100 1100 430]);
subplot(1,2,1); hold on; box on; grid on;
p1 = plot(a_base, c_base(:,1), 'r--', 'LineWidth', 1.5);
p2 = plot(a_base, c_base(:,2), 'b--', 'LineWidth', 1.5);
p3 = plot(a_covid, c_covid(:,1), 'r-',  'LineWidth', 1.5);
p4 = plot(a_covid, c_covid(:,2), 'b-',  'LineWidth', 1.5);
xlabel('Assets (a)'); ylabel('Consumption');
title(sprintf('Consumption | Typical agent (RRA_{inf}=%.3f)', sI_vector1(j_star)));
xlim([1 5]);
legend([p1 p2 p3 p4], {'Inf Base','For Base','Inf COVID+T↑','For COVID+T↑'}, 'Location','SouthEast');

subplot(1,2,2); hold on; box on; grid on;
q1 = plot(a_base,  g_base_plot(:,1),  'r--', 'LineWidth', 1.5);
q2 = plot(a_base,  g_base_plot(:,2),  'b--', 'LineWidth', 1.5);
q3 = plot(a_covid, g_covid_plot(:,1), 'r-',  'LineWidth', 1.5);
q4 = plot(a_covid, g_covid_plot(:,2), 'b-',  'LineWidth', 1.5);
xlabel('Assets (a)'); ylabel('Wealth distribution');
title('Wealth distribution (Typical agent)');
xlim([0.01 1.60]);
xline(a_base(1), ':k', 'a_{min}', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);
legend([q1 q2 q3 q4], {'Inf Base','For Base','Inf COVID+T↑','For COVID+T↑'}, 'Location','NorthEast');

%% --------------------------------------------------------------
%  r común agregado: bisección externa que vacía S_total(r)=0
%% --------------------------------------------------------------
[r_common_base, hist_base] = find_r_common_CAL(eta_vector, sI_vector1, sF_vector1, cfg_base, -0.01, 0.06);
[r_common_covid, hist_cvd] = find_r_common_CAL(eta_vector, sI_vector1, sF_vector1, cfg_covid, -0.01, 0.06);

fprintf('\n=== r*_common (BASE): %.6f ===\n', r_common_base);
fprintf('=== r*_common (COVID+T↑): %.6f ===\n', r_common_covid);

% Curva S_total(r) (Base) para ver el cruce
rgrid = linspace(-0.01, 0.06, 25);
Stot  = zeros(size(rgrid));
for k=1:numel(rgrid)
    [~, Stot(k)] = huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg_base, rgrid(k));
end
figure('Name','Exceso de ahorro agregado vs r (Base)','Color','w');
plot(rgrid, Stot, 'LineWidth', 1.8); hold on; grid on;
yline(0,'k:'); xline(r_common_base,'r--');
xlabel('r'); ylabel('S_{total}(r)');
title('Clearing agregado (Base)');

% r*_j vs RRA informal, incluyendo r_common
figure('Name','r*_j vs RRA informal','Color','w');
plot(sI_vector1, r_b, 'o-', 'LineWidth',1.5); hold on; grid on;
plot(sI_vector1, r_c, 's-', 'LineWidth',1.5);
yline(r_common_base,'--');
yline(r_common_covid,'--');
xlabel('RRA informal (s_I)'); ylabel('r*');
legend({'r*_j (Base)','r*_j (COVID+T↑)','r*_common Base','r*_common COVID+T↑'}, 'Location','best');

%% --------------------------------------------------------------
%  Exportes a Excel (opcional)
%% --------------------------------------------------------------
export_excel = true;
if export_excel
    filename = 'summary_covid_uptransfer.xlsx';

    % r*_j y r*_common
    T_r = table((1:n_agents)', r_b(:), r_c(:), ...
        'VariableNames', {'agent','r_base','r_covid_Tup'});
    writetable(T_r, filename, 'Sheet', 'r_star_per_agent');

    T_rc = table(r_common_base, r_common_covid, ...
        'VariableNames', {'r_common_base','r_common_covid'});
    writetable(T_rc, filename, 'Sheet', 'r_star_common');

    % Fracción en restricción condicional
    T_fc = table((1:n_agents)', ...
        frac_const_base(:,1), frac_const_cvd(:,1), ...
        frac_const_base(:,2), frac_const_cvd(:,2), ...
        'VariableNames', {'agent','inf_base','inf_covid_Tup','for_base','for_covid_Tup'});
    writetable(T_fc, filename, 'Sheet', 'constraint_frac');

    fprintf('Exportado a: %s\n', filename);
end

disp('Listo: r*_j, r*_common, gráficos y exportes generados.');
