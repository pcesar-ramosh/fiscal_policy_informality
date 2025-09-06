%% ==============================================================
%  main_covid_scenarios_transfer.m
%  Benchmark vs. COVID shock (incomes down; transfers go UP)
%  Autor: (tu nombre)
%  Fecha: (hoy)
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Configuración general
%% -----------------------
n_agents  = 20;                 % # de puntos para RRA informal
s_min     = 0.15;               % RRA informal min
s_max     = 0.30;               % RRA informal max

eta_vector = 0.75 * ones(1, n_agents);            % tamaño sector informal
sI_vector1 = linspace(s_min, s_max, n_agents);     % RRA informal heterogénea
sF_vector1 = 0.30 * ones(1, n_agents);             % RRA formal constante

%% -----------------------
%  Escenario 1: Benchmark
%% -----------------------
cfg_base = struct();
cfg_base.scenario              = "baseline";
cfg_base.psi_I                 = 1.00;        % sin shock
cfg_base.psi_F                 = 1.00;        % sin shock
cfg_base.transfer_multiplier   = 1.00;        % κ=1 (sin aumento)
cfg_base.keep_transfers_level  = true;        % nivel (irrelevante en baseline)
cfg_base.amin_mode             = "baseline";  % límite inferior con z1_base
cfg_base.gov_nonneg            = true;        % bien público truncado a >=0
cfg_base.phi                   = 0.13;        % regla de transferencias base
cfg_base.taxF                  = 0.10;        % impuesto renta formal
cfg_base.taxI                  = 0.00;        % impuesto renta informal
cfg_base.tauc                  = 0.00;        % impuesto al consumo (apagado)
cfg_base.theta                 = 0.02;        % prima de endeudamiento informal
cfg_base.r0                    = 0.03;        % bracket de r*
cfg_base.rmin                  = 0.01;
cfg_base.rmax                  = 0.04;

[r_b, ir_b, pop_b, a_b, g_b, c_b] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vector, sI_vector1, sF_vector1, cfg_base);

%% -----------------------
%  Escenario 2: COVID shock (ingresos ↓, transferencias ↑)
%% -----------------------
cfg_covid = cfg_base;
cfg_covid.scenario             = "covid_uptransfer";
cfg_covid.psi_I                = 0.75;        % caída 25% ingreso informal
cfg_covid.psi_F                = 0.85;        % caída 15% ingreso formal
cfg_covid.transfer_multiplier  = 1.50;        % κ>1: transfer aumenta 50% vs benchmark
cfg_covid.keep_transfers_level = true;        % mantener nivel (no tasa)
cfg_covid.amin_mode            = "baseline";  % o "shocked" si quieres más contractivo

[r_c, ir_c, pop_c, a_c, g_c, c_c] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vector, sI_vector1, sF_vector1, cfg_covid);

%% --------------------------------------------------------------
%  Comparaciones clave
%% --------------------------------------------------------------
fprintf('\n=== Equilibrium interest rate (mean across agents) ===\n');
fprintf('Benchmark r*: %.6f\n', mean(r_b));
fprintf('COVID+Transfers↑ r*: %.6f\n', mean(r_c));
fprintf('Δr (covid_up - base): %.6f\n', mean(r_c - r_b));

% Masa en la restricción (aprox: peso en a_min)
frac_const_base = zeros(n_agents, 2);
frac_const_cvd  = zeros(n_agents, 2);
for j = 1:n_agents
    Gb = g_b{j}; Gc = g_c{j};
    wb  = Gb ./ sum(Gb, 'all');
    wc  = Gc ./ sum(Gc, 'all');
    frac_const_base(j,:) = wb(1,:);   % [informal, formal]
    frac_const_cvd(j,:)  = wc(1,:);
end
fprintf('\n=== Fraction at borrowing limit (mean) ===\n');
fprintf('Informal: base=%.4f | covid+T↑=%.4f | Δ=%.4f\n', ...
    mean(frac_const_base(:,1)), mean(frac_const_cvd(:,1)), ...
    mean(frac_const_cvd(:,1)-frac_const_base(:,1)));
fprintf('Formal  : base=%.4f | covid+T↑=%.4f | Δ=%.4f\n', ...
    mean(frac_const_base(:,2)), mean(frac_const_cvd(:,2)), ...
    mean(frac_const_cvd(:,2)-frac_const_base(:,2)));

%% --------------------------------------------------------------
%  Figuras (bandas sobre todos los agentes)
%% --------------------------------------------------------------
% (1) Políticas de consumo
figure('Name','Consumption Policies: Base vs COVID+Transfers↑','Color','w');
subplot(1,2,1); hold on;
for j = 1:n_agents
    plot(a_b, c_b{j}(:,1), 'Color', [1 0 0 0.35], 'LineWidth', 1.0); % informal
    plot(a_b, c_b{j}(:,2), 'Color', [0 0 1 0.35], 'LineWidth', 1.0); % formal
end
title('Benchmark'); xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);

subplot(1,2,2); hold on;
for j = 1:n_agents
    plot(a_c, c_c{j}(:,1), 'r-', 'LineWidth', 1.2);
    plot(a_c, c_c{j}(:,2), 'b-', 'LineWidth', 1.2);
end
title('COVID shock (incomes↓, transfers↑)'); xlabel('Assets (a)'); grid on; xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast');

% (2) Distribuciones de riqueza
figure('Name','Wealth Distributions: Base vs COVID+Transfers↑','Color','w');
subplot(1,2,1); hold on;
for j = 1:n_agents
    g1 = g_b{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_b{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_b, g1, 'Color', [1 0 0 0.35], 'LineWidth', 1.0);
    plot(a_b, g2, 'Color', [0 0 1 0.35], 'LineWidth', 1.0);
end
title('Benchmark'); xlabel('Assets (a)'); ylabel('Wealth distribution'); grid on; xlim([0.01 0.60]);

subplot(1,2,2); hold on;
for j = 1:n_agents
    g1 = g_c{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_c{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_c, g1, 'r-', 'LineWidth', 1.2);
    plot(a_c, g2, 'b-', 'LineWidth', 1.2);
end
title('COVID shock (incomes↓, transfers↑)'); xlabel('Assets (a)'); grid on; xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast');

disp('Listo: escenarios ejecutados y comparados.');
