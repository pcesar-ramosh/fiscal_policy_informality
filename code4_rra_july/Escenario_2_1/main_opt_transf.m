%% ==============================================================
%  main.m
%  Objetivo: Encontrar el multiplicador de transferencias kappa*
%  que mantiene el consumo medio de los informales igual al benchmark,
%  ante un shock negativo de ingresos, y evaluar proximidad de g(a).
%  Autor: (tu nombre)
%  Fecha: (hoy)
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Configuración general
%% -----------------------
n_agents  = 20;                 % # de puntos RRA informal
s_min     = 0.15;               % RRA informal min
s_max     = 0.30;               % RRA informal max

eta_vector = 0.75 * ones(1, n_agents);            % tamaño sector informal
sI_vector1 = linspace(s_min, s_max, n_agents);     % RRA informal heterogénea
sF_vector1 = 0.30 * ones(1, n_agents);             % RRA formal constante

% Config base (política fiscal y numérica)
cfg_base = struct();
cfg_base.scenario              = "baseline";
cfg_base.psi_I                 = 1.00;        % sin shock
cfg_base.psi_F                 = 1.00;        % sin shock
cfg_base.keep_transfers_level  = true;        % transfer = phi*z1_base
cfg_base.transfer_multiplier   = 1.00;        % kappa
cfg_base.amin_mode             = "baseline";  % límite inferior con z1_base
cfg_base.gov_nonneg            = true;        % trunca bien público a >=0
cfg_base.phi                   = 0.13;        % regla de transferencias base
cfg_base.taxF                  = 0.10;        % impuesto renta formal
cfg_base.taxI                  = 0.00;        % impuesto renta informal
cfg_base.tauc                  = 0.00;        % impuesto al consumo (no usado)
cfg_base.theta                 = 0.02;        % prima endeudamiento informal
cfg_base.r0                    = 0.03;        % bracket r*
cfg_base.rmin                  = 0.01;
cfg_base.rmax                  = 0.04;

%% -----------------------
%  (A) Resolver Benchmark
%% -----------------------
res_base = huggett_Equi_RRA_transfer_solver(eta_vector, sI_vector1, sF_vector1, cfg_base);

% Consumo medio informal en benchmark (promedio sobre agentes y activos)
cmean_inf_base = mean(res_base.statsC_mean_inf);

fprintf('\n=== BENCHMARK ===\n');
fprintf('r* (mean over agents): %.6f\n', mean(res_base.r_opt));
fprintf('c_mean_inf (target): %.6f\n', cmean_inf_base);

%% -----------------------
%  (B) Resolver Shock SIN política
%% -----------------------
cfg_shock = cfg_base;
cfg_shock.scenario             = "shock";
cfg_shock.psi_I                = 0.75;   % caída 25% ingreso informal
cfg_shock.psi_F                = 0.85;   % caída 15% ingreso formal
cfg_shock.transfer_multiplier  = 1.00;   % sin aumento de transfer
cfg_shock.keep_transfers_level = true;   % nivel basado en z1_base
cfg_shock.amin_mode            = "baseline"; % o "shocked" si quieres endurecer límite

res_shock = huggett_Equi_RRA_transfer_solver(eta_vector, sI_vector1, sF_vector1, cfg_shock);
cmean_inf_shock = mean(res_shock.statsC_mean_inf);

fprintf('\n=== SHOCK (NO POLICY) ===\n');
fprintf('r* (mean): %.6f\n', mean(res_shock.r_opt));
fprintf('c_mean_inf: %.6f (Δ vs base = %.6f)\n', cmean_inf_shock, cmean_inf_shock - cmean_inf_base);

%% --------------------------------------------------------------
%  (C) Buscar kappa* por bisección para igualar c_mean_inf al benchmark
%% --------------------------------------------------------------
tol    = 1e-4;     % tolerancia sobre c_mean_inf
maxitK = 20;       % iter máx. de kappa
K_low  = 1.00;     % sin política
K_high = 3.00;     % cota superior (ajústala si es necesario)

c_low  = cmean_inf_shock;  % con K_low
% Evaluar K_high una vez
cfg_try = cfg_shock; cfg_try.transfer_multiplier = K_high;
res_try = huggett_Equi_RRA_transfer_solver(eta_vector, sI_vector1, sF_vector1, cfg_try);
c_high  = mean(res_try.statsC_mean_inf);

% Chequeo de que c_mean_inf sube con kappa
if c_high < cmean_inf_base
    warning('Con K_high=%.2f no alcanzas el nivel de c_mean_inf del benchmark. Sube K_high.', K_high);
end

for it = 1:maxitK
    K_mid = 0.5*(K_low + K_high);
    cfg_try = cfg_shock; cfg_try.transfer_multiplier = K_mid;

    res_mid = huggett_Equi_RRA_transfer_solver(eta_vector, sI_vector1, sF_vector1, cfg_try);
    c_mid   = mean(res_mid.statsC_mean_inf);

    fprintf('Iter %02d | K_low=%.4f K_high=%.4f K_mid=%.4f | c_mid=%.6f | gap=%.6e\n', ...
        it, K_low, K_high, K_mid, c_mid, c_mid - cmean_inf_base);

    if abs(c_mid - cmean_inf_base) < tol
        break
    end

    if c_mid < cmean_inf_base
        K_low = K_mid;     % subir kappa
    else
        K_high = K_mid;    % bajar kappa
    end
end

kappa_star = K_mid;
res_opt    = res_mid;  % resultado con kappa*

fprintf('\n=== OPTIMAL TRANSFER (kappa*) ===\n');
fprintf('kappa* = %.6f\n', kappa_star);
fprintf('r* (mean): %.6f\n', mean(res_opt.r_opt));
fprintf('c_mean_inf (target matched): %.6f (gap=%.3e)\n', mean(res_opt.statsC_mean_inf), ...
        mean(res_opt.statsC_mean_inf)-cmean_inf_base);

%% --------------------------------------------------------------
%  Distancia L2 de g_inf(a) vs benchmark (informativa)
%% --------------------------------------------------------------
% Definimos una métrica simple: promedio (sobre agentes) de ||g_inf^scenario - g_inf^base||_2
L2_shock = zeros(n_agents,1);
L2_opt   = zeros(n_agents,1);

for j = 1:n_agents
    gB = res_base.g_opt{j}(:,1);      % informal benchmark
    gS = res_shock.g_opt{j}(:,1);     % informal shock
    gO = res_opt.g_opt{j}(:,1);       % informal óptima

    % Igualar mallas (son iguales por construcción) y limpiar ruido
    gB(gB<1e-5)=0; gS(gS<1e-5)=0; gO(gO<1e-5)=0;

    L2_shock(j) = sqrt(trapz(res_base.a, (gS - gB).^2));
    L2_opt(j)   = sqrt(trapz(res_base.a, (gO - gB).^2));
end

fprintf('\n=== L2 distance of g_inf(a) vs Benchmark (mean across agents) ===\n');
fprintf('Shock (no policy): %.6e\n', mean(L2_shock));
fprintf('Shock + kappa*:   %.6e\n', mean(L2_opt));

%% --------------------------------------------------------------
%  Figura: Typical agent (Benchmark vs Shock vs Shock+kappa*)
%% --------------------------------------------------------------
[~, j_star] = min(abs(sI_vector1 - median(sI_vector1)));  % agente típico

a    = res_base.a;
c_b  = res_base.c_opt{j_star};  g_b = res_base.g_opt{j_star};
c_s  = res_shock.c_opt{j_star};  g_s = res_shock.g_opt{j_star};
c_o  = res_opt.c_opt{j_star};    g_o = res_opt.g_opt{j_star};

g_b_plot = g_b; g_b_plot(g_b_plot<1e-5)=NaN;
g_s_plot = g_s; g_s_plot(g_s_plot<1e-5)=NaN;
g_o_plot = g_o; g_o_plot(g_o_plot<1e-5)=NaN;

figure('Name','Typical agent: Base vs Shock vs Shock+kappa*','Color','w','Position',[80 80 1150 440]);

% Consumo
subplot(1,2,1); hold on; box on; grid on;
plot(a, c_b(:,1), 'r:',  'LineWidth', 1.8);  % inf base
plot(a, c_b(:,2), 'b:',  'LineWidth', 1.8);  % for base
plot(a, c_s(:,1), 'r--', 'LineWidth', 1.8);  % inf shock
plot(a, c_s(:,2), 'b--', 'LineWidth', 1.8);  % for shock
plot(a, c_o(:,1), 'r-',  'LineWidth', 1.8);  % inf opt
plot(a, c_o(:,2), 'b-',  'LineWidth', 1.8);  % for opt
xlabel('Assets (a)'); ylabel('Consumption'); xlim([1 5]);
title(sprintf('Consumption | Typical agent (RRA_{inf}=%.3f)', sI_vector1(j_star)));
legend({'Inf Base','For Base','Inf Shock','For Shock','Inf Opt','For Opt'}, 'Location','SouthEast'); legend boxoff;

% Distribución g(a)
subplot(1,2,2); hold on; box on; grid on;
plot(a, g_b_plot(:,1), 'r:',  'LineWidth', 1.8);
plot(a, g_b_plot(:,2), 'b:',  'LineWidth', 1.8);
plot(a, g_s_plot(:,1), 'r--', 'LineWidth', 1.8);
plot(a, g_s_plot(:,2), 'b--', 'LineWidth', 1.8);
plot(a, g_o_plot(:,1), 'r-',  'LineWidth', 1.8);
plot(a, g_o_plot(:,2), 'b-',  'LineWidth', 1.8);
xlabel('Assets (a)'); ylabel('Wealth distribution'); xlim([0.01 0.60]);
title('Wealth distribution (Typical agent)');
legend({'Inf Base','For Base','Inf Shock','For Shock','Inf Opt','For Opt'}, 'Location','NorthEast'); legend boxoff;

xline(a(1), ':k', 'a_{min}', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);

%% --------------------------------------------------------------
%  (Opcional) Exportes a Excel — APAGADOS por defecto
%% --------------------------------------------------------------
export_excel = false;   % <--- cambia a true si quieres exportar
if export_excel
    filename = 'optimal_transfer_results.xlsx';

    % r* por agente
    T_r = table((1:n_agents)', res_base.r_opt(:), res_shock.r_opt(:), res_opt.r_opt(:), ...
        'VariableNames', {'agent','r_base','r_shock','r_opt'});
    writetable(T_r, filename, 'Sheet', 'r_star');

    % c_mean_inf por agente
    T_cmean = table((1:n_agents)', res_base.statsC_mean_inf(:), res_shock.statsC_mean_inf(:), res_opt.statsC_mean_inf(:), ...
        'VariableNames', {'agent','cmean_inf_base','cmean_inf_shock','cmean_inf_opt'});
    writetable(T_cmean, filename, 'Sheet', 'cmean_inf');

    % L2 por agente
    T_L2 = table((1:n_agents)', L2_shock(:), L2_opt(:), ...
        'VariableNames', {'agent','L2_shock','L2_opt'});
    writetable(T_L2, filename, 'Sheet', 'L2_dist');

    % Series del agente típico
    T_a = table(a, 'VariableNames', {'a'});
    T_c = table(c_b(:,1), c_b(:,2), c_s(:,1), c_s(:,2), c_o(:,1), c_o(:,2), ...
       'VariableNames', {'c_inf_base','c_for_base','c_inf_shock','c_for_shock','c_inf_opt','c_for_opt'});
    T_g = table(g_b(:,1), g_b(:,2), g_s(:,1), g_s(:,2), g_o(:,1), g_o(:,2), ...
       'VariableNames', {'g_inf_base','g_for_base','g_inf_shock','g_for_shock','g_inf_opt','g_for_opt'});

    writetable(T_a, filename, 'Sheet', 'typical_a');
    writetable(T_c, filename, 'Sheet', 'typical_c');
    writetable(T_g, filename, 'Sheet', 'typical_g');

    fprintf('Exportado a: %s\n', filename);
end

disp('Listo: kappa* encontrado, resultados graficados y exportes en OFF por defecto.');
