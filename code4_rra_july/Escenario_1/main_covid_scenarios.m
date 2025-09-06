%% ==============================================================
%  main_covid_scenarios.m
%  Benchmark vs. COVID shock (ingreso cae, transferencias fijas)
%  Autor: Hamilton Galindo, Ian Carrasco y Cesar Ramos
%  Fecha: Agosto 2025
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  General Setting
%% -----------------------
% Heterogeneidad en aversión al riesgo (informales)
n_agents = 20;                 % Num of Agents
s_min    = 3.15;               % RRA informal mínima
s_max    = 5.30;               % RRA informal máxima

eta_vector = 0.75 * ones(1, n_agents);     % tamaño del sector informal (exógeno)
sI_vector1 = linspace(s_min, s_max, n_agents);   % RRA informal heterogénea
sF_vector1 = 5.30 * ones(1, n_agents);           % RRA formal (constante)

%% -----------------------
%  Escenario 1: Benchmark
%% -----------------------
cfg_base = struct();
cfg_base.scenario              = "baseline"; % sin shock
cfg_base.psi_I                 = 1.00;       % multiplicador ingreso informal
cfg_base.psi_F                 = 1.00;       % multiplicador ingreso formal
cfg_base.keep_transfers_level  = true;       % (irrelevante en baseline)
cfg_base.amin_mode             = "baseline"; % usa a_min del benchmark (z1_base)
cfg_base.gov_nonneg            = true;       % trunca bien público a >= 0
cfg_base.phi                   = 0.13;       % transferencia como % del ingreso informal BASE
cfg_base.taxF                  = 0.10;       % tasa de impuesto renta formal
cfg_base.taxI                  = 0.00;       % renta informal
cfg_base.tauc                  = 0.00;       % activar impuesto al consumo (e igual en Function)
cfg_base.theta                 = 0.02;       % prima de endeudamiento informal
cfg_base.r0                    = 0.03;       % guess de r
cfg_base.rmin                  = 0.01;       % cota inferior
cfg_base.rmax                  = 0.04;       % cota superior

[r_b, ir_b, pop_b, a_b, g_b, c_b] = ...
    huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1, cfg_base);

%% -----------------------
%  Escenario 2: COVID shock
%  (cae ingreso; transferencias se mantienen en NIVEL)
%% -----------------------
cfg_covid = cfg_base;
cfg_covid.scenario             = "covid";
cfg_covid.psi_I                = 0.75;       % caída 25% ingreso informal
cfg_covid.psi_F                = 0.85;       % caída 15% ingreso formal
cfg_covid.keep_transfers_level = true;       % mantener monto pre-shock
cfg_covid.amin_mode            = "baseline"; % O usa "shocked" si quieres que a_min caiga con z1

[r_c, ir_c, pop_c, a_c, g_c, c_c] = ...
    huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1, cfg_covid);

%% --------------------------------------------------------------
%  Comparaciones clave
%% --------------------------------------------------------------
fprintf('\n=== Tasas de interés de equilibrio (promedio sobre agentes) ===\n');
fprintf('Benchmark r*: %.6f\n', mean(r_b));
fprintf('COVID     r*: %.6f\n', mean(r_c));
fprintf('Δr (covid - base): %.6f\n', mean(r_c - r_b));

% Fracción promedio pegada a la restricción (aprox: masa en el primer punto)
frac_const_base = zeros(n_agents, 2);
frac_const_cvd  = zeros(n_agents, 2);

for j = 1:n_agents
    G  = g_b{j};  % I x 2
    Gc = g_c{j};
    w_base = G ./ sum(G, 'all');
    w_cvd  = Gc ./ sum(Gc, 'all');
    frac_const_base(j, :) = w_base(1, :); % [inf, form]
    frac_const_cvd(j, :)  = w_cvd(1, :);
end

fprintf('\n=== Fracción en restricción (media sobre agentes) ===\n');
fprintf('Informal (base): %.4f | (covid): %.4f | Δ: %.4f\n', ...
    mean(frac_const_base(:,1)), mean(frac_const_cvd(:,1)), ...
    mean(frac_const_cvd(:,1) - frac_const_base(:,1)));
fprintf('Formal   (base): %.4f | (covid): %.4f | Δ: %.4f\n', ...
    mean(frac_const_base(:,2)), mean(frac_const_cvd(:,2)), ...
    mean(frac_const_cvd(:,2) - frac_const_base(:,2)));

%% --------------------------------------------------------------
%  Figuras: bandas (todos los agentes) Benchmark vs COVID
%% --------------------------------------------------------------
% Consumo - todos los agentes
figure('Name','Consumption Policies: Benchmark vs COVID','Color','w');
subplot(1,2,1); hold on;
for j = 1:n_agents
    plot(a_b, c_b{j}(:,1), 'Color', [1 0 0 0.35], 'LineWidth', 1.0); % base inf
    plot(a_b, c_b{j}(:,2), 'Color', [0 0 1 0.35], 'LineWidth', 1.0); % base form
end
title('Benchmark'); xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);

subplot(1,2,2); hold on;
for j = 1:n_agents
    plot(a_c, c_c{j}(:,1), 'r-', 'LineWidth', 1.2); % covid inf
    plot(a_c, c_c{j}(:,2), 'b-', 'LineWidth', 1.2); % covid form
end
title('COVID shock (incomes ↓, transfers level fixed)'); xlabel('Assets (a)'); grid on; xlim([1 5]);

lg = legend({'Informal','Formal'}, 'Location','SouthEast'); lg.Box = 'off';

% Distribución de riqueza - todos los agentes
figure('Name','Wealth Distributions: Benchmark vs COVID','Color','w');
subplot(1,2,1); hold on;
for j = 1:n_agents
    g1 = g_b{j}(:,1); g1(g1<1e-5) = NaN;
    g2 = g_b{j}(:,2); g2(g2<1e-5) = NaN;
    plot(a_b, g1, 'Color', [1 0 0 0.35], 'LineWidth', 1.0);
    plot(a_b, g2, 'Color', [0 0 1 0.35], 'LineWidth', 1.0);
end
title('Benchmark'); xlabel('Assets (a)'); ylabel('Wealth distribution'); grid on; xlim([0.01 0.60]);

subplot(1,2,2); hold on;
for j = 1:n_agents
    g1 = g_c{j}(:,1); g1(g1<1e-5) = NaN;
    g2 = g_c{j}(:,2); g2(g2<1e-5) = NaN;
    plot(a_c, g1, 'r-', 'LineWidth', 1.2);
    plot(a_c, g2, 'b-', 'LineWidth', 1.2);
end
title('COVID shock (incomes ↓, transfers level fixed)'); xlabel('Assets (a)'); grid on; xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast');

%% --------------------------------------------------------------
%  Ejemplo: un agente “típico” (mediana de sI)
%% --------------------------------------------------------------
j_star = ceil(n_agents/2);
figure('Name','Typical agent policies','Color','w');
subplot(1,2,1); hold on;
plot(a_b, c_b{j_star}(:,1), 'r--', 'LineWidth', 1.6);
plot(a_b, c_b{j_star}(:,2), 'b--', 'LineWidth', 1.6);
plot(a_c, c_c{j_star}(:,1), 'r-',  'LineWidth', 1.6);
plot(a_c, c_c{j_star}(:,2), 'b-',  'LineWidth', 1.6);
xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);
title(sprintf('Consumption | RRA_i = %.3f', sI_vector1(j_star)));
legend({'Inf Base','For Base','Inf COVID','For COVID'}, 'Location','SouthEast');

subplot(1,2,2); hold on;
plot(a_b, g_b{j_star}(:,1), 'r--', 'LineWidth', 1.6);
plot(a_b, g_b{j_star}(:,2), 'b--', 'LineWidth', 1.6);
plot(a_c, g_c{j_star}(:,1), 'r-',  'LineWidth', 1.6);
plot(a_c, g_c{j_star}(:,2), 'b-',  'LineWidth', 1.6);
xlabel('Assets (a)'); ylabel('Wealth distribution'); 
grid on; 
% xlim([0.01 0.60]);
xlim([0.1 1.60]);
title('Wealth distribution (typical agent)');
legend({'Inf Base','For Base','Inf COVID','For COVID'}, 'Location','NorthEast');

%% --------------------------------------------------------------
%  (Opcional) Exportes rápidos a Excel
%% --------------------------------------------------------------
export_excel = true;
if export_excel
    filename = 'summary_r_star.xlsx';
    T = table((1:n_agents)', r_b(:), r_c(:), 'VariableNames', {'agent','r_base','r_covid'});
    writetable(T, filename, 'Sheet', 'r_star');

    % fracción en restricción
    T2 = table((1:n_agents)', frac_const_base(:,1), frac_const_cvd(:,1), ...
        frac_const_base(:,2), frac_const_cvd(:,2), ...
        'VariableNames', {'agent','inf_base','inf_covid','for_base','for_covid'});
    writetable(T2, filename, 'Sheet', 'constraint_frac');
end

disp('Listo: Benchmark vs. COVID resueltos y graficados.');


%%------------------------------------------------------------------
% 
%%------------------------------------------------------------------
%% --------------------------------------------------------------
%  Estadísticos y Gini: riqueza y consumo (Base vs COVID)
%% --------------------------------------------------------------
%% --------------------------------------------------------------
%  Estadísticos + Gini por ocupación (Informal/Formal): Base vs COVID
%% --------------------------------------------------------------
filename = 'summary_r_star.xlsx';

% Tabla riqueza (activos) por ocupación
W_occ = table('Size',[n_agents 24], ...
    'VariableTypes', repmat("double",1,24), ...
    'VariableNames',{ ...
      'mean_a_inf_base','std_a_inf_base','p10_a_inf_base','p50_a_inf_base','p90_a_inf_base','gini_a_inf_base', ...
      'mean_a_for_base','std_a_for_base','p10_a_for_base','p50_a_for_base','p90_a_for_base','gini_a_for_base', ...
      'mean_a_inf_cvd','std_a_inf_cvd','p10_a_inf_cvd','p50_a_inf_cvd','p90_a_inf_cvd','gini_a_inf_cvd', ...
      'mean_a_for_cvd','std_a_for_cvd','p10_a_for_cvd','p50_a_for_cvd','p90_a_for_cvd','gini_a_for_cvd'});

% Tabla consumo por ocupación
C_occ = table('Size',[n_agents 20], ...
    'VariableTypes', repmat("double",1,20), ...
    'VariableNames',{ ...
      'mean_c_inf_base','p10_c_inf_base','p50_c_inf_base','p90_c_inf_base','gini_c_inf_base', ...
      'mean_c_for_base','p10_c_for_base','p50_c_for_base','p90_c_for_base','gini_c_for_base', ...
      'mean_c_inf_cvd','p10_c_inf_cvd','p50_c_inf_cvd','p90_c_inf_cvd','gini_c_inf_cvd', ...
      'mean_c_for_cvd','p10_c_for_cvd','p50_c_for_cvd','p90_c_for_cvd','gini_c_for_cvd'});

for j = 1:n_agents
    % Malla y paso
    a_grid = a_b; 
    da_loc = a_grid(2) - a_grid(1);

    %% ====== BASE: Densidades por ocupación (normalizadas por ocupación) ======
    G_base = g_b{j};                         % I x 2   (col1=Informal, col2=Formal)
    % Normalización condicional por ocupación
    w_inf_b = G_base(:,1); s_inf_b = sum(w_inf_b)*da_loc; if s_inf_b>0, w_inf_b = w_inf_b/(s_inf_b); else, w_inf_b(:)=0; end
    w_for_b = G_base(:,2); s_for_b = sum(w_for_b)*da_loc; if s_for_b>0, w_for_b = w_for_b/(s_for_b); else, w_for_b(:)=0; end

    % Riqueza por ocupación (Base)
    mean_a_inf_b = sum(a_grid .* w_inf_b) * da_loc;
    std_a_inf_b  = sqrt( max(0, sum(((a_grid-mean_a_inf_b).^2).*w_inf_b)*da_loc) );
    p10_a_inf_b  = wquantile(a_grid, w_inf_b, 0.10, da_loc);
    p50_a_inf_b  = wquantile(a_grid, w_inf_b, 0.50, da_loc);
    p90_a_inf_b  = wquantile(a_grid, w_inf_b, 0.90, da_loc);
    gini_a_inf_b = gini_weighted(a_grid, w_inf_b, da_loc);

    mean_a_for_b = sum(a_grid .* w_for_b) * da_loc;
    std_a_for_b  = sqrt( max(0, sum(((a_grid-mean_a_for_b).^2).*w_for_b)*da_loc) );
    p10_a_for_b  = wquantile(a_grid, w_for_b, 0.10, da_loc);
    p50_a_for_b  = wquantile(a_grid, w_for_b, 0.50, da_loc);
    p90_a_for_b  = wquantile(a_grid, w_for_b, 0.90, da_loc);
    gini_a_for_b = gini_weighted(a_grid, w_for_b, da_loc);

    % Consumo por ocupación (Base): pesos condicionales por ocupación
    C_base = c_b{j};                          % I x 2
    % Informales
    x_c_inf_b = C_base(:,1);
    w_c_inf_b = G_base(:,1) * da_loc; 
    if sum(w_c_inf_b)>0, w_c_inf_b = w_c_inf_b/sum(w_c_inf_b); end
    vb = isfinite(x_c_inf_b) & isfinite(w_c_inf_b) & (w_c_inf_b>0);
    x_c_inf_b = x_c_inf_b(vb); w_c_inf_b = w_c_inf_b(vb);
    mean_c_inf_b = sum(x_c_inf_b .* w_c_inf_b);
    p10_c_inf_b  = wquantile_discrete(x_c_inf_b, w_c_inf_b, 0.10);
    p50_c_inf_b  = wquantile_discrete(x_c_inf_b, w_c_inf_b, 0.50);
    p90_c_inf_b  = wquantile_discrete(x_c_inf_b, w_c_inf_b, 0.90);
    gini_c_inf_b = gini_weighted_discrete(x_c_inf_b, w_c_inf_b);

    % Formales
    x_c_for_b = C_base(:,2);
    w_c_for_b = G_base(:,2) * da_loc; 
    if sum(w_c_for_b)>0, w_c_for_b = w_c_for_b/sum(w_c_for_b); end
    vb2 = isfinite(x_c_for_b) & isfinite(w_c_for_b) & (w_c_for_b>0);
    x_c_for_b = x_c_for_b(vb2); w_c_for_b = w_c_for_b(vb2);
    mean_c_for_b = sum(x_c_for_b .* w_c_for_b);
    p10_c_for_b  = wquantile_discrete(x_c_for_b, w_c_for_b, 0.10);
    p50_c_for_b  = wquantile_discrete(x_c_for_b, w_c_for_b, 0.50);
    p90_c_for_b  = wquantile_discrete(x_c_for_b, w_c_for_b, 0.90);
    gini_c_for_b = gini_weighted_discrete(x_c_for_b, w_c_for_b);

    %% ====== COVID: Densidades por ocupación (normalizadas por ocupación) ======
    G_cvd = g_c{j};
    w_inf_c = G_cvd(:,1); s_inf_c = sum(w_inf_c)*da_loc; if s_inf_c>0, w_inf_c = w_inf_c/(s_inf_c); else, w_inf_c(:)=0; end
    w_for_c = G_cvd(:,2); s_for_c = sum(w_for_c)*da_loc; if s_for_c>0, w_for_c = w_for_c/(s_for_c); else, w_for_c(:)=0; end

    % Riqueza por ocupación (COVID)
    mean_a_inf_c = sum(a_grid .* w_inf_c) * da_loc;
    std_a_inf_c  = sqrt( max(0, sum(((a_grid-mean_a_inf_c).^2).*w_inf_c)*da_loc) );
    p10_a_inf_c  = wquantile(a_grid, w_inf_c, 0.10, da_loc);
    p50_a_inf_c  = wquantile(a_grid, w_inf_c, 0.50, da_loc);
    p90_a_inf_c  = wquantile(a_grid, w_inf_c, 0.90, da_loc);
    gini_a_inf_c = gini_weighted(a_grid, w_inf_c, da_loc);

    mean_a_for_c = sum(a_grid .* w_for_c) * da_loc;
    std_a_for_c  = sqrt( max(0, sum(((a_grid-mean_a_for_c).^2).*w_for_c)*da_loc) );
    p10_a_for_c  = wquantile(a_grid, w_for_c, 0.10, da_loc);
    p50_a_for_c  = wquantile(a_grid, w_for_c, 0.50, da_loc);
    p90_a_for_c  = wquantile(a_grid, w_for_c, 0.90, da_loc);
    gini_a_for_c = gini_weighted(a_grid, w_for_c, da_loc);

    % Consumo por ocupación (COVID)
    C_cvd = c_c{j};
    % Informales
    x_c_inf_c = C_cvd(:,1);
    w_c_inf_c = G_cvd(:,1) * da_loc; 
    if sum(w_c_inf_c)>0, w_c_inf_c = w_c_inf_c/sum(w_c_inf_c); end
    vc = isfinite(x_c_inf_c) & isfinite(w_c_inf_c) & (w_c_inf_c>0);
    x_c_inf_c = x_c_inf_c(vc); w_c_inf_c = w_c_inf_c(vc);
    mean_c_inf_c = sum(x_c_inf_c .* w_c_inf_c);
    p10_c_inf_c  = wquantile_discrete(x_c_inf_c, w_c_inf_c, 0.10);
    p50_c_inf_c  = wquantile_discrete(x_c_inf_c, w_c_inf_c, 0.50);
    p90_c_inf_c  = wquantile_discrete(x_c_inf_c, w_c_inf_c, 0.90);
    gini_c_inf_c = gini_weighted_discrete(x_c_inf_c, w_c_inf_c);

    % Formales
    x_c_for_c = C_cvd(:,2);
    w_c_for_c = G_cvd(:,2) * da_loc; 
    if sum(w_c_for_c)>0, w_c_for_c = w_c_for_c/sum(w_c_for_c); end
    vc2 = isfinite(x_c_for_c) & isfinite(w_c_for_c) & (w_c_for_c>0);
    x_c_for_c = x_c_for_c(vc2); w_c_for_c = w_c_for_c(vc2);
    mean_c_for_c = sum(x_c_for_c .* w_c_for_c);
    p10_c_for_c  = wquantile_discrete(x_c_for_c, w_c_for_c, 0.10);
    p50_c_for_c  = wquantile_discrete(x_c_for_c, w_c_for_c, 0.50);
    p90_c_for_c  = wquantile_discrete(x_c_for_c, w_c_for_c, 0.90);
    gini_c_for_c = gini_weighted_discrete(x_c_for_c, w_c_for_c);

    % Guardar filas
    W_occ{j,:} = [ ...
        mean_a_inf_b, std_a_inf_b, p10_a_inf_b, p50_a_inf_b, p90_a_inf_b, gini_a_inf_b, ...
        mean_a_for_b, std_a_for_b, p10_a_for_b, p50_a_for_b, p90_a_for_b, gini_a_for_b, ...
        mean_a_inf_c, std_a_inf_c, p10_a_inf_c, p50_a_inf_c, p90_a_inf_c, gini_a_inf_c, ...
        mean_a_for_c, std_a_for_c, p10_a_for_c, p50_a_for_c, p90_a_for_c, gini_a_for_c];

    C_occ{j,:} = [ ...
        mean_c_inf_b, p10_c_inf_b, p50_c_inf_b, p90_c_inf_b, gini_c_inf_b, ...
        mean_c_for_b, p10_c_for_b, p50_c_for_b, p90_c_for_b, gini_c_for_b, ...
        mean_c_inf_c, p10_c_inf_c, p50_c_inf_c, p90_c_inf_c, gini_c_inf_c, ...
        mean_c_for_c, p10_c_for_c, p50_c_for_c, p90_c_for_c, gini_c_for_c];
end

% Añadir identificadores (agente y RRA informal)
W_occ = addvars(W_occ, (1:n_agents)', sI_vector1(:), 'Before', 1, ...
                'NewVariableNames', {'agent','RRA_informal'});
C_occ = addvars(C_occ, (1:n_agents)', sI_vector1(:), 'Before', 1, ...
                'NewVariableNames', {'agent','RRA_informal'});

% Exportar
writetable(W_occ, filename, 'Sheet', 'stats_wealth_occ');
writetable(C_occ, filename, 'Sheet', 'stats_consumption_occ');

% --------- Resumen promedios por ocupación/escenario ---------
summary_occ = table;
summary_occ.metric = { ...
    'mean_a_inf','std_a_inf','p10_a_inf','p50_a_inf','p90_a_inf','gini_a_inf', ...
    'mean_a_for','std_a_for','p10_a_for','p50_a_for','p90_a_for','gini_a_for', ...
    'mean_c_inf','p10_c_inf','p50_c_inf','p90_c_inf','gini_c_inf', ...
    'mean_c_for','p10_c_for','p50_c_for','p90_c_for','gini_c_for' }';

summary_occ.base  = [ ...
    mean(W_occ.mean_a_inf_base); mean(W_occ.std_a_inf_base);  mean(W_occ.p10_a_inf_base); mean(W_occ.p50_a_inf_base); mean(W_occ.p90_a_inf_base); mean(W_occ.gini_a_inf_base); ...
    mean(W_occ.mean_a_for_base); mean(W_occ.std_a_for_base);  mean(W_occ.p10_a_for_base); mean(W_occ.p50_a_for_base); mean(W_occ.p90_a_for_base); mean(W_occ.gini_a_for_base); ...
    mean(C_occ.mean_c_inf_base); mean(C_occ.p10_c_inf_base);  mean(C_occ.p50_c_inf_base); mean(C_occ.p90_c_inf_base); mean(C_occ.gini_c_inf_base); ...
    mean(C_occ.mean_c_for_base); mean(C_occ.p10_c_for_base);  mean(C_occ.p50_c_for_base); mean(C_occ.p90_c_for_base); mean(C_occ.gini_c_for_base) ...
    ];

summary_occ.covid = [ ...
    mean(W_occ.mean_a_inf_cvd); mean(W_occ.std_a_inf_cvd);  mean(W_occ.p10_a_inf_cvd); mean(W_occ.p50_a_inf_cvd); mean(W_occ.p90_a_inf_cvd); mean(W_occ.gini_a_inf_cvd); ...
    mean(W_occ.mean_a_for_cvd); mean(W_occ.std_a_for_cvd);  mean(W_occ.p10_a_for_cvd); mean(W_occ.p50_a_for_cvd); mean(W_occ.p90_a_for_cvd); mean(W_occ.gini_a_for_cvd); ...
    mean(C_occ.mean_c_inf_cvd); mean(C_occ.p10_c_inf_cvd);  mean(C_occ.p50_c_inf_cvd); mean(C_occ.p90_c_inf_cvd); mean(C_occ.gini_c_inf_cvd); ...
    mean(C_occ.mean_c_for_cvd); mean(C_occ.p10_c_for_cvd);  mean(C_occ.p50_c_for_cvd); mean(C_occ.p90_c_for_cvd); mean(C_occ.gini_c_for_cvd) ...
    ];

summary_occ.delta = summary_occ.covid - summary_occ.base;

writetable(summary_occ, filename, 'Sheet', 'summary_overall_occ');

disp('Exportado por ocupación: stats_wealth_occ, stats_consumption_occ y summary_overall_occ.');
