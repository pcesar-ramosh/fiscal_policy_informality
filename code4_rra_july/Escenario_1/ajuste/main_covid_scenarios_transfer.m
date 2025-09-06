%% ==============================================================
%  main_covid_scenarios_transfer.m
%  Benchmark vs. COVID shock (incomes down; transfers go UP)
%  Autor: PhD. Hamilton Galindo Gil
%  Contribuidores: Ian Carrasco and Cesar Ramos
%  Fecha: Septiembre 2025 (versión corregida)
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
%  Escenario 1: Benchmark
%% -----------------------
cfg_base = struct();
cfg_base.scenario              = "baseline";
cfg_base.psi_I                 = 1.00;        % sin shock
cfg_base.psi_F                 = 1.00;        % sin shock
cfg_base.transfer_multiplier   = 1.00;        % κ=1 (sin aumento)
cfg_base.keep_transfers_level  = true;        % nivel (irrelevante en baseline)
cfg_base.amin_mode             = "baseline";  % límite inferior con z1_base
cfg_base.gov_nonneg            = true;        % (no usado en utilidad ahora)
cfg_base.phi                   = 0.13;        % regla de transferencias base
cfg_base.taxF                  = 0.10;        % impuesto renta formal
cfg_base.taxI                  = 0.00;        % impuesto renta informal
cfg_base.tauc                  = 0.00;        % impuesto al consumo
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
cfg_covid.keep_transfers_level = true;        % mantener NIVEL basado en z1_base
cfg_covid.amin_mode            = "baseline";  % comparativa pura del límite

[r_c, ir_c, pop_c, a_c, g_c, c_c] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vector, sI_vector1, sF_vector1, cfg_covid);

%% --------------------------------------------------------------
%  Comparaciones clave (r* y fracción en restricción condicional)
%% --------------------------------------------------------------
fprintf('\n=== Equilibrium interest rate (mean across agents) ===\n');
fprintf('Benchmark r*: %.6f\n', mean(r_b));
fprintf('COVID+Transfers↑ r*: %.6f\n', mean(r_c));
fprintf('Δr (covid_up - base): %.6f\n', mean(r_c - r_b));

% Fracción en a_min condicional por ocupación (densidad cond. integra a 1)
frac_const_base = zeros(n_agents, 2);
frac_const_cvd  = zeros(n_agents, 2);

da_b = a_b(2) - a_b(1);
da_c = a_c(2) - a_c(1);

for j = 1:n_agents
    Gb = g_b{j}; Gc = g_c{j};

    mass_inf_b = sum(Gb(:,1))*da_b;  mass_for_b = sum(Gb(:,2))*da_b;
    mass_inf_c = sum(Gc(:,1))*da_c;  mass_for_c = sum(Gc(:,2))*da_c;

    frac_const_base(j,1) = (Gb(1,1)*da_b)/max(mass_inf_b, eps);
    frac_const_base(j,2) = (Gb(1,2)*da_b)/max(mass_for_b, eps);
    frac_const_cvd(j,1)  = (Gc(1,1)*da_c)/max(mass_inf_c, eps);
    frac_const_cvd(j,2)  = (Gc(1,2)*da_c)/max(mass_for_c, eps);
end

fprintf('\n=== Fraction at borrowing limit (mean, conditional by occupation) ===\n');
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

%% --------------------------------------------------------------
%  ONE FIGURE: Typical agent (before vs after) — consumption & wealth
%% --------------------------------------------------------------
[~, j_star] = min(abs(sI_vector1 - median(sI_vector1)));
a_base = a_b;  a_covid = a_c;
c_base = c_b{j_star};  c_covid = c_c{j_star};
g_base = g_b{j_star};  g_covid = g_c{j_star};

g_base_plot  = g_base;   g_base_plot(g_base_plot < 1e-5)   = NaN;
g_covid_plot = g_covid;  g_covid_plot(g_covid_plot < 1e-5) = NaN;

xlim_c  = [1 5];
xlim_g  = [0.01 1.60];

figure('Name','Typical agent: Before vs After','Color','w','Position',[100 100 1100 430]);
subplot(1,2,1); hold on; box on; grid on;
p1 = plot(a_base, c_base(:,1), 'r--', 'LineWidth', 1.8);
p2 = plot(a_base, c_base(:,2), 'b--', 'LineWidth', 1.8);
p3 = plot(a_covid, c_covid(:,1), 'r-',  'LineWidth', 1.8);
p4 = plot(a_covid, c_covid(:,2), 'b-',  'LineWidth', 1.8);
xlabel('Assets (a)'); ylabel('Consumption');
title(sprintf('Consumption | Typical agent (RRA_{inf}=%.3f)', sI_vector1(j_star)));
xlim(xlim_c);
leg1 = legend([p1 p2 p3 p4], {'Inf Base','For Base','Inf COVID+T↑','For COVID+T↑'}, 'Location','SouthEast');
set(leg1,'Box','off');

subplot(1,2,2); hold on; box on; grid on;
q1 = plot(a_base,  g_base_plot(:,1),  'r--', 'LineWidth', 1.8);
q2 = plot(a_base,  g_base_plot(:,2),  'b--', 'LineWidth', 1.8);
q3 = plot(a_covid, g_covid_plot(:,1), 'r-',  'LineWidth', 1.8);
q4 = plot(a_covid, g_covid_plot(:,2), 'b-',  'LineWidth', 1.8);
xlabel('Assets (a)'); ylabel('Wealth distribution');
title('Wealth distribution (Typical agent)');
xlim(xlim_g);
leg2 = legend([q1 q2 q3 q4], {'Inf Base','For Base','Inf COVID+T↑','For COVID+T↑'}, 'Location','NorthEast');
set(leg2,'Box','off');
xline(a_base(1), ':k', 'a_{min}', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);

%% --------------------------------------------------------------
%  Exportes a Excel — ESTADÍSTICOS & GINI por ocupación (Inf/For)
%% --------------------------------------------------------------
export_excel = true;   % <- activa exportación
if export_excel
    filename = 'summary_covid_uptransfer.xlsx';

    % ---------- r* por agente ----------
    T_r = table((1:n_agents)', r_b(:), r_c(:), 'VariableNames', ...
        {'agent','r_base','r_covid_Tup'});
    writetable(T_r, filename, 'Sheet', 'r_star');

    % ---------- fracción en restricción ----------
    T_fc = table((1:n_agents)', ...
        frac_const_base(:,1), frac_const_cvd(:,1), ...
        frac_const_base(:,2), frac_const_cvd(:,2), ...
        'VariableNames', {'agent','inf_base','inf_covid_Tup','for_base','for_covid_Tup'});
    writetable(T_fc, filename, 'Sheet', 'constraint_frac');

    % ---------- Estadísticos por ocupación (activos/consumo) ----------
    W_occ = table('Size',[n_agents 24], ...
        'VariableTypes', repmat("double",1,24), ...
        'VariableNames',{ ...
          'mean_a_inf_base','std_a_inf_base','p10_a_inf_base','p50_a_inf_base','p90_a_inf_base','gini_a_inf_base', ...
          'mean_a_for_base','std_a_for_base','p10_a_for_base','p50_a_for_base','p90_a_for_base','gini_a_for_base', ...
          'mean_a_inf_cvd','std_a_inf_cvd','p10_a_inf_cvd','p50_a_inf_cvd','p90_a_inf_cvd','gini_a_inf_cvd', ...
          'mean_a_for_cvd','std_a_for_cvd','p10_a_for_cvd','p50_a_for_cvd','p90_a_for_cvd','gini_a_for_cvd'});

    C_occ = table('Size',[n_agents 20], ...
        'VariableTypes', repmat("double",1,20), ...
        'VariableNames',{ ...
          'mean_c_inf_base','p10_c_inf_base','p50_c_inf_base','p90_c_inf_base','gini_c_inf_base', ...
          'mean_c_for_base','p10_c_for_base','p50_c_for_base','p90_c_for_base','gini_c_for_base', ...
          'mean_c_inf_cvd','p10_c_inf_cvd','p50_c_inf_cvd','p90_c_inf_cvd','gini_c_inf_cvd', ...
          'mean_c_for_cvd','p10_c_for_cvd','p50_c_for_cvd','p90_c_for_cvd','gini_c_for_cvd'});

    a_grid = a_b;  da_loc = a_grid(2) - a_grid(1);

    for j = 1:n_agents
        % ===== BASE =====
        G_base = g_b{j};   C_base = c_b{j};

        % (i) Riqueza por ocupación (normalización condicional)
        w_inf_b = G_base(:,1); s_inf_b = sum(w_inf_b)*da_loc;
        if s_inf_b>0, w_inf_b = w_inf_b/s_inf_b; else, w_inf_b(:)=0; end

        w_for_b = G_base(:,2); s_for_b = sum(w_for_b)*da_loc;
        if s_for_b>0, w_for_b = w_for_b/s_for_b; else, w_for_b(:)=0; end

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

        % (ii) Consumo por ocupación — pesos condicionales
        % Informales
        x_c_inf_b = C_base(:,1);
        w_c_inf_b = G_base(:,1) * da_loc;
        if sum(w_c_inf_b)>0, w_c_inf_b = w_c_inf_b / sum(w_c_inf_b); end
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
        if sum(w_c_for_b)>0, w_c_for_b = w_c_for_b / sum(w_c_for_b); end
        vb2 = isfinite(x_c_for_b) & isfinite(w_c_for_b) & (w_c_for_b>0);
        x_c_for_b = x_c_for_b(vb2); w_c_for_b = w_c_for_b(vb2);
        mean_c_for_b = sum(x_c_for_b .* w_c_for_b);
        p10_c_for_b  = wquantile_discrete(x_c_for_b, w_c_for_b, 0.10);
        p50_c_for_b  = wquantile_discrete(x_c_for_b, w_c_for_b, 0.50);
        p90_c_for_b  = wquantile_discrete(x_c_for_b, w_c_for_b, 0.90);
        gini_c_for_b = gini_weighted_discrete(x_c_for_b, w_c_for_b);

        % ===== COVID + TRANSFERS UP =====
        G_cvd = g_c{j};   C_cvd = c_c{j};

        w_inf_c = G_cvd(:,1); s_inf_c = sum(w_inf_c)*da_loc;
        if s_inf_c>0, w_inf_c = w_inf_c/s_inf_c; else, w_inf_c(:)=0; end

        w_for_c = G_cvd(:,2); s_for_c = sum(w_for_c)*da_loc;
        if s_for_c>0, w_for_c = w_for_c/s_for_c; else, w_for_c(:)=0; end

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

        % Consumo
        x_c_inf_c = C_cvd(:,1);
        w_c_inf_c = G_cvd(:,1) * da_loc;
        if sum(w_c_inf_c)>0, w_c_inf_c = w_c_inf_c / sum(w_c_inf_c); end
        vc = isfinite(x_c_inf_c) & isfinite(w_c_inf_c) & (w_c_inf_c>0);
        x_c_inf_c = x_c_inf_c(vc); w_c_inf_c = w_c_inf_c(vc);
        mean_c_inf_c = sum(x_c_inf_c .* w_c_inf_c);
        p10_c_inf_c  = wquantile_discrete(x_c_inf_c, w_c_inf_c, 0.10);
        p50_c_inf_c  = wquantile_discrete(x_c_inf_c, w_c_inf_c, 0.50);
        p90_c_inf_c  = wquantile_discrete(x_c_inf_c, w_c_inf_c, 0.90);
        gini_c_inf_c = gini_weighted_discrete(x_c_inf_c, w_c_inf_c);

        x_c_for_c = C_cvd(:,2);
        w_c_for_c = G_cvd(:,2) * da_loc;
        if sum(w_c_for_c)>0, w_c_for_c = w_c_for_c / sum(w_c_for_c); end
        vc2 = isfinite(x_c_for_c) & isfinite(w_c_for_c) & (w_c_for_c>0);
        x_c_for_c = x_c_for_c(vc2); w_c_for_c = w_c_for_c(vc2);
        mean_c_for_c = sum(x_c_for_c .* w_c_for_c);
        p10_c_for_c  = wquantile_discrete(x_c_for_c, w_c_for_c, 0.10);
        p50_c_for_c  = wquantile_discrete(x_c_for_c, w_c_for_c, 0.50);
        p90_c_for_c  = wquantile_discrete(x_c_for_c, w_c_for_c, 0.90);
        gini_c_for_c = gini_weighted_discrete(x_c_for_c, w_c_for_c);

        % ----- Guardar filas -----
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

    W_occ = addvars(W_occ, (1:n_agents)', sI_vector1(:), 'Before', 1, ...
                    'NewVariableNames', {'agent','RRA_informal'});
    C_occ = addvars(C_occ, (1:n_agents)', sI_vector1(:), 'Before', 1, ...
                    'NewVariableNames', {'agent','RRA_informal'});

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

    summary_occ.covidTup = [ ...
        mean(W_occ.mean_a_inf_cvd); mean(W_occ.std_a_inf_cvd);  mean(W_occ.p10_a_inf_cvd); mean(W_occ.p50_a_inf_cvd); mean(W_occ.p90_a_inf_cvd); mean(W_occ.gini_a_inf_cvd); ...
        mean(W_occ.mean_a_for_cvd); mean(W_occ.std_a_for_cvd);  mean(W_occ.p10_a_for_cvd); mean(W_occ.p50_a_for_cvd); mean(W_occ.p90_a_for_cvd); mean(W_occ.gini_a_for_cvd); ...
        mean(C_occ.mean_c_inf_cvd); mean(C_occ.p10_c_inf_cvd);  mean(C_occ.p50_c_inf_cvd); mean(C_occ.p90_c_inf_cvd); mean(C_occ.gini_c_inf_cvd); ...
        mean(C_occ.mean_c_for_cvd); mean(C_occ.p10_c_for_cvd);  mean(C_occ.p50_c_for_cvd); mean(C_occ.p90_c_for_cvd); mean(C_occ.gini_c_for_cvd) ...
        ];

    summary_occ.delta = summary_occ.covidTup - summary_occ.base;

    writetable(summary_occ, filename, 'Sheet', 'summary_overall_occ');

    fprintf('Exportado a: %s\n', filename);
end

disp('Estadísticos + Gini por ocupación (Base vs COVID+T↑) exportados al Excel.');
