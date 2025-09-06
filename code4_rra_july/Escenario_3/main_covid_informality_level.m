%% ==============================================================
%  main_covid_informality_level.m
%  COVID shock bajo baja vs. alta informalidad (eta = 0.20 vs 0.80)
%  Compara Benchmark vs. COVID para cada nivel de informalidad
%  Autor: PhD. Hamilton Galindo Gil
% Contribuidores: Ian Carrasco and Cesar Ramos
%  Fecha: Agosto 2025
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Configuración general
%% -----------------------
n_agents  = 20;                 % # de puntos para RRA informal
s_min     = 3.15;               % RRA informal min
s_max     = 5.30;               % RRA informal max

sI_vector1 = linspace(s_min, s_max, n_agents);     % RRA informal heterogénea
sF_vector1 = 5.30 * ones(1, n_agents);             % RRA formal constante

% Niveles de informalidad a comparar
eta_low  = 0.20;   % baja informalidad
eta_high = 0.80;   % alta informalidad

% Construye vectores eta para todos los agentes (cada solve usa un único eta)
eta_vec_low  = eta_low  * ones(1, n_agents);
eta_vec_high = eta_high * ones(1, n_agents);

% Config base (política fiscal y numérica)
cfg_base = struct();
cfg_base.scenario              = "baseline";
cfg_base.psi_I                 = 1.00;        % sin shock
cfg_base.psi_F                 = 1.00;        % sin shock
cfg_base.keep_transfers_level  = true;        % transfer = phi*z1_base
cfg_base.transfer_multiplier   = 1.00;        % kappa (no se usa en baseline)
cfg_base.amin_mode             = "baseline";  % límite inferior con z1_base
cfg_base.gov_nonneg            = true;        % trunca bien público a >=0
cfg_base.phi                   = 0.13;        % regla de transferencias base
cfg_base.taxF                  = 0.10;        % impuesto renta formal
cfg_base.taxI                  = 0.00;        % impuesto renta informal
cfg_base.tauc                  = 0.00;        % impuesto al consumo (no usado)
cfg_base.theta                 = 0.02;        % prima de endeudamiento informal
cfg_base.r0                    = 0.03;        % bracket r*
cfg_base.rmin                  = 0.01;
cfg_base.rmax                  = 0.04;

% Config COVID (caída de ingresos; transferencias mantienen NIVEL base)
cfg_covid = cfg_base;
cfg_covid.scenario             = "covid";   % etiqueta informativa
cfg_covid.psi_I                = 0.75;      % caída 25% en ingreso informal
cfg_covid.psi_F                = 0.85;      % caída 15% en ingreso formal
cfg_covid.transfer_multiplier  = 1.00;      % nivel base (no se aumenta aquí)
cfg_covid.keep_transfers_level = true;      % monto basado en z1_base
cfg_covid.amin_mode            = "baseline";% (o "shocked" si quieres endurecer límite)

%% -----------------------
%  (A) Resolver Benchmark: baja y alta informalidad
%% -----------------------
[r_b_L, ir_b_L, pop_b_L, a_L, g_b_L, c_b_L] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vec_low,  sI_vector1, sF_vector1, cfg_base);

[r_b_H, ir_b_H, pop_b_H, a_H, g_b_H, c_b_H] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vec_high, sI_vector1, sF_vector1, cfg_base);

%% -----------------------
%  (B) Resolver COVID: baja y alta informalidad
%% -----------------------
[r_c_L, ir_c_L, pop_c_L, a_cL, g_c_L, c_c_L] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vec_low,  sI_vector1, sF_vector1, cfg_covid);

[r_c_H, ir_c_H, pop_c_H, a_cH, g_c_H, c_c_H] = ...
    huggett_Equi_RRA_function_transfer_on(eta_vec_high, sI_vector1, sF_vector1, cfg_covid);

%% -----------------------
%  Utilidades auxiliares para estadísticos
%% -----------------------
da_L = (a_L(end) - a_L(1)) / (numel(a_L)-1);
da_H = (a_H(end) - a_H(1)) / (numel(a_H)-1);

% Consumo medio (ponderado por g) por agente y grupo
cmean_inf_b_L = zeros(n_agents,1); cmean_for_b_L = zeros(n_agents,1);
cmean_inf_c_L = zeros(n_agents,1); cmean_for_c_L = zeros(n_agents,1);
cmean_inf_b_H = zeros(n_agents,1); cmean_for_b_H = zeros(n_agents,1);
cmean_inf_c_H = zeros(n_agents,1); cmean_for_c_H = zeros(n_agents,1);

% Fracción en restricción (aprox peso en a_min)
frac_bound_b_L = zeros(n_agents,2); frac_bound_c_L = zeros(n_agents,2);
frac_bound_b_H = zeros(n_agents,2); frac_bound_c_H = zeros(n_agents,2);

for j = 1:n_agents
    % LOW informality
    gBL = g_b_L{j}; gCL = g_c_L{j}; cBL = c_b_L{j}; cCL = c_c_L{j};
    g1 = gBL(:,1); g2 = gBL(:,2); g1(g1<0)=0; g2(g2<0)=0;
    w1 = g1 / (sum(g1)*da_L); w2 = g2 / (sum(g2)*da_L);
    cmean_inf_b_L(j) = sum(w1 .* cBL(:,1)) * da_L;
    cmean_for_b_L(j) = sum(w2 .* cBL(:,2)) * da_L;
    g1 = gCL(:,1); g2 = gCL(:,2); g1(g1<0)=0; g2(g2<0)=0;
    w1 = g1 / (sum(g1)*da_L); w2 = g2 / (sum(g2)*da_L);
    cmean_inf_c_L(j) = sum(w1 .* cCL(:,1)) * da_L;
    cmean_for_c_L(j) = sum(w2 .* cCL(:,2)) * da_L;
    % Fracción en restricción (primera fila de g / total)
    frac_bound_b_L(j,:) = (g_b_L{j}(1,:) ./ sum(g_b_L{j}(:)))';
    frac_bound_c_L(j,:) = (g_c_L{j}(1,:) ./ sum(g_c_L{j}(:)))';

    % HIGH informality
    gBH = g_b_H{j}; gCH = g_c_H{j}; cBH = c_b_H{j}; cCH = c_c_H{j};
    g1 = gBH(:,1); g2 = gBH(:,2); g1(g1<0)=0; g2(g2<0)=0;
    w1 = g1 / (sum(g1)*da_H); w2 = g2 / (sum(g2)*da_H);
    cmean_inf_b_H(j) = sum(w1 .* cBH(:,1)) * da_H;
    cmean_for_b_H(j) = sum(w2 .* cBH(:,2)) * da_H;
    g1 = gCH(:,1); g2 = gCH(:,2); g1(g1<0)=0; g2(g2<0)=0;
    w1 = g1 / (sum(g1)*da_H); w2 = g2 / (sum(g2)*da_H);
    cmean_inf_c_H(j) = sum(w1 .* cCH(:,1)) * da_H;
    cmean_for_c_H(j) = sum(w2 .* cCH(:,2)) * da_H;
    % Fracción en restricción
    frac_bound_b_H(j,:) = (g_b_H{j}(1,:) ./ sum(g_b_H{j}(:)))';
    frac_bound_c_H(j,:) = (g_c_H{j}(1,:) ./ sum(g_c_H{j}(:)))';
end

%% -----------------------
%  Distancia L2 de g_inf(a) vs Benchmark (por agente)
%% -----------------------
L2_L = zeros(n_agents,1); L2_H = zeros(n_agents,1);
for j = 1:n_agents
    gB = g_b_L{j}(:,1); gC = g_c_L{j}(:,1);
    gB(gB<1e-5)=0; gC(gC<1e-5)=0;
    L2_L(j) = sqrt(trapz(a_L, (gC - gB).^2));
    gB = g_b_H{j}(:,1); gC = g_c_H{j}(:,1);
    gB(gB<1e-5)=0; gC(gC<1e-5)=0;
    L2_H(j) = sqrt(trapz(a_H, (gC - gB).^2));
end

%% -----------------------
%  Reporte de estadísticos agregados
%% -----------------------
fprintf('\n=== r* (mean across agents) ===\n');
fprintf('Low informality  | base: %.6f | covid: %.6f | Δ: %.6f\n', mean(r_b_L), mean(r_c_L), mean(r_c_L - r_b_L));
fprintf('High informality | base: %.6f | covid: %.6f | Δ: %.6f\n', mean(r_b_H), mean(r_c_H), mean(r_c_H - r_b_H));

fprintf('\n=== c_mean (mean across agents) Informal / Formal ===\n');
fprintf('Low  base: inf=%.6f for=%.6f | covid: inf=%.6f for=%.6f\n', ...
    mean(cmean_inf_b_L), mean(cmean_for_b_L), mean(cmean_inf_c_L), mean(cmean_for_c_L));
fprintf('High base: inf=%.6f for=%.6f | covid: inf=%.6f for=%.6f\n', ...
    mean(cmean_inf_b_H), mean(cmean_for_b_H), mean(cmean_inf_c_H), mean(cmean_for_c_H));

fprintf('\n=== Fraction at borrowing limit (mean) ===\n');
fprintf('Low  base: inf=%.4f for=%.4f | covid: inf=%.4f for=%.4f | Δ_inf=%.4f\n', ...
    mean(frac_bound_b_L(:,1)), mean(frac_bound_b_L(:,2)), mean(frac_bound_c_L(:,1)), mean(frac_bound_c_L(:,2)), ...
    mean(frac_bound_c_L(:,1) - frac_bound_b_L(:,1)));
fprintf('High base: inf=%.4f for=%.4f | covid: inf=%.4f for=%.4f | Δ_inf=%.4f\n', ...
    mean(frac_bound_b_H(:,1)), mean(frac_bound_b_H(:,2)), mean(frac_bound_c_H(:,1)), mean(frac_bound_c_H(:,2)), ...
    mean(frac_bound_c_H(:,1) - frac_bound_b_H(:,1)));

fprintf('\n=== L2 distance of g_inf(a) vs Benchmark (mean across agents) ===\n');
fprintf('Low informality  : %.6e\n', mean(L2_L));
fprintf('High informality : %.6e\n', mean(L2_H));

%% --------------------------------------------------------------
%  Figuras (bandas) para cada informalidad
%% --------------------------------------------------------------
% LOW: Consumo
figure('Name','LOW informality: Consumption Policies','Color','w');
subplot(1,2,1); hold on;
for j=1:n_agents
    plot(a_L, c_b_L{j}(:,1), 'Color', [1 0 0 0.35]); % inf base
    plot(a_L, c_b_L{j}(:,2), 'Color', [0 0 1 0.35]); % for base
end
title('Low informality - Benchmark'); xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);
subplot(1,2,2); hold on;
for j=1:n_agents
    plot(a_cL, c_c_L{j}(:,1), 'r-', 'LineWidth', 1.2);
    plot(a_cL, c_c_L{j}(:,2), 'b-', 'LineWidth', 1.2);
end
title('Low informality - COVID'); xlabel('Assets (a)'); grid on; xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast');

% LOW: Distribución
figure('Name','LOW informality: Wealth Distributions','Color','w');
subplot(1,2,1); hold on;
for j=1:n_agents
    g1 = g_b_L{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_b_L{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_L, g1, 'Color', [1 0 0 0.35]);
    plot(a_L, g2, 'Color', [0 0 1 0.35]);
end
title('Low informality - Benchmark'); xlabel('Assets (a)'); ylabel('Wealth distribution'); grid on; xlim([0.01 0.60]);
subplot(1,2,2); hold on;
for j=1:n_agents
    g1 = g_c_L{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_c_L{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_cL, g1, 'r-', 'LineWidth', 1.2);
    plot(a_cL, g2, 'b-', 'LineWidth', 1.2);
end
title('Low informality - COVID'); xlabel('Assets (a)'); grid on; xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast');

% HIGH: Consumo
figure('Name','HIGH informality: Consumption Policies','Color','w');
subplot(1,2,1); hold on;
for j=1:n_agents
    plot(a_H, c_b_H{j}(:,1), 'Color', [1 0 0 0.35]);
    plot(a_H, c_b_H{j}(:,2), 'Color', [0 0 1 0.35]);
end
title('High informality - Benchmark'); xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);
subplot(1,2,2); hold on;
for j=1:n_agents
    plot(a_cH, c_c_H{j}(:,1), 'r-', 'LineWidth', 1.2);
    plot(a_cH, c_c_H{j}(:,2), 'b-', 'LineWidth', 1.2);
end
title('High informality - COVID'); xlabel('Assets (a)'); grid on; xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast');

% HIGH: Distribución
figure('Name','HIGH informality: Wealth Distributions','Color','w');
subplot(1,2,1); hold on;
for j=1:n_agents
    g1 = g_b_H{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_b_H{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_H, g1, 'Color', [1 0 0 0.35]);
    plot(a_H, g2, 'Color', [0 0 1 0.35]);
end
title('High informality - Benchmark'); xlabel('Assets (a)'); ylabel('Wealth distribution'); grid on; xlim([0.01 0.60]);
subplot(1,2,2); hold on;
for j=1:n_agents
    g1 = g_c_H{j}(:,1); g1(g1<1e-5)=NaN;
    g2 = g_c_H{j}(:,2); g2(g2<1e-5)=NaN;
    plot(a_cH, g1, 'r-', 'LineWidth', 1.2);
    plot(a_cH, g2, 'b-', 'LineWidth', 1.2);
end
title('High informality - COVID'); xlabel('Assets (a)'); grid on; xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast');

%% --------------------------------------------------------------
%  Figura 2x2: Agente típico (Low vs High) Benchmark vs COVID
%% --------------------------------------------------------------
[~, j_star] = min(abs(sI_vector1 - median(sI_vector1)));  % agente típico

% Low
cLb = c_b_L{j_star}; gLb = g_b_L{j_star}; cLc = c_c_L{j_star}; gLc = g_c_L{j_star};
gLb_plot = gLb; gLb_plot(gLb_plot<1e-5)=NaN;
gLc_plot = gLc; gLc_plot(gLc_plot<1e-5)=NaN;

% High
cHb = c_b_H{j_star}; gHb = g_b_H{j_star}; cHc = c_c_H{j_star}; gHc = g_c_H{j_star};
gHb_plot = gHb; gHb_plot(gHb_plot<1e-5)=NaN;
gHc_plot = gHc; gHc_plot(gHc_plot<1e-5)=NaN;

figure('Name','Typical agent: Low vs High informality','Color','w','Position',[60 60 1200 720]);

% (1) Low - Consumo
subplot(2,2,1); hold on; box on; grid on;
plot(a_L, cLb(:,1), 'r--','LineWidth',1.5); plot(a_L, cLb(:,2), 'b--','LineWidth',1.5);
plot(a_cL, cLc(:,1), 'r-','LineWidth',1.8); plot(a_cL, cLc(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Consumption'); title(sprintf('LOW: Cons | RRA_{inf}=%.3f', sI_vector1(j_star))); xlim([1 5]);
legend({'Inf Base','For Base','Inf COVID','For COVID'}, 'Location','SouthEast'); legend boxoff;

% (2) Low - Distribución
subplot(2,2,2); hold on; box on; grid on;
plot(a_L, gLb_plot(:,1), 'r--','LineWidth',1.5); plot(a_L, gLb_plot(:,2), 'b--','LineWidth',1.5);
plot(a_cL, gLc_plot(:,1), 'r-','LineWidth',1.8);  plot(a_cL, gLc_plot(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Wealth distribution'); title('LOW: Wealth dist'); xlim([0.01 0.60]);
legend({'Inf Base','For Base','Inf COVID','For COVID'}, 'Location','NorthEast'); legend boxoff;
xline(a_L(1), ':k', 'a_{min}', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);

% (3) High - Consumo
subplot(2,2,3); hold on; box on; grid on;
plot(a_H, cHb(:,1), 'r--','LineWidth',1.5); plot(a_H, cHb(:,2), 'b--','LineWidth',1.5);
plot(a_cH, cHc(:,1), 'r-','LineWidth',1.8); plot(a_cH, cHc(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Consumption'); title(sprintf('HIGH: Cons | RRA_{inf}=%.3f', sI_vector1(j_star))); xlim([1 5]);
legend({'Inf Base','For Base','Inf COVID','For COVID'}, 'Location','SouthEast'); legend boxoff;

% (4) High - Distribución
subplot(2,2,4); hold on; box on; grid on;
plot(a_H, gHb_plot(:,1), 'r--','LineWidth',1.5); plot(a_H, gHb_plot(:,2), 'b--','LineWidth',1.5);
plot(a_cH, gHc_plot(:,1), 'r-','LineWidth',1.8);  plot(a_cH, gHc_plot(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Wealth distribution'); title('HIGH: Wealth dist'); xlim([0.01 0.60]);
legend({'Inf Base','For Base','Inf COVID','For COVID'}, 'Location','NorthEast'); legend boxoff;
xline(a_H(1), ':k', 'a_{min}', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);

%% --------------------------------------------------------------
%  (Opcional) Exportes a Excel — APAGADOS por defecto
%% --------------------------------------------------------------
export_excel = false;   % <--- cambia a true si quieres exportar
if export_excel
    filename = 'covid_informality_levels.xlsx';

    % r*
    T_r = table((1:n_agents)', r_b_L(:), r_c_L(:), r_b_H(:), r_c_H(:), ...
        'VariableNames', {'agent','r_base_low','r_covid_low','r_base_high','r_covid_high'});
    writetable(T_r, filename, 'Sheet', 'r_star');

    % c_mean
    T_c = table((1:n_agents)', ...
        cmean_inf_b_L, cmean_for_b_L, cmean_inf_c_L, cmean_for_c_L, ...
        cmean_inf_b_H, cmean_for_b_H, cmean_inf_c_H, cmean_for_c_H, ...
        'VariableNames', {'agent', ...
        'cmean_inf_base_low','cmean_for_base_low','cmean_inf_covid_low','cmean_for_covid_low', ...
        'cmean_inf_base_high','cmean_for_base_high','cmean_inf_covid_high','cmean_for_covid_high'});
    writetable(T_c, filename, 'Sheet', 'c_means');

    % frac bound
    T_fb = table((1:n_agents)', ...
        frac_bound_b_L(:,1), frac_bound_b_L(:,2), frac_bound_c_L(:,1), frac_bound_c_L(:,2), ...
        frac_bound_b_H(:,1), frac_bound_b_H(:,2), frac_bound_c_H(:,1), frac_bound_c_H(:,2), ...
        'VariableNames', {'agent', ...
        'fb_inf_base_low','fb_for_base_low','fb_inf_covid_low','fb_for_covid_low', ...
        'fb_inf_base_high','fb_for_base_high','fb_inf_covid_high','fb_for_covid_high'});
    writetable(T_fb, filename, 'Sheet', 'frac_bound');

    % L2 de g_inf(a)
    T_L2 = table((1:n_agents)', L2_L(:), L2_H(:), ...
        'VariableNames', {'agent','L2_low','L2_high'});
    writetable(T_L2, filename, 'Sheet', 'L2_dist');

    fprintf('Exportado a: %s\n', filename);
end

disp('Listo: Benchmark vs COVID bajo baja y alta informalidad. Figuras y estadísticos generados.');
