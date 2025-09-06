%% ==============================================================
%  main_high_low_informality_benchmark.m
%  Benchmark (no COVID): comparar baja vs alta informalidad
%  Autor: your_name
%  Fecha: today
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Configuracion
%% -----------------------
n_agents  = 20;           % numero de puntos RRA informal
s_min     = 3.15;         % RRA informal min
s_max     = 5.30;         % RRA informal max

sI_vector1 = linspace(s_min, s_max, n_agents);   % RRA informal heterogenea
sF_vector1 = 5.30 * ones(1, n_agents);           % RRA formal constante

% Niveles de informalidad a comparar
eta_low  = 0.20;  % baja informalidad
eta_high = 0.80;  % alta informalidad
eta_vec_low  = eta_low  * ones(1, n_agents);
eta_vec_high = eta_high * ones(1, n_agents);

% Configuracion baseline (sin shock)
cfg_base = struct();
cfg_base.scenario              = 'baseline';
cfg_base.psi_I                 = 1.00;      % sin shock
cfg_base.psi_F                 = 1.00;      % sin shock
cfg_base.keep_transfers_level  = true;      % Transfer = phi * z1_base
cfg_base.transfer_multiplier   = 1.00;      % kappa
cfg_base.amin_mode             = 'baseline';
cfg_base.gov_nonneg            = true;
cfg_base.phi                   = 0.13;
cfg_base.taxF                  = 0.10;
cfg_base.taxI                  = 0.00;
cfg_base.tauc                  = 0.00;      % no usado en presupuesto aqui
cfg_base.theta                 = 0.02;
cfg_base.r0                    = 0.03;
cfg_base.rmin                  = 0.01;
cfg_base.rmax                  = 0.04;

%% -----------------------
%  Resolver: Benchmark baja y alta informalidad
%% -----------------------
[r_b_L, ir_b_L, pop_b_L, a_L, g_b_L, c_b_L] = ...
    huggett_Equi_RRA_function_transfer_base(eta_vec_low,  sI_vector1, sF_vector1, cfg_base);

[r_b_H, ir_b_H, pop_b_H, a_H, g_b_H, c_b_H] = ...
    huggett_Equi_RRA_function_transfer_base(eta_vec_high, sI_vector1, sF_vector1, cfg_base);

%% -----------------------
%  Estadisticos
%% -----------------------
da_L = (a_L(end) - a_L(1)) / (numel(a_L)-1);
da_H = (a_H(end) - a_H(1)) / (numel(a_H)-1);

cmean_inf_b_L = zeros(n_agents,1); cmean_for_b_L = zeros(n_agents,1);
cmean_inf_b_H = zeros(n_agents,1); cmean_for_b_H = zeros(n_agents,1);

frac_bound_b_L = zeros(n_agents,2);
frac_bound_b_H = zeros(n_agents,2);

% L2 distancia entre distribuciones informales High vs Low (benchmark)
L2_HvsL = zeros(n_agents,1);

for j = 1:n_agents
    % Low
    gBL = g_b_L{j}; cBL = c_b_L{j};
    g1 = gBL(:,1); g2 = gBL(:,2); g1(g1<0)=0; g2(g2<0)=0;
    w1 = g1 / (sum(g1)*da_L); w2 = g2 / (sum(g2)*da_L);
    cmean_inf_b_L(j) = sum(w1 .* cBL(:,1)) * da_L;
    cmean_for_b_L(j) = sum(w2 .* cBL(:,2)) * da_L;
    frac_bound_b_L(j,:) = (gBL(1,:) ./ sum(gBL(:)))';

    % High
    gBH = g_b_H{j}; cBH = c_b_H{j};
    g1 = gBH(:,1); g2 = gBH(:,2); g1(g1<0)=0; g2(g2<0)=0;
    w1 = g1 / (sum(g1)*da_H); w2 = g2 / (sum(g2)*da_H);
    cmean_inf_b_H(j) = sum(w1 .* cBH(:,1)) * da_H;
    cmean_for_b_H(j) = sum(w2 .* cBH(:,2)) * da_H;
    frac_bound_b_H(j,:) = (gBH(1,:) ./ sum(gBH(:)))';

    % L2 informal High vs Low (sobre la misma malla; asumimos a_H==a_L)
    gL = g_b_L{j}(:,1); gH = g_b_H{j}(:,1);
    gL(gL<1e-5)=0; gH(gH<1e-5)=0;
    % si las mallas difieren, interpolar gH en a_L
    if numel(a_H) == numel(a_L) && max(abs(a_H - a_L)) < 1e-12
        L2_HvsL(j) = sqrt(trapz(a_L, (gH - gL).^2));
    else
        gH_on_L = interp1(a_H, gH, a_L, 'linear', 'extrap');
        L2_HvsL(j) = sqrt(trapz(a_L, (gH_on_L - gL).^2));
    end
end

%% -----------------------
%  Reporte agregado
%% -----------------------
fprintf('\n=== r* (promedio entre agentes) ===\n');
fprintf('Low informality  | r*: %.6f\n', mean(r_b_L));
fprintf('High informality | r*: %.6f\n', mean(r_b_H));
fprintf('Delta (High - Low): %.6f\n', mean(r_b_H - r_b_L));

fprintf('\n=== c_mean (promedio entre agentes) Informal / Formal ===\n');
fprintf('Low  base: inf=%.6f  for=%.6f\n', mean(cmean_inf_b_L), mean(cmean_for_b_L));
fprintf('High base: inf=%.6f  for=%.6f\n', mean(cmean_inf_b_H), mean(cmean_for_b_H));

fprintf('\n=== Fraction at borrowing limit (promedio) ===\n');
fprintf('Low  base: inf=%.4f  for=%.4f\n', mean(frac_bound_b_L(:,1)), mean(frac_bound_b_L(:,2)));
fprintf('High base: inf=%.4f  for=%.4f\n', mean(frac_bound_b_H(:,1)), mean(frac_bound_b_H(:,2)));

fprintf('\n=== L2 distance of g_inf(a): High vs Low (benchmark) ===\n');
fprintf('Mean L2 across agents: %.6e\n', mean(L2_HvsL));

%% -----------------------
%  Graficos: bandas Low vs High (Benchmark)
%% -----------------------
% Consumo: Low vs High (dos paneles)
figure('Name','Benchmark: Consumption Policies (Low vs High)','Color','w');
subplot(1,2,1); hold on;
for j=1:n_agents
    plot(a_L, c_b_L{j}(:,1), 'Color', [1 0 0 0.35]); % informal
    plot(a_L, c_b_L{j}(:,2), 'Color', [0 0 1 0.35]); % formal
end
title('Low informality'); xlabel('Assets (a)'); ylabel('Consumption'); grid on; xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast');

subplot(1,2,2); hold on;
for j=1:n_agents
    plot(a_H, c_b_H{j}(:,1), 'r-', 'LineWidth', 1.2);
    plot(a_H, c_b_H{j}(:,2), 'b-', 'LineWidth', 1.2);
end
title('High informality'); xlabel('Assets (a)'); grid on; xlim([1 5]);

% Distribucion de riqueza: Low vs High (dos paneles)
figure('Name','Benchmark: Wealth Distributions (Low vs High)','Color','w');
subplot(1,2,1); hold on;
for j=1:n_agents
    g1 = g_b_L{j}(:,1); g2 = g_b_L{j}(:,2); g1(g1<1e-5)=NaN; g2(g2<1e-5)=NaN;
    plot(a_L, g1, 'Color', [1 0 0 0.35]);
    plot(a_L, g2, 'Color', [0 0 1 0.35]);
end
title('Low informality'); xlabel('Assets (a)'); ylabel('Wealth distribution'); grid on; xlim([0.01 0.60]);

subplot(1,2,2); hold on;
for j=1:n_agents
    g1 = g_b_H{j}(:,1); g2 = g_b_H{j}(:,2); g1(g1<1e-5)=NaN; g2(g2<1e-5)=NaN;
    plot(a_H, g1, 'r-', 'LineWidth', 1.2);
    plot(a_H, g2, 'b-', 'LineWidth', 1.2);
end
title('High informality'); xlabel('Assets (a)'); grid on; xlim([0.01 0.60]);

%% -----------------------
%  Figura 2x2: agente tipico (Low vs High)
%% -----------------------
[~, j_star] = min(abs(sI_vector1 - median(sI_vector1)));  % agente tipico

% Low
cLb = c_b_L{j_star}; gLb = g_b_L{j_star};
gLb_plot = gLb; gLb_plot(gLb_plot<1e-5)=NaN;

% High
cHb = c_b_H{j_star}; gHb = g_b_H{j_star};
gHb_plot = gHb; gHb_plot(gHb_plot<1e-5)=NaN;

figure('Name','Typical agent: Low vs High informality (Benchmark)','Color','w','Position',[60 60 1200 720]);

% (1) Low - Consumo
subplot(2,2,1); hold on; box on; grid on;
plot(a_L, cLb(:,1), 'r-','LineWidth',1.8); plot(a_L, cLb(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Consumption'); title(sprintf('LOW: Cons | RRA_inf=%.3f', sI_vector1(j_star))); xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast'); legend boxoff;

% (2) Low - Distribucion
subplot(2,2,2); hold on; box on; grid on;
plot(a_L, gLb_plot(:,1), 'r-','LineWidth',1.8); plot(a_L, gLb_plot(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Wealth distribution'); title('LOW: Wealth dist'); xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast'); legend boxoff;
xline(a_L(1), ':k', 'a_min', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);

% (3) High - Consumo
subplot(2,2,3); hold on; box on; grid on;
plot(a_H, cHb(:,1), 'r-','LineWidth',1.8); plot(a_H, cHb(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Consumption'); title(sprintf('HIGH: Cons | RRA_inf=%.3f', sI_vector1(j_star))); xlim([1 5]);
legend({'Informal','Formal'}, 'Location','SouthEast'); legend boxoff;

% (4) High - Distribucion
subplot(2,2,4); hold on; box on; grid on;
plot(a_H, gHb_plot(:,1), 'r-','LineWidth',1.8); plot(a_H, gHb_plot(:,2), 'b-','LineWidth',1.8);
xlabel('Assets (a)'); ylabel('Wealth distribution'); title('HIGH: Wealth dist'); xlim([0.01 0.60]);
legend({'Informal','Formal'}, 'Location','NorthEast'); legend boxoff;
xline(a_H(1), ':k', 'a_min', 'LabelVerticalAlignment','bottom', 'Alpha',0.6);

%% -----------------------
%  Exportes a Excel (APAGADOS por defecto)
%% -----------------------
export_excel = false;   % cambiar a true si deseas exportar
if export_excel
    filename = 'benchmark_high_low_informality.xlsx';

    % r*
    T_r = table((1:n_agents)', r_b_L(:), r_b_H(:), ...
        'VariableNames', {'agent','r_base_low','r_base_high'});
    writetable(T_r, filename, 'Sheet', 'r_star');

    % c_mean
    T_c = table((1:n_agents)', ...
        cmean_inf_b_L, cmean_for_b_L, cmean_inf_b_H, cmean_for_b_H, ...
        'VariableNames', {'agent', ...
        'cmean_inf_base_low','cmean_for_base_low', ...
        'cmean_inf_base_high','cmean_for_base_high'});
    writetable(T_c, filename, 'Sheet', 'c_means');

    % fraction at constraint
    T_fb = table((1:n_agents)', ...
        frac_bound_b_L(:,1), frac_bound_b_L(:,2), ...
        frac_bound_b_H(:,1), frac_bound_b_H(:,2), ...
        'VariableNames', {'agent', ...
        'fb_inf_base_low','fb_for_base_low', ...
        'fb_inf_base_high','fb_for_base_high'});
    writetable(T_fb, filename, 'Sheet', 'frac_bound');

    % L2 High vs Low
    T_L2 = table((1:n_agents)', L2_HvsL(:), ...
        'VariableNames', {'agent','L2_high_vs_low'});
    writetable(T_L2, filename, 'Sheet', 'L2_dist');

    fprintf('Exportado a: %s\n', filename);
end

disp('Listo: Benchmark sin COVID, comparacion baja vs alta informalidad.');
