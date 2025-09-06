%% ==============================================================
%  main_scenarios_transfer_taxes.m
%  Escenarios con transferencias + impuestos:
%   (S0) Base: transferencias a informales, impuesto a la renta (solo formales),
%        e impuesto al consumo (ambos)
%   (S1) Transferencias aumentan
%   (S2) Impuesto a la renta (formales) aumenta
%   (S3) Impuesto al consumo (ambos) aumenta
%  Graficos + estadisticos (incluye Gini). Exportes OFF por defecto.
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Parametrizacion comun
%% -----------------------
eta = 0.75;                 % informalidad (fija)
n_agents  = 20;             % heterogeneidad en RRA informal
s_min     = 0.15;
s_max     = 0.30;
sI_vector1 = linspace(s_min, s_max, n_agents);   % RRA informal heterogenea
sF_vector1 = 0.30 * ones(1, n_agents);           % RRA formal constante
eta_vector = eta * ones(1, n_agents);

% Config base del solucionador
base_cfg = struct();
base_cfg.psi_I = 1.00;                 % sin choque de ingreso
base_cfg.psi_F = 1.00;
base_cfg.taxI  = 0.00;                 % ingreso (informal) - 0 por defecto
base_cfg.taxF  = 0.10;                 % ingreso (formal) 10% en BASE
base_cfg.tauc  = 0.08;                 % IVA/consumo 8% en BASE
base_cfg.phi   = 0.10;                 % transferencia a informales = phi*z1_base
base_cfg.transfer_mode = 'level_base'; % 'level_base' | 'proportional'
base_cfg.keep_transfers_level = true;  % transferencia en nivel sobre z1_base
base_cfg.gshare_public = 0.00;         % aqui no usamos bien publico
base_cfg.gov_nonneg = true;            % trunca gasto publico >= 0
base_cfg.theta = 0.02;                 % prima de endeudamiento informal
base_cfg.r0 = 0.03; base_cfg.rmin=0.01; base_cfg.rmax=0.04;
base_cfg.amin_mode = 'baseline';       % limite de deuda segun benchmark

%% -----------------------
%  Escenarios
%% -----------------------
% (S0) Base
cfg0 = base_cfg;

% (S1) Transferencias aumentan
cfg1 = base_cfg; 
cfg1.phi  = 0.18;    % subir transferencias (nivel sobre z1_base)

% (S2) Impuesto a la renta (formales) aumenta
cfg2 = base_cfg;
cfg2.taxF = 0.20;    % subir impuesto a renta formal

% (S3) Impuesto al consumo (ambos) aumenta
cfg3 = base_cfg;
cfg3.tauc = 0.15;    % subir IVA/consumo

%% -----------------------
%  Resolver cada escenario
%% -----------------------
[r0, ~, pop0, a0, g0, c0, stats0] = ...
    huggett_Equi_RRA_function_transfer_taxes(eta_vector, sI_vector1, sF_vector1, cfg0);

[r1, ~, pop1, a1, g1, c1, stats1] = ...
    huggett_Equi_RRA_function_transfer_taxes(eta_vector, sI_vector1, sF_vector1, cfg1);

[r2, ~, pop2, a2, g2, c2, stats2] = ...
    huggett_Equi_RRA_function_transfer_taxes(eta_vector, sI_vector1, sF_vector1, cfg2);

[r3, ~, pop3, a3, g3, c3, stats3] = ...
    huggett_Equi_RRA_function_transfer_taxes(eta_vector, sI_vector1, sF_vector1, cfg3);

% malla (misma por construccion)
a = a0;

%% -----------------------
%  Reporte sintetico en consola
%% -----------------------
fprintf('\n==== Informalidad fija = 0.75 | Resumen por escenario ====\n');
fprintf('r* mean:     S0=%.6f  S1=%.6f  S2=%.6f  S3=%.6f\n', ...
    mean(r0), mean(r1), mean(r2), mean(r3));
fprintf('c_mean inf:  S0=%.6f  S1=%.6f  S2=%.6f  S3=%.6f\n', ...
    mean(stats0.cmean_inf), mean(stats1.cmean_inf), mean(stats2.cmean_inf), mean(stats3.cmean_inf));
fprintf('c_mean for:  S0=%.6f  S1=%.6f  S2=%.6f  S3=%.6f\n', ...
    mean(stats0.cmean_for), mean(stats1.cmean_for), mean(stats2.cmean_for), mean(stats3.cmean_for));
fprintf('Frac bound inf (mean): S0=%.4f  S1=%.4f  S2=%.4f  S3=%.4f\n', ...
    mean(stats0.frac_bound(:,1)), mean(stats1.frac_bound(:,1)), ...
    mean(stats2.frac_bound(:,1)), mean(stats3.frac_bound(:,1)));
fprintf('Gini wealth (inf):     S0=%.4f  S1=%.4f  S2=%.4f  S3=%.4f\n', ...
    mean(stats0.gini_inf), mean(stats1.gini_inf), mean(stats2.gini_inf), mean(stats3.gini_inf));
fprintf('Gini wealth (total):   S0=%.4f  S1=%.4f  S2=%.4f  S3=%.4f\n\n', ...
    stats0.gini_total, stats1.gini_total, stats2.gini_total, stats3.gini_total);

%% -----------------------
%  Figura: bandas por escenario (consumo)
%% -----------------------
figure('Name','Policies: Consumption across scenarios','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% S0
nexttile; hold on; title('S0: Base'); grid on; xlim([1 5]);
for j=1:n_agents
    plot(a, c0{j}(:,1), 'r-'); % informal
    plot(a, c0{j}(:,2), 'b-'); % formal
end
xlabel('Assets (a)'); ylabel('Consumption');

% S1
nexttile; hold on; title('S1: Transfers \uparrow'); grid on; xlim([1 5]);
for j=1:n_agents
    plot(a, c1{j}(:,1), 'r-'); plot(a, c1{j}(:,2), 'b-');
end
xlabel('Assets (a)');

% S2
nexttile; hold on; title('S2: Income tax (formal) \uparrow'); grid on; xlim([1 5]);
for j=1:n_agents
    plot(a, c2{j}(:,1), 'r-'); plot(a, c2{j}(:,2), 'b-');
end
xlabel('Assets (a)'); ylabel('Consumption');

% S3
nexttile; hold on; title('S3: Consumption tax \uparrow'); grid on; xlim([1 5]);
for j=1:n_agents
    plot(a, c3{j}(:,1), 'r-'); plot(a, c3{j}(:,2), 'b-');
end
xlabel('Assets (a)');

legend({'Informal','Formal'},'Location','SouthOutside','NumColumns',2);

%% -----------------------
%  Figura: bandas por escenario (distribucion de riqueza)
%% -----------------------
figure('Name','Distributions: Wealth across scenarios','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

mask_small = @(x) (x .* (x>=1e-5) + NaN*(x<1e-5));

% S0
nexttile; hold on; title('S0: Base'); grid on; xlim([0.01 0.60]);
for j=1:n_agents
    gi = mask_small(g0{j}(:,1)); gf = mask_small(g0{j}(:,2));
    plot(a, gi, 'r-'); plot(a, gf, 'b-');
end
xlabel('Assets (a)'); ylabel('Wealth dist');

% S1
nexttile; hold on; title('S1: Transfers \uparrow'); grid on; xlim([0.01 0.60]);
for j=1:n_agents
    gi = mask_small(g1{j}(:,1)); gf = mask_small(g1{j}(:,2));
    plot(a, gi, 'r-'); plot(a, gf, 'b-');
end
xlabel('Assets (a)');

% S2
nexttile; hold on; title('S2: Income tax (formal) \uparrow'); grid on; xlim([0.01 0.60]);
for j=1:n_agents
    gi = mask_small(g2{j}(:,1)); gf = mask_small(g2{j}(:,2));
    plot(a, gi, 'r-'); plot(a, gf, 'b-');
end
xlabel('Assets (a)'); ylabel('Wealth dist');

% S3
nexttile; hold on; title('S3: Consumption tax \uparrow'); grid on; xlim([0.01 0.60]);
for j=1:n_agents
    gi = mask_small(g3{j}(:,1)); gf = mask_small(g3{j}(:,2));
    plot(a, gi, 'r-'); plot(a, gf, 'b-');
end
xlabel('Assets (a)');

legend({'Informal','Formal'},'Location','SouthOutside','NumColumns',2);

%% -----------------------
%  Figura: Agente tipico (Base vs cambios)
%% -----------------------
[~, j_star] = min(abs(sI_vector1 - median(sI_vector1)));

figure('Name','Typical agent: Base vs policy changes','Color','w','Position',[60 60 1200 900]);

% Base vs Transfers up
subplot(3,2,1); hold on; box on; grid on; xlim([1 5]);
plot(a, c0{j_star}(:,1), 'r--','LineWidth',1.4);
plot(a, c1{j_star}(:,1), 'r-','LineWidth',1.6);
plot(a, c0{j_star}(:,2), 'b--','LineWidth',1.4);
plot(a, c1{j_star}(:,2), 'b-','LineWidth',1.6);
title('Consumption: Base vs Transfers \uparrow'); xlabel('Assets (a)'); ylabel('c(a)');
legend({'Inf Base','Inf S1','For Base','For S1'},'Location','SouthEast');

subplot(3,2,2); hold on; box on; grid on; xlim([0.01 0.60]);
plot(a, mask_small(g0{j_star}(:,1)), 'r--','LineWidth',1.4);
plot(a, mask_small(g1{j_star}(:,1)), 'r-','LineWidth',1.6);
plot(a, mask_small(g0{j_star}(:,2)), 'b--','LineWidth',1.4);
plot(a, mask_small(g1{j_star}(:,2)), 'b-','LineWidth',1.6);
title('Wealth: Base vs Transfers \uparrow'); xlabel('Assets (a)'); ylabel('g(a)');

% Base vs Income tax up
subplot(3,2,3); hold on; box on; grid on; xlim([1 5]);
plot(a, c0{j_star}(:,1), 'r--','LineWidth',1.4);
plot(a, c2{j_star}(:,1), 'r-','LineWidth',1.6);
plot(a, c0{j_star}(:,2), 'b--','LineWidth',1.4);
plot(a, c2{j_star}(:,2), 'b-','LineWidth',1.6);
title('Consumption: Base vs Income tax \uparrow'); xlabel('Assets (a)'); ylabel('c(a)');
legend({'Inf Base','Inf S2','For Base','For S2'},'Location','SouthEast');

subplot(3,2,4); hold on; box on; grid on; xlim([0.01 0.60]);
plot(a, mask_small(g0{j_star}(:,1)), 'r--','LineWidth',1.4);
plot(a, mask_small(g2{j_star}(:,1)), 'r-','LineWidth',1.6);
plot(a, mask_small(g0{j_star}(:,2)), 'b--','LineWidth',1.4);
plot(a, mask_small(g2{j_star}(:,2)), 'b-','LineWidth',1.6);
title('Wealth: Base vs Income tax \uparrow'); xlabel('Assets (a)'); ylabel('g(a)');

% Base vs Consumption tax up
subplot(3,2,5); hold on; box on; grid on; xlim([1 5]);
plot(a, c0{j_star}(:,1), 'r--','LineWidth',1.4);
plot(a, c3{j_star}(:,1), 'r-','LineWidth',1.6);
plot(a, c0{j_star}(:,2), 'b--','LineWidth',1.4);
plot(a, c3{j_star}(:,2), 'b-','LineWidth',1.6);
title('Consumption: Base vs VAT \uparrow'); xlabel('Assets (a)'); ylabel('c(a)');
legend({'Inf Base','Inf S3','For Base','For S3'},'Location','SouthEast');

subplot(3,2,6); hold on; box on; grid on; xlim([0.01 0.60]);
plot(a, mask_small(g0{j_star}(:,1)), 'r--','LineWidth',1.4);
plot(a, mask_small(g3{j_star}(:,1)), 'r-','LineWidth',1.6);
plot(a, mask_small(g0{j_star}(:,2)), 'b--','LineWidth',1.4);
plot(a, mask_small(g3{j_star}(:,2)), 'b-','LineWidth',1.6);
title('Wealth: Base vs VAT \uparrow'); xlabel('Assets (a)'); ylabel('g(a)');

%% -----------------------
%  Tabla de estadisticos y exportes (OFF)
%% -----------------------
Summary = table( (1:n_agents)', ...
    r0(:), r1(:), r2(:), r3(:), ...
    stats0.cmean_inf, stats1.cmean_inf, stats2.cmean_inf, stats3.cmean_inf, ...
    stats0.cmean_for, stats1.cmean_for, stats2.cmean_for, stats3.cmean_for, ...
    stats0.frac_bound(:,1), stats1.frac_bound(:,1), stats2.frac_bound(:,1), stats3.frac_bound(:,1), ...
    stats0.frac_bound(:,2), stats1.frac_bound(:,2), stats2.frac_bound(:,2), stats3.frac_bound(:,2), ...
    stats0.gini_inf, stats1.gini_inf, stats2.gini_inf, stats3.gini_inf, ...
    'VariableNames', { ...
      'agent', ...
      'r_S0','r_S1','r_S2','r_S3', ...
      'cmean_inf_S0','cmean_inf_S1','cmean_inf_S2','cmean_inf_S3', ...
      'cmean_for_S0','cmean_for_S1','cmean_for_S2','cmean_for_S3', ...
      'fbound_inf_S0','fbound_inf_S1','fbound_inf_S2','fbound_inf_S3', ...
      'fbound_for_S0','fbound_for_S1','fbound_for_S2','fbound_for_S3', ...
      'gini_inf_S0','gini_inf_S1','gini_inf_S2','gini_inf_S3' ...
    });

Gini_total = table( ...
    stats0.gini_total, stats1.gini_total, stats2.gini_total, stats3.gini_total, ...
    'VariableNames', {'GiniTot_S0','GiniTot_S1','GiniTot_S2','GiniTot_S3'});

Fiscal = table( ...
    stats0.rev_income, stats1.rev_income, stats2.rev_income, stats3.rev_income, ...
    stats0.rev_tauc,   stats1.rev_tauc,   stats2.rev_tauc,   stats3.rev_tauc, ...
    stats0.transfers,  stats1.transfers,  stats2.transfers,  stats3.transfers, ...
    stats0.gov_pub,    stats1.gov_pub,    stats2.gov_pub,    stats3.gov_pub, ...
    stats0.gov_surplus,stats1.gov_surplus,stats2.gov_surplus,stats3.gov_surplus, ...
    'VariableNames', {'IncomeTax_S0','IncomeTax_S1','IncomeTax_S2','IncomeTax_S3', ...
                      'VAT_S0','VAT_S1','VAT_S2','VAT_S3', ...
                      'Transfers_S0','Transfers_S1','Transfers_S2','Transfers_S3', ...
                      'PubGood_S0','PubGood_S1','PubGood_S2','PubGood_S3', ...
                      'Surplus_S0','Surplus_S1','Surplus_S2','Surplus_S3'});

disp('--- Gini total por escenario ---'); disp(Gini_total);

export_excel = false; % cambiar a true si deseas exportar
if export_excel
    fname = 'scenarios_transfer_taxes.xlsx';
    writetable(Summary,     fname, 'Sheet', 'by_agent');
    writetable(Gini_total,  fname, 'Sheet', 'gini_total');
    writetable(Fiscal,      fname, 'Sheet', 'fiscal');
    fprintf('Exportado a: %s\n', fname);
end

disp('Listo: S0..S3 resueltos, graficas y estadisticos generados.');
