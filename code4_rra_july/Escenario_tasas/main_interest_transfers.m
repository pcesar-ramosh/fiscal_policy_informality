%% ==============================================================
%  main_interest_transfers.m  (versión corregida)
%  Analiza r* con transferencias bajo 4 experimentos y grafica:
%  (A) posiciones medias de activos (informal/formal) y (B) r*
%% ==============================================================

clear; clc; close all;

%% Heterogeneidad en RRA
n_agents  = 20;
s_min     = 3.15;
s_max     = 5.30;
sI_vec0   = linspace(s_min, s_max, n_agents);   % RRA informal
sF_vec0   = 5.30 * ones(1, n_agents);           % RRA formal

%% Configuración base
cfg_base = struct();
cfg_base.scenario             = "baseline";
cfg_base.psi_I                = 1.00;
cfg_base.psi_F                = 1.00;
cfg_base.keep_transfers_level = true;
cfg_base.transfer_multiplier  = 1.00;
cfg_base.amin_mode            = "baseline";
cfg_base.gov_nonneg           = true;
cfg_base.phi                  = 0.13;
cfg_base.taxF                 = 0.10;
cfg_base.taxI                 = 0.00;
cfg_base.tauc                 = 0.00;
cfg_base.theta                = 0.02;
cfg_base.r0                   = 0.03;
cfg_base.rmin                 = 0.01;
cfg_base.rmax                 = 0.04;

% (opcional) configuración COVID
cfg_covid = cfg_base;
cfg_covid.scenario            = "covid";
cfg_covid.psi_I               = 0.75;
cfg_covid.psi_F               = 0.85;
cfg_covid.keep_transfers_level= true;

%% ==============================================================
%  EXPERIMENTO 1: r* vs tamaño de informalidad (eta)
%% ==============================================================

eta_grid = 0.20:0.10:0.90;
r_eta = zeros(size(eta_grid));
Ainf_eta = zeros(size(eta_grid));
Afor_eta = zeros(size(eta_grid));

for k = 1:numel(eta_grid)
    eta_k = eta_grid(k) * ones(1, n_agents);
    [r_opt,~,~,a, g_opt, ~] = ...
        huggett_Equi_RRA_function_transfer_on(eta_k, sI_vec0, sF_vec0, cfg_base);
    [rbar, Ain, Afor] = summarize(a, g_opt, r_opt);
    r_eta(k)    = rbar;
    Ainf_eta(k) = Ain;
    Afor_eta(k) = Afor;
end

plot_three(eta_grid, Ainf_eta, Afor_eta, r_eta, 'Informality Size (\eta)', ...
    'EXP 1: r* vs \eta (con transferencias)');

%% ==============================================================
%  EXPERIMENTO 2: r* vs impuesto informal (tau_I), eta=0.64
%% ==============================================================

eta_fix = 0.64 * ones(1, n_agents);
tauI_grid = linspace(0.00, 0.18, 10);

r_tau = zeros(size(tauI_grid));
Ainf_tau = zeros(size(tauI_grid));
Afor_tau = zeros(size(tauI_grid));

for k = 1:numel(tauI_grid)
    cfg = cfg_base; 
    cfg.taxI = tauI_grid(k);
    [r_opt,~,~,a, g_opt, ~] = ...
        huggett_Equi_RRA_function_transfer_on(eta_fix, sI_vec0, sF_vec0, cfg);
    [rbar, Ain, Afor] = summarize(a, g_opt, r_opt);
    r_tau(k)    = rbar;
    Ainf_tau(k) = Ain;
    Afor_tau(k) = Afor;
end

plot_three(tauI_grid, Ainf_tau, Afor_tau, r_tau, 'Tax Rate Informal (\tau_I)', ...
    'EXP 2: r* vs \tau_I (con transferencias)');

%% ==============================================================
%  EXPERIMENTO 3: r* vs prima de endeudamiento (theta), eta=0.64
%% ==============================================================

theta_grid = linspace(0.00, 0.02, 11);

r_th = zeros(size(theta_grid));
Ainf_th = zeros(size(theta_grid));
Afor_th = zeros(size(theta_grid));

for k = 1:numel(theta_grid)
    cfg = cfg_base; 
    cfg.taxI  = 0.00;
    cfg.theta = theta_grid(k);
    [r_opt,~,~,a, g_opt, ~] = ...
        huggett_Equi_RRA_function_transfer_on(eta_fix, sI_vec0, sF_vec0, cfg);
    [rbar, Ain, Afor] = summarize(a, g_opt, r_opt);
    r_th(k)    = rbar;
    Ainf_th(k) = Ain;
    Afor_th(k) = Afor;
end

% r y r+theta
figure('Color','w','Name','EXP 3: r* vs \theta');
subplot(1,3,1); hold on; box on; grid on;
plot(theta_grid, Ainf_th,'--','LineWidth',1.8);
plot(theta_grid, Afor_th,'-','LineWidth',1.8);
xlabel('Interest Rate Premium (\theta)','Interpreter','tex');
ylabel('Mean assets (a)'); title('Bond Market (A)');
legend({'Informal','Formal'}, 'Location','best');

subplot(1,3,2); hold on; box on; grid on;
plot(theta_grid, Afor_th,'--','LineWidth',1.8);
plot(theta_grid, Ainf_th,'-','LineWidth',1.8);
xlabel('\theta','Interpreter','tex'); ylabel('Mean assets (a)');
title('Bond Market (B)'); legend({'Formal','Informal'}, 'Location','best');

subplot(1,3,3); hold on; box on; grid on;
plot(theta_grid, 100*r_th,'-','LineWidth',2.0, 'DisplayName','Formal (r)');
plot(theta_grid, 100*(r_th + theta_grid),':','LineWidth',2.0, 'DisplayName','Informal (r+\theta)');
xlabel('Interest Rate Premium (\theta)','Interpreter','tex');
ylabel('Equilibrium Interest Rate (%)'); title('r*'); legend('Location','best');

%% ==============================================================
%  EXPERIMENTO 4: r* vs RRA_{inf} (más riesgo en informales)
%% ==============================================================

gamma_mult_grid = [1.0 1.25 1.5 2.0];

r_gm = zeros(size(gamma_mult_grid));
Ainf_gm = zeros(size(gamma_mult_grid));
Afor_gm = zeros(size(gamma_mult_grid));

for k = 1:numel(gamma_mult_grid)
    sI_vec = min(8.0, gamma_mult_grid(k) * sI_vec0);  % cap por estabilidad
    cfg = cfg_base; 
    cfg.taxI  = 0.00; 
    cfg.theta = 0.02;
    [r_opt,~,~,a, g_opt, ~] = ...
        huggett_Equi_RRA_function_transfer_on(eta_fix, sI_vec, sF_vec0, cfg);
    [rbar, Ain, Afor] = summarize(a, g_opt, r_opt);
    r_gm(k)    = rbar;
    Ainf_gm(k) = Ain;
    Afor_gm(k) = Afor;
end

plot_three(gamma_mult_grid, Ainf_gm, Afor_gm, r_gm, ...
    'RRA multiplier for informal (\gamma_I/\gamma_{I,0})', ...
    'EXP 4: r* vs RRA_{inf} (\theta=0.02, \tau_I=0, \eta=0.64)');

disp('Listo. Cuatro experimentos generados (A/B markets y r*).');

%% ===================== FUNCIONES LOCALES ======================

function [rbar, Ain, Afor] = summarize(a, gcell, rvec)
    rbar = mean(rvec);
    Ain  = mean(cellfun(@(g) trapz(a, g(:,1).*a), gcell)); % ā informal
    Afor = mean(cellfun(@(g) trapz(a, g(:,2).*a), gcell)); % ā formal
end

function plot_three(x, Ain, Afor, rstar, xlab, ttl)
    figure('Color','w','Name',ttl);

    % (i) Bond Market A
    subplot(1,3,1); hold on; box on; grid on;
    plot(x, Ain,'--','LineWidth',1.8);
    plot(x, Afor,'-','LineWidth',1.8);
    xlabel(xlab,'Interpreter','tex'); ylabel('Mean assets (a)');
    title('Bond Market (A)');
    legend({'Asset Position (Informal)','Asset Position (Formal)'}, 'Location','best');

    % (ii) Bond Market B
    subplot(1,3,2); hold on; box on; grid on;
    plot(x, Afor,'--','LineWidth',1.8);
    plot(x, Ain,'-','LineWidth',1.8);
    xlabel(xlab,'Interpreter','tex'); ylabel('Mean assets (a)');
    title('Bond Market (B)');
    legend({'Formal','Informal'}, 'Location','best');

    % (iii) r*
    subplot(1,3,3); hold on; box on; grid on;
    plot(x, 100*rstar,'-','LineWidth',2.0);
    xlabel(xlab,'Interpreter','tex'); ylabel('Equilibrium Interest Rate (%)');
    title('r*');
end
