%% lambda_dynamics_demo.m
% Demostración de dinámicas en intensidades de transición:
% - Factor cíclico x_t ~ AR(1)
% - λ2_t = λ2_bar * exp(sigma2 * x_t)
% - λ1_t = λ1_bar * exp(sigma1 * x_t)
% - Dependencia en riqueza (estado): mayor prob. pasar a informal si a<0
%
% En cada t resolvemos el equilibrio y graficamos series temporales.
% ---------------------------------------------------------------
clear; clc; close all;

T   = 8;            % nº de periodos (ajusta si quieres)
RRA = 4.20;         % preferencias
eta_seed = 0.64;    % solo como "semilla" (solver devolverá popI endógena)

% --- AR(1) para ciclo
rho_x = 0.75; sigma_eps = 0.8;
x = zeros(T,1);
for t=2:T
    x(t) = rho_x*x(t-1) + sigma_eps*randn(1);
end

% --- Intensidades base
p22_bar = 0.8155;
la2_bar = -log(p22_bar);
la1_bar = (1-eta_seed)*la2_bar/eta_seed;

% Sensibilidades a ciclo
sigma1 = 0.35;   % para λ1
sigma2 = 0.25;   % para λ2

% Dependencia en riqueza (estado) — parámetros
k_neg  = 0.8;    % multiplica la( a<0 )
k_pos  = 0.3;    % reduce la( a>0 ), suavemente

% Contenedores
r_t      = zeros(T,1);
eta_t    = zeros(T,1);
CmeanI_t = zeros(T,1); CmeanF_t = zeros(T,1);
SmeanI_t = zeros(T,1); SmeanF_t = zeros(T,1);
Gini_t   = zeros(T,1);
Tl_t     = zeros(T,1); Tc_t = zeros(T,1); Tr_t = zeros(T,1); Gx_t = zeros(T,1);
B_t      = zeros(T,1);

% Para distribuir gráficos de riqueza en tres cortes (t1, t2, t3)
snap_idx = unique(round(linspace(1,T,3)));

% Guarda g(a) y c(a) en snaps
snap = struct();

for t=1:T
    % --- λ_t con ciclo
    la1_t = la1_bar * exp(sigma1 * x(t));
    la2_t = la2_bar * exp(sigma2 * x(t));

    % --- definimos funciones estado-dependientes (vectoriza sobre 'a')
    lambdaSpec.mode = 'custom';
    lambdaSpec.la1_fun = @(a,eta_guess) ...
        max( la1_t .* ( 1 + k_neg*(a<0).*min(1, -a/abs(min(a))) ...
                          - k_pos*(a>0).*min(0.5, a/max(max(a),1e-6)) ), 1e-8);

    lambdaSpec.la2_fun = @(a,eta_guess) ...
        max( la2_t .* ( 1 + 0.5*k_neg*(a<0).*min(1, -a/abs(min(a))) ...
                          - 0.2*k_pos*(a>0).*min(0.5, a/max(max(a),1e-6)) ), 1e-8);

    % --- resolver equilibrio (eta_seed se ignora para lambda custom)
    [r, ~, pop1_vec, statsM, statsCM, GDist, a, Dist, ~, ~, C, Sav, ~, ~, I, amin] = ...
        huggett_Equi_INFO_function(RRA, eta_seed, lambdaSpec);

    da = (a(end)-a(1))/(numel(a)-1);

    % --- extraer variables
    r_t(t)   = r(1);
    eta_t(t) = pop1_vec(1);
    Gini_t(t)= statsM(1,9);

    % Consumo medio por tipo
    CmeanI_t(t) = statsCM(1,1);
    CmeanF_t(t) = statsCM(1,2);

    % Ahorro medio por tipo (promedio de s(a) ponderado por g(a))
    g = Dist{1};
    sI = Sav{1}(:,1); sF = Sav{1}(:,2);
    SmeanI_t(t) = sum(sI .* g(:,1))*da / max(sum(g(:,1))*da, 1e-12);
    SmeanF_t(t) = sum(sF .* g(:,2))*da / max(sum(g(:,2))*da, 1e-12);

    % Cuentas fiscales (mismos params del solver)
    tau_l = 0.15; tau_c = 0.10; phi = 0.15; Gov = 0.07; z1=0.33; z2=1.00;
    PopI = sum(g(:,1))*da; PopF = sum(g(:,2))*da;
    Ctot = sum(C{1}(:,1).*g(:,1))*da + sum(C{1}(:,2).*g(:,2))*da;
    Tl_t(t) = tau_l * z2 * PopF;
    Tc_t(t) = tau_c * Ctot;
    Tr_t(t) = phi   * z1  * PopI;
    Gx_t(t) = Gov * (PopI + PopF);
    B_t(t)  = (Tl_t(t) + Tc_t(t) - (Gx_t(t) + Tr_t(t))) / max(r_t(t),1e-9);

    % Guardar snapshots
    if ismember(t, snap_idx)
        k = find(snap_idx==t,1);
        snap(k).t = t;
        snap(k).a = a;
        snap(k).g1 = g(:,1);
        snap(k).g2 = g(:,2);
    end
end

%% --- GRÁFICOS
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% 1) Lambda ciclo y composición/ r
figure('Name','Dinamica: composicion y tasa');
subplot(2,2,1); plot(1:T, eta_t,'-o','LineWidth',1.6); grid on;
xlabel('t'); ylabel('Informalidad \eta_t'); title('Composición informal');

subplot(2,2,2); plot(1:T, r_t,'-s','LineWidth',1.6); grid on;
xlabel('t'); ylabel('r_t'); title('Tasa de interés');

subplot(2,2,3); plot(1:T, CmeanI_t,'-o', 1:T, CmeanF_t,'-s','LineWidth',1.6);
grid on; xlabel('t'); ylabel('Consumo medio'); legend('Informal','Formal','Location','best');
title('Consumo medio por tipo');

subplot(2,2,4); plot(1:T, SmeanI_t,'-o', 1:T, SmeanF_t,'-s','LineWidth',1.6);
grid on; xlabel('t'); ylabel('Ahorro medio'); legend('Informal','Formal','Location','best');
title('Ahorro medio por tipo');

set(gcf,'Color','w');

% 2) Cuentas fiscales y deuda
figure('Name','Cuentas fiscales y deuda');
subplot(2,2,1); plot(1:T,Tl_t,'-o',1:T,Tc_t,'-s','LineWidth',1.6);
grid on; xlabel('t'); ylabel('Recaudación'); legend('T_l','T_c','Location','best'); title('Recaudación');

subplot(2,2,2); plot(1:T,Tr_t,'-o',1:T,Gx_t,'-s','LineWidth',1.6);
grid on; xlabel('t'); ylabel('Gasto primario'); legend('Tr','G','Location','best'); title('Gasto primario');

subplot(2,2,3); plot(1:T,B_t,'-^','LineWidth',1.6);
grid on; xlabel('t'); ylabel('Deuda B_t'); title('Deuda pública');

subplot(2,2,4); plot(1:T,Gini_t,'-^','LineWidth',1.6);
grid on; xlabel('t'); ylabel('Gini riqueza'); title('Desigualdad (Gini)');

set(gcf,'Color','w');

% 3) Cambios en la distribución de riqueza por tipo (3 cortes)
figure('Name','Distribucion de riqueza (cortes en t)');
for k=1:numel(snap_idx)
    subplot(2, numel(snap_idx), k)
    bar(snap(k).a, snap(k).g1, 'FaceAlpha',0.6, 'EdgeColor','none'); 
    xlim([amin, 0.5]); title(sprintf('Informal (t=%d)', snap(k).t));
    xlabel('a'); ylabel('g_1(a)'); grid on;
    subplot(2, numel(snap_idx), numel(snap_idx)+k)
    bar(snap(k).a, snap(k).g2, 'FaceAlpha',0.6, 'EdgeColor','none'); 
    xlim([amin, 0.5]); title(sprintf('Formal (t=%d)', snap(k).t));
    xlabel('a'); ylabel('g_2(a)'); grid on;
end
set(gcf,'Color','w');

disp('Listo: lambda_dynamics_demo ejecutado y figuras generadas.');
