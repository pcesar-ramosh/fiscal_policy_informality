%% =============================================================
%  main_base_model.m  — GE + simulación fiscal con transiciones
%  - r*(sI)
%  - Fiscal vs sI y vs eta
%  - Simulación temporal: Tl_t, Tc_t, Tr_t, G_t, rB_t, Gap_t, B_t
%  - Barras: distribución de riqueza por agente (informal/formal)
% =============================================================
clear; clc; close all; rng(123);

% --- usar TEX por defecto (para tildes/ñ); fórmulas se ponen con 'Interpreter','latex'
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

%% ----- Barrido en sI (RRA informal)
n_agents  = 19;
s_min     = 3.15; 
s_max     = 5.30;

% CLAVE: dejar eta libre (NaN) para que dependa de lambdas(sI)
eta_vector = NaN(1, n_agents);                 % <---- composición endógena
sI_vector1 = linspace(s_min, s_max, n_agents);  % barrido sI
sF_vector1 = 5.30 * ones(1, n_agents);          % RRA formal fija

%% ----- Resolver GE en el barrido sI
[r_vec, ~, pop1_vector, statsMatrix, statsCMatrix, GDistribution, a, ...
 Distribution, Fiscal, Cpolicies, Spolicies, Params] = ...
    huggett_Equi_RRA_function_transfer(eta_vector, sI_vector1, sF_vector1);

% Baseline (último jj)
jj0 = n_agents;
g   = Distribution{jj0};
aa  = a(:);
da  = (aa(end)-aa(1))/(numel(aa)-1);

theta = Params.theta; 
tau_c = Params.tau_c; 
tau_l = Params.tau_l; 
phi   = Params.phi; 
Gov   = Params.Gov;   
p11   = Params.p11_vec(jj0); 
p22   = Params.p22_vec(jj0);
r_bas = r_vec(jj0);

%% ===================== DISTRIBUCIÓN DE RIQUEZA (barras) =====================
% bins y conteos ponderados por la densidad * da
edges = linspace(min(aa), max(aa), 41);      % 40 bins
centers = 0.5*(edges(1:end-1)+edges(2:end));
idx = discretize(aa, edges);

mass_inf = g(:,1)*da; 
mass_for = g(:,2)*da;
valid = ~isnan(idx);
counts_inf = accumarray(idx(valid), mass_inf(valid), [numel(edges)-1,1], @sum, 0);
counts_for = accumarray(idx(valid), mass_for(valid), [numel(edges)-1,1], @sum, 0);
counts_tot = counts_inf + counts_for;

figure;
bar(centers, [counts_inf counts_for], 'grouped');
xlabel('Activos a'); ylabel('Masa');
legend('Informal','Formal','Location','best');
title('Distribución de riqueza por tipo (barras)'); grid on; set(gcf,'Color','w');

% (opcional) total en barras apiladas
figure;
bar(centers, [counts_inf counts_for], 'stacked');
xlabel('Activos a'); ylabel('Masa'); grid on;
legend('Informal','Formal','Location','best');
title('Distribución de riqueza total (apilada)'); set(gcf,'Color','w');

%% ===================== POLÍTICAS =====================
c_baseline   = Cpolicies{jj0};
sflow_baseln = Spolicies{jj0};

figure; 
plot(aa,c_baseline(:,1),'-','LineWidth',1.6); hold on;
plot(aa,c_baseline(:,2),'--','LineWidth',1.6);
xlabel('Activos a'); ylabel('Consumo c(a)'); grid on;
legend('Informal','Formal','Location','best'); 
title('Política de consumo por tipo'); set(gcf,'Color','w');

figure; 
plot(aa,sflow_baseln(:,1),'-','LineWidth',1.6); hold on;
plot(aa,sflow_baseln(:,2),'--','LineWidth',1.6); yline(0,'k:');
xlabel('Activos a'); 
ylabel('$\dot{a}(a)$','Interpreter','latex');   % <<-- solo aquí usamos LaTeX
grid on;
legend('Informal','Formal','Location','best');
title('Política de ahorro por tipo'); set(gcf,'Color','w');

%% ===================== r*(sI) =====================
figure; 
plot(sI_vector1, r_vec,'-o','LineWidth',1.6);
xlabel('$s_I$','Interpreter','latex'); 
ylabel('$r^\*(s_I)$','Interpreter','latex');
grid on;
title('$r^\*$ vs. $s_I$','Interpreter','latex'); 
set(gcf,'Color','w');

%% ===================== FISCAL vs sI =====================
% Fiscal{jj} = [Tc Tl Tr G rB Gap popI popF Y Btarget]
F = cell2mat(Fiscal(:));
Tc = F(:,1); Tl = F(:,2); Tr = F(:,3); Gg = F(:,4); rB = F(:,5); Gap = F(:,6);

figure; 
plot(sI_vector1,Tl,'-o', sI_vector1,Tc,'-s','LineWidth',1.6);
legend('$T_{\ell}$ (laboral)','$T_c$ (IVA)','Location','best','Interpreter','latex');
xlabel('$s_I$','Interpreter','latex'); ylabel('Ingresos'); grid on;
title('Ingresos fiscales vs. $s_I$','Interpreter','latex'); set(gcf,'Color','w');

figure; 
plot(sI_vector1,Tr,'-o', sI_vector1,Gg,'-s', sI_vector1,rB,'-^','LineWidth',1.6);
legend('Transferencias','$G$ (proxy)','$rB$','Location','best','Interpreter','latex');
xlabel('$s_I$','Interpreter','latex'); ylabel('Gastos'); grid on;
title('Gastos fiscales vs. $s_I$','Interpreter','latex'); set(gcf,'Color','w');

figure; 
plot(sI_vector1,Gap,'-o','LineWidth',1.6); yline(0,'k:');
xlabel('$s_I$','Interpreter','latex'); ylabel('Gap'); grid on;
title('Balance fiscal (Gap) vs. $s_I$','Interpreter','latex'); set(gcf,'Color','w');

%% ===================== Barrido en \eta (para comparar) =====================
eta_axis = linspace(0.40,0.80,13);
sI_mid   = median(sI_vector1); 
sF_mid   = median(sF_vector1);

[r_eta, ~, eta_out, ~, ~, ~, ~, ~, Fiscal_eta, ~, ~, Params2] = ...
    huggett_Equi_RRA_function_transfer(eta_axis, sI_mid*ones(size(eta_axis)), sF_mid*ones(size(eta_axis)));

figure; 
plot(eta_out, r_eta, '-o','LineWidth',1.6);
xlabel('Informalidad \eta'); 
ylabel('$r$','Interpreter','latex'); grid on;
title('$r(\eta)$ con $B(r)$ endógeno','Interpreter','latex'); 
set(gcf,'Color','w');

%% ===================== SIMULACIÓN DINÁMICA (transiciones por periodo) ===
% Consumo medio por tipo en el equilibrio jj0
popI_ss = sum(g(:,1))*da; popF_ss = sum(g(:,2))*da;
CbarI   = (sum(c_baseline(:,1).*g(:,1))*da) / max(popI_ss,eps);
CbarF   = (sum(c_baseline(:,2).*g(:,2))*da) / max(popF_ss,eps);

% Probabilidades un-periodo
p_stay_I = p11; p_stay_F = p22; 
p_ItoF   = 1 - p_stay_I; 
p_FtoI   = 1 - p_stay_F;

% Tamaño y horizonte
N = 5000; T = 200;

% Estados iniciales ≈ estacionario
S = rand(N,1) > popF_ss; % 1=Informal, 0=Formal

shareI = nan(T,1); Tl_t=shareI; Tc_t=shareI; Tr_t=shareI; G_t=shareI; rB_t=shareI; Gap_t=shareI; B_t=shareI;

% punto de partida para deuda (usa Btarget del barrido en eta si existe)
if exist('Fiscal_eta','var') && ~isempty(Fiscal_eta)
    B_t(1) = Fiscal_eta{end}(10);
else
    B_t(1) = F(end,10); % fallback: último Btarget del barrido sI
end

for t=1:T
    shareI(t) = mean(S==1);
    shareF_t  = 1 - shareI(t);

    % Impuestos/transferencias periodo t
    Tl_t(t) = tau_l * 1.0 * shareF_t;                        % z2=1
    Tc_t(t) = tau_c * ( CbarI*shareI(t) + CbarF*shareF_t );  % IVA
    Tr_t(t) = phi   * 0.33 * shareI(t);                      % z1=0.33
    G_t(t)  = Gov * (shareI(t) + shareF_t);

    % Intereses y deuda (manejo t=1 sin t-1)
    if t==1
        rB_t(t) = r_bas * B_t(1);
    else
        rB_t(t) = r_bas * B_t(t-1);
    end
    Gap_t(t) = G_t(t) + Tr_t(t) + rB_t(t) - (Tl_t(t) + Tc_t(t));
    if t < T
        B_t(t+1) = B_t(t) + Gap_t(t);
    end

    % Transiciones al periodo t+1
    U = rand(N,1);
    I_now = (S==1); F_now = ~I_now;
    S(I_now & (U<=p_ItoF)) = 0;   % I->F
    S(F_now & (U<=p_FtoI)) = 1;   % F->I
end

% Gráficos dinámicos
figure;
subplot(3,1,1); plot(shareI,'LineWidth',1.2); ylim([0 1]); grid on;
ylabel('share informal'); title('Transiciones: composición por periodo'); set(gcf,'Color','w');
subplot(3,1,2); plot(Tl_t,'-','LineWidth',1.1); hold on; plot(Tc_t,'-','LineWidth',1.1);
legend('$T_{\ell}$','$T_c$','Interpreter','latex'); grid on; ylabel('Ingresos');
subplot(3,1,3); plot(Tr_t,'-','LineWidth',1.1); hold on; plot(G_t,'-','LineWidth',1.1); plot(rB_t,'-','LineWidth',1.1);
legend('Tr','$G$','$rB$','Interpreter','latex'); grid on; xlabel('t'); ylabel('Gastos');

figure;
subplot(2,1,1); plot(Gap_t,'-','LineWidth',1.2); yline(0,'k:'); grid on;
title('Gap_t por periodo'); ylabel('Gap');
subplot(2,1,2); plot(B_t,'-','LineWidth',1.2); grid on; xlabel('t'); ylabel('Deuda B_t');

%% ===================== EXPORTES CSV =====================
% Wealth stats (incluye Gini por tipo y total)
T_wealth = table(pop1_vector(:), sI_vector1(:), sF_vector1(:), r_vec(:), ...
    statsMatrix(:,1),statsMatrix(:,2),statsMatrix(:,3),statsMatrix(:,4), ...
    statsMatrix(:,5),statsMatrix(:,6),statsMatrix(:,7),statsMatrix(:,8), ...
    statsMatrix(:,9),statsMatrix(:,10), ...
    'VariableNames', {'eta_out','sI','sF','r','gmean_inf','gmean_for','gmed_inf','gmed_for','gmean_tot','gmed_tot','gini_inf','gini_for','gini_tot','p11'});
writetable(T_wealth, 'stats_wealth.csv');

% Consumo
T_cons = table(pop1_vector(:), sI_vector1(:), sF_vector1(:), r_vec(:), ...
    statsCMatrix(:,1),statsCMatrix(:,2),statsCMatrix(:,3), ...
    statsCMatrix(:,4),statsCMatrix(:,5),statsCMatrix(:,6), ...
    'VariableNames', {'eta_out','sI','sF','r','Cmean_inf','Cmean_for','Cmed_inf','Cmed_for','Cmean_tot','Cmed_tot'});
writetable(T_cons, 'stats_consumption.csv');

% Fiscal vs sI
T_fiscal = table(pop1_vector(:), sI_vector1(:), sF_vector1(:), r_vec(:), ...
    F(:,1),F(:,2),F(:,3),F(:,4),F(:,5),F(:,6),F(:,7),F(:,8),F(:,9),F(:,10), ...
    'VariableNames', {'eta_out','sI','sF','r','Tc','Tl','Tr','G','rB','Gap','popI','popF','Y','Btarget'});
writetable(T_fiscal, 'fiscal_summary_vs_sI.csv');

% r vs sI y vs eta
writetable(table(sI_vector1(:), r_vec(:), 'VariableNames', {'sI','r'}), 'r_vs_sI.csv');
writetable(table(eta_out(:), r_eta(:), 'VariableNames', {'eta','r'}), 'r_vs_eta.csv');

% Fiscal vs eta
Fe = cell2mat(Fiscal_eta(:));
T_fiscal_eta = table(eta_out(:), Fe(:,1), Fe(:,2), Fe(:,3), Fe(:,4), Fe(:,5), Fe(:,6), Fe(:,7), Fe(:,8), Fe(:,9), Fe(:,10), ...
    'VariableNames', {'eta','Tc','Tl','Tr','G','rB','Gap','popI','popF','Y','Btarget'});
writetable(T_fiscal_eta, 'fiscal_summary_vs_eta.csv');

% Dinámica simulada
T_dyn = table((1:T)', shareI, Tl_t, Tc_t, Tr_t, G_t, rB_t, Gap_t, B_t, ...
    'VariableNames', {'t','shareI','Tl','Tc','Tr','G','rB','Gap','B'});
writetable(T_dyn, 'fiscal_dynamics_simulation.csv');

disp('CSVs exportados (carpeta actual).');
