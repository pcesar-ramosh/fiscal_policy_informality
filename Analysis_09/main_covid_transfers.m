%% main_covid_transfers.m
% Compara escenario BASE vs. COVID (caída de ingresos y +transferencias)
% - Formales: z2 ↓ 15%
% - Informales: z1 ↓ 25%
% - Transferencias: phi ↑ (por defecto 0.25)
% - Bien público: se mantiene
% Graficos y CSV con estadísticas y cuentas fiscales.

clear; clc; close all;

% ----------- PARÁMETROS COMUNES ------------
RRA = 4.20;                 % RRA misma para I/F (puedes diferenciar si quieres)
eta_target = 0.64;          % Tamaño objetivo de informalidad (legacy)
theta = 0.02; rho=0.05;

tau_l = 0.15;               % impuesto laboral (solo formales)
tau_c = 0.15;               % IVA
Gov   = 0.05;               % bien público (por persona)
phi_base  = 0.10;           % transferencias base (informales)
phi_covid = 0.15;           % transferencias en COVID (ajústalo si deseas)

z1_base = 0.33; z2_base = 1.00;        % ingresos pre-COVID
z1_cvd  = z1_base*(1-0.25);            % -25% informales
z2_cvd  = z2_base*(1-0.15);            % -15% formales

% Discretización y r-bracketing
Igrid = 700; amax=5;
paramsCommon = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov, ...
    'I',Igrid,'amax',amax,'r_guess',0.03,'rmin',0.005,'rmax',0.08, ...
    'p22_bar',0.8155,'eta_target',eta_target);

% ----------- ESCENARIO BASE ---------------
parBase = paramsCommon; 
parBase.z1  = z1_base; 
parBase.z2  = z2_base; 
parBase.phi = phi_base; 
parBase.amin = -0.30*parBase.z1;

base = huggett_base_covid_function(parBase);

% ----------- ESCENARIO COVID --------------
parCVD = paramsCommon; 
parCVD.z1  = z1_cvd; 
parCVD.z2  = z2_cvd; 
parCVD.phi = phi_covid; 
parCVD.amin = -0.30*parCVD.z1;

covid = huggett_base_covid_function(parCVD);

% ----------- EXPORT CSV -------------------
T = table;
T.scenario = ["BASE"; "COVID"];
T.r        = [base.r; covid.r];
T.popI     = [base.popI; covid.popI];
T.popF     = [base.popF; covid.popF];
T.Y        = [base.Y;   covid.Y];
T.Ctot     = [base.Ctot; covid.Ctot];

% stats riqueza/consumo
wmeanI = [base.stats.wealth_mean(1); covid.stats.wealth_mean(1)];
wmeanF = [base.stats.wealth_mean(2); covid.stats.wealth_mean(2)];
wmeanT = [base.stats.wealth_mean(3); covid.stats.wealth_mean(3)];
giniI  = [base.stats.gini(1); covid.stats.gini(1)];
giniF  = [base.stats.gini(2); covid.stats.gini(2)];
giniT  = [base.stats.gini(3); covid.stats.gini(3)];
cmeanI = [base.stats.cons_mean(1); covid.stats.cons_mean(1)];
cmeanF = [base.stats.cons_mean(2); covid.stats.cons_mean(2)];
cmeanT = [base.stats.cons_mean(3); covid.stats.cons_mean(3)];

T.w_mean_I = wmeanI; T.w_mean_F = wmeanF; T.w_mean_T = wmeanT;
T.gini_I   = giniI;  T.gini_F   = giniF;  T.gini_T   = giniT;
T.c_mean_I = cmeanI; T.c_mean_F = cmeanF; T.c_mean_T = cmeanT;

% prestatarios / prestamistas
T.fracBorrow_I = [base.borrowers.fracBorrow(1); covid.borrowers.fracBorrow(1)];
T.fracBorrow_F = [base.borrowers.fracBorrow(2); covid.borrowers.fracBorrow(2)];
T.volBorrow_I  = [base.borrowers.volBorrow(1);  covid.borrowers.volBorrow(1)];
T.volBorrow_F  = [base.borrowers.volBorrow(2);  covid.borrowers.volBorrow(2)];
T.volLend_I    = [base.borrowers.volLend(1);    covid.borrowers.volLend(1)];
T.volLend_F    = [base.borrowers.volLend(2);    covid.borrowers.volLend(2)];

% fiscal
fb = base.fiscal; fc = covid.fiscal;
T.Tl   = [fb.Tl;  fc.Tl];
T.Tc   = [fb.Tc;  fc.Tc];
T.Tr   = [fb.Tr;  fc.Tr];
T.G    = [fb.G;   fc.G];
T.B    = [fb.B;   fc.B];
T.rB   = [fb.rB;  fc.rB];
T.PB   = [fb.PB;  fc.PB];
T.BB   = [fb.BB;  fc.BB];

writetable(T, 'compare_base_covid.csv');
disp('Exportado: compare_base_covid.csv');

% ----------- GRAFICOS ---------------------

% Evitar problemas de intérpretes
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

a = base.a; a2 = covid.a;  % mismas grillas por construcción

% 1) Consumo por tipo
figure('Name','Consumo por tipo (BASE vs COVID)');
subplot(1,2,1);
plot(a, base.c(:,1), 'LineWidth',1.5); hold on;
plot(a2, covid.c(:,1),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('c informal');
legend('Base','COVID','Location','best'); title('Informales'); set(gcf,'Color','w');

subplot(1,2,2);
plot(a, base.c(:,2), 'LineWidth',1.5); hold on;
plot(a2, covid.c(:,2),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('c formal');
legend('Base','COVID','Location','best'); title('Formales');

% 2) Ahorro s(a) por tipo
figure('Name','Ahorro por tipo (BASE vs COVID)');
subplot(1,2,1);
plot(a, base.s(:,1), 'LineWidth',1.5); hold on;
plot(a2, covid.s(:,1),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('s informal');
legend('Base','COVID','Location','best'); title('Informales'); set(gcf,'Color','w');

subplot(1,2,2);
plot(a, base.s(:,2), 'LineWidth',1.5); hold on;
plot(a2, covid.s(:,2),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('s formal');
legend('Base','COVID','Location','best'); title('Formales');

% 3) Distribución de riqueza g(a) por tipo (barras suaves)
figure('Name','Distribucion de riqueza (por tipo)');
subplot(2,1,1);
bar(a, base.g(:,1), 'FaceAlpha',0.5, 'EdgeColor','none'); hold on;
bar(a2, covid.g(:,1),'FaceAlpha',0.5, 'EdgeColor','none'); 
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g informal');
legend('Base','COVID','Location','northeast'); title('Informales'); set(gcf,'Color','w');

subplot(2,1,2);
bar(a, base.g(:,2), 'FaceAlpha',0.5, 'EdgeColor','none'); hold on;
bar(a2, covid.g(:,2),'FaceAlpha',0.5, 'EdgeColor','none'); 
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g formal');
legend('Base','COVID','Location','northeast'); title('Formales');

% 4) Prestatarios / prestamistas: fracciones
figure('Name','Fraccion prestatarios/prestamistas');
catX = categorical({'Informal','Formal'});
catX = reordercats(catX,{'Informal','Formal'});
subplot(1,2,1);
bar(catX, [base.borrowers.fracBorrow(:), covid.borrowers.fracBorrow(:)]);
legend('Base','COVID','Location','best'); ylabel('Frac. prestatarios'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [base.borrowers.fracLend(:), covid.borrowers.fracLend(:)]);
legend('Base','COVID','Location','best'); ylabel('Frac. prestamistas'); grid on;

% 5) Volumenes de deuda/ahorro (promedios ponderados de a)
figure('Name','Volumenes de deuda/ahorro');
subplot(1,2,1);
bar(catX, [abs(base.borrowers.volBorrow(:)), abs(covid.borrowers.volBorrow(:))]);
legend('Base','COVID','Location','best'); ylabel('|Deuda agregada|'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [base.borrowers.volLend(:), covid.borrowers.volLend(:)]);
legend('Base','COVID','Location','best'); ylabel('Ahorro agregado'); grid on;

% 6) Tasa de interés y composición
figure('Name','Tasa de interes y composicion');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.r, covid.r]); grid on; ylabel('r');
title('Tasa de interes de equilibrio'); set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.popI, covid.popI]);
grid on; ylabel('Informalidad'); title('Composicion informal');

% 7) Cuentas fiscales (niveles)
figure('Name','Cuentas fiscales: ingresos y gastos');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.fiscal.Tl+base.fiscal.Tc, covid.fiscal.Tl+covid.fiscal.Tc]);
ylabel('Ingresos'); title('Ingresos totales (Tl+Tc)'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.fiscal.G+base.fiscal.Tr+base.fiscal.rB, covid.fiscal.G+covid.fiscal.Tr+covid.fiscal.rB]);
ylabel('Gastos'); title('Gastos: G + Tr + rB'); grid on;

% 8) Saldo primario (PB) y balance total (BB)
figure('Name','Saldo primario y balance');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.fiscal.PB, covid.fiscal.PB]);
ylabel('Saldo primario'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.fiscal.BB, covid.fiscal.BB]);
ylabel('Balance total'); grid on; title('Deberia ser ~0 por construcción');

% 9) Gini (riqueza)
figure('Name','Gini de riqueza');
bar(categorical({'Base','COVID'}), [base.stats.gini(3), covid.stats.gini(3)]);
ylabel('Gini total'); grid on; set(gcf,'Color','w');
