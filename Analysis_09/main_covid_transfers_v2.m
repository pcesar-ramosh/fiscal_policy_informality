%% main_covid_transfers_v2.m
% Comparación BASE vs COVID con B exógena (no se fuerza Ingresos = Gastos)

clear; clc; close all;

% Parámetros comunes
RRA = 4.2;  eta_target = 0.64;  theta=0.02;  rho=0.05;
tau_l=0.15; tau_c=0.10; Gov=0.07;
phi_base=0.15;  phi_covid=0.25;
z1_base=0.33;   z2_base=1.00;
z1_cvd = z1_base*(1-0.25);
z2_cvd = z2_base*(1-0.15);

Igrid=700; amax=5;

paramsCommon = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov, ...
    'I',Igrid,'amax',amax,'r_guess',0.03,'rmin',0.005,'rmax',0.08, ...
    'p22_bar',0.8155,'eta_target',eta_target);

% ---- Deuda exógena: porcentaje del PIB (misma regla en ambos) ----
debtRule.mode = 'exog_ratio';
debtRule.bY_ratio = 0.40;   % 40% del PIB (ajústalo)
% (alternativa: debtRule.mode='exog_level'; debtRule.B_level=0.40; )

% BASE
parBase = paramsCommon; parBase.z1=z1_base; parBase.z2=z2_base; parBase.phi=phi_base;
parBase.amin = -0.30*parBase.z1; parBase.debt = debtRule;
base = huggett_base_function_v2(parBase);

% COVID
parCVD = paramsCommon; parCVD.z1=z1_cvd; parCVD.z2=z2_cvd; parCVD.phi=phi_covid;
parCVD.amin = -0.30*parCVD.z1; parCVD.debt = debtRule;
covid = huggett_base_function_v2(parCVD);

% ---- Export CSV "de verdad" (sin forzar cierre) ----
T = table;
T.scenario = ["BASE"; "COVID"];
T.r      = [base.r; covid.r];
T.popI   = [base.popI; covid.popI];
T.Y      = [base.Y;   covid.Y];
T.Ctot   = [base.Ctot; covid.Ctot];
T.Tl     = [base.fiscal.Tl;  covid.fiscal.Tl];
T.Tc     = [base.fiscal.Tc;  covid.fiscal.Tc];
T.Tr     = [base.fiscal.Tr;  covid.fiscal.Tr];
T.G      = [base.fiscal.G;   covid.fiscal.G];
T.B      = [base.fiscal.B;   covid.fiscal.B];
T.rB     = [base.fiscal.rB;  covid.fiscal.rB];
T.PB     = [base.fiscal.PB;  covid.fiscal.PB];
T.BB     = [base.fiscal.BB;  covid.fiscal.BB];
T.gapBudget = [base.fiscal.gapBudget; covid.fiscal.gapBudget];
T.giniTot   = [base.stats.gini(3); covid.stats.gini(3)];
writetable(T,'compare_base_covid_v2.csv');
disp('Exportado: compare_base_covid_v2.csv');

% ---- Gráficos (sin LaTeX) ----
a = base.a; a2 = covid.a;

% 1) Ingresos vs Gastos (no forzados)
figure('Name','Ingresos y Gastos (no forzados)');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.fiscal.Tl+base.fiscal.Tc, covid.fiscal.Tl+covid.fiscal.Tc]);
ylabel('Ingresos'); title('Ingresos: Tl+Tc'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.fiscal.G+base.fiscal.Tr+base.fiscal.rB, covid.fiscal.G+covid.fiscal.Tr+covid.fiscal.rB]);
ylabel('Gastos'); title('Gastos: G + Tr + rB'); grid on;

% 2) Saldo primario y balance total
figure('Name','Saldo primario y balance total');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.fiscal.PB, covid.fiscal.PB]);
ylabel('Saldo primario'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.fiscal.BB, covid.fiscal.BB]);
ylabel('Balance total'); title('BB = PB - rB'); grid on;

% 3) Tasa r y composición informal
figure('Name','Tasa r y composicion informal');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.r, covid.r]); ylabel('r'); title('Tasa de interes'); grid on;
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.popI, covid.popI]); ylabel('Informalidad'); title('Composicion informal'); grid on;

% 4) Consumo por tipo
figure('Name','Consumo por tipo');
subplot(1,2,1);
plot(a, base.c(:,1),'LineWidth',1.5); hold on; plot(a2,covid.c(:,1),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('Consumo informal'); legend('Base','COVID','Location','best');
subplot(1,2,2);
plot(a, base.c(:,2),'LineWidth',1.5); hold on; plot(a2,covid.c(:,2),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('Consumo formal');  legend('Base','COVID','Location','best');

% 5) Ahorro por tipo
figure('Name','Ahorro por tipo');
subplot(1,2,1);
plot(a, base.s(:,1),'LineWidth',1.5); hold on; plot(a2,covid.s(:,1),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('s informal'); legend('Base','COVID','Location','best');
subplot(1,2,2);
plot(a, base.s(:,2),'LineWidth',1.5); hold on; plot(a2,covid.s(:,2),'--','LineWidth',1.5);
yline(0,'k:'); grid on; xlabel('a'); ylabel('s formal');   legend('Base','COVID','Location','best');

% 6) Distribucion de riqueza (barras)
figure('Name','Distribucion de riqueza');
subplot(2,1,1);
bar(a, base.g(:,1), 'FaceAlpha',0.5, 'EdgeColor','none'); hold on;
bar(a2, covid.g(:,1),'FaceAlpha',0.5,'EdgeColor','none'); xlim([min(a) 1.0]);
xlabel('a'); ylabel('g informal'); legend('Base','COVID','Location','northeast'); grid on;
subplot(2,1,2);
bar(a, base.g(:,2), 'FaceAlpha',0.5, 'EdgeColor','none'); hold on;
bar(a2, covid.g(:,2),'FaceAlpha',0.5,'EdgeColor','none'); xlim([min(a) 1.0]);
xlabel('a'); ylabel('g formal');   legend('Base','COVID','Location','northeast'); grid on;

% 7) Prestatarios/prestamistas (fracciones)
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});
figure('Name','Prestatarios y prestamistas');
subplot(1,2,1);
bar(catX, [base.borrowers.fracBorrow(:), covid.borrowers.fracBorrow(:)]);
ylabel('Frac. prestatarios'); legend('Base','COVID','Location','best'); grid on;
subplot(1,2,2);
bar(catX, [base.borrowers.fracLend(:),   covid.borrowers.fracLend(:)]);
ylabel('Frac. prestamistas'); legend('Base','COVID','Location','best'); grid on;
