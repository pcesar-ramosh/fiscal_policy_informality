%% main_covid_sintransf.m
% Comparación BASE vs COVID (solo caída de ingresos; SIN cambiar phi, tau_l,
% tau_c, Gov). Produce:
% - CSV: (i) cuentas fiscales, (ii) mercado de activos, (iii) estadísticas,
%        (iv) prestatarios/prestamistas
% - Gráficos: c(a), s(a), g(a), mercado de activos, cuentas fiscales,
%             r y composición, prestatarios/prestamistas, Gini/Lorenz.

clear; clc; close all;

%% ===== 1) PARÁMETROS BASE =====
RRA   = 3.40;         % RRA_I = RRA_F
rho   = 0.05;
theta = 0.02;
tau_l = 0.15;         % impuesto laboral (formales)
tau_c = 0.18;         % IVA
Gov   = 0.05;         % bien público
phi   = 0.09;         % transferencias a informales (NO cambia en COVID)

z1_base = 0.33;       % ingreso informal (baseline)
z2_base = 1.00;       % ingreso formal (baseline)

% Choque COVID: SOLO ingresos
z1_cvd  = z1_base*(1-0.25);   % -25% informales
z2_cvd  = z2_base*(1-0.15);   % -15% formales

eta_target = 0.654;           % informalidad de referencia (legacy)
p22_bar    = 0.8155;          % persistencia formal (para λ2)
Igrid = 700; amax = 5.0;

% Límite inferior ligado al ingreso informal de cada escenario
amin_base = -0.30*z1_base;
amin_cvd  = -0.30*z1_cvd;

% r: bracket y guess
r_guess = 0.03; rmin = 0.005; rmax = 0.08;

paramsCommon = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov,'phi',phi, ...
    'I',Igrid,'amax',amax,'r_guess',r_guess,'rmin',rmin,'rmax',rmax, ...
    'p22_bar',p22_bar,'eta_target',eta_target);

%% ===== 2) Resolver BASE =====
parBase = paramsCommon;
parBase.z1  = z1_base; parBase.z2  = z2_base; parBase.amin = amin_base;

base = hugguet_covid_sintrasf(parBase);

%% ===== 3) Resolver COVID (sin cambio de política) =====
parCVD = paramsCommon;
parCVD.z1  = z1_cvd;  parCVD.z2  = z2_cvd;  parCVD.amin = amin_cvd;

covid = hugguet_covid_sintrasf(parCVD);

%% ===== 4) EXPORTAR CSV COMPARATIVOS =====
if ~exist('./tables','dir'), mkdir('./tables'); end

% 4.1 Cuentas fiscales (niveles)
fb = base.fiscal; fc = covid.fiscal;
T_fiscal = table;
T_fiscal.scenario = ["BASE"; "COVID"];
T_fiscal.Tl       = [fb.Tl;  fc.Tl];
T_fiscal.Tc       = [fb.Tc;  fc.Tc];
T_fiscal.Ingresos = [fb.Tl+fb.Tc; fc.Tl+fc.Tc];
T_fiscal.Tr       = [fb.Tr;  fc.Tr];
T_fiscal.G        = [fb.G;   fc.G];
T_fiscal.rB       = [fb.rB;  fc.rB];
T_fiscal.Gastos   = [fb.Tr+fb.G+fb.rB; fc.Tr+fc.G+fc.rB];
T_fiscal.PB       = [fb.PB;  fc.PB];
T_fiscal.B        = [fb.B;   fc.B];
T_fiscal.BB       = [fb.BB;  fc.BB];
writetable(T_fiscal,'./tables/compare_fiscal_base_covid_noTr.csv');

% 4.2 Mercado de activos (demanda privada vs B pública)
da = base.a(2)-base.a(1);
A_priv_base = sum( (base.g(:,1)+base.g(:,2)).*base.a ) * da;
A_priv_cvd  = sum( (covid.g(:,1)+covid.g(:,2)).*covid.a ) * (covid.a(2)-covid.a(1));
T_assets = table(["BASE";"COVID"], [A_priv_base;A_priv_cvd], [fb.B;fc.B], ...
    'VariableNames', {'scenario','A_private','B_public'});
writetable(T_assets,'./tables/compare_asset_market_base_covid_noTr.csv');

% 4.3 Estadísticas de hogares
S = base.stats; C = covid.stats;
T_stats = table;
T_stats.scenario     = ["BASE";"COVID"];
T_stats.r            = [base.r; covid.r];
T_stats.popI         = [base.popI; covid.popI];
T_stats.popF         = [base.popF; covid.popF];
T_stats.Y            = [base.Y;   covid.Y];
T_stats.Ctot         = [base.Ctot; covid.Ctot];
T_stats.wealth_mean_I= [S.wealth_mean(1); C.wealth_mean(1)];
T_stats.wealth_mean_F= [S.wealth_mean(2); C.wealth_mean(2)];
T_stats.wealth_mean_T= [S.wealth_mean(3); C.wealth_mean(3)];
T_stats.gini_I       = [S.gini(1); C.gini(1)];
T_stats.gini_F       = [S.gini(2); C.gini(2)];
T_stats.gini_T       = [S.gini(3); C.gini(3)];
T_stats.cons_mean_I  = [S.cons_mean(1); C.cons_mean(1)];
T_stats.cons_mean_F  = [S.cons_mean(2); C.cons_mean(2)];
T_stats.cons_mean_T  = [S.cons_mean(3); C.cons_mean(3)];
T_stats.p11_rep      = [S.p11; C.p11];
writetable(T_stats,'./tables/compare_household_stats_base_covid_noTr.csv');

% 4.4 Prestatarios / prestamistas
B = base.borrowers; K = covid.borrowers;
T_borr = table;
T_borr.scenario       = ["BASE"; "COVID"];
T_borr.fracBorrow_I   = [B.fracBorrow(1); K.fracBorrow(1)];
T_borr.fracBorrow_F   = [B.fracBorrow(2); K.fracBorrow(2)];
T_borr.fracLend_I     = [B.fracLend(1);   K.fracLend(1)];
T_borr.fracLend_F     = [B.fracLend(2);   K.fracLend(2)];
T_borr.volBorrow_I    = [B.volBorrow(1);  K.volBorrow(1)];
T_borr.volBorrow_F    = [B.volBorrow(2);  K.volBorrow(2)];
T_borr.volLend_I      = [B.volLend(1);    K.volLend(1)];
T_borr.volLend_F      = [B.volLend(2);    K.volLend(2)];
writetable(T_borr,'./tables/compare_borrowers_lenders_base_covid_noTr.csv');

disp('CSV exportados (compare_*_noTr).');

%% ===== 5) GRÁFICOS =====
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

a = base.a; a2 = covid.a;

% 5.1 Consumo por tipo
figure('Name','Consumo por tipo: BASE vs COVID (sin cambios de politica)');
subplot(1,2,1);
plot(a, base.c(:,1), 'LineWidth',1.6); hold on;
plot(a2,covid.c(:,1),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('a'); ylabel('c informal');
legend('Base','COVID','Location','best'); title('Informales'); set(gcf,'Color','w');
subplot(1,2,2);
plot(a, base.c(:,2), 'LineWidth',1.6); hold on;
plot(a2,covid.c(:,2),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('a'); ylabel('c formal');
legend('Base','COVID','Location','best'); title('Formales');

% 5.2 Ahorro por tipo
figure('Name','Ahorro por tipo: BASE vs COVID');
subplot(1,2,1);
plot(a, base.s(:,1), 'LineWidth',1.6); hold on;
plot(a2,covid.s(:,1),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('a'); ylabel('s informal');
legend('Base','COVID','Location','best'); title('Informales'); set(gcf,'Color','w');
subplot(1,2,2);
plot(a, base.s(:,2), 'LineWidth',1.6); hold on;
plot(a2,covid.s(:,2),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('a'); ylabel('s formal');
legend('Base','COVID','Location','best'); title('Formales');

% 5.3 Distribución g(a)
figure('Name','Distribucion de riqueza g(a): BASE vs COVID');
subplot(2,1,1);
bar(a, base.g(:,1), 'FaceAlpha',0.55, 'EdgeColor','none'); hold on;
bar(a2,covid.g(:,1),'FaceAlpha',0.55,'EdgeColor','none');
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g_1(a)'); title('Informales'); set(gcf,'Color','w');
subplot(2,1,2);
bar(a, base.g(:,2), 'FaceAlpha',0.55, 'EdgeColor','none'); hold on;
bar(a2,covid.g(:,2),'FaceAlpha',0.55,'EdgeColor','none');
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g_2(a)'); title('Formales');

% 5.4 Prestatarios / prestamistas (fracciones y volúmenes)
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});

figure('Name','Fraccion prestatarios/prestamistas');
subplot(1,2,1);
bar(catX, [base.borrowers.fracBorrow(:), covid.borrowers.fracBorrow(:)]);
legend('Base','COVID','Location','best'); ylabel('Frac. prestatarios'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [base.borrowers.fracLend(:), covid.borrowers.fracLend(:)]);
legend('Base','COVID','Location','best'); ylabel('Frac. prestamistas'); grid on;

figure('Name','Volumenes de deuda/ahorro');
subplot(1,2,1);
bar(catX, [abs(base.borrowers.volBorrow(:)), abs(covid.borrowers.volBorrow(:))]);
legend('Base','COVID','Location','best'); ylabel('|Deuda agregada|'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [base.borrowers.volLend(:), covid.borrowers.volLend(:)]);
legend('Base','COVID','Location','best'); ylabel('Ahorro agregado'); grid on;

% 5.5 Tasa de interés y composición
figure('Name','Tasa de interes y composicion');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [base.r, covid.r]); grid on; ylabel('r');
title('Tasa de interes de equilibrio'); set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [base.popI, covid.popI]);
grid on; ylabel('Informalidad'); title('Composicion informal');

% 5.6 Cuentas fiscales (niveles)
figure('Name','Cuentas fiscales: ingresos y gastos');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [fb.Tl+fb.Tc, fc.Tl+fc.Tc]);
ylabel('Ingresos'); title('Ingresos totales (Tl+Tc)'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [fb.Tr+fb.G+fb.rB, fc.Tr+fc.G+fc.rB]);
ylabel('Gastos'); title('Gastos: Tr + G + rB'); grid on;

% 5.7 Saldo primario y balance
figure('Name','Saldo primario y balance');
subplot(1,2,1);
bar(categorical({'Base','COVID'}), [fb.PB, fc.PB]);
ylabel('Saldo primario'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Base','COVID'}), [fb.BB, fc.BB]);
ylabel('Balance total'); grid on; title('≈0 por construcción');

% 5.8 Curva de Lorenz (riqueza total)
daB  = a(2)-a(1); gTB = base.g(:,1)+base.g(:,2); WB=sum(gTB)*daB;
[asB,ixB] = sort(a); gsB = gTB(ixB)*daB; cumPopB = cumsum(gsB)/WB;
wealth_posB = asB - min(0,min(asB)) + 1e-12; cumWealthB = cumsum(wealth_posB.*gsB);

daC  = a2(2)-a2(1); gTC = covid.g(:,1)+covid.g(:,2); WC=sum(gTC)*daC;
[asC,ixC] = sort(a2); gsC = gTC(ixC)*daC; cumPopC = cumsum(gsC)/WC;
wealth_posC = asC - min(0,min(asC)) + 1e-12; cumWealthC = cumsum(wealth_posC.*gsC);

if cumWealthB(end)>0 && cumWealthC(end)>0
    LB = cumWealthB/cumWealthB(end);
    LC = cumWealthC/cumWealthC(end);
    figure('Name','Curva de Lorenz (riqueza total): Base vs COVID');
    plot(cumPopB, LB, 'LineWidth',1.8); hold on;
    plot(cumPopC, LC, '--','LineWidth',1.8);
    plot([0,1],[0,1],'k:'); grid on; axis square;
    xlabel('Fracción población'); ylabel('Fracción riqueza acumulada');
    title(sprintf('Lorenz: Gini Base=%.3f, COVID=%.3f', base.stats.gini(3), covid.stats.gini(3)));
    legend('Base','COVID','Location','southeast'); set(gcf,'Color','w');
end

disp('Listo. Revisa los CSVs y las figuras (COVID sin cambios de política).');
