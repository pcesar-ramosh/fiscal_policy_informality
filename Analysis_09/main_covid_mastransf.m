%% main_covid_mastransf.m
% Escenario COVID: caída de ingresos (z1 ↓25%, z2 ↓15%) y
% comparación de transferencias: phi = 0.09 (base COVID) vs 0.14 (más transferencias).
% Produce:
% - CSV: comparación de (i) cuentas fiscales, (ii) mercado de activos,
%        (iii) estadísticas, (iv) prestatarios/prestamistas, (v) puntos de Lorenz
% - Gráficos: c(a), s(a), g(a), mercado de activos, cuentas fiscales,
%             prestatarios/prestamistas, tasa de interés, gini y Curva de Lorenz.

clear; clc; close all;

addpath(pwd);  % aseguramos path local

%% ===== 1) PARÁMETROS COMUNES COVID =====
RRA   = 3.40;          % RRA_I = RRA_F
rho   = 0.05;
theta = 0.02;
tau_l = 0.15;          % impuesto laboral (solo formales)
tau_c = 0.18;          % IVA
Gov   = 0.05;          % bien público (flujo aditivo)
phi_covid_base = 0.09; % transferencias base (informales) EN COVID
phi_covid_high = 0.14; % transferencias altas EN COVID

% Ingresos pre-COVID
z1_pre = 0.33;
z2_pre = 1.00;

% Caída por COVID
z1_cvd = z1_pre*(1-0.25);   % informales -25%
z2_cvd = z2_pre*(1-0.15);   % formales  -15%

% Informalidad objetivo (legacy) y persistencia
eta_target = 0.654;
p22_bar    = 0.8155;

% Discretización y r-bracketing
Igrid = 700; amax=5.0;
amin_base = -0.30*z1_cvd;   % límite de crédito atado al ingreso informal COVID

r_guess = 0.03; rmin = 0.005; rmax = 0.08;

paramsCommon = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov, ...
    'I',Igrid,'amax',amax,'amin',amin_base, ...
    'r_guess',r_guess,'rmin',rmin,'rmax',rmax, ...
    'p22_bar',p22_bar,'eta_target',eta_target);

%% ===== 2) Resolver: COVID con phi=0.09 vs phi=0.14 =====
parCVD_low = paramsCommon;  parCVD_low.z1 = z1_cvd; parCVD_low.z2 = z2_cvd; parCVD_low.phi = phi_covid_base;
parCVD_hi  = paramsCommon;  parCVD_hi.z1  = z1_cvd; parCVD_hi.z2  = z2_cvd; parCVD_hi.phi  = phi_covid_high;

covid_low = hugguet_covid_mastransf(parCVD_low);
covid_hi  = hugguet_covid_mastransf(parCVD_hi);

% Atajos
a  = covid_low.a;   a2 = covid_hi.a;   % misma grilla por construcción
g1 = covid_low.g;   g2 = covid_hi.g;
c1 = covid_low.c;   c2 = covid_hi.c;
s1 = covid_low.s;   s2 = covid_hi.s;

fb1 = covid_low.fiscal; fb2 = covid_hi.fiscal;

%% ===== 3) EXPORTAR CSVs =====
if ~exist('./tables','dir'); mkdir('./tables'); end

% 3.1 Cuentas fiscales
T_fiscal = table;
T_fiscal.scenario  = ["COVID_phi=0.09"; "COVID_phi=0.14"];
T_fiscal.Tl        = [fb1.Tl;  fb2.Tl];
T_fiscal.Tc        = [fb1.Tc;  fb2.Tc];
T_fiscal.Ingresos  = [fb1.Tl+fb1.Tc; fb2.Tl+fb2.Tc];
T_fiscal.Tr        = [fb1.Tr;  fb2.Tr];
T_fiscal.G         = [fb1.G;   fb2.G];
T_fiscal.rB        = [fb1.rB;  fb2.rB];
T_fiscal.Gastos    = [fb1.Tr+fb1.G+fb1.rB; fb2.Tr+fb2.G+fb2.rB];
T_fiscal.PB        = [fb1.PB;  fb2.PB];
T_fiscal.B         = [fb1.B;   fb2.B];
T_fiscal.BB        = [fb1.BB;  fb2.BB];
writetable(T_fiscal, './tables/covid_mastransf_fiscal.csv');

% 3.2 Mercado de activos
da = a(2)-a(1);
Apriv1 = sum((g1(:,1)+g1(:,2)).*a)*da;
Apriv2 = sum((g2(:,1)+g2(:,2)).*a)*da;
T_assets = table(["COVID_phi=0.09"; "COVID_phi=0.14"], [Apriv1; Apriv2], [fb1.B; fb2.B], ...
    'VariableNames', {'scenario','A_private','B_public'});
writetable(T_assets, './tables/covid_mastransf_assets.csv');

% 3.3 Estadísticas de hogar
S1 = covid_low.stats; S2 = covid_hi.stats;
T_stats = table;
T_stats.scenario = ["COVID_phi=0.09"; "COVID_phi=0.14"];
T_stats.r        = [covid_low.r; covid_hi.r];
T_stats.popI     = [covid_low.popI; covid_hi.popI];
T_stats.popF     = [covid_low.popF; covid_hi.popF];
T_stats.Y        = [covid_low.Y;   covid_hi.Y];
T_stats.Ctot     = [covid_low.Ctot; covid_hi.Ctot];

T_stats.w_mean_I = [S1.wealth_mean(1); S2.wealth_mean(1)];
T_stats.w_mean_F = [S1.wealth_mean(2); S2.wealth_mean(2)];
T_stats.w_mean_T = [S1.wealth_mean(3); S2.wealth_mean(3)];
T_stats.w_med_I  = [S1.wealth_median(1); S2.wealth_median(1)];
T_stats.w_med_F  = [S1.wealth_median(2); S2.wealth_median(2)];
T_stats.w_med_T  = [S1.wealth_median(3); S2.wealth_median(3)];
T_stats.gini_I   = [S1.gini(1); S2.gini(1)];
T_stats.gini_F   = [S1.gini(2); S2.gini(2)];
T_stats.gini_T   = [S1.gini(3); S2.gini(3)];
T_stats.c_mean_I = [S1.cons_mean(1); S2.cons_mean(1)];
T_stats.c_mean_F = [S1.cons_mean(2); S2.cons_mean(2)];
T_stats.c_mean_T = [S1.cons_mean(3); S2.cons_mean(3)];
T_stats.c_med_I  = [S1.cons_median(1); S2.cons_median(1)];
T_stats.c_med_F  = [S1.cons_median(2); S2.cons_median(2)];
T_stats.c_med_T  = [S1.cons_median(3); S2.cons_median(3)];
writetable(T_stats, './tables/covid_mastransf_stats.csv');

% 3.4 Prestatarios / prestamistas
B1 = covid_low.borrowers; B2 = covid_hi.borrowers;
T_borr = table;
T_borr.scenario     = ["COVID_phi=0.09"; "COVID_phi=0.14"];
T_borr.fracBorrow_I = [B1.fracBorrow(1); B2.fracBorrow(1)];
T_borr.fracBorrow_F = [B1.fracBorrow(2); B2.fracBorrow(2)];
T_borr.fracLend_I   = [B1.fracLend(1);   B2.fracLend(1)];
T_borr.fracLend_F   = [B1.fracLend(2);   B2.fracLend(2)];
T_borr.volBorrow_I  = [B1.volBorrow(1);  B2.volBorrow(1)];
T_borr.volBorrow_F  = [B1.volBorrow(2);  B2.volBorrow(2)];
T_borr.volLend_I    = [B1.volLend(1);    B2.volLend(1)];
T_borr.volLend_F    = [B1.volLend(2);    B2.volLend(2)];
writetable(T_borr, './tables/covid_mastransf_borrowers.csv');

%% ===== 4) GRÁFICOS =====
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% 4.1 Consumo c(a) (por tipo)
figure('Name','Consumo por tipo: COVID (phi=0.09 vs 0.14)');
subplot(1,2,1);
plot(a,  c1(:,1),'LineWidth',1.6); hold on;
plot(a2, c2(:,1),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('c_I(a)');
legend('\phi=0.09','\phi=0.14','Location','best'); title('Informales'); set(gcf,'Color','w');

subplot(1,2,2);
plot(a,  c1(:,2),'LineWidth',1.6); hold on;
plot(a2, c2(:,2),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('c_F(a)');
legend('\phi=0.09','\phi=0.14','Location','best'); title('Formales');

% 4.2 Ahorro s(a)
figure('Name','Ahorro por tipo: COVID (phi=0.09 vs 0.14)');
subplot(1,2,1);
plot(a,  s1(:,1),'LineWidth',1.6); hold on;
plot(a2, s2(:,1),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('s_I(a)');
legend('\phi=0.09','\phi=0.14','Location','best'); title('Informales'); set(gcf,'Color','w');

subplot(1,2,2);
plot(a,  s1(:,2),'LineWidth',1.6); hold on;
plot(a2, s2(:,2),'--','LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('s_F(a)');
legend('\phi=0.09','\phi=0.14','Location','best'); title('Formales');

% 4.3 Distribución de riqueza g(a)
figure('Name','Distribución de riqueza: COVID (phi=0.09 vs 0.14)');
subplot(2,1,1);
bar(a, g1(:,1),'FaceAlpha',0.55,'EdgeColor','none'); hold on;
bar(a2,g2(:,1),'FaceAlpha',0.40,'EdgeColor','none'); 
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g_1(a)');
legend('\phi=0.09','\phi=0.14','Location','northeast'); title('Informales'); set(gcf,'Color','w');

subplot(2,1,2);
bar(a, g1(:,2),'FaceAlpha',0.55,'EdgeColor','none'); hold on;
bar(a2,g2(:,2),'FaceAlpha',0.40,'EdgeColor','none'); 
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g_2(a)');
legend('\phi=0.09','\phi=0.14','Location','northeast'); title('Formales');

% 4.4 Mercado de activos: Demanda privada vs Oferta pública
figure('Name','Mercado de activos: COVID (phi=0.09 vs 0.14)');
bar(categorical({'Apriv(\phi=0.09)','Apriv(\phi=0.14)','B(\phi=0.09)','B(\phi=0.14)'}), ...
    [Apriv1, Apriv2, fb1.B, fb2.B]);
ylabel('Nivel'); grid on; title(sprintf('Cierre a r_{0.09}=%.4f, r_{0.14}=%.4f', covid_low.r, covid_hi.r));
set(gcf,'Color','w');

% 4.5 Cuentas fiscales: DESAGREGADO (ingresos y gastos) + deltas
figure('Name','Cuentas fiscales (COVID) - desagregado');

% --- Ingresos desagregados: Tl y Tc ---
subplot(1,3,1);
Xing = categorical({'Tl','Tc'}); Xing = reordercats(Xing, {'Tl','Tc'});
Ying = [fb1.Tl fb2.Tl; fb1.Tc fb2.Tc];
bar(Xing, Ying, 'grouped');
ylabel('Ingresos'); title('Ingresos desagregados');
legend('\phi=0.09','\phi=0.14','Location','best'); grid on; set(gcf,'Color','w');

% --- Gastos desagregados: Tr, G, rB ---
subplot(1,3,2);
Xgas = categorical({'Tr','G','rB'}); Xgas = reordercats(Xgas, {'Tr','G','rB'});
Ygas = [fb1.Tr fb2.Tr; fb1.G fb2.G; fb1.rB fb2.rB];
bar(Xgas, Ygas, 'grouped');
ylabel('Gastos'); title('Gastos desagregados');
legend('\phi=0.09','\phi=0.14','Location','best'); grid on;

% --- Cambios absolutos (phi 0.14 − phi 0.09) por componente ---
subplot(1,3,3);
Xall = categorical({'Tl','Tc','Tr','G','rB','PB','B'});
Xall = reordercats(Xall, {'Tl','Tc','Tr','G','rB','PB','B'});
Delta = [fb2.Tl - fb1.Tl, ...
         fb2.Tc - fb1.Tc, ...
         fb2.Tr - fb1.Tr, ...
         fb2.G  - fb1.G,  ...
         fb2.rB - fb1.rB, ...
         fb2.PB - fb1.PB, ...
         fb2.B  - fb1.B];
bar(Xall, Delta);
yline(0,'k:'); grid on;
ylabel('\Delta nivel (\phi=0.14 - \phi=0.09)');
title('Cambios por componente');

% (Opcional) Etiquetas de valor encima de cada barra de Delta
yl = ylim;
for i = 1:numel(Delta)
    text(Xall(i), Delta(i) + 0.02*(yl(2)-yl(1)), sprintf('%.4f', Delta(i)), ...
        'HorizontalAlignment','center', 'Rotation', 90, 'FontSize', 8);
end


% 4.6 Prestatarios / prestamistas (fracciones y volúmenes)
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});

figure('Name','Fracciones prestatarios/prestamistas (COVID)');
subplot(1,2,1);
bar(catX, [covid_low.borrowers.fracBorrow(:), covid_hi.borrowers.fracBorrow(:)]);
ylabel('Frac. prestatarios (a<0)'); legend('\phi=0.09','\phi=0.14','Location','best'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [covid_low.borrowers.fracLend(:), covid_hi.borrowers.fracLend(:)]);
ylabel('Frac. prestamistas (a>0)'); legend('\phi=0.09','\phi=0.14','Location','best'); grid on;

figure('Name','Volúmenes deuda/ahorro (COVID)');
subplot(1,2,1);
bar(catX, [abs(covid_low.borrowers.volBorrow(:)), abs(covid_hi.borrowers.volBorrow(:))]);
ylabel('|Deuda agregada|'); legend('\phi=0.09','\phi=0.14','Location','best'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [covid_low.borrowers.volLend(:), covid_hi.borrowers.volLend(:)]);
ylabel('Ahorro agregado'); legend('\phi=0.09','\phi=0.14','Location','best'); grid on;

% 4.7 Tasa de interés y Gini total
figure('Name','r y Gini (COVID)');
subplot(1,2,1);
bar(categorical({'\phi=0.09','\phi=0.14'}), [covid_low.r, covid_hi.r]);
ylabel('r'); grid on; title('Tasa de interés de equilibrio'); set(gcf,'Color','w');

subplot(1,2,2);
bar(categorical({'\phi=0.09','\phi=0.14'}), [S1.gini(3), S2.gini(3)]);
ylabel('Gini riqueza total'); grid on; title('Desigualdad');

% 4.8 Curva de Lorenz (riqueza total) para ambos escenarios
gT1 = g1(:,1)+g1(:,2);
gT2 = g2(:,1)+g2(:,2);
[as, ix] = sort(a);
da  = a(2)-a(1);

% Escenario φ=0.09
gs1 = gT1(ix)*da;
W1  = sum(gs1);
cumPop1 = cumsum(gs1)/W1;
wealth_pos1 = as - min(0,min(as)) + 1e-12; % desplazar a>=0
cumWealth1 = cumsum(wealth_pos1.*gs1);
L1 = cumWealth1 / max(cumWealth1(end),1e-16);

% Escenario φ=0.14
gs2 = gT2(ix)*da;
W2  = sum(gs2);
cumPop2 = cumsum(gs2)/W2;
wealth_pos2 = as - min(0,min(as)) + 1e-12;
cumWealth2 = cumsum(wealth_pos2.*gs2);
L2 = cumWealth2 / max(cumWealth2(end),1e-16);

% Gráfico Lorenz
figure('Name','Curva de Lorenz de riqueza (COVID: \phi=0.09 vs \phi=0.14)');
plot(cumPop1, L1, 'LineWidth',1.8); hold on;
plot(cumPop2, L2, '--','LineWidth',1.8);
plot([0,1],[0,1],'k:'); grid on; axis square;
xlabel('Fracción población'); ylabel('Fracción riqueza acumulada');
title(sprintf('Lorenz total (Gini: %.3f vs %.3f)', S1.gini(3), S2.gini(3)));
legend('\phi=0.09','\phi=0.14','Location','southeast');
set(gcf,'Color','w');

% (Opcional) Exportar puntos de la Lorenz a CSV
T_lorenz = table;
T_lorenz.cumPop_phi009 = cumPop1(:);
T_lorenz.L_phi009      = L1(:);
T_lorenz.cumPop_phi014 = cumPop2(:);
T_lorenz.L_phi014      = L2(:);
writetable(T_lorenz, './tables/covid_mastransf_lorenz.csv');

disp('Listo. Revisa los CSVs en ./tables y las figuras.');
