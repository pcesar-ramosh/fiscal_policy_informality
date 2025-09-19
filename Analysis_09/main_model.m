%% main_base_only.m
% Ejecuta el modelo BASE (sin shocks) y produce:
% - CSVs con: (i) cuentas fiscales, (ii) mercado de activos, (iii) estadisticas hogar,
%             (iv) prestatarios/prestamistas
% - Graficos: politicas c(a), s(a), densidades g, mercado de activos,
%             cuentas fiscales (ingresos/gastos), Lorenz de riqueza, prestatarios/prestamistas.
clear; clc; close all;

%% ===== 1) PARÁMETROS BASE =====
RRA   = 2.30;         % RRA_I = RRA_F (puedes diferenciar)
rho   = 0.05;
theta = 0.02;         % prima por endeudarse informal
tau_l = 0.15;         % impuesto laboral
tau_c = 0.15;         % IVA
Gov   = 0.05;         % bien público (flujo aditivo en utilidad)
phi   = 0.10;         % transferencias a informales (proporcionales a z1)
z1    = 0.33;         % ingreso informal
z2    = 1.00;         % ingreso formal

eta_target = 0.54;    % estacionaria objetivo de informalidad
p22_bar    = 0.8155;  % persistencia formal (mapea a lambda_2)
Igrid = 700;
amax  = 5.0;
amin  = -0.30*z1;     % LIMITE INFERIOR (verifica que no "ahogue" el consumo)

% bracket y guess para r
r_guess = 0.03; rmin = 0.005; rmax = 0.08;

paramsBase = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov,'phi',phi, ...
    'z1',z1,'z2',z2,'I',Igrid,'amax',amax,'amin',amin, ...
    'r_guess',r_guess,'rmin',rmin,'rmax',rmax, ...
    'p22_bar',p22_bar,'eta_target',eta_target);

%% ===== 2) RESOLVER ECONOMÍA BASE =====
base = huggett_base_covid_function(paramsBase);

% Alias cortos
a   = base.a;     g  = base.g;     c  = base.c;     s = base.s;
popI = base.popI; popF = base.popF;
Y    = base.Y;    Ctot= base.Ctot; r  = base.r;

fb = base.fiscal; % Tl, Tc, Tr, G, B, rB, PB, BB

fprintf('\n== BASE ==\n');
fprintf('r         = %.6f\n', r);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y         = %.6f,   Ctot = %.6f\n', Y, Ctot);
fprintf('PB        = %.6f,   rB   = %.6f,   BB = %.6f (≈0 en estacionario)\n', fb.PB, fb.rB, fb.BB);

%% ===== 3) EXPORTAR CSVs =====
% 3.1 Cuentas fiscales (desagregado ingresos/gastos)
T_fiscal = table;
T_fiscal.scenario = "BASE";
T_fiscal.Tl = fb.Tl;                 % impuesto laboral
T_fiscal.Tc = fb.Tc;                 % impuesto al consumo
T_fiscal.IngresosTot = fb.Tl + fb.Tc;

T_fiscal.Tr = fb.Tr;                 % transferencias
T_fiscal.G  = fb.G;                  % bien público
T_fiscal.rB = fb.rB;                 % servicio de deuda
T_fiscal.GastosTot = fb.Tr + fb.G + fb.rB;

T_fiscal.PB = fb.PB;                 % saldo primario
T_fiscal.B  = fb.B;                  % stock de deuda
T_fiscal.BB = fb.BB;                 % balance global (≈0)

writetable(T_fiscal, 'base_fiscal_breakdown.csv');

% 3.2 Mercado de activos (demanda privada vs oferta pública)
A_priv = sum( (g(:,1)+g(:,2)).*a ) * (a(2)-a(1));  % ∫ a(g1+g2) da
T_assets = table("BASE", A_priv, fb.B, 'VariableNames', {'scenario','A_private','B_public'});
writetable(T_assets, 'base_asset_market.csv');

% 3.3 Estadísticas de hogar (del struct .stats)
S = base.stats;
T_stats = table;
T_stats.scenario = "BASE";
T_stats.r = r; T_stats.popI = popI; T_stats.popF = popF; T_stats.Y = Y; T_stats.Ctot = Ctot;
T_stats.wealth_mean_I = S.wealth_mean(1);
T_stats.wealth_mean_F = S.wealth_mean(2);
T_stats.wealth_mean_T = S.wealth_mean(3);
T_stats.wealth_med_I  = S.wealth_median(1);
T_stats.wealth_med_F  = S.wealth_median(2);
T_stats.wealth_med_T  = S.wealth_median(3);
T_stats.gini_I = S.gini(1);
T_stats.gini_F = S.gini(2);
T_stats.gini_T = S.gini(3);
T_stats.cons_mean_I = S.cons_mean(1);
T_stats.cons_mean_F = S.cons_mean(2);
T_stats.cons_mean_T = S.cons_mean(3);
T_stats.cons_med_I  = S.cons_median(1);
T_stats.cons_med_F  = S.cons_median(2);
T_stats.cons_med_T  = S.cons_median(3);
T_stats.p11_rep     = S.p11;  % proxy de persistencia informal
writetable(T_stats, 'base_household_stats.csv');

% 3.4 Prestatarios / prestamistas
Borr = base.borrowers;
T_borr = table;
T_borr.scenario = "BASE";
T_borr.fracBorrow_I = Borr.fracBorrow(1);
T_borr.fracBorrow_F = Borr.fracBorrow(2);
T_borr.fracLend_I   = Borr.fracLend(1);
T_borr.fracLend_F   = Borr.fracLend(2);
T_borr.volBorrow_I  = Borr.volBorrow(1);  % <0
T_borr.volBorrow_F  = Borr.volBorrow(2);
T_borr.volLend_I    = Borr.volLend(1);    % >0
T_borr.volLend_F    = Borr.volLend(2);
writetable(T_borr, 'base_borrowers_lenders.csv');

disp('Exportados: base_fiscal_breakdown.csv, base_asset_market.csv, base_household_stats.csv, base_borrowers_lenders.csv');

%% ===== 4) GRÁFICOS =====
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% 4.1 Políticas de consumo
figure('Name','Consumo por tipo (BASE)');
plot(a, c(:,1), 'LineWidth',1.6); hold on;
plot(a, c(:,2), 'LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('Consumo c(a)');
legend('Informal','Formal','Location','best'); title('Políticas c(a)');
set(gcf,'Color','w');

% 4.2 Políticas de ahorro
figure('Name','Ahorro por tipo (BASE)');
plot(a, s(:,1), 'LineWidth',1.6); hold on;
plot(a, s(:,2), 'LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('Ahorro s(a)');
legend('Informal','Formal','Location','best'); title('Políticas s(a)');
set(gcf,'Color','w');

% 4.3 Distribución de riqueza g(a) por tipo
figure('Name','Distribucion de riqueza g(a)');
subplot(2,1,1);
bar(a, g(:,1), 'FaceAlpha',0.6, 'EdgeColor','none'); grid on;
xlim([min(a) 1.0]); xlabel('a'); ylabel('g_1(a)'); title('Informales'); set(gcf,'Color','w');
subplot(2,1,2);
bar(a, g(:,2), 'FaceAlpha',0.6, 'EdgeColor','none'); grid on;
xlim([min(a) 1.0]); xlabel('a'); ylabel('g_2(a)'); title('Formales');

% 4.4 Mercado de activos: Demanda privada vs Oferta pública
figure('Name','Mercado de activos (BASE)');
bar(categorical({'Demanda Privada','Oferta B(r)'}), [A_priv, fb.B]);
ylabel('Nivel de activos'); grid on; title(sprintf('Cierre de mercado a r=%.4f', r));
set(gcf,'Color','w');

% 4.5 Cuentas fiscales: Ingresos (Tl,Tc) vs Gastos (Tr,G,rB)
figure('Name','Cuentas fiscales (BASE)');
subplot(1,2,1);
bar(categorical({'Tl','Tc'}), [fb.Tl, fb.Tc]);
ylabel('Ingresos'); title('Desagregado de Ingresos'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'Tr','G','rB'}), [fb.Tr, fb.G, fb.rB]);
ylabel('Gastos'); title('Desagregado de Gastos'); grid on;

% 4.6 Prestatarios / prestamistas: fracciones y volúmenes
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});

figure('Name','Fracciones prestatarios/prestamistas (BASE)');
subplot(1,2,1);
bar(catX, [Borr.fracBorrow(:)]);
ylabel('Frac. prestatarios'); title('a<0'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [Borr.fracLend(:)]);
ylabel('Frac. prestamistas'); title('a>0'); grid on;

figure('Name','Volumenes de deuda/ahorro (BASE)');
subplot(1,2,1);
bar(catX, abs([Borr.volBorrow(:)]));
ylabel('|Deuda agregada|'); title('Suma ponderada en a<0'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [Borr.volLend(:)]);
ylabel('Ahorro agregado'); title('Suma ponderada en a>0'); grid on;

% 4.7 Curva de Lorenz de riqueza (total)
da  = a(2)-a(1);
gT  = g(:,1)+g(:,2);
W   = sum(gT)*da;
[as, ix]   = sort(a);
gs         = gT(ix)*da;
cumPop     = cumsum(gs)/W;
wealth_pos = as - min(0,min(as)) + 1e-12;     % desplazar a >= 0
cumWealth  = cumsum(wealth_pos.*gs);
if cumWealth(end) > 0
    L = cumWealth / cumWealth(end);
    figure('Name','Curva de Lorenz (riqueza total)');
    plot(cumPop, L, 'LineWidth',1.8); hold on;
    plot([0,1],[0,1],'k--'); grid on; axis square;
    xlabel('Fracción población'); ylabel('Fracción riqueza acumulada');
    title(sprintf('Lorenz de riqueza (Gini=%.3f)', base.stats.gini(3)));
    set(gcf,'Color','w');
end

% 4.8 Composición: consumo medio por tipo y participación poblacional
figure('Name','Composicion poblacional y consumo medio');
subplot(1,2,1);
bar(categorical({'Informal','Formal'}), [popI, popF]);
ylabel('Población'); title('Composición poblacional'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'c_I','c_F'}), [base.stats.cons_mean(1), base.stats.cons_mean(2)]);
ylabel('Consumo medio'); title('Consumo medio por tipo'); grid on;

disp('Listo. Revisa los CSVs y las figuras abiertas.');
