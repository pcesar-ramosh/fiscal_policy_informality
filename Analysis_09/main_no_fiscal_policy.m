%% main_no_fiscal_policy.m
% Modelo SIN política fiscal:
%  - tau_l=tau_c=phi=Gov=0, B=0
%  - r cierra ∫ a g(a) da = 0
% Exporta CSV y grafica lo mismo que el base (ingresos/gastos serán 0).

clear; clc; close all;

%% ===== 1) PARÁMETROS =====
RRA   = 3.40;
rho   = 0.05;
theta = 0.02;
z1    = 0.33;
z2    = 1.00;

eta_target = 0.654;
p22_bar    = 0.8155;

Igrid = 700;
amax  = 5.0;
amin  = -0.30*z1;

r_guess = 0.03; rmin = 0.005; rmax = 0.08;

params = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'z1',z1,'z2',z2,'I',Igrid,'amax',amax,'amin',amin, ...
    'r_guess',r_guess,'rmin',rmin,'rmax',rmax, ...
    'p22_bar',p22_bar,'eta_target',eta_target);

%% ===== 2) Resolver =====
res = hugguet_no_fiscal_policy(params);

% Aliases
a   = res.a;  g = res.g;  c = res.c;  s = res.s;
popI = res.popI; popF = res.popF; Y = res.Y; Ctot = res.Ctot; r = res.r;
fb = res.fiscal;   % serán ceros (Tl,Tc,Tr,G,B,rB,PB,BB)

fprintf('\n== NO-FISCAL ==\n');
fprintf('r         = %.6f\n', r);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y         = %.6f,   Ctot = %.6f\n', Y, Ctot);

%% ===== 3) CSV =====
if ~exist('./tables','dir'), mkdir('./tables'); end

% 3.1 "Cuentas" fiscales (serán 0)
T_fiscal = table;
T_fiscal.scenario = "NO_FISCAL";
T_fiscal.Tl = fb.Tl;  T_fiscal.Tc = fb.Tc;  T_fiscal.IngresosTot = fb.Tl + fb.Tc;
T_fiscal.Tr = fb.Tr;  T_fiscal.G  = fb.G;   T_fiscal.rB = fb.rB;
T_fiscal.GastosTot = fb.Tr + fb.G + fb.rB;
T_fiscal.PB = fb.PB;  T_fiscal.B  = fb.B;   T_fiscal.BB = fb.BB;
writetable(T_fiscal, './tables/nofiscal_fiscal_breakdown.csv');

% 3.2 Mercado de activos (privado vs público=0)
da = a(2)-a(1);
A_priv = sum( (g(:,1)+g(:,2)).*a ) * da;
T_assets = table("NO_FISCAL", A_priv, 0, 'VariableNames', {'scenario','A_private','B_public'});
writetable(T_assets, './tables/nofiscal_asset_market.csv');

% 3.3 Estadísticas
S = res.stats;
T_stats = table;
T_stats.scenario = "NO_FISCAL";
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
T_stats.p11_rep     = S.p11;
writetable(T_stats, './tables/nofiscal_household_stats.csv');

% 3.4 Prestatarios / prestamistas
Borr = res.borrowers;
T_borr = table;
T_borr.scenario = "NO_FISCAL";
T_borr.fracBorrow_I = Borr.fracBorrow(1);
T_borr.fracBorrow_F = Borr.fracBorrow(2);
T_borr.fracLend_I   = Borr.fracLend(1);
T_borr.fracLend_F   = Borr.fracLend(2);
T_borr.volBorrow_I  = Borr.volBorrow(1);
T_borr.volBorrow_F  = Borr.volBorrow(2);
T_borr.volLend_I    = Borr.volLend(1);
T_borr.volLend_F    = Borr.volLend(2);
writetable(T_borr, './tables/nofiscal_borrowers_lenders.csv');

disp('Exportados: nofiscal_* CSVs');

%% ===== 4) GRÁFICOS =====
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% 4.1 c(a)
figure('Name','Consumo por tipo (SIN fiscal)');
plot(a, c(:,1), 'LineWidth',1.6); hold on;
plot(a, c(:,2), 'LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('Consumo c(a)');
legend('Informal','Formal','Location','best'); title('Políticas c(a)'); set(gcf,'Color','w');

% 4.2 s(a)
figure('Name','Ahorro por tipo (SIN fiscal)');
plot(a, s(:,1), 'LineWidth',1.6); hold on;
plot(a, s(:,2), 'LineWidth',1.6);
yline(0,'k:'); grid on; xlabel('Activos a'); ylabel('Ahorro s(a)');
legend('Informal','Formal','Location','best'); title('Políticas s(a)'); set(gcf,'Color','w');

% 4.3 g(a)
figure('Name','Distribucion de riqueza g(a) (SIN fiscal)');
subplot(2,1,1);
bar(a, g(:,1), 'FaceAlpha',0.6, 'EdgeColor','none'); grid on;
xlim([min(a) 1.0]); xlabel('a'); ylabel('g_1(a)'); title('Informales'); set(gcf,'Color','w');
subplot(2,1,2);
bar(a, g(:,2), 'FaceAlpha',0.6, 'EdgeColor','none'); grid on;
xlim([min(a) 1.0]); xlabel('a'); ylabel('g_2(a)'); title('Formales');

% 4.4 Mercado de activos (público=0)
figure('Name','Mercado de activos (SIN fiscal)');
bar(categorical({'Demanda Privada','Oferta Pública'}), [A_priv, 0]);
ylabel('Nivel de activos'); grid on; title(sprintf('Cierre privado a r=%.4f', r)); set(gcf,'Color','w');

% 4.5 Prestatarios/prestamistas
catX = categorical({'Informal','Formal'}); catX = reordercats(catX,{'Informal','Formal'});
figure('Name','Fracciones prestatarios/prestamistas (SIN fiscal)');
subplot(1,2,1);
bar(catX, [Borr.fracBorrow(:)]); ylabel('Frac. prestatarios'); title('a<0'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [Borr.fracLend(:)]);   ylabel('Frac. prestamistas'); title('a>0'); grid on;

figure('Name','Volumenes de deuda/ahorro (SIN fiscal)');
subplot(1,2,1);
bar(catX, abs([Borr.volBorrow(:)])); ylabel('|Deuda agregada|'); title('Suma ponderada en a<0'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(catX, [Borr.volLend(:)]);    ylabel('Ahorro agregado'); title('Suma ponderada en a>0'); grid on;

% 4.6 Curva de Lorenz (riqueza total)
gT  = g(:,1)+g(:,2); W = sum(gT)*da;
[as, ix]   = sort(a); gs = gT(ix)*da; cumPop = cumsum(gs)/W;
wealth_pos = as - min(0,min(as)) + 1e-12;
cumWealth  = cumsum(wealth_pos.*gs);
if cumWealth(end)>0
    L = cumWealth / cumWealth(end);
    figure('Name','Curva de Lorenz (SIN fiscal)');
    plot(cumPop, L, 'LineWidth',1.8); hold on;
    plot([0,1],[0,1],'k--'); grid on; axis square;
    xlabel('Fracción población'); ylabel('Fracción riqueza acumulada');
    title(sprintf('Lorenz riqueza (Gini=%.3f)', res.stats.gini(3))); set(gcf,'Color','w');
end

% 4.7 Composición y consumo medio
figure('Name','Composicion y consumo medio (SIN fiscal)');
subplot(1,2,1);
bar(categorical({'Informal','Formal'}), [popI, popF]);
ylabel('Población'); title('Composición poblacional'); grid on; set(gcf,'Color','w');
subplot(1,2,2);
bar(categorical({'c_I','c_F'}), [res.stats.cons_mean(1), res.stats.cons_mean(2)]);
ylabel('Consumo medio'); title('Consumo medio por tipo'); grid on;

disp('Listo (SIN política fiscal). Revisa los CSVs y figuras.');
