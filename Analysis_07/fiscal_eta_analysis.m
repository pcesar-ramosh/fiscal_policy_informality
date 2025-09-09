%% fiscal_eta_analysis.m
% Analiza la política fiscal por nivel de informalidad (eta):
% - Recaudación:    Tl (ingreso laboral), Tc (consumo)
% - Gasto:          Transferencias (Tr), Bien público (G), Intereses (r*B)
% - Saldo Primario: PB = Tl + Tc - (G + Tr)
% - Deuda endógena: B = PB / r   (por construcción, BB = PB - rB ≈ 0)
% Exporta todo a CSV y grafica componentes vs eta.
%
% NOTA DE CONSISTENCIA: Los parámetros deben coincidir con los del solver
% (huggett_Equi_INFO_function). Si cambias allí, cambia aquí también.
% -----------------------------------------------------------------------
clear; clc; close all;

% --- Grilla de informalidad (puedes editar)
eta_grid = [0.20 0.30 0.40 0.50 0.60 0.64 0.70 0.80 0.90];

% --- Preferencias base (mismo valor usado en los ejemplos anteriores)
RRA = 4.20;

% -----------------------------------------------------------------------
% PARAMETROS DE POLITICA (DEBEN COINCIDIR CON LOS DEL SOLVER)
% -----------------------------------------------------------------------
tau_l  = 0.15;   % impuesto laboral (solo formales)
tau_c  = 0.10;   % IVA (ambos)
phi    = 0.15;   % transferencias (solo informales, proporcionales a y1)
Gov    = 0.07;   % bien público (entra a utilidad)
z1     = 0.33;   % ingreso informal
z2     = 1.00;   % ingreso formal

% -----------------------------------------------------------------------
% CORRER EL SOLVER (una sola llamada con toda la grilla eta_grid)
% -----------------------------------------------------------------------
[r_vec, ~, pop1_vector, ~, ~, GDistribution, a, Distribution, ~, ~, ...
 C, ~, ~, ~, I, amin] = huggett_Equi_INFO_function(RRA, eta_grid);

% Chequeos básicos
if any(pop1_vector(:) ~= eta_grid(:))
    warning('pop1_vector devuelto no coincide con eta_grid.');
end

% Paso de grilla
da = (a(end)-a(1))/(numel(a)-1);

% -----------------------------------------------------------------------
% AGREGACION FISCAL POR CADA ETA
% -----------------------------------------------------------------------
nE  = numel(eta_grid);
eta = pop1_vector(:);
r   = r_vec(:);

% Pre-alocación
PopI = zeros(nE,1);  % masa informal
PopF = zeros(nE,1);  % masa formal
Ctot = zeros(nE,1);  % consumo agregado (base IVA)
Tl   = zeros(nE,1);  % impuesto laboral
Tc   = zeros(nE,1);  % impuesto al consumo
Tr   = zeros(nE,1);  % transferencias
Gx   = zeros(nE,1);  % gasto en bien público
B    = zeros(nE,1);  % deuda pública endógena
rB   = zeros(nE,1);  % servicio de deuda
PB   = zeros(nE,1);  % saldo primario
BB   = zeros(nE,1);  % balance total (~0 por construcción)
Y    = zeros(nE,1);  % “PIB” proxy: ingreso laboral exógeno

for j = 1:nE
    g = Distribution{j};           % [I x 2] densidades por a (informal, formal)
    c = C{j};                      % [I x 2] consumo por tipo (política)

    % Poblaciones por tipo
    PopI(j) = sum(g(:,1))*da;
    PopF(j) = sum(g(:,2))*da;

    % Consumo agregado (base del IVA)
    Ctot(j) = sum(c(:,1).*g(:,1))*da + sum(c(:,2).*g(:,2))*da;

    % Recaudación
    Tl(j) = tau_l * z2 * PopF(j);
    Tc(j) = tau_c * Ctot(j);

    % Gasto
    Tr(j) = phi  * z1 * PopI(j);
    Gx(j) = Gov  * (PopI(j)+PopF(j));

    % PIB proxy (ingreso exógeno, sin contar impuestos/transferencias)
    Y(j)  = z1*PopI(j) + z2*PopF(j);

    % Saldo primario y deuda endógena (cierra con r)
    PB(j) = Tl(j) + Tc(j) - (Gx(j) + Tr(j));
    B(j)  = PB(j) / max(r(j), 1e-9);
    rB(j) = r(j) * B(j);

    % Balance total (≈ 0 numéricamente)
    BB(j) = PB(j) - rB(j);
end

% --- Razones sobre PIB (útiles para análisis)
Tl_Y  = Tl ./ max(Y,1e-12);
Tc_Y  = Tc ./ max(Y,1e-12);
Tr_Y  = Tr ./ max(Y,1e-12);
Gx_Y  = Gx ./ max(Y,1e-12);
PB_Y  = PB ./ max(Y,1e-12);
B_Y   = B  ./ max(Y,1e-12);
rB_Y  = rB ./ max(Y,1e-12);

% -----------------------------------------------------------------------
% EXPORT A CSV
% -----------------------------------------------------------------------
T = table(eta, r, PopI, PopF, Y, Ctot, Tl, Tc, Tr, Gx, PB, B, rB, BB, ...
          Tl_Y, Tc_Y, Tr_Y, Gx_Y, PB_Y, B_Y, rB_Y, ...
    'VariableNames', { ...
        'eta','r','PopI','PopF','Y','Ctot','Tl','Tc','Tr','G','PB','B','rB','BB', ...
        'Tl_Y','Tc_Y','Tr_Y','G_Y','PB_Y','B_Y','rB_Y'});

writetable(T, 'fiscal_by_eta.csv');
disp('Exportado: fiscal_by_eta.csv');

% -----------------------------------------------------------------------
% GRAFICOS
% -----------------------------------------------------------------------
% Evitar problemas de interpretes con acentos en LaTeX:
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% 1) Recaudación vs eta
figure('Name','Recaudacion vs eta');
plot(eta, Tl,'-o','LineWidth',1.4); hold on;
plot(eta, Tc,'-s','LineWidth',1.4);
grid on; xlabel('\eta'); ylabel('Recaudacion');
legend('T_l (ingreso laboral)','T_c (consumo)','Location','best');
title('Recaudacion: T_l y T_c'); set(gcf,'Color','w');

% 2) Gasto vs eta
figure('Name','Gasto vs eta');
plot(eta, Tr,'-o','LineWidth',1.4); hold on;
plot(eta, Gx,'-s','LineWidth',1.4);
plot(eta, rB,'-^','LineWidth',1.4);
grid on; xlabel('\eta'); ylabel('Gasto');
legend('Transferencias','G (bien publico)','Intereses rB','Location','best');
title('Gasto: Tr, G y rB'); set(gcf,'Color','w');

% 3) Saldo primario y balance total
figure('Name','Saldos fiscales');
yyaxis left
plot(eta, PB,'-o','LineWidth',1.6); ylabel('Saldo primario (PB)'); grid on;
yyaxis right
plot(eta, BB,'-s','LineWidth',1.6); ylabel('Balance total (BB)'); 
xlabel('\eta'); title('PB y BB'); set(gcf,'Color','w');

% 4) Deuda y tasa de interes
figure('Name','Deuda y tasa');
yyaxis left
plot(eta, B,'-o','LineWidth',1.6); ylabel('Deuda B');
yyaxis right
plot(eta, r,'-s','LineWidth',1.6); ylabel('r(\eta)');
xlabel('\eta'); grid on; title('B(\eta) y r(\eta)'); set(gcf,'Color','w');

% 5) Stacked bars: Recaudacion vs Gasto (en niveles)
figure('Name','Stacked: Recaudacion y Gasto');
subplot(1,2,1)
bar(eta, [Tl Tc], 'stacked'); 
xlabel('\eta'); ylabel('Recaudacion'); grid on;
legend('T_l','T_c','Location','best'); title('Recaudacion (stacked)'); set(gca,'Box','off');

subplot(1,2,2)
bar(eta, [Gx Tr rB], 'stacked');
xlabel('\eta'); ylabel('Gasto'); grid on;
legend('G','Tr','rB','Location','best'); title('Gasto (stacked)'); set(gca,'Box','off');
set(gcf,'Color','w');

% 6) Version en % del PIB
figure('Name','Cuentas fiscales / PIB');
subplot(2,2,1)
plot(eta, Tl_Y,'-o', eta, Tc_Y,'-s','LineWidth',1.4); grid on;
xlabel('\eta'); ylabel('Recaudacion / Y'); legend('T_l/Y','T_c/Y','Location','best');
title('Recaudacion relativa'); set(gca,'Box','off');

subplot(2,2,2)
plot(eta, Tr_Y,'-o', eta, Gx_Y,'-s','LineWidth',1.4); grid on;
xlabel('\eta'); ylabel('Gasto primario / Y'); legend('Tr/Y','G/Y','Location','best');
title('Gasto primario relativo'); set(gca,'Box','off');

subplot(2,2,3)
plot(eta, PB_Y,'-^','LineWidth',1.4); grid on;
xlabel('\eta'); ylabel('PB / Y'); title('Saldo primario relativo'); set(gca,'Box','off');

subplot(2,2,4)
plot(eta, B_Y,'-o', eta, rB_Y,'-s','LineWidth',1.4); grid on;
xlabel('\eta'); ylabel('Deuda e intereses / Y'); legend('B/Y','rB/Y','Location','best');
title('Deuda e intereses relativos'); set(gca,'Box','off');
set(gcf,'Color','w');

disp('Listo: gráficos y fiscal_by_eta.csv generados.');
