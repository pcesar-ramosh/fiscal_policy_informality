%% ==============================================================
%  main_viz_formal_informal_CAL.m
%  Visualización integral (formal vs informal):
%   - c(a), s(a) del agente típico
%   - g(a) por ocupación
%   - Oferta/Demanda de bonos por r (agregado y por ocupación)
%   - S_total(r) y r*_common
%  Sept 2025
%% ==============================================================

clear; clc; close all;

%% -----------------------
%  Heterogeneidad
%% -----------------------
n_agents  = 20;
s_min     = 3.15; s_max = 5.30;
eta_vector = 0.75 * ones(1, n_agents);
sI_vector1 = linspace(s_min, s_max, n_agents);
sF_vector1 = 5.30 * ones(1, n_agents);

[~, j_star] = min(abs(sI_vector1 - median(sI_vector1))); % agente típico

%% -----------------------
%  Configuración: BASE y COVID
%% -----------------------
cfg_base = struct( ...
  'scenario',"baseline", 'psi_I',1.0, 'psi_F',1.0, ...
  'transfer_multiplier',1.0, 'keep_transfers_level',true, ...
  'amin_policy',"absolute", 'amin_abs',-1.00, ...   % crédito más laxo
  'theta',0.005, 'phi',0.05, ...
  'taxF',0.10, 'taxI',0.00, 'tauc',0.00 ...
);

cfg_covid = cfg_base;
cfg_covid.scenario            = "covid_uptransfer";
cfg_covid.psi_I               = 0.75;
cfg_covid.psi_F               = 0.85;
cfg_covid.transfer_multiplier = 1.25;

%% -----------------------
%  r* común agregado por bisección
%% -----------------------
[r_common_base, hist_base] = find_r_common_CAL(eta_vector, sI_vector1, sF_vector1, cfg_base, -0.01, 0.06);
[r_common_covid, hist_cvd] = find_r_common_CAL(eta_vector, sI_vector1, sF_vector1, cfg_covid, -0.01, 0.06);

fprintf('\n=== r*_common (Base)=%.4f | r*_common (COVID)=%.4f ===\n', r_common_base, r_common_covid);

%% -----------------------
%  Objeto a r* (para políticas y g del típico)
%% -----------------------
[~, ~, gB, cB, a_grid] = huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg_base, r_common_base);
[~, ~, gC, cC, ~     ] = huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg_covid, r_common_covid);

I   = numel(a_grid); da = a_grid(2)-a_grid(1);
amin= a_grid(1); amax = a_grid(end);

% Recursos para reconstruir ahorro s(a) del agente típico
% BASE
theta = cfg_base.theta; taxI = cfg_base.taxI; taxF = cfg_base.taxF;
z1b = 0.33; z2b = 1.00;
phi = cfg_base.phi; kappaB = 1.0;
TransferB = (cfg_base.keep_transfers_level)*(kappaB*phi*z1b) + (~cfg_base.keep_transfers_level)*0; %#ok<NASGU>

rrB = r_common_base + theta*(a_grid<0);
resI_B = (1-taxI)*z1b + rrB.*a_grid + kappaB*phi*z1b;
resF_B = (1-taxF)*z2b + r_common_base*a_grid;

% COVID
kappaC = cfg_covid.transfer_multiplier; z1c = cfg_covid.psi_I*z1b; z2c = cfg_covid.psi_F*z2b;
TransferC = (cfg_covid.keep_transfers_level)*(kappaC*phi*z1b) + (~cfg_covid.keep_transfers_level)*(kappaC*phi*z1c); %#ok<NASGU>

rrC = r_common_covid + theta*(a_grid<0);
resI_C = (1-taxI)*z1c + rrC.*a_grid + ((cfg_covid.keep_transfers_level)*(kappaC*phi*z1b) + (~cfg_covid.keep_transfers_level)*(kappaC*phi*z1c));
resF_C = (1-taxF)*z2c + r_common_covid*a_grid;

% Políticas del agente típico
c_base_typ  = cB{j_star};     % I x 2
c_covid_typ = cC{j_star};     % I x 2

s_base_inf  = resI_B - c_base_typ(:,1);
s_base_for  = resF_B - c_base_typ(:,2);
s_covid_inf = resI_C - c_covid_typ(:,1);
s_covid_for = resF_C - c_covid_typ(:,2);

%% -----------------------
%  Gráfico 1: c(a) y s(a) — Agente típico (formal/inf)
%% -----------------------
figure('Name','Agente típico: consumo y ahorro','Color','w','Position',[100 100 1100 430]);

subplot(1,2,1); hold on; grid on; box on;
plot(a_grid, c_base_typ(:,1), 'r--', 'LineWidth',1.6);
plot(a_grid, c_base_typ(:,2), 'b--', 'LineWidth',1.6);
plot(a_grid, c_covid_typ(:,1),'r-',  'LineWidth',1.6);
plot(a_grid, c_covid_typ(:,2),'b-',  'LineWidth',1.6);
xline(amin, ':k', 'a_{min}', 'LabelVerticalAlignment','bottom');
title('Consumo c(a) | Típico'); xlabel('Activos a'); ylabel('c');
legend('Inf Base','For Base','Inf COVID','For COVID','Location','SouthEast');

subplot(1,2,2); hold on; grid on; box on;
plot(a_grid, s_base_inf,  'r--', 'LineWidth',1.6);
plot(a_grid, s_base_for,  'b--', 'LineWidth',1.6);
plot(a_grid, s_covid_inf, 'r-',  'LineWidth',1.6);
plot(a_grid, s_covid_for, 'b-',  'LineWidth',1.6);
yline(0,'k:');
xline(amin, ':k', 'a_{min}', 'LabelVerticalAlignment','bottom');
title('Ahorro s(a)=Ingreso+r·a−c | Típico'); xlabel('Activos a'); ylabel('s(a)');
legend('Inf Base','For Base','Inf COVID','For COVID','Location','SouthEast');

%% -----------------------
%  Gráfico 2: g(a) por ocupación a r*_common
%% -----------------------
gB_inf = gB{j_star}(:,1); gB_for = gB{j_star}(:,2);
gC_inf = gC{j_star}(:,1); gC_for = gC{j_star}(:,2);

figure('Name','Distribución de riqueza g(a) por ocupación (típico)','Color','w','Position',[120 120 1050 430]);
subplot(1,2,1); hold on; grid on; box on;
plot(a_grid, gB_inf, 'r--','LineWidth',1.6);
plot(a_grid, gB_for, 'b--','LineWidth',1.6);
plot(a_grid, gC_inf, 'r-','LineWidth',1.6);
plot(a_grid, gC_for, 'b-','LineWidth',1.6);
xline(amin, ':k', 'a_{min}');
title('g(a) | Típico'); xlabel('Activos a'); ylabel('densidad');
legend('Inf Base','For Base','Inf COVID','For COVID','Location','NorthEast');

subplot(1,2,2); hold on; grid on; box on;
bar_data = [sum(max(a_grid,0).*gB_inf)*da, sum(max(a_grid,0).*gB_for)*da; ...
            sum(-min(a_grid,0).*gB_inf)*da, sum(-min(a_grid,0).*gB_for)*da];
bar(bar_data,'stacked'); xticklabels({'Oferta (a>0)','Demanda (-a<0)'}); 
title('Oferta/Demanda por ocupación (@ r^*_{Base})'); ylabel('masa de bonos');

%% -----------------------
%  Gráfico 3: Oferta/Demanda de bonos vs r (agregado y por ocupación)
%% -----------------------
rgrid = linspace(-0.01,0.06,25);
SupplyTot_B = zeros(size(rgrid)); DemandTot_B = zeros(size(rgrid));
SupplyI_B   = zeros(size(rgrid)); DemandI_B   = zeros(size(rgrid));
SupplyF_B   = zeros(size(rgrid)); DemandF_B   = zeros(size(rgrid));

SupplyTot_C = zeros(size(rgrid)); DemandTot_C = zeros(size(rgrid));
SupplyI_C   = zeros(size(rgrid)); DemandI_C   = zeros(size(rgrid));
SupplyF_C   = zeros(size(rgrid)); DemandF_C   = zeros(size(rgrid));

for k=1:numel(rgrid)
    r_ = rgrid(k);

    % BASE
    [~, ~, gallB, ~, a_g] = huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg_base, r_);
    da_ = a_g(2)-a_g(1);

    supI=0; demI=0; supF=0; demF=0;
    for j=1:n_agents
        gI = gallB{j}(:,1); gF = gallB{j}(:,2);
        supI = supI + sum(max(a_g,0).*gI)*da_;
        demI = demI + sum(-min(a_g,0).*gI)*da_;
        supF = supF + sum(max(a_g,0).*gF)*da_;
        demF = demF + sum(-min(a_g,0).*gF)*da_;
    end
    SupplyI_B(k)=supI; DemandI_B(k)=demI; SupplyF_B(k)=supF; DemandF_B(k)=demF;
    SupplyTot_B(k)=supI+supF; DemandTot_B(k)=demI+demF;

    % COVID
    [~, ~, gallC, ~, a_g] = huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg_covid, r_);
    da_ = a_g(2)-a_g(1);

    supI=0; demI=0; supF=0; demF=0;
    for j=1:n_agents
        gI = gallC{j}(:,1); gF = gallC{j}(:,2);
        supI = supI + sum(max(a_g,0).*gI)*da_;
        demI = demI + sum(-min(a_g,0).*gI)*da_;
        supF = supF + sum(max(a_g,0).*gF)*da_;
        demF = demF + sum(-min(a_g,0).*gF)*da_;
    end
    SupplyI_C(k)=supI; DemandI_C(k)=demI; SupplyF_C(k)=supF; DemandF_C(k)=demF;
    SupplyTot_C(k)=supI+supF; DemandTot_C(k)=demI+demF;
end

figure('Name','Mercado de bonos: oferta y demanda vs r','Color','w','Position',[140 140 1200 680]);

% Panel (a): Agregado
subplot(2,2,1); hold on; grid on; box on;
plot(rgrid, SupplyTot_B, 'b-', 'LineWidth',1.8);
plot(rgrid, DemandTot_B, 'r-', 'LineWidth',1.8);
plot(rgrid, SupplyTot_C, 'b--', 'LineWidth',1.8);
plot(rgrid, DemandTot_C, 'r--', 'LineWidth',1.8);
xline(r_common_base, 'k:','Base r^*'); xline(r_common_covid,'k--','COVID r^*');
title('(a) Agregado'); xlabel('r'); ylabel('Bonos (masa)');
legend('Oferta Base','Demanda Base','Oferta COVID','Demanda COVID','Location','NorthWest');

% Panel (b): Oferta por ocupación
subplot(2,2,2); hold on; grid on; box on;
plot(rgrid, SupplyI_B, 'm-', 'LineWidth',1.6);
plot(rgrid, SupplyF_B, 'c-', 'LineWidth',1.6);
plot(rgrid, SupplyI_C, 'm--','LineWidth',1.6);
plot(rgrid, SupplyF_C, 'c--','LineWidth',1.6);
xline(r_common_base, 'k:'); xline(r_common_covid,'k--');
title('(b) Oferta por ocupación'); xlabel('r'); ylabel('a>0');
legend('Inf Base','For Base','Inf COVID','For COVID','Location','NorthWest');

% Panel (c): Demanda por ocupación
subplot(2,2,3); hold on; grid on; box on;
plot(rgrid, DemandI_B, 'm-', 'LineWidth',1.6);
plot(rgrid, DemandF_B, 'c-', 'LineWidth',1.6);
plot(rgrid, DemandI_C, 'm--','LineWidth',1.6);
plot(rgrid, DemandF_C, 'c--','LineWidth',1.6);
xline(r_common_base, 'k:'); xline(r_common_covid,'k--');
title('(c) Demanda por ocupación'); xlabel('r'); ylabel('-a<0');

% Panel (d): Exceso agregado S_total(r)
subplot(2,2,4); hold on; grid on; box on;
StotB = SupplyTot_B - DemandTot_B;
StotC = SupplyTot_C - DemandTot_C;
plot(rgrid, StotB, 'k-', 'LineWidth',1.8);
plot(rgrid, StotC, 'k--','LineWidth',1.8);
yline(0,'k:');
xline(r_common_base, 'r:','Base r^*'); xline(r_common_covid,'r--','COVID r^*');
title('(d) Exceso S_{total}(r)'); xlabel('r'); ylabel('S_{total}');
legend('Base','COVID','Location','NorthWest');

disp('Listo: figuras generadas (consumo/ahorro, g(a), oferta/demanda por r, S_total).');
