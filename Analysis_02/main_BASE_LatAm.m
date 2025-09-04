%% ======================== BASE (LatAm) =========================
%  Objetivo: r* = 3% con política fiscal BASE
%  Fallback: si no hay sign-change en phi, auto-calibra {a_min,theta,tau_c,tau_l,phi_max}
% ================================================================

clear; clc; close all;

%% Heterogeneidad y transición ocupacional
n_agents  = 20;
sI_lo     = 3.00;  sI_hi    = 5.00;
sI_vector = linspace(sI_lo, sI_hi, n_agents);
sF_vector = 4.00 * ones(1, n_agents);

share_informal = 0.55;   % π_I (LatAm emergente)
p22  = 0.75;             % persistencia formal
la2  = -log(p22);        % F->F
la1  = (1-share_informal)*la2/share_informal; % I->I
transitions.la1 = la1; transitions.la2 = la2;

%% Calibración BASE (LatAm) — punto de partida
cfg = struct();
cfg.r_target  = 0.03;      % objetivo r*
cfg.tau_l     = 0.15;      % impuesto laboral solo FORMAL
cfg.tau_c     = 0.08;      % IVA efectivo
cfg.theta     = 0.01;      % spread por endeudamiento
cfg.amin_abs  = -1.50;     % límite de deuda
cfg.amin_policy = "absolute";
cfg.z1        = 0.33;      % ingreso informal
cfg.z2        = 1.00;      % ingreso formal
cfg.rho       = 0.05;      % descuento
cfg.I         = 700;       % puntos de malla
cfg.amax      = 5.0;
cfg.xiG       = 0.05;      % peso del bien público en utilidad (aditivo)
cfg.epsG      = 1e-8;      % para log(G+eps)

% Intervalo base para phi
phi_lo = 0.00;  phi_hi = 0.60;  tol_phi = 1e-5;

%% === 1) Intento directo: calibrar phi para r_target ===
try
    [phi_star, S_at_target] = calibrate_phi_for_target_r_BASE( ...
        share_informal, sI_vector, sF_vector, transitions, cfg, phi_lo, phi_hi, tol_phi);

    cfg.phi = phi_star;
    fprintf('\n[BASE] phi* calibrado = %.6f  (S_total@r=%.4f = %.3e)\n', phi_star, cfg.r_target, S_at_target);

catch ME
    fprintf('\n[WARN] Calibración directa de phi no encontró cambio de signo.\n%s\n', ME.message);

    %% === 2) Fallback automático: ajustes mínimos para habilitar el cruce ===
    auto_opts = struct();
    % Listas de candidatos (orden de “preferencia” para cambios mínimos)
    auto_opts.amin_list  = [cfg.amin_abs, -2.00, -2.50];
    auto_opts.theta_list = [cfg.theta, 0.005, 0.002];
    auto_opts.tauc_list  = [cfg.tau_c, 0.04, 0.00];
    auto_opts.taul_list  = [cfg.tau_l, 0.20];   % opcional
    auto_opts.phimax_list= [phi_hi, 0.80, 1.00];

    [cfg, phi_star, S_at_target] = calibrate_baseline_to_target_r_AUTO_BASE( ...
        share_informal, sI_vector, sF_vector, transitions, cfg, auto_opts, tol_phi);

    fprintf('\n[AUTO] Calibración mínima encontrada:\n');
    fprintf('      amin_abs = %.2f | theta = %.3f | tau_c = %.3f | tau_l = %.3f | phi* = %.6f\n', ...
        cfg.amin_abs, cfg.theta, cfg.tau_c, cfg.tau_l, phi_star);
    fprintf('      S_total@r=%.4f = %.3e (≈ 0)\n', cfg.r_target, S_at_target);
end

%% === Resolver a r_target con phi final ===
r_fixed = cfg.r_target;
[S_by_agent, S_total, g_out, c_out, a_grid, G_val, Ec_total] = ...
    huggett_S_given_r_BASE(share_informal, sI_vector, sF_vector, transitions, cfg, r_fixed);

da = a_grid(2)-a_grid(1); amin=a_grid(1); amax=a_grid(end);
[~, j_star] = min(abs(sI_vector - median(sI_vector)));

%% === Reportes fiscales y de equilibrio ===
fprintf('\n=== Reporte fiscal (BASE @ r=%.4f) ===\n', r_fixed);
pi_I = la2/(la1+la2); pi_F = la1/(la1+la2);
fprintf('E[c] = %.4f | G = %.4f | tau_l*z2*pi_F = %.4f | phi*z1*pi_I = %.4f\n', ...
    Ec_total, G_val, cfg.tau_l*cfg.z2*pi_F, cfg.phi*cfg.z1*pi_I);

fprintf('\n=== Clearing agregado ===\n');
fprintf('S_total = %.3e (debe ser ~0)\n', S_total);

%% ============ GRÁFICOS: comparación formal vs informal ============
% (1) c(a) y s(a) del agente típico (@ r_target)
figure('Name','BASE: Consumo y Ahorro (Típico)','Color','w','Position',[100 100 1100 430]);

% Reconstruir recursos para s(a)
r = r_fixed; th = cfg.theta; tau_c = cfg.tau_c; tau_l = cfg.tau_l;
z1=cfg.z1; z2=cfg.z2; phi=cfg.phi;

rr_vec = r + th*(a_grid<0);
resI = (1-0)*z1 + rr_vec.*a_grid + phi*z1;   % ingreso neto + transfer
resF = (1-tau_l)*z2 + r*a_grid;              % ingreso formal neto

c_typ = c_out{j_star}; % I x 2
sI = resI - (1+tau_c).*c_typ(:,1);
sF = resF - (1+tau_c).*c_typ(:,2);

subplot(1,2,1); hold on; grid on; box on;
plot(a_grid, c_typ(:,1), 'r-', 'LineWidth',1.6);
plot(a_grid, c_typ(:,2), 'b-', 'LineWidth',1.6);
xline(amin, ':k'); title('Consumo c(a) | Agente típico');
xlabel('Activos a'); ylabel('c'); legend('Informal','Formal','Location','SouthEast');

subplot(1,2,2); hold on; grid on; box on;
plot(a_grid, sI, 'r-', 'LineWidth',1.6);
plot(a_grid, sF, 'b-', 'LineWidth',1.6);
yline(0,'k:'); xline(amin, ':k');
title('Ahorro s(a) | Agente típico'); xlabel('Activos a'); ylabel('s(a)');
legend('Informal','Formal','Location','SouthEast');

% (2) g(a) por ocupación en el típico
g_typ = g_out{j_star};
figure('Name','BASE: g(a) por ocupación (Típico)','Color','w','Position',[120 120 1050 430]);
subplot(1,2,1); hold on; grid on; box on;
plot(a_grid, g_typ(:,1), 'r-', 'LineWidth',1.6);
plot(a_grid, g_typ(:,2), 'b-', 'LineWidth',1.6);
xline(amin, ':k'); title('g(a) | Típico'); xlabel('a'); ylabel('densidad');
legend('Informal','Formal','Location','NorthEast');

% (3) Oferta/Demanda por r y S_total(r)
rgrid = linspace(-0.005, 0.06, 25);
[SupplyI, DemandI, SupplyF, DemandF, SupplyTot, DemandTot, Stot] = ...
    market_curves_by_r_BASE(share_informal, sI_vector, sF_vector, transitions, cfg, rgrid);

figure('Name','BASE: Mercado de bonos por r','Color','w','Position',[140 140 1200 680]);
subplot(2,2,1); hold on; grid on; box on;
plot(rgrid, SupplyTot, 'b-', 'LineWidth',1.8);
plot(rgrid, DemandTot, 'r-', 'LineWidth',1.8);
xline(cfg.r_target,'k--','r^* target');
title('(a) Agregado'); xlabel('r'); ylabel('Bonos'); legend('Oferta','Demanda','Location','NorthWest');

subplot(2,2,2); hold on; grid on; box on;
plot(rgrid, SupplyI, 'm-', 'LineWidth',1.6);
plot(rgrid, SupplyF, 'c-', 'LineWidth',1.6);
xline(cfg.r_target,'k--');
title('(b) Oferta por ocupación'); xlabel('r'); ylabel('a>0'); legend('Informal','Formal');

subplot(2,2,3); hold on; grid on; box on;
plot(rgrid, DemandI, 'm-', 'LineWidth',1.6);
plot(rgrid, DemandF, 'c-', 'LineWidth',1.6);
xline(cfg.r_target,'k--');
title('(c) Demanda por ocupación'); xlabel('r'); ylabel('-a<0');

subplot(2,2,4); hold on; grid on; box on;
plot(rgrid, Stot, 'k-', 'LineWidth',1.8); yline(0,'k:');
xline(cfg.r_target,'r--','r^*'); title('(d) Exceso S_{total}(r)');
xlabel('r'); ylabel('S_{total}');

disp('Listo: BASE calibrado a r*=3% (con fallback automático si fue necesario).');
