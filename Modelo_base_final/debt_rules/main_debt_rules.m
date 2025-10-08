% ======================= main_debt_rules_target.m =======================
clear; clc; close all;
outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------- Parámetros base ----------
cfg = struct();
cfg.RRA_I = 3.40; cfg.RRA_F = 3.40;
cfg.rho   = 0.08;
cfg.tau_l = 0.15; cfg.tau_c = 0.18; cfg.phi = 0.09;
cfg.z1 = 0.33;  cfg.z2 = 1.00;
cfg.theta_I = 0.06; cfg.theta_F = 0.01;

cfg.eta_target = 0.654; cfg.p22_bar = 0.8155;
cfg.I = 600;  cfg.amax = 3.0; cfg.amin = -2.0*cfg.z1;

cfg.r_guess = 0.03; cfg.rmin = 0.000; cfg.rmax = 0.12;  % deja rmin=0 si no quieres r<0
cfg.maxit_V = 120; cfg.crit_V = 1e-6; cfg.Delta = 1400;
cfg.maxit_r = 1500; cfg.crit_S = 1e-5; cfg.fix_r = 0;

% Bien público
cfg.psi_G = 0.08; cfg.omegaG = 0.50; cfg.sigma_a = 0.010;
cfg.alphaG = 0.50; cfg.clamp_G_to_zero = true; cfg.G_cap_ratio = 0.08;
cfg.report_G_effects = 1;

% ---------- Escenarios de deuda ----------
sc = struct([]);

% (1) RATIO_to_Y: B = bbar * Y
sc(1).name = "RATIO_to_Y";
sc(1).cfg  = cfg; sc(1).cfg.DebtRule = 'ratio_to_Y'; sc(1).cfg.Bbar = 0.35;

% (2) LEVEL_fixed: B = Bbar (nivel)
sc(2).name = "LEVEL_fixed";
sc(2).cfg  = cfg; sc(2).cfg.DebtRule = 'level';      sc(2).cfg.Bbar = 0.20;

% (3) FREE_nonneg: opción (c) asesor (estacionario: B = 0)
sc(3).name = "FREE_nonneg_Beq0";
sc(3).cfg  = cfg; sc(3).cfg.DebtRule = 'free_nonneg';
% Nota: si quieres permitir r<0 para ver si A_priv(r)=0 cae a la izquierda:
% sc(3).cfg.rmin = -0.02;

% ---------- (OPCIONAL) Subir r* a 0.0235 moviendo deuda ----------
tune_r = true; r_target = 0.0235;
for k=1:numel(sc)
    if tune_r && ~strcmpi(sc(k).cfg.DebtRule,'free_nonneg')
        [sc(k).cfg, sc(k).sol] = calibrate_B_to_r(sc(k).cfg, r_target, 30, 2.5e-4);
    else
        sc(k).sol = solve_two_type_huggett_fiscal_B_debtrules(sc(k).cfg);
    end
    fprintf('[%s] r*=%.4f | B=%.3f | Y=%.3f | C=%.3f | G=%.3f | PB=%.3f | BB=%.3f | no_root=%d\n',...
        sc(k).name, sc(k).sol.r, sc(k).sol.fiscal.B, sc(k).sol.Y, sc(k).sol.Ctot,...
        sc(k).sol.fiscal.G, sc(k).sol.fiscal.PB, sc(k).sol.fiscal.BB, sc(k).sol.no_root);
end

% ---------- Tabla resumen ----------
tmp = [sc.sol];
T = table( string({sc.name})', [tmp.r]', [tmp.Y]', [tmp.Ctot]', ...
           arrayfun(@(x)x.sol.fiscal.B, sc)', arrayfun(@(x)x.sol.fiscal.G, sc)', ...
           arrayfun(@(x)x.sol.fiscal.PB, sc)', arrayfun(@(x)x.sol.fiscal.BB, sc)', ...
           'VariableNames',{'scenario','r','Y','C','B','G','PB','BB'});
disp(T); writetable(T, fullfile(outdir_tabs,'compare_debt_rules_target.csv'));

% ---------- Curvas A_priv(r) vs B(r) centradas en r* ----------
paper_style();
for k=1:numel(sc)
    base = sc(k).sol; cfgk = sc(k).cfg;
    r_star = base.r;
    left   = max(cfgk.rmin, max(0, r_star - 0.015));
    right  = min(cfgk.rmax, r_star + 0.035);
    if left>=right, left=max(cfgk.rmin,r_star*0.5); right=min(cfgk.rmax,r_star*1.5); end
    r_span = linspace(left, right, 51);

    Apriv = nan(size(r_span)); Bspan = Apriv;
    alt = cfgk; alt.fix_r=1; alt.maxit_r=1; alt.maxit_V = min(80, alt.maxit_V);
    for i=1:numel(r_span)
        alt.r_guess = r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmpSol = solve_two_type_huggett_fiscal_B_debtrules(alt);
        da = tmpSol.a(2)-tmpSol.a(1);
        Apriv(i) = sum((tmpSol.g(:,1)+tmpSol.g(:,2)).*tmpSol.a)*da;
        Bspan(i) = tmpSol.fiscal.B;
    end
    % punto en el cruce
    da = base.a(2)-base.a(1);
    Astar = sum((base.g(:,1)+base.g(:,2)).*base.a)*da;
    Bstar = base.fiscal.B;

    fig=figure('Name',['Apriv_vs_B_' char(sc(k).name)]);
    plot(r_span,Apriv,'LineWidth',2); hold on; plot(r_span,Bspan,'k--','LineWidth',2);
    xline(r_star,'r:','LineWidth',1.5); plot(r_star, Astar,'ro','MarkerSize',6,'MarkerFaceColor','r');
    grid on; xlabel('Interest rate r'); ylabel('Assets / Bonds');
    legend({'A_{priv}(r)','B(r)','r^*','Equilibrium point'},'Location','northwest');
    title(['Private demand vs Public supply — ' char(sc(k).name)]);
    print(fig, fullfile(outdir_figs, ['Apriv_vs_B_' char(sc(k).name) '.png']), '-dpng','-r220');
end

fprintf('Outputs: %s (tables) | %s (figures)\n', outdir_tabs, outdir_figs);

% -------- helpers --------
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultFigureColor','w');
end

% ---- calibrador: mueve Bbar (ratio o nivel) para alcanzar r_target ----
function [cfg_out, sol_best] = calibrate_B_to_r(cfg_in, r_target, maxit, tol)
    if nargin<3, maxit=25; end, if nargin<4, tol=2.5e-4; end
    cfg_out = cfg_in;

    switch lower(cfg_in.DebtRule)
        case 'ratio_to_y'
            BLo = 0.08; BHi = 1.80;
        case 'level'
            BLo = 0.02; BHi = 1.00;   % niveles razonables para tu escala de Y
        otherwise
            error('calibrate_B_to_r: DebtRule must be ratio_to_Y or level.');
    end

    [rLo,solLo] = r_given_Bparam(cfg_out,BLo);
    [rHi,solHi] = r_given_Bparam(cfg_out,BHi);
    it=0;
    while ~((rLo<=r_target && r_target<=rHi) || (rHi<=r_target && r_target<=rLo))
        it=it+1; if it>8, break; end
        % expande el bracket por arriba
        BHi = BHi*1.5; [rHi,solHi] = r_given_Bparam(cfg_out,BHi);
        % y por abajo si hace falta
        if r_target < min(rLo,rHi), BLo = max(0.01,BLo*0.6); [rLo,solLo] = r_given_Bparam(cfg_out,BLo); end
    end

    for it=1:maxit
        Bmid = 0.5*(BLo+BHi);
        [rmid,solmid] = r_given_Bparam(cfg_out,Bmid);
        if abs(rmid - r_target) < tol
            cfg_out.Bbar = Bmid; sol_best = solmid; return;
        end
        if (rLo <= rHi && rmid < r_target) || (rHi < rLo && rmid > r_target)
            BLo = Bmid; rLo = rmid; solLo = solmid;
        else
            BHi = Bmid; rHi = rmid; solHi = solmid;
        end
    end
    % salida por el lado más cercano
    if abs(rLo - r_target) < abs(rHi - r_target)
        cfg_out.Bbar = BLo; sol_best = solLo;
    else
        cfg_out.Bbar = BHi; sol_best = solHi;
    end
end

function [r_star, sol] = r_given_Bparam(cfg_in, Bparam)
    cfg_tmp = cfg_in; cfg_tmp.Bbar = Bparam;
    sol = solve_two_type_huggett_fiscal_B_debtrules(cfg_tmp);
    r_star = sol.r;
end
