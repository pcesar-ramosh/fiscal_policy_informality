% ====================== main_all_scenarios_tables.m =======================
clear; clc; close all;

outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end

% ================== Parámetros base del modelo (limpio) ===================
cfg = struct();

% Preferencias e impuestos (base)
cfg.RRA_I = 2.40;
cfg.RRA_F = 3.40;
cfg.rho   = 0.08;
cfg.tau_l = 0.15;
cfg.tau_c = 0.18;
cfg.phi   = 0.09;

% Ingresos por tipo (base)
cfg.z1 = 0.33;   % informal
cfg.z2 = 1.00;   % formal

% Spreads (costo de endeudamiento)
cfg.theta_I = 0.06;
cfg.theta_F = 0.01;

% Markov de tipos (base)
cfg.eta_target = 0.654;
cfg.p22_bar    = 0.8155;

% Grid de activos
cfg.I    = 700;
cfg.amin = -2.0*cfg.z1;
cfg.amax = 3.0;

% Interés (bisección)
cfg.r_guess = 0.03;
cfg.rmin    = 0.005;
cfg.rmax    = 0.10;
cfg.maxit_r = 100;
cfg.crit_S  = 1e-5;

% HJB / numérico
cfg.maxit_V = 160;
cfg.crit_V  = 1e-6;
cfg.Delta   = 1400;
cfg.sigma_a = 0.007;

% Bien público en utilidad (no esencial para tablas, pero consistente)
cfg.psi_G  = 0.08;
cfg.omegaG = 0.50;
cfg.report_G_effects = 0;

% Gobierno y deuda
cfg.B_mode = 'level';   % mantener B en nivel (base)
cfg.Bbar   = 0.25;      % provisional (afinamos abajo con pre-pase)

% =================== Pre-pase para fijar Bbar razonable ===================
% Igual que en tus scripts limpios: elegimos Bbar en medio de A_priv(rmin) y
% A_priv(rmax) con Bbar=0 temporal en el pre-pase para no contaminar.
cfg_lo = cfg; cfg_lo.fix_r=1; cfg_lo.r_guess=cfg.rmin; cfg_lo.rmin=cfg.rmin; cfg_lo.rmax=cfg.rmin; cfg_lo.Bbar=0;
cfg_hi = cfg; cfg_hi.fix_r=1; cfg_hi.r_guess=cfg.rmax; cfg_hi.rmin=cfg.rmax; cfg_hi.rmax=cfg.rmax; cfg_hi.Bbar=0;

sol_lo = solve_two_type_huggett_fiscal_Bfixed(cfg_lo);
sol_hi = solve_two_type_huggett_fiscal_Bfixed(cfg_hi);

da = sol_lo.a(2)-sol_lo.a(1);
Apriv_lo = sum(sol_lo.g(:,1).*sol_lo.a)*da + sum(sol_lo.g(:,2).*sol_lo.a)*da;
Apriv_hi = sum(sol_hi.g(:,1).*sol_hi.a)*da + sum(sol_hi.g(:,2).*sol_hi.a)*da;

B_pre = 0.5*(Apriv_lo + Apriv_hi);
B_pre = max(0.05, min(1.50, B_pre));
cfg.Bbar = B_pre;

fprintf('Pre-pase Bbar: Apriv(rmin)=%.4f, Apriv(rmax)=%.4f => Bbar=%.4f\n', Apriv_lo, Apriv_hi, cfg.Bbar);

% ============================ Definir escenarios ==========================
% Ajustes (puedes cambiarlos fácilmente)
d_tau_c     = +0.05;     % i) aumento IVA
d_tau_l     = +0.05;     % ii) aumento impuesto laboral
phi_up_pct  = +0.15;     % iii) +15% transferencias (factor 1.15)
z1_drop     = 0.80;      % iv-v) caída 20% informal
z2_drop     = 0.90;      % iv-v) caída 10% formal
eta_high    = 0.75;      % vii) mayor informalidad
eta_low     = 0.45;      % viii) menor informalidad
d_theta_I   = +0.04;     % x) costo endeudamiento (informal)
d_theta_F   = +0.02;     % x) costo endeudamiento (formal)

SC = {}; k = 0;

% 0) Base
k=k+1; SC{k}.label='BASE';            SC{k}.f=@(c) c;

% i) IVA up
k=k+1; SC{k}.label='VAT_UP';          SC{k}.f=@(c) setf(c,'tau_c',c.tau_c + d_tau_c);

% ii) Impuesto laboral up
k=k+1; SC{k}.label='LABTAX_UP';       SC{k}.f=@(c) setf(c,'tau_l',c.tau_l + d_tau_l);

% iii) Transferencias up (+15%)
k=k+1; SC{k}.label='TRANSF_UP';       SC{k}.f=@(c) setf(c,'phi',c.phi*(1+phi_up_pct));

% iv) Caída ingresos sin transferencias
k=k+1; SC{k}.label='INCDROP_NO_TR';   SC{k}.f=@(c) setf(setf(c,'z1',c.z1*z1_drop),'z2',c.z2*z2_drop);

% v) Caída ingresos + transferencias up
k=k+1; SC{k}.label='INCDROP_TR_UP';   SC{k}.f=@(c) setf(setf(setf(c,'z1',c.z1*z1_drop),'z2',c.z2*z2_drop),'phi',c.phi*(1+phi_up_pct));

% vi) Informalidad base (explícito)
k=k+1; SC{k}.label='INF_BASE';        SC{k}.f=@(c) setf(c,'eta_target',cfg.eta_target);

% vii) Informalidad alta
k=k+1; SC{k}.label='INF_HIGH';        SC{k}.f=@(c) setf(c,'eta_target',eta_high);

% viii) Informalidad baja
k=k+1; SC{k}.label='INF_LOW';         SC{k}.f=@(c) setf(c,'eta_target',eta_low);

% ix) Sin política fiscal: sin impuestos, sin transferencias y sin deuda
% (G=0 endógeno con nuestra regla exacta si B=0)
k=k+1; SC{k}.label='NO_FISC';         SC{k}.f=@(c) set_no_fisc(c);

% x) Costo endeudamiento más alto
k=k+1; SC{k}.label='SPREAD_UP';       SC{k}.f=@(c) setf(setf(c,'theta_I',c.theta_I + d_theta_I),'theta_F',c.theta_F + d_theta_F);

% ====================== Correr escenarios y armar tabla ===================
Rows = [];
for s=1:numel(SC)
    cfg_s = SC{s}.f(cfg);

    % Mantener Bbar de referencia en todos, excepto NO_FISC:
    if ~strcmp(SC{s}.label,'NO_FISC')
        cfg_s.B_mode = cfg.B_mode;
        cfg_s.Bbar   = cfg.Bbar;
    end

    sol = solve_two_type_huggett_fiscal_Bfixed(cfg_s);
    Row = collect_stats(sol, SC{s}.label);
    Rows = [Rows; Row]; %#ok<AGROW>
    fprintf('Escenario %-14s | r*=%.5f | Y=%.4f | C=%.4f | B=%.4f | Tl=%.4f Tc=%.4f Tr=%.4f G=%.4f\n', ...
        SC{s}.label, sol.r, sol.Y, sol.Ctot, sol.fiscal.B, sol.fiscal.Tl, sol.fiscal.Tc, sol.fiscal.Tr, sol.fiscal.G);
end

% Exportar
writetable(Rows, fullfile(outdir_tabs,'all_scenarios_stats.csv'));
disp('CSV exportado: ./tables/all_scenarios_stats.csv');

% ========================= Helpers (solo del main) ========================
function c2 = setf(c, field, val)
    c2 = c; c2.(field) = val;
end

function c2 = set_no_fisc(c)
    c2 = c;
    c2.tau_c = 0.0; c2.tau_l = 0.0; c2.phi = 0.0;
    c2.psi_G = 0.0;            % utilidad no afecta (opcional)
    c2.B_mode = 'level';
    c2.Bbar   = 0.0;           % sin deuda => G=0 con cierre exacto
end

function T = collect_stats(sol, label)
    % Empaqueta todas las estadísticas en una sola fila (tabla).
    a = sol.a; g = sol.g; c = sol.c; s = sol.s; da = a(2)-a(1);
    wI = g(:,1)*da; wF = g(:,2)*da; wT = (g(:,1)+g(:,2))*da;

    % Atajos (fiscal y agregados)
    fb = sol.fiscal;
    popI = sum(wI); popF = sum(wF);

    % ---------- Estadísticos ponderados ----------
    % Wealth (a)
    [WmI,WmedI,WsdI,WvarI,Wp90I,Wp10I,WratI,GW_I] = wstats_all(a,wI,true);
    [WmF,WmedF,WsdF,WvarF,Wp90F,Wp10F,WratF,GW_F] = wstats_all(a,wF,true);
    [WmT,WmedT,WsdT,WvarT,Wp90T,Wp10T,WRatT,GW_T] = wstats_all(a,wT,true);

    % Consumption c(:,i)
    [CmI,CmedI,CsdI,CvarI,Cp90I,Cp10I,CRatI,GC_I] = wstats_all(c(:,1),wI,true);
    [CmF,CmedF,CsdF,CvarF,Cp90F,Cp10F,CRatF,GC_F] = wstats_all(c(:,2),wF,true);
    [CmT,CmedT,CsdT,CvarT,Cp90T,Cp10T,CRatT,GC_T] = wstats_all(c(:,1)+c(:,2),wT,true);

    % Savings policy s(:,i)
    [SmI,SmedI,SsdI,SvarI,Sp90I,Sp10I,SRatI,GS_I] = wstats_all(s(:,1),wI,true);
    [SmF,SmedF,SsdF,SvarF,Sp90F,Sp10F,SRatF,GS_F] = wstats_all(s(:,2),wF,true);
    [SmT,SmedT,SsdT,SvarT,Sp90T,Sp10T,SRatT,GS_T] = wstats_all(s(:,1)+s(:,2),wT,true);

    % Endeudamiento / Prestamistas (shares y stats condicionales)
    % Shares
    fracB_I = sum(wI(a<0))/max(popI,1e-12);
    fracL_I = sum(wI(a>0))/max(popI,1e-12);
    fracB_F = sum(wF(a<0))/max(popF,1e-12);
    fracL_F = sum(wF(a>0))/max(popF,1e-12);
    fracB_T = sum(wT(a<0))/max(sum(wT),1e-12);
    fracL_T = sum(wT(a>0))/max(sum(wT),1e-12);

    % Deuda condicional (entre deudores): x = -a, a<0
    [DmeanI,DmedI,DsdI,DvarI,Dp90I,Dp10I,DRatI] = cond_debt_stats(a,wI);
    [DmeanF,DmedF,DsdF,DvarF,Dp90F,Dp10F,DRatF] = cond_debt_stats(a,wF);
    [DmeanT,DmedT,DsdT,DvarT,Dp90T,Dp10T,DRatT] = cond_debt_stats(a,wT);

    % Prestamistas condicional (entre a>0): activos positivos
    [LmeanI,LmedI,LsdI,LvarI,Lp90I,Lp10I,LRatI] = cond_pos_stats(a,wI);
    [LmeanF,LmedF,LsdF,LvarF,Lp90F,Lp10F,LRatF] = cond_pos_stats(a,wF);
    [LmeanT,LmedT,LsdT,LvarT,Lp90T,Lp10T,LRatT] = cond_pos_stats(a,wT);

    % ---------- Armar tabla de una fila ----------
    T = table(string(label), sol.r, fb.B, popI, popF, sol.Y, sol.Ctot, ...
        fb.Tl, fb.Tc, fb.Tl+fb.Tc, fb.Tr, fb.G, fb.rB, (fb.G+fb.Tr+fb.rB), fb.PB, fb.BB, ...
        ... % Wealth
        WmI,WmedI,WsdI,WvarI,GW_I,Wp90I,Wp10I,WratI, ...
        WmF,WmedF,WsdF,WvarF,GW_F,Wp90F,Wp10F,WratF, ...
        WmT,WmedT,WsdT,WvarT,GW_T,Wp90T,Wp10T,WRatT, ...
        ... % Consumption
        CmI,CmedI,CsdI,CvarI,GC_I,Cp90I,Cp10I,CRatI, ...
        CmF,CmedF,CsdF,CvarF,GC_F,Cp90F,Cp10F,CRatF, ...
        CmT,CmedT,CsdT,CvarT,GC_T,Cp90T,Cp10T,CRatT, ...
        ... % Savings
        SmI,SmedI,SsdI,SvarI,GS_I,Sp90I,Sp10I,SRatI, ...
        SmF,SmedF,SsdF,SvarF,GS_F,Sp90F,Sp10F,SRatF, ...
        SmT,SmedT,SsdT,SvarT,GS_T,Sp90T,Sp10T,SRatT, ...
        ... % Debt / Lenders: shares + cond stats
        fracB_I,fracL_I, DmeanI,DmedI,DsdI,DvarI,Dp90I,Dp10I,DRatI,  LmeanI,LmedI,LsdI,LvarI,Lp90I,Lp10I,LRatI, ...
        fracB_F,fracL_F, DmeanF,DmedF,DsdF,DvarF,Dp90F,Dp10F,DRatF,  LmeanF,LmedF,LsdF,LvarF,Lp90F,Lp10F,LRatF, ...
        fracB_T,fracL_T, DmeanT,DmedT,DsdT,DvarT,Dp90T,Dp10T,DRatT,  LmeanT,LmedT,LsdT,LvarT,Lp90T,Lp10T,LRatT, ...
        'VariableNames', { ...
        'scenario','r','B','popI','popF','Y','Ctot', ...
        'Tl','Tc','Rev_total','Tr','G','rB','Exp_total','PB','BB', ...
        ... % Wealth names
        'wealth_mean_I','wealth_med_I','wealth_sd_I','wealth_var_I','giniW_I','wealth_p90_I','wealth_p10_I','wealth_p90p10_I', ...
        'wealth_mean_F','wealth_med_F','wealth_sd_F','wealth_var_F','giniW_F','wealth_p90_F','wealth_p10_F','wealth_p90p10_F', ...
        'wealth_mean_T','wealth_med_T','wealth_sd_T','wealth_var_T','giniW_T','wealth_p90_T','wealth_p10_T','wealth_p90p10_T', ...
        ... % Consumption names
        'cons_mean_I','cons_med_I','cons_sd_I','cons_var_I','giniC_I','cons_p90_I','cons_p10_I','cons_p90p10_I', ...
        'cons_mean_F','cons_med_F','cons_sd_F','cons_var_F','giniC_F','cons_p90_F','cons_p10_F','cons_p90p10_F', ...
        'cons_mean_T','cons_med_T','cons_sd_T','cons_var_T','giniC_T','cons_p90_T','cons_p10_T','cons_p90p10_T', ...
        ... % Savings names
        'sav_mean_I','sav_med_I','sav_sd_I','sav_var_I','giniS_I','sav_p90_I','sav_p10_I','sav_p90p10_I', ...
        'sav_mean_F','sav_med_F','sav_sd_F','sav_var_F','giniS_F','sav_p90_F','sav_p10_F','sav_p90p10_F', ...
        'sav_mean_T','sav_med_T','sav_sd_T','sav_var_T','giniS_T','sav_p90_T','sav_p10_T','sav_p90p10_T', ...
        ... % Debt/Lender names
        'share_borrow_I','share_lend_I','debt_mean_I','debt_med_I','debt_sd_I','debt_var_I','debt_p90_I','debt_p10_I','debt_p90p10_I', ...
        'lendpos_mean_I','lendpos_med_I','lendpos_sd_I','lendpos_var_I','lendpos_p90_I','lendpos_p10_I','lendpos_p90p10_I', ...
        'share_borrow_F','share_lend_F','debt_mean_F','debt_med_F','debt_sd_F','debt_var_F','debt_p90_F','debt_p10_F','debt_p90p10_F', ...
        'lendpos_mean_F','lendpos_med_F','lendpos_sd_F','lendpos_var_F','lendpos_p90_F','lendpos_p10_F','lendpos_p90p10_F', ...
        'share_borrow_T','share_lend_T','debt_mean_T','debt_med_T','debt_sd_T','debt_var_T','debt_p90_T','debt_p10_T','debt_p90p10_T', ...
        'lendpos_mean_T','lendpos_med_T','lendpos_sd_T','lendpos_var_T','lendpos_p90_T','lendpos_p10_T','lendpos_p90p10_T' ...
        });
end

function [m,med,sd,varx,p90,p10,rat,gini] = wstats_all(x,w,do_gini)
    W = sum(w); if W<=0, m=NaN; med=NaN; sd=NaN; varx=NaN; p90=NaN; p10=NaN; rat=NaN; gini=NaN; return; end
    m   = sum(x(:).*w(:))/W;
    varx= sum(((x(:)-m).^2).*w(:))/W;
    sd  = sqrt(max(varx,0));
    med = wpercentile(x,w,0.5);
    p90 = wpercentile(x,w,0.9);
    p10 = wpercentile(x,w,0.1);
    rat = p90 / max(p10,1e-12);
    if do_gini
        % Gini con desplazamiento para evitar negativos
        xmin = min(x(:));
        xs   = x(:) - min(0,xmin) + 1e-12;
        gini = wgini(xs,w);
    else
        gini = NaN;
    end
end

function q = wpercentile(x,w,p)
    x = x(:); w = w(:);
    [xs,ix] = sort(x); ws = w(ix);
    cw = cumsum(ws); W = sum(ws);
    if W<=0, q = NaN; return; end
    idx = find(cw>=p*W,1,'first'); 
    if isempty(idx), q = xs(end); else, q = xs(idx); end
end

function g = wgini(x,w)
    x = x(:); w = w(:); W = sum(w);
    if W<=0, g=NaN; return; end
    [xs,ix] = sort(x); ws = w(ix);
    cw = cumsum(ws)/W;
    cx = cumsum(xs.*ws);
    if cx(end)<=0, g=NaN; return; end
    L = cx / cx(end);           % curva de Lorenz
    g = 1 - 2*trapz(cw, L);
end

function [m,med,sd,varx,p90,p10,rat] = cond_debt_stats(a,w)
    % Estadísticos sobre la "deuda" condicional a a<0, definiendo x = -a (positiva)
    mask = (a<0);
    if ~any(mask)
        [m,med,sd,varx,p90,p10,rat] = deal(NaN);
        return;
    end
    x = -a(mask); ww = w(mask);
    [m,med,sd,varx,p90,p10,rat,~] = wstats_all(x,ww,false);
end

function [m,med,sd,varx,p90,p10,rat] = cond_pos_stats(a,w)
    % Estadísticos sobre activos positivos condicionales a a>0
    mask = (a>0);
    if ~any(mask)
        [m,med,sd,varx,p90,p10,rat] = deal(NaN);
        return;
    end
    x = a(mask); ww = w(mask);
    [m,med,sd,varx,p90,p10,rat,~] = wstats_all(x,ww,false);
end
