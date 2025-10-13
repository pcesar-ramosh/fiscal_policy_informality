% ==================== main_TRANSF_tax_financed_compare.m ====================
clear; clc; close all;
outdir_tabs = './tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end
paper_style();

% ---------- Parámetros base (idénticos a tus runs previos) ----------
cfg = struct();
cfg.RRA_I = 3.40; cfg.RRA_F = 3.40;
cfg.rho   = 0.08;

cfg.tau_l = 0.15;                 % impuesto laboral (sobre formales)
cfg.tau_c = 0.18;                 % IVA
cfg.phi   = 0.09;                 % transferencias a informales (como fracción de z1)
cfg.z1    = 0.33; cfg.z2 = 1.00;

cfg.theta_I = 0.06; cfg.theta_F = 0.01;

cfg.eta_target = 0.654; cfg.p22_bar = 0.8155;

cfg.I = 700; cfg.amax = 3.0;
cfg.amin = -2.0*cfg.z1;

cfg.r_guess = 0.03; cfg.rmin = 0.005; cfg.rmax = 0.10;
cfg.maxit_V = 160; cfg.crit_V = 1e-6; cfg.Delta = 1400;
cfg.maxit_r = 80;  cfg.crit_S = 1e-5; cfg.fix_r = 0;

% Bien público (sólo multiplicador suave, para no distorsionar FOCs)
cfg.psi_G   = 0.08;
cfg.omegaG  = 0.50;
cfg.report_G_effects = 1;

cfg.sigma_a = 0.007;

% Deuda pública: fija y exógena
cfg.B_mode  = 'level';            % usa 'level' (más estable para comparar)
cfg.Bbar    = 0.20;

% --------- Paso 1: resolver base con G residual para medir G/Y ----------
cfg.fiscal_mode = 'G_resid_B_fixed';
base0 = solve_two_type_huggett_fiscal_Bfixed(cfg);
Gshare_base = base0.fiscal.G / max(base0.Y,1e-12);

% --------- Paso 2: “BASE” limpio: fijar G y ajustar IVA para BB=0 ---------
cfg_base = cfg;
cfg_base.fiscal_mode   = 'G_fixed_B_fixed_adjust_tax';
cfg_base.G_target_ratio= Gshare_base;   % mantener mismo G/Y del base
cfg_base.adjust_tax    = 'vat';         % << o 'labor' si quieres ajustar el impuesto laboral
base = solve_two_type_huggett_fiscal_Bfixed(cfg_base);

% --------- Paso 3: Shock de transferencias (financiado con impuestos) -----
cfg_tr  = cfg_base;
cfg_tr.phi = cfg_base.phi * 1.50;   % +50% de transferencias a informales
% (el solver subirá automáticamente τ_c (o τ_l) para que BB=0 con el mismo G)
tran = solve_two_type_huggett_fiscal_Bfixed(cfg_tr);

% --------------------- Reporte por consola ---------------------
print_summary('BASE', base);
print_summary('TRANSF↑', tran);

% --------------------- CSV resumen ------------------------------
T = table( ...
    ["BASE";"TRANSF_UP"], ...
    [base.r; tran.r], ...
    [base.fiscal.Tl; tran.fiscal.Tl], ...
    [base.fiscal.Tc; tran.fiscal.Tc], ...
    [base.fiscal.Tr; tran.fiscal.Tr], ...
    [base.fiscal.G;  tran.fiscal.G], ...
    [base.fiscal.PB; tran.fiscal.PB], ...
    [base.fiscal.BB; tran.fiscal.BB], ...
    [base.popI; tran.popI], [base.popF; tran.popF], ...
    [base.Ctot; tran.Ctot], [base.Y; tran.Y], ...
    'VariableNames', {'scenario','r','Tl','Tc','Tr','G','PB','BB','popI','popF','Ctot','Y'});
writetable(T, fullfile(outdir_tabs,'compare_TRANSF_taxfin.csv'));
fprintf('CSV exportado: %s\n', fullfile(outdir_tabs,'compare_TRANSF_taxfin.csv'));

% --------------------- Gráficos comparativos -------------------
a = base.a;
lab = {'BASE','TRANS↑'};

% c(a)
fig=figure('Name','c(a): BASE vs TRANS↑ (tax adj)');
plot(a,base.c(:,1),'LineWidth',2); hold on; plot(a,base.c(:,2),'LineWidth',2);
plot(a,tran.c(:,1),'--','LineWidth',2);      plot(a,tran.c(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('c(a)');
legend('Inf BASE','For BASE','Inf TR↑','For TR↑','Location','best');
export_fig(fig, fullfile(outdir_figs,'transf_c_of_a'));

% s(a)
fig=figure('Name','s(a): BASE vs TRANS↑ (tax adj)');
plot(a,base.s(:,1),'LineWidth',2); hold on; plot(a,base.s(:,2),'LineWidth',2);
plot(a,tran.s(:,1),'--','LineWidth',2);      plot(a,tran.s(:,2),'--','LineWidth',2);
grid on; xlabel('a'); ylabel('s(a)');
legend('Inf BASE','For BASE','Inf TR↑','For TR↑','Location','best');
export_fig(fig, fullfile(outdir_figs,'transf_s_of_a'));

% Densidades de riqueza por tipo
fig=figure('Name','Wealth densities by type');
subplot(1,2,1);
bar(a,base.g(:,1),'FaceAlpha',0.8,'EdgeColor','none'); hold on;
bar(a,tran.g(:,1),'FaceAlpha',0.4,'EdgeColor','none');
grid on; title('Informal g_1(a)'); xlabel('a'); xlim([min(a) max(a)]);
legend(lab,'Location','best');

subplot(1,2,2);
bar(a,base.g(:,2),'FaceAlpha',0.8,'EdgeColor','none'); hold on;
bar(a,tran.g(:,2),'FaceAlpha',0.4,'EdgeColor','none');
grid on; title('Formal g_2(a)'); xlabel('a'); xlim([min(a) max(a)]);
legend(lab,'Location','best');
export_fig(fig, fullfile(outdir_figs,'transf_densities'));

% Cierre financiero (barras)
da = a(2)-a(1);
A_I_b = sum(base.g(:,1).*a)*da;  A_F_b = sum(base.g(:,2).*a)*da;  A_priv_b = A_I_b + A_F_b;
A_I_t = sum(tran.g(:,1).*a)*da;  A_F_t = sum(tran.g(:,2).*a)*da;  A_priv_t = A_I_t + A_F_t;

fig=figure('Name','Asset market closure (aggregate)');
X = categorical({'A_priv','B (public)'}); X = reordercats(X,{'A_priv','B (public)'});
bar(X,[A_priv_b base.fiscal.B; A_priv_t tran.fiscal.B]); grid on;
ylabel('level'); legend(lab,'Location','best');
title(sprintf('Closes at r*: BASE %.4f | TR↑ %.4f', base.r, tran.r));
export_fig(fig, fullfile(outdir_figs,'transf_asset_closure'));

% Curvas S(r) por escenario
[r_span_b,S_b,~,~,rstar_b] = financial_curves(cfg_base, base.r);
[r_span_t,S_t,~,~,rstar_t] = financial_curves(cfg_tr,   tran.r);
fig=figure('Name','Excess asset supply S(r) comparison');
plot(r_span_b,S_b,'LineWidth',2); hold on; plot(r_span_t,S_t,'--','LineWidth',2);
yline(0,'k:'); xline(rstar_b,':','Color',[0 .4 .8]); xline(rstar_t,':','Color',[.85 .33 .10]);
grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B'); legend('BASE','TRANS↑','S=0','r* BASE','r* TR↑','Location','northwest');
export_fig(fig, fullfile(outdir_figs,'transf_S_of_r'));

% Fiscal (ingresos/egresos/balances)
fb = base.fiscal; ft = tran.fiscal;
fig=figure('Name','Fiscal accounts: BASE vs TRANS↑');
subplot(1,3,1); 
cats=categorical({'Labor','VAT','Total'}); cats=reordercats(cats,{'Labor','VAT','Total'});
bar(cats,[fb.Tl fb.Tc fb.Tl+fb.Tc; ft.Tl ft.Tc ft.Tl+ft.Tc]); title('Revenues'); ylabel('amount'); grid on; legend(lab);
subplot(1,3,2);
cats=categorical({'Debt serv','G','Transfers'}); cats=reordercats(cats,{'Debt serv','G','Transfers'});
bar(cats,[fb.rB fb.G fb.Tr; ft.rB ft.G ft.Tr]); title('Expenditures'); grid on;
subplot(1,3,3);
cats=categorical({'Primary','Overall'}); cats=reordercats(cats,{'Primary','Overall'});
bar(cats,[fb.PB fb.BB; ft.PB ft.BB]); yline(0,'k--'); title('Balances'); grid on;
export_fig(fig, fullfile(outdir_figs,'transf_fiscal'));

% Borrowers/Lenders
fig=figure('Name','Borrowers & Lenders');
cat2 = categorical({'Informal','Formal'}); cat2=reordercats(cat2,{'Informal','Formal'});
subplot(1,2,1);
bar(cat2,[base.borrowers.fracBorrow(:) tran.borrowers.fracBorrow(:)]); title('Borrowers share'); legend(lab); grid on;
subplot(1,2,2);
bar(cat2,[base.borrowers.fracLend(:)   tran.borrowers.fracLend(:)]);   title('Lenders share');   legend(lab); grid on;
export_fig(fig, fullfile(outdir_figs,'transf_borrowers_lenders_shares'));

fig=figure('Name','Borrowers & Lenders volumes');
subplot(1,2,1);
bar(cat2,[abs(base.borrowers.volBorrow(:)) abs(tran.borrowers.volBorrow(:))]); title('|Debt|'); legend(lab); grid on;
subplot(1,2,2);
bar(cat2,[base.borrowers.volLend(:) tran.borrowers.volLend(:)]);       title('Savings'); legend(lab); grid on;
export_fig(fig, fullfile(outdir_figs,'transf_borrowers_lenders_vols'));

% Lorenz riqueza (por tipo)
fig=figure('Name','Lorenz wealth by type: BASE vs TRANS↑');
[LwI_b,cI_b] = lorenz_from_assets_single(a, base.g(:,1));
[LwI_t,cI_t] = lorenz_from_assets_single(a, tran.g(:,1));
[LwF_b,cF_b] = lorenz_from_assets_single(a, base.g(:,2));
[LwF_t,cF_t] = lorenz_from_assets_single(a, tran.g(:,2));
subplot(1,2,1);
plot(cI_b,LwI_b,'LineWidth',2); hold on; plot(cI_t,LwI_t,'--','LineWidth',2); plot([0 1],[0 1],'k:');
axis square; grid on; title(sprintf('Informal: Gini %.3f | %.3f', base.stats.giniW(1), tran.stats.giniW(1)));
subplot(1,2,2);
plot(cF_b,LwF_b,'LineWidth',2); hold on; plot(cF_t,LwF_t,'--','LineWidth',2); plot([0 1],[0 1],'k:');
axis square; grid on; title(sprintf('Formal: Gini %.3f | %.3f', base.stats.giniW(2), tran.stats.giniW(2)));
export_fig(fig, fullfile(outdir_figs,'transf_lorenz_bytype'));

% Lorenz consumo (total)
fig=figure('Name','Lorenz consumption total');
[~,Lc_b,cPop_b]=lorenz_from_values(a, base.g, base.c(:,1)+base.c(:,2));
[~,Lc_t,cPop_t]=lorenz_from_values(a, tran.g, tran.c(:,1)+tran.c(:,2));
plot(cPop_b,Lc_b,'LineWidth',2); hold on; plot(cPop_t,Lc_t,'--','LineWidth',2); plot([0 1],[0 1],'k:');
axis square; grid on; xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Gini_C total: %.3f | %.3f', base.stats.giniC(3), tran.stats.giniC(3)));
legend(lab,'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'transf_lorenz_consumption_total'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ================== Helpers locales ==================
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r250'); print(fig,[basepath '.pdf'],'-dpdf');
end
function print_summary(tag, S)
    fb=S.fiscal;
    fprintf('\n== %s ==\n', tag);
    fprintf('r*=%.4f | Y=%.4f | Ctot=%.4f | pop(I,F)=(%.3f, %.3f)\n', S.r, S.Y, S.Ctot, S.popI, S.popF);
    fprintf('Taxes: (Tl=%.4f, Tc=%.4f)  Transfers Tr=%.4f  G=%.4f  B=%.4f  rB=%.4f\n', fb.Tl, fb.Tc, fb.Tr, fb.G, fb.B, fb.rB);
    fprintf('Balances: PB=%.4e  BB=%.4e\n', fb.PB, fb.BB);
    fprintf('Ginis W (I,F,T)=(%.3f, %.3f, %.3f) | Ginis C (I,F,T)=(%.3f, %.3f, %.3f)\n', ...
        S.stats.giniW, S.stats.giniC);
end
function [r_span,Sgrid,Apriv,Bg,r_star] = financial_curves(cfg_in, r_star_guess)
    r_lo = max(cfg_in.rmin, max(0.6*r_star_guess, 0.007));
    r_hi = min(cfg_in.rmax, max(1.6*r_star_guess, r_star_guess+0.02));
    r_span = linspace(r_lo, r_hi, 35);
    cfg = cfg_in; cfg.fix_r=1;
    Sgrid=nan(size(r_span)); Apriv=Sgrid; Bg=Sgrid;
    for k=1:numel(r_span)
        cfg.r_guess=r_span(k); cfg.rmin=cfg.r_guess; cfg.rmax=cfg.r_guess;
        sol= solve_two_type_huggett_fiscal_Bfixed(cfg);
        da=sol.a(2)-sol.a(1);
        Apriv(k)=sum(sol.g(:,1).*sol.a)*da + sum(sol.g(:,2).*sol.a)*da;
        Bg(k)    =sol.fiscal.B; Sgrid(k)=Apriv(k)-Bg(k);
    end
    idx = find(diff(sign(Sgrid))~=0,1,'first');
    if ~isempty(idx)
        r1=r_span(idx); r2=r_span(idx+1); s1=Sgrid(idx); s2=Sgrid(idx+1);
        r_star = r1 - s1*(r2-r1)/max(s2-s1,1e-12);
    else
        r_star = r_star_guess;
    end
end
function [gT,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w);
    vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix);
    cumPop=cumsum(w_s)/W; cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end
function [L,cumPop]=lorenz_from_assets_single(a, gk)
    da=a(2)-a(1); W=sum(gk)*da; [as,ix]=sort(a); gs=gk(ix)*da; cumPop=cumsum(gs)/max(W,1e-12);
    wealth_pos = as - min(0,min(as)) + 1e-12;  cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
