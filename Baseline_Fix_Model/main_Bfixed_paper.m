% ============================ main_Bfixed_paper.m ============================
clear; clc; close all;
outdir_tabs='./tables'; if ~exist(outdir_tabs,'dir'), mkdir(outdir_tabs); end
outdir_figs='./figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------- Parámetros base (Perú-like y r objetivo ≈ 3.5%) ----------
cfg=struct();
cfg.RRA_I=3.40; cfg.RRA_F=3.40; cfg.rho=0.08;
cfg.tau_l=0.15; cfg.tau_c=0.18; cfg.phi=0.09;
cfg.z1=0.33; cfg.z2=1.00;

cfg.theta_I=0.08;        % prima mayor informal
cfg.theta_F=0.02;

cfg.eta_target=0.654; cfg.p22_bar=0.8155;

cfg.I=400; cfg.amax=5.0; cfg.amin=-0.30*cfg.z1;

cfg.r_guess=0.03; cfg.rmin=0.005; cfg.rmax=0.10;
cfg.maxit_V=140; cfg.crit_V=1e-6; cfg.Delta=1400;
cfg.maxit_r=1500; cfg.crit_S=1e-5; cfg.fix_r=0;

cfg.sigma_a=0.012;     % suaviza FP
cfg.psi_G=0.30; cfg.omegaG=1.0; % G multiplicativo
cfg.B_mode='ratio_to_Y'; cfg.Bbar=0.35; % punto de partida Perú-like
cfg.alphaG=0.50; cfg.clamp_G_to_zero=true; cfg.G_cap_ratio=0.08;

% ---- Calibrar opcionalmente B/Y para lograr r ~ target_r ----
CALIBRATE_R_TARGET = true; target_r = 0.035;
if CALIBRATE_R_TARGET
    lo=0.20; hi=0.70;  % rango plausible de B/Y
    for it=1:18
        mid=0.5*(lo+hi); alt=cfg; alt.Bbar=mid;
        sol=solve_two_type_huggett_fiscal_Bfixed(alt);
        if sol.r < target_r, lo=mid; else, hi=mid; end
    end
    cfg.Bbar=0.5*(lo+hi);
end

% ------------------------------ Resolver ------------------------------
base=solve_two_type_huggett_fiscal_Bfixed(cfg);

a=base.a; g=base.g; c=base.c; s=base.s; r=base.r;
popI=base.popI; popF=base.popF; Y=base.Y; Ctot=base.Ctot; Gpc=base.Gpc;
fb=base.fiscal; S_res=base.S_residual; hjb_res=base.hjb_residual;

fprintf('\n== BASE (B fixed, G multiplicativo) ==\n');
fprintf('r = %.4f  | S_excess = %.3e\n', r, S_res);
fprintf('popI/popF = %.4f / %.4f (eta=%.4f)\n', popI, popF, popI/(popI+popF));
fprintf('Y = %.6f,   Ctot = %.6f,  Gpc = %.6f  |  B/Y = %.3f\n', Y, Ctot, Gpc, fb.B/max(Y,1e-12));
fprintf('Tl=%.6f  Tc=%.6f  Tr=%.6f  G=%.6f  rB=%.6f  -> PB=%.6f  BB=%.6e\n', ...
        fb.Tl, fb.Tc, fb.Tr, fb.G, fb.rB, fb.PB, fb.BB);
fprintf('||HJB residual||_∞ ≈ %.3e\n', hjb_res);

% ------------------------------ CSVs ---------------------------------
T_fiscal=table("BASE_Bfixed",fb.Tl,fb.Tc,fb.Tl+fb.Tc,fb.Tr,fb.G,fb.rB,fb.PB,fb.B,fb.BB, ...
'VariableNames',{'scenario','labor_tax','vat_tax','rev_total','transfers','public_good','debt_serv','primary_bal','debt_stock','global_bal'});
writetable(T_fiscal,fullfile(outdir_tabs,'base_fiscal_breakdown.csv'));

A_priv = base.Apriv;  % del solver
T_assets=table("BASE_Bfixed",A_priv,fb.B,'VariableNames',{'scenario','A_private','B_public'});
writetable(T_assets,fullfile(outdir_tabs,'base_asset_market.csv'));

S=base.stats;
T_stats=table("BASE_Bfixed",r,popI,popF,Y,Ctot, ...
    S.wealth_mean(1),S.wealth_mean(2),S.wealth_mean(3), ...
    S.wealth_median(1),S.wealth_median(2),S.wealth_median(3), ...
    S.giniW(1),S.giniW(2),S.giniW(3), ...
    S.cons_mean(1),S.cons_mean(2),S.cons_mean(3), ...
    S.cons_median(1),S.cons_median(2),S.cons_median(3), ...
    S.giniC(1),S.giniC(2),S.giniC(3), S.p11, ...
'VariableNames',{'scenario','r','popI','popF','Y','Ctot', ...
'wealth_mean_I','wealth_mean_F','wealth_mean_T', ...
'wealth_med_I','wealth_med_F','wealth_med_T', ...
'giniW_I','giniW_F','giniW_T', ...
'cons_mean_I','cons_mean_F','cons_mean_T', ...
'cons_med_I','cons_med_F','cons_med_T', ...
'giniC_I','giniC_F','giniC_T','p11_rep'});
writetable(T_stats,fullfile(outdir_tabs,'base_household_stats.csv'));

Borr=base.borrowers;
T_borr=table("BASE_Bfixed",Borr.fracBorrow(1),Borr.fracBorrow(2),Borr.fracLend(1),Borr.fracLend(2), ...
             Borr.volBorrow(1),Borr.volBorrow(2),Borr.volLend(1),Borr.volLend(2), ...
'VariableNames',{'scenario','fracBorrow_I','fracBorrow_F','fracLend_I','fracLend_F','volBorrow_I','volBorrow_F','volLend_I','volLend_F'});
writetable(T_borr,fullfile(outdir_tabs,'base_borrowers_lenders.csv'));

fprintf('CSV exportados en %s\n', outdir_tabs);

% ------------------------------ Figuras ------------------------------
paper_style(); da=a(2)-a(1); %#ok<NASGU>

% c(a)
fig=figure('Name','Policy: Consumption'); plot(a,c(:,1),'LineWidth',2); hold on; plot(a,c(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('Consumption c(a)'); legend({'Informal','Formal'},'Location','best'); title('Consumption Policies');
export_fig(fig,fullfile(outdir_figs,'policy_consumption'));

% s(a)
fig=figure('Name','Policy: Savings'); plot(a,s(:,1),'LineWidth',2); hold on; plot(a,s(:,2),'LineWidth',2);
grid on; xlabel('Assets a'); ylabel('Savings s(a)'); legend({'Informal','Formal'},'Location','best'); title('Savings Policies');
export_fig(fig,fullfile(outdir_figs,'policy_savings'));

% Densidades
fig=figure('Name','Wealth distributions');
subplot(2,1,1); bar(a,g(:,1),'FaceAlpha',0.9,'EdgeColor','none'); grid on; xlim([min(a) min(1.2,max(a))]);
xlabel('Assets a'); ylabel('g_1(a)'); title('Informal density');
subplot(2,1,2); bar(a,g(:,2),'FaceAlpha',0.9,'EdgeColor','none'); grid on; xlim([min(a) min(5.0,max(a))]);
xlabel('Assets a'); ylabel('g_2(a)'); title('Formal density');
export_fig(fig,fullfile(outdir_figs,'wealth_distributions'));

% Cierre financiero (barras)
fig=figure('Name','Asset market closure');
bar(categorical({'Private demand','Public debt B'}),[A_priv,fb.B]); ylabel('Assets level');
grid on; title(sprintf('Asset market closes at r = %.4f',r)); export_fig(fig,fullfile(outdir_figs,'asset_market_closure'));

% S(r) y A_priv(r) vs B + descomposición por tipo
r_span=linspace(max(0.005,cfg.r_guess*0.4), cfg.r_guess*2.2, 31);
Sgrid=nan(size(r_span)); AT=nan(size(r_span)); AI=nan(size(r_span)); AF=nan(size(r_span));
alt=cfg; alt.fix_r=1; alt.maxit_r=1;
for i=1:numel(r_span)
    alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
    tmp=solve_two_type_huggett_fiscal_Bfixed(alt);
    Sgrid(i)=tmp.S_residual; AT(i)=tmp.Apriv; AI(i)=tmp.Apriv_I; AF(i)=tmp.Apriv_F;
end
fig=figure('Name','Excess asset supply S(r)'); plot(r_span,Sgrid,'LineWidth',2); yline(0,'k--');
grid on; xlabel('Interest rate r'); ylabel('Excess assets S(r)'); title('Asset market diagnostic (B fixed)'); xline(r,'r:','LineWidth',1.5);
export_fig(fig,fullfile(outdir_figs,'excess_supply_curve'));

fig=figure('Name','A_priv(r) vs B (by type)');
plot(r_span,AT,'LineWidth',2); hold on; plot(r_span,AI,'--','LineWidth',1.6); plot(r_span,AF,':','LineWidth',1.6);
yline(fb.B,'k--','LineWidth',1.5); xline(r,'r:','LineWidth',1.5);
grid on; xlabel('Interest rate r'); ylabel('Levels');
legend({'A_{priv}(r)','A_I(r)','A_F(r)','B (supply)','r^*'},'Location','northwest');
title('Private demand (total y por tipo) vs Public supply'); export_fig(fig,fullfile(outdir_figs,'asset_demand_vs_B_bytype'));

% Fiscal
fig=figure('Name','Fiscal accounts and balances');
subplot(1,3,1); bar(categorical({'Labor tax','VAT'}),[fb.Tl,fb.Tc]); ylabel('Revenues'); title('Revenues'); grid on;
subplot(1,3,2); bar(categorical({'Debt service','Public good','Transfers'}),[fb.rB,fb.G,fb.Tr]); ylabel('Expenditures'); title('Expenditures'); grid on;
subplot(1,3,3); bar(categorical({'Primary balance','Overall (global)'}),[fb.PB,fb.BB]); yline(0,'k--'); ylabel('Balance'); title('Balances'); grid on;
export_fig(fig,fullfile(outdir_figs,'fiscal_accounts_and_balances'));

% Borrowers/Lenders
catX=categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
fig=figure('Name','Borrowers-Lenders shares & volumes');
subplot(1,2,1); bar(catX,[base.borrowers.fracBorrow(:) base.borrowers.fracLend(:)]); grid on;
ylabel('Share'); legend({'Borrowers (a<0)','Lenders (a>0)'},'Location','best'); title('Shares by type');
subplot(1,2,2); bar(catX,[abs(base.borrowers.volBorrow(:)) base.borrowers.volLend(:)]); grid on;
ylabel('Volume'); legend({'|Debt|','Savings'},'Location','best'); title('Volumes by type');
export_fig(fig,fullfile(outdir_figs,'borrowers_lenders'));

% Lorenz total y por tipo (riqueza)
fig=figure('Name','Lorenz wealth (total & by type)');
[LwT,cumPopT]=lorenz_from_density_single(a,g(:,1)+g(:,2));
[LwI,cumPopI]=lorenz_from_density_single(a,g(:,1));
[LwF,cumPopF]=lorenz_from_density_single(a,g(:,2));
plot(cumPopT,LwT,'k','LineWidth',2); hold on;
plot(cumPopI,LwI,'LineWidth',1.8); plot(cumPopF,LwF,'LineWidth',1.8);
plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Wealth share');
legend({'Total','Informal','Formal','45°'},'Location','southeast');
title(sprintf('Lorenz wealth  (Gini_T=%.3f, Gini_I=%.3f, Gini_F=%.3f)', S.giniW(3), S.giniW(1), S.giniW(2)));
export_fig(fig,fullfile(outdir_figs,'lorenz_wealth_by_type'));

% Lorenz por tipo (consumo)
fig=figure('Name','Lorenz consumption (by type)');
[LcI,cPopI]=lorenz_from_values_single(a,g(:,1),c(:,1));
[LcF,cPopF]=lorenz_from_values_single(a,g(:,2),c(:,2));
plot(cPopI,LcI,'LineWidth',1.8); hold on; plot(cPopF,LcF,'LineWidth',1.8);
plot([0,1],[0,1],'k--'); grid on; axis square;
xlabel('Population share'); ylabel('Consumption share');
title(sprintf('Lorenz consumption  (Gini_I=%.3f, Gini_F=%.3f)', S.giniC(1), S.giniC(2)));
legend({'Informal','Formal','45°'},'Location','southeast');
export_fig(fig,fullfile(outdir_figs,'lorenz_consumption_by_type'));

% MPC(a)
eps_z=0.01; alt=cfg; alt.fix_r=1; alt.r_guess=r; alt.rmin=r; alt.rmax=r;
alt.z1=cfg.z1*(1+eps_z); alt.z2=cfg.z2*(1+eps_z);
solP=solve_two_type_huggett_fiscal_Bfixed(alt);
dres1=(cfg.z1*(1+eps_z)-cfg.z1)*(1+cfg.phi); dres2=(cfg.z2*(1+eps_z)-cfg.z2)*(1-cfg.tau_l);
MPC1=(solP.c(:,1)-c(:,1))/max(dres1,1e-12); MPC2=(solP.c(:,2)-c(:,2))/max(dres2,1e-12);
fig=figure('Name','MPC by assets'); plot(a,MPC1,'LineWidth',2); hold on; plot(a,MPC2,'LineWidth',2);
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)'); legend({'Informal','Formal'},'Location','best');
title('MPC along asset grid (finite-diff)'); export_fig(fig,fullfile(outdir_figs,'mpc_by_assets'));

fprintf('Figuras guardadas en %s\n', outdir_figs);

% ---------------- Helpers (gráficos) ----------------
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
end
function [L,cumPop]=lorenz_from_density_single(a,gcol)
    da=a(2)-a(1); W=sum(gcol)*da; W=max(W,1e-12);
    [as,ix]=sort(a); gs=gcol(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as - min(0,min(as)) + 1e-12; cumWealth=cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
function [L,cumPop]=lorenz_from_values_single(a,gcol,x)
    da=a(2)-a(1); w=gcol*da; W=sum(w);
    vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix); cumPop=cumsum(w_s)/W;
    cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end
