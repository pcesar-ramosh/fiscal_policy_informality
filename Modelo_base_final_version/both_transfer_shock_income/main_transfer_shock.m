% ======================== main_transfer_shock.m ==========================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------------------- BASE configuration -------------------------------
cfg0 = struct();
% Preferencias e impuestos
cfg0.RRA_I = 3.40; cfg0.RRA_F = 3.40; cfg0.rho = 0.08;
cfg0.tau_l = 0.15; cfg0.tau_c = 0.18; cfg0.phi = 0.09;
% Ingresos base
cfg0.z1 = 0.33; cfg0.z2 = 1.00;
% Primas
cfg0.theta_I = 0.06; cfg0.theta_F = 0.01;
% Ocupacional
cfg0.eta_target = 0.654; cfg0.p22_bar = 0.8155;
% Grid y numérico
cfg0.I = 700; cfg0.amax = 3.0; cfg0.amin = -2.0*cfg0.z1;
cfg0.r_guess=0.03; cfg0.rmin=0.005; cfg0.rmax=0.10;
cfg0.maxit_V=140; cfg0.crit_V=1e-6; cfg0.Delta=1400;
cfg0.maxit_r=1500; cfg0.crit_S=1e-5; cfg0.fix_r=0;
% Bien público
cfg0.psi_G = 0.08; cfg0.omegaG = 0.50; cfg0.report_G_effects = 0;
% Difusión
cfg0.sigma_a = 0.010;
% Gobierno / Deuda
cfg0.B_mode = 'ratio_to_Y'; cfg0.Bbar=0.80; cfg0.alphaG=0.50;
cfg0.clamp_G_to_zero = true; cfg0.G_cap_ratio = 1;

paper_style();

% --------------------- Escenarios de choque y transfer -------------------
pol0 = struct('shockI',0.00,'shockF',0.00,'xi_I',0.00,'xi_F',0.00);     % BASE
pol1 = struct('shockI',0.10,'shockF',0.05,'xi_I',0.00,'xi_F',0.00);     % Choque sin TR
pol2 = struct('shockI',0.10,'shockF',0.05,'xi_I',0.70,'xi_F',0.30);     % Choque + top-up

SC = struct('name',{},'cfg',{},'pol',{},'sol',{});
SC(1).name='BASE';        SC(1).cfg=cfg0; SC(1).pol=pol0;
SC(2).name='SHOCK_noTR';  SC(2).cfg=cfg0; SC(2).pol=pol1;
SC(3).name='SHOCK_TR';    SC(3).cfg=cfg0; SC(3).pol=pol2;

% ----------------------------- Resolver ----------------------------------
for k=1:numel(SC)
    tic;
    SC(k).sol = solve_two_type_huggett_fiscal_transfer_shock(SC(k).cfg, SC(k).pol);
    s=SC(k).sol;
    fprintf('[%s] r*=%.4f | eta=%.3f | B=%.3f | G=%.3f | Y=%.3f | C=%.3f | Tr=%.3f (%.2fs)\n', ...
        SC(k).name, s.r, s.popI/(s.popI+s.popF), s.fiscal.B, s.fiscal.G, s.Y, s.Ctot, s.fiscal.Tr, toc);

    % ---- CSVs por escenario ----
    a=s.a; g=s.g; da=a(2)-a(1);
    A_I=sum(g(:,1).*a)*da; A_F=sum(g(:,2).*a)*da; A_priv=A_I+A_F;
    scen=string(SC(k).name);

    T_fiscal = table(scen, s.fiscal.Tl, s.fiscal.Tc, s.fiscal.Tr, s.fiscal.G, s.fiscal.rB, ...
        s.fiscal.PB, s.fiscal.B, s.fiscal.BB, s.fiscal.tI_pc, s.fiscal.tF_pc, ...
        'VariableNames', {'scenario','Tl','Tc','Tr','G','rB','PB','B','BB','tI_pc','tF_pc'});
    writetable(T_fiscal, fullfile(outdir_tabs, sprintf('fiscal_%s.csv', SC(k).name)));

    T_assets = table(scen, A_I, A_F, A_priv, s.fiscal.B, ...
        'VariableNames', {'scenario','A_I','A_F','A_priv','B_public'});
    writetable(T_assets, fullfile(outdir_tabs, sprintf('assets_%s.csv', SC(k).name)));

    S = s.stats;
    T_stats = table(scen, s.r, s.popI, s.popF, s.Y, s.Ctot, ...
        S.wealth_mean(1),S.wealth_mean(2),S.wealth_mean(3), ...
        S.giniW(1),S.giniW(2),S.giniW(3), ...
        S.cons_mean(1),S.cons_mean(2),S.cons_mean(3), ...
        S.giniC(1),S.giniC(2),S.giniC(3), S.p11, ...
        'VariableNames', {'scenario','r','popI','popF','Y','C', ...
         'wmean_I','wmean_F','wmean_T','giniW_I','giniW_F','giniW_T', ...
         'cmean_I','cmean_F','cmean_T','giniC_I','giniC_F','giniC_T','p11_rep'});
    writetable(T_stats, fullfile(outdir_tabs, sprintf('hh_stats_%s.csv', SC(k).name)));
end

% ------------- Figuras solicitadas (solo las requeridas) -----------------
names = {SC.name};
cols = lines(numel(SC));

% (1) Consumo por tipo
fig=figure('Name','Consumption by type - shock & transfers'); hold on;
for k=1:numel(SC)
    a=SC(k).sol.a; c=SC(k).sol.c;
    plot(a,c(:,1),'LineWidth',2,'Color',cols(k,:),'LineStyle','-');   % I
    plot(a,c(:,2),'LineWidth',2,'Color',cols(k,:),'LineStyle','--');  % F
end
grid on; xlabel('Assets a'); ylabel('c(a)'); title('Consumo: I sólido / F punteado');
legend(names,'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'shockTR_policy_consumption'));

% (2) Ahorro por tipo
fig=figure('Name','Savings by type - shock & transfers'); hold on;
for k=1:numel(SC)
    a=SC(k).sol.a; s=SC(k).sol.s;
    plot(a,s(:,1),'LineWidth',2,'Color',cols(k,:),'LineStyle','-');
    plot(a,s(:,2),'LineWidth',2,'Color',cols(k,:),'LineStyle','--');
end
grid on; xlabel('Assets a'); ylabel('s(a)'); title('Ahorro: I sólido / F punteado');
legend(names,'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'shockTR_policy_savings'));

% (3) Lorenz riqueza (total)
fig=figure('Name','Lorenz wealth (total) - shock & transfers'); hold on;
for k=1:numel(SC)
    [~,Lw,cumPop]=lorenz_from_density(SC(k).sol.a, SC(k).sol.g);
    plot(cumPop,Lw,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
legend(names,'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'shockTR_lorenz_wealth_total'));

% (4) Lorenz consumo (total)
fig=figure('Name','Lorenz consumption (total) - shock & transfers'); hold on;
for k=1:numel(SC)
    [~,Lc,cumPopC]=lorenz_from_values(SC(k).sol.a, SC(k).sol.g, SC(k).sol.c(:,1)+SC(k).sol.c(:,2));
    plot(cumPopC,Lc,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
legend(names,'Location','southeast');
export_fig(fig, fullfile(outdir_figs,'shockTR_lorenz_consumption_total'));

% (5) Fiscal: cuentas y balances
fig=figure('Name','Fiscal accounts - shock & transfers');
% Ingresos
subplot(1,3,1);
vals=zeros(2,numel(SC));
for k=1:numel(SC), vals(:,k)=[SC(k).sol.fiscal.Tl; SC(k).sol.fiscal.Tc]; end
h=bar(categorical({'Labor','VAT'}), vals, 'grouped'); grid on; ylabel('Revenues'); title('Revenues'); legend(names,'Location','best');
add_grouped_labels(h,'%.3f');
% Gastos
subplot(1,3,2);
vals=zeros(3,numel(SC));
for k=1:numel(SC), vals(:,k)=[SC(k).sol.fiscal.rB; SC(k).sol.fiscal.G; SC(k).sol.fiscal.Tr]; end
h=bar(categorical({'Debt serv','Public G','Transfers'}), vals, 'grouped'); grid on; ylabel('Expenditures'); title('Expenditures'); legend(names,'Location','best');
add_grouped_labels(h,'%.3f');
% Balances
subplot(1,3,3);
vals=zeros(2,numel(SC));
for k=1:numel(SC), vals(:,k)=[SC(k).sol.fiscal.PB; SC(k).sol.fiscal.BB]; end
h=bar(categorical({'Primary','Overall'}), vals, 'grouped'); yline(0,'k--');
grid on; ylabel('Balance'); title('Balances'); legend(names,'Location','best');
add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'shockTR_fiscal_accounts'));

% (6) Borrowers/Lenders
fig=figure('Name','Borrowers/Lenders - shock & transfers');
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
subplot(1,2,1); % shares
Borrow = zeros(2,numel(SC));
for k=1:numel(SC), Borrow(:,k) = SC(k).sol.borrowers.fracBorrow(:); end
h=bar(catX, Borrow, 'grouped'); grid on; ylabel('Share'); title('Borrowers (a<0)'); legend(names,'Location','best');
add_grouped_labels(h,'%.3f');
subplot(1,2,2); % |Debt| volumes
VolB = zeros(2,numel(SC));
for k=1:numel(SC), VolB(:,k) = abs(SC(k).sol.borrowers.volBorrow(:)); end
h=bar(catX, VolB, 'grouped'); grid on; ylabel('Volume'); title('|Debt| volumes'); legend(names,'Location','best');
add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'shockTR_borrowers_lenders'));

% (7) MPC(a) por tipo (impulso 1% a z1 y z2, sosteniendo r)
fig=figure('Name','MPC(a) - shock & transfers'); hold on;
for k=1:numel(SC)
    cfg=SC(k).cfg; r=SC(k).sol.r; a=SC(k).sol.a;
    alt=cfg; alt.fix_r=1; alt.r_guess=r; alt.rmin=r; alt.rmax=r;
    eps=0.01; pol=SC(k).pol; pol.z1_base=cfg.z1; pol.z2_base=cfg.z2;
    pol.shockI=pol.shockI; pol.shockF=pol.shockF; % mismos shocks
    % impulso sobre z1 y z2:
    alt.z1 = cfg.z1*(1+eps); alt.z2 = cfg.z2*(1+eps);
    solP = solve_two_type_huggett_fiscal_transfer_shock(alt, pol);
    dres1 = (cfg.z1*(1+eps)-cfg.z1)*(1 + cfg.phi);
    dres2 = (cfg.z2*(1+eps)-cfg.z2)*(1 - cfg.tau_l);
    MPC1=(solP.c(:,1)-SC(k).sol.c(:,1))/max(dres1,1e-12);
    MPC2=(solP.c(:,2)-SC(k).sol.c(:,2))/max(dres2,1e-12);
    plot(a,MPC1,'Color',cols(k,:),'LineWidth',2,'LineStyle','-');
    plot(a,MPC2,'Color',cols(k,:),'LineWidth',2,'LineStyle','--');
end
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)');
legend(names,'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'shockTR_mpc_by_assets'));

% (8) A_priv(r) & B(r) y (9) S(r)
r_span = linspace(0.006,0.07,41);
fig=figure('Name','A_priv & B vs r - shock & transfers'); hold on;
hAp=gobjects(numel(SC),1); hB=hAp; hr=hAp;
for k=1:numel(SC)
    alt=SC(k).cfg; pol=SC(k).pol; alt.fix_r=1; alt.maxit_r=1;
    Apriv=nan(size(r_span)); Bgrid=Apriv;
    for i=1:numel(r_span)
        alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal_transfer_shock(alt, pol);
        da_i=tmp.a(2)-tmp.a(1);
        Apriv(i)=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
        Bgrid(i)=tmp.fiscal.B;
    end
    hAp(k)=plot(r_span,Apriv,'LineWidth',2,'Color',cols(k,:));
    hB(k)=plot(r_span,Bgrid,'--','LineWidth',1.6,'Color',cols(k,:));
    hr(k)=xline(SC(k).sol.r,':','Color',cols(k,:),'LineWidth',1.6);
end
grid on; xlabel('Interest rate r'); ylabel('Assets / Bonds level');
legend([hAp;hB;hr], [strcat(names,' - A_{priv}'), strcat(names,' - B'), strcat(names,' - r^*')], 'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'shockTR_asset_demand_vs_B'));

fig=figure('Name','Excess S(r) - shock & transfers'); hold on;
for k=1:numel(SC)
    alt=SC(k).cfg; pol=SC(k).pol; alt.fix_r=1; alt.maxit_r=1;
    Sgrid=nan(size(r_span));
    for i=1:numel(r_span)
        alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_fiscal_transfer_shock(alt, pol);
        da_i=tmp.a(2)-tmp.a(1);
        Apriv_i=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
        Sgrid(i)=Apriv_i - tmp.fiscal.B;
    end
    plot(r_span,Sgrid,'LineWidth',2,'Color',cols(k,:));
    xline(SC(k).sol.r,':','Color',cols(k,:),'LineWidth',1.3);
end
yline(0,'k--'); grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
legend(names,'Location','best');
export_fig(fig, fullfile(outdir_figs,'shockTR_excess_supply_curve'));

fprintf('Tablas en %s | Figuras en %s\n', outdir_tabs, outdir_figs);

% =============================== Helpers =================================
function paper_style()
    set(groot,'defaulttextinterpreter','tex'); set(groot,'defaultAxesTickLabelInterpreter','tex');
    set(groot,'defaultLegendInterpreter','tex'); set(groot,'DefaultAxesFontSize',12);
    set(groot,'DefaultLineLineWidth',1.9); set(groot,'DefaultFigureColor','w');
end
function export_fig(fig,basepath)
    if nargin<1||isempty(fig), fig=gcf; end
    print(fig,[basepath '.png'],'-dpng','-r300'); print(fig,[basepath '.pdf'],'-dpdf');
end
function [gT,L,cumPop]=lorenz_from_density(a,g)
    da=a(2)-a(1); gT=g(:,1)+g(:,2); W=sum(gT)*da; [as,ix]=sort(a); gs=gT(ix)*da; cumPop=cumsum(gs)/W;
    wealth_pos=as-min(0,min(as))+1e-12; cumWealth=cumsum(wealth_pos.*gs); L=cumWealth/max(cumWealth(end),1e-12);
end
function [vals,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w); vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix);
    cumPop=cumsum(w_s)/W; cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end
function add_grouped_labels(h, fmt)
    if nargin<2, fmt = '%.3f'; end
    for i=1:numel(h)
        x = h(i).XEndPoints; y = h(i).YEndPoints;
        for j=1:numel(x)
            text(x(j), y(j), sprintf([' ',fmt], y(j)), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
        end
    end
end
    