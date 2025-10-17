% ===================== main_compare_taxes_slim.m =========================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------------------- Configuración base común -------------------------
cfg0 = struct();
% Preferencias
cfg0.RRA_I = 3.40; cfg0.RRA_F = 3.40; cfg0.rho = 0.08;
% Impuestos base
cfg0.tau_l = 0.15; cfg0.tau_c = 0.18; cfg0.phi = 0.09;
% Ingresos
cfg0.z1 = 0.33; cfg0.z2 = 1.00;
% Primas por tipo
cfg0.theta_I = 0.06; cfg0.theta_F = 0.01;
% Ocupacional
cfg0.eta_target = 0.654; cfg0.p22_bar = 0.8155;
% Grid y numérico
cfg0.I=700; cfg0.amax=3.0; cfg0.amin=-2.0*cfg0.z1;
cfg0.r_guess=0.03; cfg0.rmin=0.005; cfg0.rmax=0.10;
cfg0.maxit_V=140; cfg0.crit_V=1e-6; cfg0.Delta=1400;
cfg0.maxit_r=1500; cfg0.crit_S=1e-5; cfg0.fix_r=0;
% Bien público
cfg0.psi_G=0.08; cfg0.omegaG=0.50; cfg0.report_G_effects=0;
% Difusión
cfg0.sigma_a=0.010;
% Deuda pública
cfg0.B_mode='ratio_to_Y'; cfg0.Bbar=0.80; cfg0.alphaG=0.50;
cfg0.clamp_G_to_zero=true; cfg0.G_cap_ratio=1;
% Incidencia IVA (por si se usa VAT no neutral)
cfg0.chi_I = 0.40;         % IVA efectivo informal = chi_I * tau_c
% Opciones experimentales (apagadas)
cfg0.rebate_rule='none'; cfg0.rebate_share=0.0;
cfg0.taul_swap=0.0;  cfg0.tau_c_base=0.18;      % (para VAT)
cfg0.vat_swap = 0.0; cfg0.tau_l_base=0.15;      % (para LAB)

paper_style();

% ========================== Escenarios IVA ===============================
tau_vec = [0.18 0.20 0.22 0.26];
VAT = run_tax_scenarios('VAT', cfg0, tau_vec);

% Guardar tablas por escenario (I/F/Total)
writetable(VAT.T_stats,  fullfile(outdir_tabs,'VAT_hh_stats.csv'));
writetable(VAT.T_fiscal, fullfile(outdir_tabs,'VAT_fiscal.csv'));
writetable(VAT.T_assets, fullfile(outdir_tabs,'VAT_assets.csv'));
writetable(VAT.T_summary,fullfile(outdir_tabs,'VAT_summary.csv'));

% ====================== Escenarios impuesto laboral ======================
taul_vec = [0.12 0.15 0.18 0.22];
LAB = run_tax_scenarios('LAB', cfg0, taul_vec);

writetable(LAB.T_stats,  fullfile(outdir_tabs,'LAB_hh_stats.csv'));
writetable(LAB.T_fiscal, fullfile(outdir_tabs,'LAB_fiscal.csv'));
writetable(LAB.T_assets, fullfile(outdir_tabs,'LAB_assets.csv'));
writetable(LAB.T_summary,fullfile(outdir_tabs,'LAB_summary.csv'));

% ============================= FIGURAS ===================================
% (1) Consumo y (2) Ahorro por tipo, (3) Lorenz riqueza, (4) Lorenz consumo,
% (5) Cuentas fiscales, (6) Borrowers/Lenders, (7) MPC(a) por tipo,
% (8) A_priv(r) & B(r) y (9) S(r).

% ---------- Utilidades internas para plot ----------
plot_consumption_savings_bytype(VAT, outdir_figs, 'VAT');
plot_consumption_savings_bytype(LAB, outdir_figs, 'LAB');

plot_lorenz_bytype(VAT, outdir_figs, 'VAT');
plot_lorenz_bytype(LAB, outdir_figs, 'LAB');

plot_fiscal_accounts(VAT, outdir_figs, 'VAT');
plot_fiscal_accounts(LAB, outdir_figs, 'LAB');

plot_borrowers_lenders(VAT, outdir_figs, 'VAT');
plot_borrowers_lenders(LAB, outdir_figs, 'LAB');

plot_mpc_by_assets(VAT, outdir_figs, 'VAT');
plot_mpc_by_assets(LAB, outdir_figs, 'LAB');

plot_asset_market(VAT, outdir_figs, 'VAT');
plot_asset_market(LAB, outdir_figs, 'LAB');

fprintf('Listo. Tablas en %s y figuras en %s\n', outdir_tabs, outdir_figs);

% ============================ Helpers plot ===============================
function plot_consumption_savings_bytype(R, outdir, tag)
    cols = lines(numel(R.sc));
    % Consumo
    fig=figure('Name',['Consumption by type - ' tag]); hold on;
    for k=1:numel(R.sc)
        a=R.sc(k).sol.a; c=R.sc(k).sol.c;
        plot(a,c(:,1),'LineWidth',2,'Color',cols(k,:),'LineStyle','-');   % Informal
        plot(a,c(:,2),'LineWidth',2,'Color',cols(k,:),'LineStyle','--');  % Formal
    end
    grid on; xlabel('Assets a'); ylabel('c(a)');
    legend(R.labels,'Location','bestoutside'); title(['Consumo por tipo - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_policy_consumption']));
    % Ahorro
    fig=figure('Name',['Savings by type - ' tag]); hold on;
    for k=1:numel(R.sc)
        a=R.sc(k).sol.a; s=R.sc(k).sol.s;
        plot(a,s(:,1),'LineWidth',2,'Color',cols(k,:),'LineStyle','-');
        plot(a,s(:,2),'LineWidth',2,'Color',cols(k,:),'LineStyle','--');
    end
    grid on; xlabel('Assets a'); ylabel('s(a)');
    legend(R.labels,'Location','bestoutside'); title(['Ahorro por tipo - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_policy_savings']));
end

function plot_lorenz_bytype(R, outdir, tag)
    cols = lines(numel(R.sc));
    % Lorenz riqueza (total)
    fig=figure('Name',['Lorenz wealth total - ' tag]); hold on;
    for k=1:numel(R.sc)
        [~,Lw,cumPop]=lorenz_from_density(R.sc(k).sol.a, R.sc(k).sol.g);
        plot(cumPop,Lw,'LineWidth',2,'Color',cols(k,:));
    end
    plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
    labs = arrayfun(@(k) sprintf('%s (Gini=%.3f)', R.labels{k}, R.sc(k).sol.stats.giniW(3)), 1:numel(R.sc),'uni',0);
    legend(labs,'Location','southeast'); title(['Lorenz riqueza - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_lorenz_wealth_total']));
    % Lorenz consumo (total)
    fig=figure('Name',['Lorenz consumption total - ' tag]); hold on;
    for k=1:numel(R.sc)
        [~,Lc,cumPopC]=lorenz_from_values(R.sc(k).sol.a, R.sc(k).sol.g, R.sc(k).sol.c(:,1)+R.sc(k).sol.c(:,2));
        plot(cumPopC,Lc,'LineWidth',2,'Color',cols(k,:));
    end
    plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
    labs = arrayfun(@(k) sprintf('%s (Gini=%.3f)', R.labels{k}, R.sc(k).sol.stats.giniC(3)), 1:numel(R.sc),'uni',0);
    legend(labs,'Location','southeast'); title(['Lorenz consumo - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_lorenz_consumption_total']));
end

function plot_fiscal_accounts(R, outdir, tag)
    fig=figure('Name',['Fiscal accounts - ' tag]);
    % Ingresos
    subplot(1,3,1);
    vals=zeros(2,numel(R.sc));
    for k=1:numel(R.sc), vals(:,k)=[R.sc(k).sol.fiscal.Tl; R.sc(k).sol.fiscal.Tc]; end
    h=bar(categorical({'Labor','VAT'}), vals, 'grouped'); grid on; ylabel('Revenues'); title('Revenues'); legend(R.labels,'Location','best');
    add_grouped_labels(h,'%.3f');
    % Gastos
    subplot(1,3,2);
    vals=zeros(3,numel(R.sc));
    for k=1:numel(R.sc), vals(:,k)=[R.sc(k).sol.fiscal.rB; R.sc(k).sol.fiscal.G; R.sc(k).sol.fiscal.Tr]; end
    h=bar(categorical({'Debt serv','Public G','Transfers'}), vals, 'grouped'); grid on; ylabel('Expenditures'); title('Expenditures'); legend(R.labels,'Location','best');
    add_grouped_labels(h,'%.3f');
    % Balances
    subplot(1,3,3);
    vals=zeros(2,numel(R.sc));
    for k=1:numel(R.sc), vals(:,k)=[R.sc(k).sol.fiscal.PB; R.sc(k).sol.fiscal.BB]; end
    h=bar(categorical({'Primary','Overall'}), vals, 'grouped'); yline(0,'k--');
    grid on; ylabel('Balance'); title('Balances'); legend(R.labels,'Location','best');
    add_grouped_labels(h,'%.3f');
    export_fig(fig, fullfile(outdir, [lower(tag) '_fiscal_accounts']));
end

function plot_borrowers_lenders(R, outdir, tag)
    fig=figure('Name',['Borrowers/Lenders - ' tag]);
    catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
    subplot(1,2,1); % shares
    Borrow = zeros(2,numel(R.sc));
    for k=1:numel(R.sc), Borrow(:,k) = R.sc(k).sol.borrowers.fracBorrow(:); end
    h=bar(catX, Borrow, 'grouped'); grid on; ylabel('Share'); title('Borrowers (a<0)'); legend(R.labels,'Location','best');
    add_grouped_labels(h,'%.3f');
    subplot(1,2,2); % volumes |Debt|
    VolB = zeros(2,numel(R.sc));
    for k=1:numel(R.sc), VolB(:,k) = abs(R.sc(k).sol.borrowers.volBorrow(:)); end
    h=bar(catX, VolB, 'grouped'); grid on; ylabel('Volume'); title('|Debt| volumes'); legend(R.labels,'Location','best');
    add_grouped_labels(h,'%.3f');
    export_fig(fig, fullfile(outdir, [lower(tag) '_borrowers_lenders']));
end

function plot_mpc_by_assets(R, outdir, tag)
    % Calcula MPC con incremento porcentual pequeño en z1 y z2 manteniendo r
    fig=figure('Name',['MPC(a) - ' tag]); hold on;
    cols = lines(numel(R.sc));
    for k=1:numel(R.sc)
        cfg=R.sc(k).cfg; r=R.sc(k).sol.r; a=R.sc(k).sol.a;
        alt=cfg; alt.fix_r=1; alt.r_guess=r; alt.rmin=r; alt.rmax=r;
        eps=0.01; alt.z1=cfg.z1*(1+eps); alt.z2=cfg.z2*(1+eps);
        % Resolver con el solver correspondiente al impuesto (usa handle guardado)
        solP = R.solver_handle(alt, R.sc(k).param_value, R.extra);
        dres1 = (cfg.z1*(1+eps)-cfg.z1) * (1 + cfg.phi);
        dres2 = (cfg.z2*(1+eps)-cfg.z2) * (1 - cfg.tau_l);
        MPC1=(solP.c(:,1)-R.sc(k).sol.c(:,1))/max(dres1,1e-12);
        MPC2=(solP.c(:,2)-R.sc(k).sol.c(:,2))/max(dres2,1e-12);
        plot(a,MPC1,'Color',cols(k,:),'LineWidth',2,'LineStyle','-');
        plot(a,MPC2,'Color',cols(k,:),'LineWidth',2,'LineStyle','--');
    end
    yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)');
    legend(R.labels,'Location','bestoutside'); title(['MPC por activos - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_mpc_by_assets']));
end

function plot_asset_market(R, outdir, tag)
    r_span = linspace(0.006,0.07,41);
    cols = lines(numel(R.sc));
    % A_priv(r) y B(r)
    fig=figure('Name',['A_priv & B vs r - ' tag]); hold on;
    hAp=gobjects(numel(R.sc),1); hB=hAp; hr=hAp;
    for k=1:numel(R.sc)
        alt = R.sc(k).cfg; alt.fix_r=1; alt.maxit_r=1;
        Apriv=nan(size(r_span)); Bgrid=Apriv;
        for i=1:numel(r_span)
            alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
            tmp = R.solver_handle(alt, R.sc(k).param_value, R.extra);
            da_i=tmp.a(2)-tmp.a(1);
            Apriv(i)=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
            Bgrid(i)=tmp.fiscal.B;
        end
        hAp(k)=plot(r_span,Apriv,'LineWidth',2,'Color',cols(k,:));
        hB(k) =plot(r_span,Bgrid,'--','LineWidth',1.6,'Color',cols(k,:));
        hr(k) =xline(R.sc(k).sol.r,':','Color',cols(k,:),'LineWidth',1.6);
    end
    grid on; xlabel('Interest rate r'); ylabel('Assets / Bonds level');
    legend([hAp;hB;hr],[strcat(R.labels,' - A_{priv}'), strcat(R.labels,' - B'), strcat(R.labels,' - r^*')],'Location','bestoutside');
    title(['A_{priv}(r) & B(r) - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_asset_demand_vs_B']));
    % S(r)
    fig=figure('Name',['Excess supply S(r) - ' tag]); hold on;
    for k=1:numel(R.sc)
        alt = R.sc(k).cfg; alt.fix_r=1; alt.maxit_r=1;
        Sgrid=nan(size(r_span));
        for i=1:numel(r_span)
            alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
            tmp=R.solver_handle(alt, R.sc(k).param_value, R.extra);
            da_i=tmp.a(2)-tmp.a(1);
            Apriv_i=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
            Sgrid(i)=Apriv_i - tmp.fiscal.B;
        end
        plot(r_span,Sgrid,'LineWidth',2,'Color',cols(k,:));
        xline(R.sc(k).sol.r,':','Color',cols(k,:),'LineWidth',1.3);
    end
    yline(0,'k--'); grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
    legend(R.labels,'Location','best'); title(['Excess-asset curve - ' tag]);
    export_fig(fig, fullfile(outdir, [lower(tag) '_excess_supply_curve']));
end

% ------------------------------ utilidades -------------------------------
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
            text(x(j), y(j), sprintf([' ',fmt], y(j)), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
        end
    end
end
