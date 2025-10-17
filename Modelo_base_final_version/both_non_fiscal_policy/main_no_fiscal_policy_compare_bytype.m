% ========== main_no_fiscal_policy_compare_bytype.m ==========
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------------------- Configuración BASE -------------------------------
cfg0 = struct();
cfg0.RRA_I=3.40; cfg0.RRA_F=3.40; cfg0.rho=0.08;
cfg0.tau_l=0.15; cfg0.tau_c=0.18; cfg0.phi=0.09;
cfg0.z1=0.33; cfg0.z2=1.00;
cfg0.theta_I=0.06; cfg0.theta_F=0.01;
cfg0.eta_target=0.654; cfg0.p22_bar=0.8155;
cfg0.I=700; cfg0.amax=3.0; cfg0.amin=-2.0*cfg0.z1;
cfg0.r_guess=0.03; cfg0.rmin=0.005; cfg0.rmax=0.10;
cfg0.maxit_V=140; cfg0.crit_V=1e-6; cfg0.Delta=1400;
cfg0.maxit_r=1500; cfg0.crit_S=1e-5; cfg0.fix_r=0;
% Bien público
cfg0.psi_G=0.08; cfg0.omegaG=0.50; cfg0.report_G_effects=0;
% Difusión
cfg0.sigma_a=0.010;
% Gobierno / Deuda
cfg0.B_mode='ratio_to_Y'; cfg0.Bbar=0.80; cfg0.alphaG=0.50;
cfg0.clamp_G_to_zero=true; cfg0.G_cap_ratio=1;

paper_style();

% -------------------- Escenarios: BASE vs NOPOL --------------------------
% Flags del “apagón fiscal”
nopol = struct('zero_taxes',true,'zero_transfers',true, ...
               'zero_public_good',true,'zero_debt',true);

SC = struct('name',{},'cfg',{},'flags',{},'sol',{});
SC(1).name  = 'BASE';  SC(1).cfg = cfg0; SC(1).flags = struct();     % nada apagado
SC(2).name  = 'NOPOL'; SC(2).cfg = cfg0; SC(2).flags = nopol;        % todo off

% ----------------------------- Resolver -----------------------------------
for k=1:numel(SC)
    tic;
    SC(k).sol = solve_two_type_huggett_no_fiscal_bytype(SC(k).cfg, SC(k).flags);
    s=SC(k).sol;
    fprintf('[%s] r*=%.4f | eta=%.3f | B=%.3f | G=%.3f | Y=%.3f | C=%.3f | Tr=%.3f (%.2fs)\n', ...
        SC(k).name, s.r, s.popI/(s.popI+s.popF), s.fiscal.B, s.fiscal.G, s.Y, s.Ctot, s.fiscal.Tr, toc);

    % --------- Tablas por escenario (I/F/Total) --------------------------
    a=s.a; g=s.g; da=a(2)-a(1); A_I=sum(g(:,1).*a)*da; A_F=sum(g(:,2).*a)*da; A_priv=A_I+A_F;
    scen=string(SC(k).name);

    f=SC(k).sol.fiscal;
    T_fiscal = table(scen, f.Tl, f.Tc, f.Tr, f.G, f.rB, f.PB, f.B, f.BB, ...
        'VariableNames', {'scenario','Tl','Tc','Tr','G','rB','PB','B','BB'});
    writetable(T_fiscal, fullfile(outdir_tabs, sprintf('fiscal_%s.csv', SC(k).name)));

    T_assets = table(scen, A_I, A_F, A_priv, f.B, ...
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

% ============================ FIGURAS =====================================
names = {SC.name};
cols  = lines(numel(SC));
mkI = @(k) {'LineStyle','-','Color',cols(k,:),'LineWidth',2};   % I sólido
mkF = @(k) {'LineStyle','--','Color',cols(k,:),'LineWidth',2};  % F punteado
labelsIF = labels_by_order(names);  % (1-I,1-F,2-I,2-F,...)

% (1) Consumo por tipo
fig=figure('Name','Consumption by type'); hold on;
for k=1:numel(SC)
    a=SC(k).sol.a; c=SC(k).sol.c;
    sty=mkI(k); plot(a,c(:,1), sty{:});
    sty=mkF(k); plot(a,c(:,2), sty{:});
end
grid on; xlabel('Assets a'); ylabel('c(a)'); title('Consumo por tipo (I sólido / F punteado)');
legend(labelsIF, 'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'nopol_policy_consumption'));

% (2) Ahorro por tipo
fig=figure('Name','Savings by type'); hold on;
for k=1:numel(SC)
    a=SC(k).sol.a; s=SC(k).sol.s;
    sty=mkI(k); plot(a,s(:,1), sty{:});
    sty=mkF(k); plot(a,s(:,2), sty{:});
end
grid on; xlabel('Assets a'); ylabel('s(a)'); title('Ahorro por tipo (I sólido / F punteado)');
legend(labelsIF, 'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'nopol_policy_savings'));

% (3) Lorenz riqueza (total) + Gini
fig=figure('Name','Lorenz wealth (total)'); hold on;
for k=1:numel(SC)
    [~,Lw,cumPop]=lorenz_from_density(SC(k).sol.a, SC(k).sol.g);
    plot(cumPop,Lw,'Color',cols(k,:),'LineWidth',2);
end
plot([0,1],[0,1],'k--'); axis square; grid on;
labs = arrayfun(@(k) sprintf('%s  (Gini=%.3f)', names{k}, SC(k).sol.stats.giniW(3)), 1:numel(SC), 'uni', 0);
legend(labs,'Location','southeast');
xlabel('Population share'); ylabel('Wealth share'); title('Lorenz riqueza (total)');
export_fig(fig, fullfile(outdir_figs,'nopol_lorenz_wealth_total'));

% (4) Lorenz consumo (total) + Gini
fig=figure('Name','Lorenz consumption (total)'); hold on;
for k=1:numel(SC)
    [~,Lc,cumPopC]=lorenz_from_values(SC(k).sol.a, SC(k).sol.g, SC(k).sol.c(:,1)+SC(k).sol.c(:,2));
    plot(cumPopC,Lc,'Color',cols(k,:),'LineWidth',2);
end
plot([0,1],[0,1],'k--'); axis square; grid on;
labs = arrayfun(@(k) sprintf('%s  (Gini=%.3f)', names{k}, SC(k).sol.stats.giniC(3)), 1:numel(SC), 'uni', 0);
legend(labs,'Location','southeast');
xlabel('Population share'); ylabel('Consumption share'); title('Lorenz consumo (total)');
export_fig(fig, fullfile(outdir_figs,'nopol_lorenz_consumption_total'));

% (5) Gini por tipo (riqueza y consumo)
fig=figure('Name','Gini by type'); 
T = array2table(zeros(numel(SC),5),'VariableNames',{'idx','gW_I','gW_F','gC_I','gC_F'});
for k=1:numel(SC)
    T.idx(k)=k; T.gW_I(k)=SC(k).sol.stats.giniW(1); T.gW_F(k)=SC(k).sol.stats.giniW(2);
    T.gC_I(k)=SC(k).sol.stats.giniC(1); T.gC_F(k)=SC(k).sol.stats.giniC(2);
end
subplot(1,2,1); hold on; title('Gini riqueza por tipo');
plot(1:numel(SC), T.gW_I,'-o','LineWidth',2); plot(1:numel(SC), T.gW_F,'--s','LineWidth',2);
set(gca,'XTick',1:numel(SC),'XTickLabel',names); grid on; ylabel('Gini');
legend({'Informal','Formal'},'Location','best'); xtickangle(20);
subplot(1,2,2); hold on; title('Gini consumo por tipo');
plot(1:numel(SC), T.gC_I,'-o','LineWidth',2); plot(1:numel(SC), T.gC_F,'--s','LineWidth',2);
set(gca,'XTick',1:numel(SC),'XTickLabel',names); grid on; ylabel('Gini');
legend({'Informal','Formal'},'Location','best'); xtickangle(20);
export_fig(fig, fullfile(outdir_figs,'nopol_ginis'));

% (6) Fiscal (debe mostrar todo en 0 para NOPOL)
fig=figure('Name','Fiscal accounts (BASE vs NOPOL)');
subplot(1,3,1); % ingresos
vals=zeros(2,numel(SC));
for k=1:numel(SC), vals(:,k)=[SC(k).sol.fiscal.Tl; SC(k).sol.fiscal.Tc]; end
h=bar(categorical(names), vals'); grid on; ylabel('Revenues'); title('Ingresos (Tl,Tc)');
legend({'Labor','VAT'},'Location','best'); add_grouped_labels(h,'%.3f');

subplot(1,3,2); % gastos
vals=zeros(3,numel(SC));
for k=1:numel(SC), vals(:,k)=[SC(k).sol.fiscal.rB; SC(k).sol.fiscal.G; SC(k).sol.fiscal.Tr]; end
h=bar(categorical(names), vals'); grid on; ylabel('Expenditures'); title('Gastos (rB,G,Tr)');
legend({'Debt serv','Public G','Transfers'},'Location','best'); add_grouped_labels(h,'%.3f');

subplot(1,3,3); % balances
vals=zeros(2,numel(SC));
for k=1:numel(SC), vals(:,k)=[SC(k).sol.fiscal.PB; SC(k).sol.fiscal.BB]; end
h=bar(categorical(names), vals'); yline(0,'k--');
grid on; ylabel('Balance'); title('Balances (PB,BB)'); legend({'Primary','Overall'},'Location','best');
add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'nopol_fiscal_accounts'));

% (7) Borrowers/Lenders por tipo
fig=figure('Name','Borrowers/Lenders by type');
catX = categorical({'Informal','Formal'}); catX=reordercats(catX,{'Informal','Formal'});
subplot(1,2,1); % shares borrowers
Bshare=zeros(2,numel(SC));
for k=1:numel(SC), Bshare(:,k)=SC(k).sol.borrowers.fracBorrow(:); end
h=bar(catX, Bshare, 'grouped'); grid on; ylabel('Share'); title('Borrowers (a<0)');
legend(names,'Location','best'); add_grouped_labels(h,'%.3f');
subplot(1,2,2); % |Debt| volumes
Bvol=zeros(2,numel(SC));
for k=1:numel(SC), Bvol(:,k)=abs(SC(k).sol.borrowers.volBorrow(:)); end
h=bar(catX, Bvol, 'grouped'); grid on; ylabel('Volume'); title('|Debt| volumes');
legend(names,'Location','best'); add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'nopol_borrowers_lenders'));

% (8) MPC(a) por tipo (r fijo en cada escenario)  *** SECCIÓN CORREGIDA ***
fig=figure('Name','MPC(a) by type'); hold on;
for k=1:numel(SC)
    cfg=SC(k).cfg; r=SC(k).sol.r; a=SC(k).sol.a; flags=SC(k).flags;
    alt=cfg; alt.fix_r=1; alt.r_guess=r; alt.rmin=r; alt.rmax=r;
    eps=0.01; alt.z1=cfg.z1*(1+eps); alt.z2=cfg.z2*(1+eps);
    solP = solve_two_type_huggett_no_fiscal_bytype(alt, flags);

    % Efectos efectivos de política para el impulso
    phi_eff = cfg.phi;
    if flags_on(flags,'zero_transfers'), phi_eff = 0; end
    tau_l_eff = cfg.tau_l;
    if flags_on(flags,'zero_taxes'), tau_l_eff = 0; end

    dres1 = (cfg.z1*(1+eps)-cfg.z1) * (1 + phi_eff);   % informal
    dres2 = (cfg.z2*(1+eps)-cfg.z2) * (1 - tau_l_eff); % formal

    MPC1=(solP.c(:,1)-SC(k).sol.c(:,1))/max(dres1,1e-12);
    MPC2=(solP.c(:,2)-SC(k).sol.c(:,2))/max(dres2,1e-12);
    sty=mkI(k); plot(a,MPC1, sty{:});
    sty=mkF(k); plot(a,MPC2, sty{:});
end
yline(0,'k:'); grid on; xlabel('Assets a'); ylabel('MPC(a)');
legend(labelsIF,'Location','bestoutside'); title('MPC por activos y tipo');
export_fig(fig, fullfile(outdir_figs,'nopol_mpc'));

% (9) Mercado de activos: A_priv(r), B(r) y S(r)
r_span = linspace(0.006,0.07,41);
fig=figure('Name','A_priv & B vs r'); hold on;
for k=1:numel(SC)
    alt=SC(k).cfg; flags=SC(k).flags; alt.fix_r=1; alt.maxit_r=1;
    Apriv=nan(size(r_span)); Bgrid=Apriv;
    for i=1:numel(r_span)
        alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_no_fiscal_bytype(alt, flags);
        da_i=tmp.a(2)-tmp.a(1);
        Apriv(i)=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
        Bgrid(i)=tmp.fiscal.B;
    end
    plot(r_span,Apriv,'-','Color',cols(k,:),'LineWidth',2);
    plot(r_span,Bgrid,'--','Color',cols(k,:),'LineWidth',1.8);
    xline(SC(k).sol.r,':','Color',cols(k,:),'LineWidth',1.6);
end
grid on; xlabel('Interest rate r'); ylabel('Levels');
legend([strcat(names,' - A_{priv}(r)'), strcat(names,' - B(r)'), strcat(names,' - r^*')], 'Location','bestoutside');
export_fig(fig, fullfile(outdir_figs,'nopol_asset_demand_vs_B'));

fig=figure('Name','Excess S(r)'); hold on;
for k=1:numel(SC)
    alt=SC(k).cfg; flags=SC(k).flags; alt.fix_r=1; alt.maxit_r=1;
    Sgrid=nan(size(r_span));
    for i=1:numel(r_span)
        alt.r_guess=r_span(i); alt.rmin=alt.r_guess; alt.rmax=alt.r_guess;
        tmp=solve_two_type_huggett_no_fiscal_bytype(alt, flags);
        da_i=tmp.a(2)-tmp.a(1);
        Apriv_i=sum((tmp.g(:,1)+tmp.g(:,2)).*tmp.a)*da_i;
        Sgrid(i)=Apriv_i - tmp.fiscal.B;
    end
    plot(r_span,Sgrid,'Color',cols(k,:),'LineWidth',2);
    xline(SC(k).sol.r,':','Color',cols(k,:),'LineWidth',1.3);
end
yline(0,'k--'); grid on; xlabel('r'); ylabel('S(r)=A_{priv}-B');
legend(names,'Location','best');
export_fig(fig, fullfile(outdir_figs,'nopol_excess_supply'));

fprintf('Listo. Tablas en %s | Figuras en %s\n', outdir_tabs, outdir_figs);

% ============================ Helpers =====================================
function on = flags_on(flags,field)
    on = isfield(flags,field) && islogical(flags.(field)) && flags.(field);
end
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
function labs = labels_by_order(names)
    n = numel(names); labs = cell(1, 2*n); j=1;
    for k=1:n
        labs{j}   = [names{k} ' - I']; 
        labs{j+1} = [names{k} ' - F']; 
        j=j+2;
    end
end
