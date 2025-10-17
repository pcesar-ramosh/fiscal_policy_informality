% =========================== main_compare_vat.m ===========================
clear; clc; close all;
outdir_tabs = './tables';  if ~exist(outdir_tabs,'dir'),  mkdir(outdir_tabs);  end
outdir_figs = './figures'; if ~exist(outdir_figs,'dir'), mkdir(outdir_figs); end

% ---------------------- BASE configuration (igual que tu base) -----------
cfg0 = struct();
cfg0.RRA_I = 3.40; cfg0.RRA_F = 3.40;
cfg0.rho   = 0.08;
cfg0.tau_l = 0.15; cfg0.tau_c = 0.18; cfg0.phi = 0.09;
cfg0.z1 = 0.33; cfg0.z2 = 1.00;

cfg0.theta_I = 0.06;                 % informal premium
cfg0.theta_F = 0.01;                 % formal premium

cfg0.eta_target = 0.654;             % informality target (fijo en este experimento)
cfg0.p22_bar   = 0.8155;

cfg0.I = 700; cfg0.amax = 3.0; cfg0.amin = -2.0*cfg0.z1;

cfg0.r_guess = 0.03; cfg0.rmin = 0.005; cfg0.rmax = 0.10;
cfg0.maxit_V = 140; cfg0.crit_V = 1e-6; cfg0.Delta = 1400;
cfg0.maxit_r = 1500; cfg0.crit_S = 1e-5; cfg0.fix_r = 0;

% Bien público (moderado, como en tu versión)
cfg0.psi_G   = 0.08;
cfg0.omegaG  = 0.50;
cfg0.report_G_effects = 0;

% Difusión
cfg0.sigma_a = 0.010;

% Gobierno / Deuda
cfg0.B_mode  = 'ratio_to_Y';     % 'ratio_to_Y' or 'level'
cfg0.Bbar    = 0.80;             % B/Y objetivo
cfg0.alphaG  = 0.50;
cfg0.clamp_G_to_zero = true;
cfg0.G_cap_ratio = 1;            % sin cap efectivo (1 * Y)

paper_style();

% --------------------------- Escenarios τ_c -------------------------------
% Primero BASE (reportado con su τ_c), luego subidas del IVA:
tau_vec = [0.18, 0.20, 0.22, 0.26];  % puedes ajustar o extender
labels  = arrayfun(@(x) sprintf('\\tau_c = %.2f', x), tau_vec, 'uni', 0);

S = struct('name',{},'tau_c',{},'cfg',{},'sol',{});
for k=1:numel(tau_vec)
    S(k).name  = sprintf('VAT_%02d', round(100*tau_vec(k)));
    S(k).tau_c = tau_vec(k);
    S(k).cfg   = cfg0;
end

% ----------------------------- Resolver ----------------------------------
cols = lines(numel(tau_vec));
for k=1:numel(S)
    tic;
    S(k).sol = solve_two_type_huggett_fiscal_tauC(S(k).cfg, S(k).tau_c);
    s = S(k).sol;
    fprintf('[%s] τ_c=%.2f | r*=%.4f | η=%.3f | B=%.3f | G=%.3f | Y=%.3f | C=%.3f  (%.2fs)\n', ...
        S(k).name, S(k).tau_c, s.r, s.popI/(s.popI+s.popF), ...
        s.fiscal.B, s.fiscal.G, s.Y, s.Ctot, toc);

    % ---- CSVs por escenario ----
    a  = s.a; g = s.g; da = a(2)-a(1);
    A_I = sum(g(:,1).*a)*da; A_F = sum(g(:,2).*a)*da; A_priv = A_I + A_F;

    scen = string(S(k).name);

    T_fiscal = table(scen, s.fiscal.Tl, s.fiscal.Tc, ...
        s.fiscal.Tl+s.fiscal.Tc, s.fiscal.Tr, ...
        s.fiscal.G, s.fiscal.rB, s.fiscal.PB, ...
        s.fiscal.B, s.fiscal.BB, ...
        'VariableNames', {'scenario','labor_tax','vat_tax','rev_total','transfers', ...
                          'public_good','debt_serv','primary_bal','debt_stock','global_bal'});
    writetable(T_fiscal, fullfile(outdir_tabs, sprintf('fiscal_%s.csv', S(k).name)));

    T_assets = table(scen, A_priv, s.fiscal.B, ...
        'VariableNames', {'scenario','A_private','B_public'});
    writetable(T_assets, fullfile(outdir_tabs, sprintf('assets_%s.csv', S(k).name)));

    st = s.stats;
    T_stats = table(scen, s.r, s.popI, s.popF, s.Y, s.Ctot, ...
        st.wealth_mean(1),st.wealth_mean(2),st.wealth_mean(3), ...
        st.wealth_median(1),st.wealth_median(2),st.wealth_median(3), ...
        st.giniW(1),st.giniW(2),st.giniW(3), ...
        st.cons_mean(1),st.cons_mean(2),st.cons_mean(3), ...
        st.cons_median(1),st.cons_median(2),st.cons_median(3), ...
        st.giniC(1),st.giniC(2),st.giniC(3), s.stats.p11, ...
        'VariableNames', {'scenario','r','popI','popF','Y','Ctot', ...
        'wealth_mean_I','wealth_mean_F','wealth_mean_T', ...
        'wealth_med_I','wealth_med_F','wealth_med_T', ...
        'giniW_I','giniW_F','giniW_T', ...
        'cons_mean_I','cons_mean_F','cons_mean_T', ...
        'cons_med_I','cons_med_F','cons_med_T', ...
        'giniC_I','giniC_F','giniC_T','p11_rep'});
    writetable(T_stats, fullfile(outdir_tabs, sprintf('hh_stats_%s.csv', S(k).name)));
end

% ------------------- Resumen comparativo y elasticidades ------------------
RES = cell(numel(S), 17);
for k=1:numel(S)
    s = S(k).sol; da = s.a(2)-s.a(1);
    A_priv = sum((s.g(:,1)+s.g(:,2)).*s.a)*da;
    RES(k,:) = { labels{k}, S(k).tau_c, s.r, s.fiscal.B, A_priv, ...
                 s.popI/(s.popI+s.popF), s.Y, s.Ctot, s.fiscal.G, ...
                 s.fiscal.Tc, s.fiscal.PB, s.fiscal.BB, ...
                 s.stats.giniW(3), s.stats.giniC(3), ...
                 s.fiscal.Tl, s.fiscal.Tr, s.fiscal.rB };
end
Tsum = cell2table(RES, 'VariableNames', ...
 {'scenario','tau_c','r','B','A_priv','eta_real','Y','C','G','Tc','PB','BB','giniW_T','giniC_T','Tl','Tr','rB'});
writetable(Tsum, fullfile(outdir_tabs,'compare_vat_summary.csv'));
disp(Tsum);

% Aproximaciones de semi-elasticidades (centrales si hay >=3 puntos)
if numel(S)>=3
    tau = Tsum.tau_c; [tauS,ix] = sort(tau);
    C  = Tsum.C(ix);   Yv = Tsum.Y(ix);  rv = Tsum.r(ix);  PBv = Tsum.PB(ix);
    Tc = Tsum.Tc(ix);  Gv = Tsum.G(ix);

    dCdt   = grad_central(tauS, C);
    dYdt   = grad_central(tauS, Yv);
    drdt   = grad_central(tauS, rv);
    dPBdt  = grad_central(tauS, PBv);
    dTc_dt = grad_central(tauS, Tc);
    dGdt   = grad_central(tauS, Gv);

    E = table(tauS, dCdt./max(C,1e-12), dYdt./max(Yv,1e-12), drdt, dPBdt, dTc_dt, dGdt, ...
        'VariableNames', {'tau_c','semi_elast_logC','semi_elast_logY','dr_dtau','dPB_dtau','dTc_dtau','dG_dtau'});
    writetable(E, fullfile(outdir_tabs,'compare_vat_semielasticities.csv'));
    disp(E);
end

% ========================= FIGURAS COMPARATIVAS ===========================
cols = lines(numel(S));

% ---------- POLICIES c(a) ----------
fig=figure('Name','Policies: consumption by τ_c'); hold on;
hI=gobjects(numel(S),1); hF=hI;
for k=1:numel(S)
    a=S(k).sol.a; c=S(k).sol.c;
    hI(k)=plot(a,c(:,1),'Color',cols(k,:),'LineWidth',2,'LineStyle','-');     % Informal
    hF(k)=plot(a,c(:,2),'Color',cols(k,:),'LineWidth',2,'LineStyle','--');    % Formal
end
grid on; xlabel('Assets a'); ylabel('Consumption c(a)');
leg = [cellfun(@(s)[s ' - I'],labels,'uni',0), cellfun(@(s)[s ' - F'],labels,'uni',0)];
legend([hI;hF], leg, 'Location','bestoutside');
title('Policies c(a): I sólido, F punteado');
export_fig(fig, fullfile(outdir_figs,'vat_policy_consumption'));

% ---------- POLICIES s(a) ----------
fig=figure('Name','Policies: savings by τ_c'); hold on;
hI=gobjects(numel(S),1); hF=hI;
for k=1:numel(S)
    a=S(k).sol.a; s=S(k).sol.s;
    hI(k)=plot(a,s(:,1),'Color',cols(k,:),'LineWidth',2,'LineStyle','-');
    hF(k)=plot(a,s(:,2),'Color',cols(k,:),'LineWidth',2,'LineStyle','--');
end
grid on; xlabel('Assets a'); ylabel('Savings s(a)');
legend([hI;hF], leg, 'Location','bestoutside');
title('Policies s(a)');
export_fig(fig, fullfile(outdir_figs,'vat_policy_savings'));

% ---------- Densidades de riqueza ----------
fig=figure('Name','Wealth densities by τ_c');
subplot(2,1,1); hold on;
for k=1:numel(S), plot(S(k).sol.a, S(k).sol.g(:,1), 'LineWidth',2,'Color',cols(k,:)); end
grid on; xlabel('Assets a'); ylabel('g_I(a)'); title('Informal density'); legend(labels,'Location','best');
subplot(2,1,2); hold on;
for k=1:numel(S), plot(S(k).sol.a, S(k).sol.g(:,2), 'LineWidth',2,'Color',cols(k,:)); end
grid on; xlabel('Assets a'); ylabel('g_F(a)'); title('Formal density'); legend(labels,'Location','best');
export_fig(fig, fullfile(outdir_figs,'vat_wealth_distributions'));

% ---------- Cierre del mercado de activos ----------
fig=figure('Name','Asset market closure by τ_c (stacked)');
tiledlayout(1,numel(S),'Padding','compact','TileSpacing','compact');
for k=1:numel(S)
    nexttile;
    a=S(k).sol.a; g=S(k).sol.g; da=a(2)-a(1);
    A_I=sum(g(:,1).*a)*da; A_F=sum(g(:,2).*a)*da; B=S(k).sol.fiscal.B;
    Ybars = [A_I A_F 0; 0 0 B];
    bar(categorical({'Private demand','Public debt B'}), Ybars, 'stacked');
    grid on; ylabel('Level'); title(sprintf('%s (r^*=%.4f)', labels{k}, S(k).sol.r));
    legend({'A_I','A_F','B'},'Location','northeast');
    xt = [1 2]; yt = [A_I+A_F, B];
    for j=1:2, text(xt(j), yt(j), sprintf(' %.3f',yt(j)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontSize',10); end
end
export_fig(fig, fullfile(outdir_figs,'vat_asset_market_closure'));

% ---------- Cuentas fiscales ----------
fig=figure('Name','Fiscal accounts and balances (τ_c compare)');
% Ingresos
subplot(1,3,1);
vals=zeros(2,numel(S));
for k=1:numel(S), vals(:,k)=[S(k).sol.fiscal.Tl; S(k).sol.fiscal.Tc]; end
h=bar(categorical({'Labor tax','VAT'}), vals, 'grouped'); grid on; ylabel('Revenues'); title('Revenues'); legend(labels,'Location','best');
add_grouped_labels(h,'%.3f');
% Gastos
subplot(1,3,2);
vals=zeros(3,numel(S));
for k=1:numel(S), vals(:,k)=[S(k).sol.fiscal.rB; S(k).sol.fiscal.G; S(k).sol.fiscal.Tr]; end
h=bar(categorical({'Debt service','Public good','Transfers'}), vals, 'grouped'); grid on; ylabel('Expenditures'); title('Expenditures'); legend(labels,'Location','best');
add_grouped_labels(h,'%.3f');
% Balances
subplot(1,3,3);
vals=zeros(2,numel(S));
for k=1:numel(S), vals(:,k)=[S(k).sol.fiscal.PB; S(k).sol.fiscal.BB]; end
h=bar(categorical({'Primary balance','Overall (global)'}), vals, 'grouped'); yline(0,'k--');
grid on; ylabel('Balance'); title('Balances'); legend(labels,'Location','best');
add_grouped_labels(h,'%.3f');
export_fig(fig, fullfile(outdir_figs,'vat_fiscal_accounts'));

% ---------- Lorenz (total, riqueza y consumo) ----------
fig=figure('Name','Lorenz wealth (total) vs τ_c'); hold on;
hW = gobjects(numel(S),1);
for k=1:numel(S)
    [~,Lw,cumPop]=lorenz_from_density(S(k).sol.a, S(k).sol.g);
    hW(k)=plot(cumPop,Lw,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Wealth share');
labs = arrayfun(@(k) sprintf('%s (Gini=%.3f)', labels{k}, S(k).sol.stats.giniW(3)), 1:numel(S),'uni',0);
legend(hW, labs, 'Location','southeast');
title('Lorenz curve (Wealth, total)');
export_fig(fig, fullfile(outdir_figs,'vat_lorenz_wealth_total'));

fig=figure('Name','Lorenz consumption (total) vs τ_c'); hold on;
hC = gobjects(numel(S),1);
for k=1:numel(S)
    [~,Lc,cumPopC]=lorenz_from_values(S(k).sol.a, S(k).sol.g, S(k).sol.c(:,1)+S(k).sol.c(:,2));
    hC(k)=plot(cumPopC,Lc,'LineWidth',2,'Color',cols(k,:));
end
plot([0,1],[0,1],'k--'); grid on; axis square; xlabel('Population share'); ylabel('Consumption share');
labs = arrayfun(@(k) sprintf('%s (Gini=%.3f)', labels{k}, S(k).sol.stats.giniC(3)), 1:numel(S),'uni',0);
legend(hC, labs, 'Location','southeast');
title('Lorenz curve (Consumption, total)');
export_fig(fig, fullfile(outdir_figs,'vat_lorenz_consumption_total'));

% ---------- Curvas r(τ_c), Y(τ_c), C(τ_c), G(τ_c), PB(τ_c) ----------
if exist('Tsum','var')
    fig=figure('Name','Aggregates vs τ_c'); hold on;
    yyaxis left
    plot(Tsum.tau_c, Tsum.r, '-o', 'LineWidth',2); ylabel('Interest rate r');
    yyaxis right
    plot(Tsum.tau_c, Tsum.Y, '-s', 'LineWidth',2); 
    plot(Tsum.tau_c, Tsum.C, '-d', 'LineWidth',2);
    plot(Tsum.tau_c, Tsum.G, '-^', 'LineWidth',2);
    plot(Tsum.tau_c, Tsum.PB,'-v', 'LineWidth',2);
    grid on; xlabel('\tau_c'); ylabel('Levels'); legend('r','Y','C','G','PB','Location','best');
    title('General equilibrium vs VAT (\tau_c)');
    export_fig(fig, fullfile(outdir_figs,'vat_aggregates_vs_tauC'));
end

fprintf('Figures saved to %s\n', outdir_figs);

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
function [L,cumPop]=lorenz_from_assets_single(a, gk)
    da=a(2)-a(1); W=sum(gk)*da; [as,ix]=sort(a); gs=gk(ix)*da; cumPop=cumsum(gs)/max(W,1e-12);
    wealth_pos = as - min(0,min(as)) + 1e-12; cumWealth = cumsum(wealth_pos.*gs);
    L = cumWealth / max(cumWealth(end),1e-12);
end
function [vals,L,cumPop]=lorenz_from_values(a,g,x)
    da=a(2)-a(1); w=(g(:,1)+g(:,2))*da; W=sum(w); vals=x(:); [vals_s,ix]=sort(vals); w_s=w(ix);
    cumPop=cumsum(w_s)/W; cumx=cumsum(vals_s.*w_s); L=cumx/max(cumx(end),1e-12);
end
function [L,cumPop]=lorenz_from_values_single(xk, gk, a)
    da=a(2)-a(1); w = gk(:)*da; W=sum(w); [xs,ix]=sort(xk(:)); ws=w(ix);
    cumPop = cumsum(ws)/max(W,1e-12); cumx = cumsum(xs.*ws); L = cumx / max(cumx(end),1e-12);
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
function g = grad_central(x, y)
    % derivada central simple (extremos con primera diferencia)
    n = numel(x); g = nan(size(y));
    if n>=3
        g(1)   = (y(2)-y(1))/max(x(2)-x(1),1e-12);
        for i=2:n-1
            g(i) = (y(i+1)-y(i-1))/max(x(i+1)-x(i-1),1e-12);
        end
        g(n)   = (y(n)-y(n-1))/max(x(n)-x(n-1),1e-12);
    end
end
