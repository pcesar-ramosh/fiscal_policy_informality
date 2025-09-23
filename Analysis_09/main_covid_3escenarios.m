%% main_covid_3escenarios.m  (FIX TITLES)
clear; clc; close all;

%% ---- 1) Parámetros comunes ----
RRA   = 3.40; rho=0.05; theta=0.02;
tau_l = 0.15; tau_c=0.18; Gov=0.05;

phi_base = 0.09;    % normal
phi_high = 0.14;    % alto
phi_zero = 0.00;    % sin transferencias

% Ingresos COVID
z1_base = 0.33;  z2_base = 1.00;
z1_cvd  = z1_base*(1-0.25);   % -25% informales
z2_cvd  = z2_base*(1-0.15);   % -15% formales

eta_target = 0.654; p22_bar = 0.8155;
Igrid=700; amax=5.0;

common = struct('RRA_I',RRA,'RRA_F',RRA,'rho',rho,'theta',theta, ...
    'tau_l',tau_l,'tau_c',tau_c,'Gov',Gov,'I',Igrid,'amax',amax, ...
    'r_guess',0.03,'rmin',0.005,'rmax',0.08, ...
    'p22_bar',p22_bar,'eta_target',eta_target);

%% ---- 2) Escenarios COVID ----
parC0 = common; parC0.z1=z1_cvd; parC0.z2=z2_cvd; parC0.phi=phi_zero; parC0.amin=-0.30*parC0.z1;
COV_noTr = huggett_base_covid_function(parC0);

parCB = common; parCB.z1=z1_cvd; parCB.z2=z2_cvd; parCB.phi=phi_base; parCB.amin=-0.30*parCB.z1;
COV_phiB = huggett_base_covid_function(parCB);

parCH = parCB; parCH.phi=phi_high;
COV_phiH = huggett_base_covid_function(parCH);

%% ---- 3) Consola rápida ----
fprintf('\n== r == COV sinTr %.5f | COV phi_base %.5f | COV phi_alto %.5f\n', ...
    COV_noTr.r, COV_phiB.r, COV_phiH.r);
fprintf('   popI (informal) : %.4f | %.4f | %.4f\n', ...
    COV_noTr.popI, COV_phiB.popI, COV_phiH.popI);
fprintf('   Y (PIB ingreso) : %.4f | %.4f | %.4f\n', ...
    COV_noTr.Y, COV_phiB.Y, COV_phiH.Y);
fprintf('   Ctot (consumo)  : %.4f | %.4f | %.4f\n', ...
    COV_noTr.Ctot, COV_phiB.Ctot, COV_phiH.Ctot);

%% ---- 4) Exportar CSVs ----
if ~exist('./tables','dir'), mkdir('./tables'); end
F1=COV_noTr.fiscal; F2=COV_phiB.fiscal; F3=COV_phiH.fiscal;

F = table(["COVID_sinTr";"COVID_phiBase";"COVID_phiAlto"], ...
    [F1.Tl;F2.Tl;F3.Tl],[F1.Tc;F2.Tc;F3.Tc], ...
    [F1.Tr;F2.Tr;F3.Tr],[F1.G;F2.G;F3.G],[F1.rB;F2.rB;F3.rB], ...
    [F1.B;F2.B;F3.B],[F1.PB;F2.PB;F3.PB],[F1.BB;F2.BB;F3.BB], ...
    'VariableNames',{'scenario','Tl','Tc','Tr','G','rB','B','PB','BB'});
F.IngresosTot = F.Tl + F.Tc;
F.GastosTot   = F.Tr + F.G + F.rB;
writetable(F,'./tables/covid3_fiscal.csv');

da = COV_noTr.a(2)-COV_noTr.a(1);
Apriv = @(out) sum((out.g(:,1)+out.g(:,2)).*out.a)*da;
A = table(["COVID_sinTr";"COVID_phiBase";"COVID_phiAlto"], ...
    [Apriv(COV_noTr); Apriv(COV_phiB); Apriv(COV_phiH)], ...
    [F1.B;F2.B;F3.B], 'VariableNames',{'scenario','A_private','B_public'});
writetable(A,'./tables/covid3_asset_market.csv');

packStats = @(S,out) table(string(S), out.r, out.popI, out.popF, out.Y, out.Ctot, ...
   out.stats.wealth_mean(1), out.stats.wealth_mean(2), out.stats.wealth_mean(3), ...
   out.stats.wealth_median(1), out.stats.wealth_median(2), out.stats.wealth_median(3), ...
   out.stats.gini(1), out.stats.gini(2), out.stats.gini(3), ...
   out.stats.cons_mean(1), out.stats.cons_mean(2), out.stats.cons_mean(3), ...
   out.stats.cons_median(1), out.stats.cons_median(2), out.stats.cons_median(3), ...
   out.stats.p11, ...
   'VariableNames',{'scenario','r','popI','popF','Y','Ctot', ...
     'w_mean_I','w_mean_F','w_mean_T','w_med_I','w_med_F','w_med_T', ...
     'gini_I','gini_F','gini_T','c_mean_I','c_mean_F','c_mean_T', ...
     'c_med_I','c_med_F','c_med_T','p11'});
TS = [packStats("COVID_sinTr",COV_noTr);
      packStats("COVID_phiBase",COV_phiB);
      packStats("COVID_phiAlto",COV_phiH)];
writetable(TS,'./tables/covid3_household_stats.csv');

B = table(["COVID_sinTr";"COVID_phiBase";"COVID_phiAlto"], ...
    [COV_noTr.borrowers.fracBorrow(1); COV_phiB.borrowers.fracBorrow(1); COV_phiH.borrowers.fracBorrow(1)], ...
    [COV_noTr.borrowers.fracBorrow(2); COV_phiB.borrowers.fracBorrow(2); COV_phiH.borrowers.fracBorrow(2)], ...
    [COV_noTr.borrowers.fracLend(1);   COV_phiB.borrowers.fracLend(1);   COV_phiH.borrowers.fracLend(1)], ...
    [COV_noTr.borrowers.fracLend(2);   COV_phiB.borrowers.fracLend(2);   COV_phiH.borrowers.fracLend(2)], ...
    [COV_noTr.borrowers.volBorrow(1);  COV_phiB.borrowers.volBorrow(1);  COV_phiH.borrowers.volBorrow(1)], ...
    [COV_noTr.borrowers.volBorrow(2);  COV_phiB.borrowers.volBorrow(2);  COV_phiH.borrowers.volBorrow(2)], ...
    [COV_noTr.borrowers.volLend(1);    COV_phiB.borrowers.volLend(1);    COV_phiH.borrowers.volLend(1)], ...
    [COV_noTr.borrowers.volLend(2);    COV_phiB.borrowers.volLend(2);    COV_phiH.borrowers.volLend(2)], ...
    'VariableNames',{'scenario','fracBorrow_I','fracBorrow_F','fracLend_I','fracLend_F', ...
                     'volBorrow_I','volBorrow_F','volLend_I','volLend_F'});
writetable(B,'./tables/covid3_borrowers.csv');

%% ---- 5) Gráficos ----
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

a = COV_noTr.a;

% (1) Consumo c(a)
figure('Name','Consumo por tipo - COVID 3 escenarios');
for j=1:2
    subplot(1,2,j)
    plot(a, COV_noTr.c(:,j),'LineWidth',1.6); hold on;
    plot(a, COV_phiB.c(:,j),'--','LineWidth',1.6);
    plot(a, COV_phiH.c(:,j),':','LineWidth',1.8);
    yline(0,'k:'); grid on; xlabel('a'); ylabel(sprintf('c_{%d}(a)',j));
    legend('sinTr','phiBase','phiAlto','Location','best');
    if j==1, title('Informales'); else, title('Formales'); end
end

% (2) Ahorro s(a)
figure('Name','Ahorro por tipo - COVID 3 escenarios');
for j=1:2
    subplot(1,2,j)
    plot(a, COV_noTr.s(:,j),'LineWidth',1.6); hold on;
    plot(a, COV_phiB.s(:,j),'--','LineWidth',1.6);
    plot(a, COV_phiH.s(:,j),':','LineWidth',1.8);
    yline(0,'k:'); grid on; xlabel('a'); ylabel(sprintf('s_{%d}(a)',j));
    legend('sinTr','phiBase','phiAlto','Location','best');
    if j==1, title('Informales'); else, title('Formales'); end
end

% (3) Distribución g(a)
figure('Name','Distribucion de riqueza - COVID 3 escenarios');
subplot(2,1,1);
bar(a, COV_noTr.g(:,1),'FaceAlpha',0.35,'EdgeColor','none'); hold on;
bar(a, COV_phiB.g(:,1),'FaceAlpha',0.35,'EdgeColor','none');
bar(a, COV_phiH.g(:,1),'FaceAlpha',0.35,'EdgeColor','none');
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g_1(a)');
legend('sinTr','phiBase','phiAlto'); title('Informales');
subplot(2,1,2);
bar(a, COV_noTr.g(:,2),'FaceAlpha',0.35,'EdgeColor','none'); hold on;
bar(a, COV_phiB.g(:,2),'FaceAlpha',0.35,'EdgeColor','none');
bar(a, COV_phiH.g(:,2),'FaceAlpha',0.35,'EdgeColor','none');
xlim([min(a) 1.0]); grid on; xlabel('a'); ylabel('g_2(a)');
legend('sinTr','phiBase','phiAlto'); title('Formales');

% (4) Mercado de activos (demanda privada)
Apr_no  = Apriv(COV_noTr);  Apr_b = Apriv(COV_phiB);  Apr_h = Apriv(COV_phiH);
figure('Name','Mercado de activos - COVID 3 escenarios');
bar(categorical({'sinTr','phiBase','phiAlto'}), [Apr_no Apr_b Apr_h]);
ylabel('Demanda privada (∫ a g da)'); grid on; title('Demanda privada');

% (5) Ingresos y gastos desagregados
X = categorical({'sinTr','phiBase','phiAlto'});

figure('Name','Ingresos desagregados - COVID');
subplot(1,2,1); bar(X,[F1.Tl F2.Tl F3.Tl]); grid on; ylabel('Ingresos'); title('T_\ell (laboral)');
subplot(1,2,2); bar(X,[F1.Tc F2.Tc F3.Tc]); grid on; ylabel('Ingresos'); title('T_c (IVA)');

figure('Name','Gastos desagregados - COVID');
subplot(1,3,1); bar(X,[F1.Tr F2.Tr F3.Tr]); title('Tr'); grid on;
subplot(1,3,2); bar(X,[F1.G  F2.G  F3.G ]); title('G');  grid on;
subplot(1,3,3); bar(X,[F1.rB F2.rB F3.rB]); title('rB'); grid on;

% (6) Curva de Lorenz
figure('Name','Curva de Lorenz - COVID 3 escenarios');
hold on; outs = {COV_noTr, COV_phiB, COV_phiH}; names={'sinTr','phiBase','phiAlto'};
for k=1:3
    out = outs{k}; gT = out.g(:,1)+out.g(:,2); da = out.a(2)-out.a(1);
    W = sum(gT)*da; [as,ix]=sort(out.a); gs=gT(ix)*da;
    cumPop=cumsum(gs)/W; wealth_pos = as - min(0,min(as)) + 1e-12;
    L=cumsum(wealth_pos.*gs); L=L/L(end);
    plot(cumPop,L,'LineWidth',1.6);
end
plot([0,1],[0,1],'k--'); axis square; grid on;
legend(names,'Location','southeast'); xlabel('Frac. población'); ylabel('Frac. riqueza');

% (7) Borrowers / Lenders
cats = categorical({'Informal','Formal'}); cats=reordercats(cats,{'Informal','Formal'});
figure('Name','Fracciones prestatarios - COVID');
subplot(1,3,1); bar(cats, COV_noTr.borrowers.fracBorrow(:)); title('Borrow sinTr'); grid on;
subplot(1,3,2); bar(cats, COV_phiB.borrowers.fracBorrow(:)); title('Borrow phiBase'); grid on;
subplot(1,3,3); bar(cats, COV_phiH.borrowers.fracBorrow(:)); title('Borrow phiAlto'); grid on;

figure('Name','Deuda y Ahorro agregados - COVID');
subplot(1,2,1);
bar(cats, [abs(COV_noTr.borrowers.volBorrow(:)), abs(COV_phiB.borrowers.volBorrow(:)), abs(COV_phiH.borrowers.volBorrow(:))]);
legend('sinTr','phiBase','phiAlto','Location','best'); ylabel('|Deuda agregada|'); grid on;
subplot(1,2,2);
bar(cats, [COV_noTr.borrowers.volLend(:), COV_phiB.borrowers.volLend(:), COV_phiH.borrowers.volLend(:)]);
legend('sinTr','phiBase','phiAlto','Location','best'); ylabel('Ahorro agregado'); grid on;

disp('Listo: corre sin error y muestra todos los gráficos/CSVs.');
