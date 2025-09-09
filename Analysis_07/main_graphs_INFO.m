%% main_graphs_INFO.m
% Construye las figuras 1..5 usando las variables del solver
clear; clc; close all;

% TEX por defecto; usar LaTeX solo en ecuaciones
set(groot,'defaulttextinterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');

% --- ejecutar el script que llena el workspace
run huggett_Equi_INFO_Jan24

amax1 = 0.4;
amin1 = amin-0.03;

%% [1] Policy functions
figure('Name','Fig1-INFO')
% CONSUMPTION
subplot(2,2,1)
h1=plot(a,C{1}(:,1),'r',a,C{1}(:,2),'b',...
        a,C{2}(:,1),'r--',a,C{2}(:,2),'b--',...
        linspace(amin1,amax1,I),zeros(1,I),'k--','LineWidth',1.4);
xline(amin,'k--','LineWidth',1.0);
lgd = legend(h1, ...
    ['$c_1$ ($\eta$=',num2str(pop1_vector(1)),')'], '$c_2$', ...
    ['$c_1$ ($\eta$=',num2str(pop1_vector(2)),')'], '$c_2$', ...
    'Location','south','Interpreter','latex'); 
lgd.NumColumns = 2;
xlabel('Wealth, $a$','Interpreter','latex'); 
ylabel('Consumption, $c_i(a)$','Interpreter','latex')
xlim([amin1 amax1]); title('Consumption'); grid on;

% SAVING
subplot(2,2,2)
h2=plot(a,Saving{1}(:,1),'r',a,Saving{1}(:,2),'b',...
        a,Saving{2}(:,1),'r--',a,Saving{2}(:,2),'b--',...
        linspace(amin1,amax1,I),zeros(1,I),'k--','LineWidth',1.4);
xline(amin,'k--','LineWidth',1.0);
lgd = legend(h2, ...
    ['$s_1$ ($\eta$=',num2str(pop1_vector(1)),')'], '$s_2$', ...
    ['$s_1$ ($\eta$=',num2str(pop1_vector(2)),')'], '$s_2$', ...
    'Location','south','Interpreter','latex'); 
lgd.NumColumns = 2;
xlabel('Wealth, $a$','Interpreter','latex'); 
ylabel('Saving, $s_i(a)$','Interpreter','latex')
xlim([amin1 amax1]); title('Saving'); grid on;

set(gcf,'Color','w'); set(gcf,'PaperPosition',[0 0 1 1]);
print(gcf,'-dpdf','Fig1_INFO.pdf');

%% [2] Distribution: wealth & consumption
figure('Name','Fig2-INFO')

% Wealth densities (eta1)
subplot(2,3,1)
plot(a,Distribution{1}(:,1),'r',a,Distribution{1}(:,2),'b:','LineWidth',1.4)
xline(amin,'k--','LineWidth',1.0)
legend('$g_1(a)$ (informal)','$g_2(a)$ (formal)','Interpreter','latex','Location','northeast');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('Densities, $g_i(a)$','Interpreter','latex');
xlim([amin1 amax1]); title('Wealth Density'); subtitle(['($\eta$=',num2str(pop1_vector(1)),')'],'Interpreter','latex'); grid on;

% Wealth densities (eta2)
subplot(2,3,2)
plot(a,Distribution{2}(:,1),'r',a,Distribution{2}(:,2),'b:','LineWidth',1.4)
xline(amin,'k--','LineWidth',1.0)
legend('$g_1(a)$ (informal)','$g_2(a)$ (formal)','Interpreter','latex','Location','northeast');
xlabel('Wealth, $a$','Interpreter','latex');
xlim([amin1 amax1]); title('Wealth Density'); subtitle(['($\eta$=',num2str(pop1_vector(2)),')'],'Interpreter','latex'); grid on;

% Wealth densities (eta2, 2 ejes)
subplot(2,3,3); colororder({'r','b'}); yyaxis left
plot(a,Distribution{2}(:,1),'r','LineWidth',1.4); ylabel('$g_1(a)$ (informal)','Interpreter','latex');
yyaxis right
plot(a,Distribution{2}(:,2),'b:','LineWidth',1.4); ylabel('$g_2(a)$ (formal)','Interpreter','latex');
xline(amin,'k--','LineWidth',1.0)
xlabel('Wealth, $a$','Interpreter','latex'); xlim([amin1 amax1]); title('Wealth Density'); subtitle(['($\eta$=',num2str(pop1_vector(2)),')'],'Interpreter','latex'); grid on;

% Consumption densities (proxy g vs c)
subplot(2,3,4)
plot(C{1}(:,1),CDistribution{1}(:,1),'r', C{1}(:,2),CDistribution{1}(:,2),'b:','LineWidth',1.4)
legend('$g_1(c)$','$g_2(c)$','Interpreter','latex','Location','northeast');
xlabel('Consumption, $c$','Interpreter','latex'); ylabel('Densities, $g_i(c)$','Interpreter','latex');
title('Consumption Density'); subtitle(['($\eta$=',num2str(pop1_vector(1)),')'],'Interpreter','latex'); grid on;

subplot(2,3,5)
plot(C{2}(:,1),CDistribution{2}(:,1),'r', C{2}(:,2),CDistribution{2}(:,2),'b:','LineWidth',1.4)
legend('$g_1(c)$','$g_2(c)$','Interpreter','latex','Location','northeast');
xlabel('Consumption, $c$','Interpreter','latex'); ylabel('Densities, $g_i(c)$','Interpreter','latex');
title('Consumption Density'); subtitle(['($\eta$=',num2str(pop1_vector(2)),')'],'Interpreter','latex'); grid on;

subplot(2,3,6); colororder({'r','b'}); yyaxis left
plot(C{2}(:,1),CDistribution{2}(:,1),'r','LineWidth',1.4); ylabel('$g_1(c)$','Interpreter','latex');
yyaxis right
plot(C{2}(:,2),CDistribution{2}(:,2),'b:','LineWidth',1.4); ylabel('$g_2(c)$','Interpreter','latex');
xlabel('Consumption, $c$','Interpreter','latex'); title('Consumption Density'); subtitle(['($\eta$=',num2str(pop1_vector(2)),')'],'Interpreter','latex'); grid on;

set(gcf,'Color','w'); set(gcf,'PaperPosition',[0 0 1 1]);
print(gcf,'-dpdf','Fig2_INFO.pdf');

%% [3] Total distribution: wealth
figure('Name','Fig3-INFO')
subplot(2,2,1)
plot(a,GDistribution{1},'k--',a,Distribution{1}(:,1),'r',a,Distribution{1}(:,2),'b:','LineWidth',1.4)
xline(amin,'k--','LineWidth',1.0)
legend('Total','$g_1(a)$','$g_2(a)$','Interpreter','latex','Location','northeast');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('Densities','Interpreter','latex')
xlim([amin1 amax1]); title(['Wealth Density ($\eta$=',num2str(pop1_vector(1)),')'],'Interpreter','latex'); grid on;

subplot(2,2,2)
plot(a,GDistribution{2},'k--',a,Distribution{2}(:,1),'r',a,Distribution{2}(:,2),'b:','LineWidth',1.4)
xline(amin,'k--','LineWidth',1.0)
legend('Total','$g_1(a)$','$g_2(a)$','Interpreter','latex','Location','northeast');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('Densities','Interpreter','latex')
xlim([amin1 amax1]); title(['Wealth Density ($\eta$=',num2str(pop1_vector(2)),')'],'Interpreter','latex'); grid on;

% Total consumption per wealth level (proxy)
subplot(2,2,3)
plot(a, GCdistribution{1},'r', a, CDistribution{1}(:,1),'b:', a, CDistribution{1}(:,2),'k--','LineWidth',1.4)
xline(amin,'k--','LineWidth',1.0)
legend('Total Cons. (proxy)','Informal (proxy)','Formal (proxy)','Location','best');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('Consumption'); xlim([amin1 amax1]); title('Total Consumption vs Wealth'); grid on;

subplot(2,2,4)
plot(a, GCdistribution{2},'r', a, CDistribution{2}(:,1),'b:', a, CDistribution{2}(:,2),'k--','LineWidth',1.4)
xline(amin,'k--','LineWidth',1.0)
legend('Total Cons. (proxy)','Informal (proxy)','Formal (proxy)','Location','best');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('Consumption'); xlim([amin1 amax1]); title('Total Consumption vs Wealth'); grid on;

set(gcf,'Color','w'); set(gcf,'PaperPosition',[0 0 1 1]);
print(gcf,'-dpdf','Fig3_INFO.pdf');

%% [4] Comparación de distribuciones por tipo entre etas
figure('Name','Fig4-INFO')
subplot(2,2,1)
plot(a,Distribution{1}(:,1),'b--',a,Distribution{2}(:,1),'r:','LineWidth',1.4)
legend(['$\eta$=',num2str(pop1_vector(1))],['$\eta$=',num2str(pop1_vector(2))],'Interpreter','latex','Location','northeast');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('$g_1(a)$','Interpreter','latex')
xlim([amin1 amax1]); title('Wealth Density — Informal'); grid on;

subplot(2,2,2)
plot(a,Distribution{1}(:,2),'b--',a,Distribution{2}(:,2),'r:','LineWidth',1.4)
legend(['$\eta$=',num2str(pop1_vector(1))],['$\eta$=',num2str(pop1_vector(2))],'Interpreter','latex','Location','northeast');
xlabel('Wealth, $a$','Interpreter','latex'); ylabel('$g_2(a)$','Interpreter','latex')
xlim([amin1 amax1]); title('Wealth Density — Formal'); grid on;

subplot(2,2,3)
plot(C{1}(:,1),CDistribution{1}(:,1),'b--',C{2}(:,1),CDistribution{2}(:,1),'r:','LineWidth',1.4)
legend(['$\eta$=',num2str(pop1_vector(1))],['$\eta$=',num2str(pop1_vector(2))],'Interpreter','latex','Location','northeast');
xlabel('Consumption, $c$','Interpreter','latex'); ylabel('$g_1(c)$','Interpreter','latex')
title('Consumption Density — Informal'); grid on;

subplot(2,2,4)
plot(C{1}(:,2),CDistribution{1}(:,2),'b--',C{2}(:,2),CDistribution{2}(:,2),'r:','LineWidth',1.4)
legend(['$\eta$=',num2str(pop1_vector(1))],['$\eta$=',num2str(pop1_vector(2))],'Interpreter','latex','Location','northeast');
xlabel('Consumption, $c$','Interpreter','latex'); ylabel('$g_2(c)$','Interpreter','latex')
title('Consumption Density — Formal'); grid on;

set(gcf,'Color','w'); set(gcf,'PaperPosition',[0 0 1 1]);
print(gcf,'-dpdf','Fig4_INFO.pdf');
