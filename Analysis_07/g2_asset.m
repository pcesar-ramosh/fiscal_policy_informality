% g2_asset.m
if ~exist('AssetDS','var') || ~exist('eta','var')
    warning('AssetDS o eta no existen. Omite g2_asset.'); return;
end

IaB = AssetDS(:,1); FaB = AssetDS(:,2);
IaL = AssetDS(:,3); FaL = AssetDS(:,4);
totalSupply = AssetDS(:,5); totalDemand = AssetDS(:,6);

figure('Name','Assets: supply/demand + r(eta)')
yyaxis left
plot(eta, IaB, '-o','LineWidth',1.4); hold on;
plot(eta, FaB, '-s','LineWidth',1.4);
plot(eta, IaL, '--o','LineWidth',1.4);
plot(eta, FaL, '--s','LineWidth',1.4);
plot(eta, totalSupply, '-.','LineWidth',1.6);
plot(eta, totalDemand, ':','LineWidth',1.6);
ylabel('Activos agregados (∫ a g(a) da)'); grid on;

yyaxis right
if exist('r_all','var')
    plot(eta, r_all, 'k-','LineWidth',1.6);
    ylabel('$r(\eta)$','Interpreter','latex');
end
xlabel('Informalidad $\eta$','Interpreter','latex');
legend('IaB (inf, a<0)','FaB (for, a<0)','IaL (inf, a>0)','FaL (for, a>0)', ...
       'Total Supply','Total Demand','r(\eta)','Location','best','Interpreter','latex');
title('Demanda/Oferta de activos por tipo y $r(\eta)$','Interpreter','latex');
set(gcf,'Color','w');
