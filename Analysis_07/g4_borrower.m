% g4_borrower.m
if ~exist('BL','var') || ~exist('eta','var')
    % Construir BL si no viene de stats_Eta_v2
    if exist('N_borrowers','var')
        gini_tot = eta_matrix(:,9);
        BL = [ eta, N_borrowers, frac_N, frac_NI, frac_NF, ...
               frac_NIT, frac_NFT, borrowing, lending, 100*gini_tot ];
    else
        warning('No hay datos para BL. Omite g4_borrower.'); return;
    end
end

figure('Name','Borrowers: fracciones y agregados')
subplot(2,2,1)
plot(eta, BL(:,3),'-o','LineWidth',1.4);
xlabel('$\eta$','Interpreter','latex'); ylabel('Frac. prestatarios (total)'); grid on; title('Prestatarios totales');

subplot(2,2,2)
plot(eta, BL(:,4),'-o', eta, BL(:,5),'-s','LineWidth',1.4);
legend('Informal','Formal','Location','best');
xlabel('$\eta$','Interpreter','latex'); ylabel('Frac. prestatarios (por tipo)'); grid on; title('Por tipo');

subplot(2,2,3)
plot(eta, BL(:,8),'-o', eta, BL(:,9),'-s','LineWidth',1.4);
legend('Borrowing (a<0)','Lending (a>0)','Location','best');
xlabel('$\eta$','Interpreter','latex'); ylabel('Agregados de activos'); grid on; title('Borrowing vs Lending');

subplot(2,2,4)
plot(eta, BL(:,10),'-^','LineWidth',1.4);
xlabel('$\eta$','Interpreter','latex'); ylabel('Gini total (%)'); grid on; title('Desigualdad (Gini)');
set(gcf,'Color','w');
