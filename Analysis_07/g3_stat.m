% g3_stat.m
if ~exist('eta_matrix','var') || ~exist('statsCMatrix','var') || ~exist('eta','var')
    warning('Faltan matrices de momentos. Omite g3_stat.'); return;
end

% Riqueza
gmean_inf = eta_matrix(:,1); gmean_for = eta_matrix(:,2);
gmed_inf  = eta_matrix(:,3); gmed_for  = eta_matrix(:,4);
gmean_tot = eta_matrix(:,5); gmed_tot  = eta_matrix(:,6);
gini_inf  = eta_matrix(:,7); gini_for  = eta_matrix(:,8); gini_tot = eta_matrix(:,9);

% Consumo
Cmean_inf = statsCMatrix(:,1); Cmean_for = statsCMatrix(:,2);
Cmed_inf  = statsCMatrix(:,3); Cmed_for  = statsCMatrix(:,4);
Cmean_tot = statsCMatrix(:,5); Cmed_tot  = statsCMatrix(:,6);

figure('Name','Stats: wealth & consumption vs eta')
subplot(2,2,1)
plot(eta, gmean_inf,'-o', eta, gmean_for,'-s', eta, gmean_tot,'-^','LineWidth',1.4);
legend('$E[a|inf]$','$E[a|for]$','$E[a]$','Interpreter','latex','Location','best');
xlabel('$\eta$','Interpreter','latex'); ylabel('Media de riqueza'); grid on; title('Medias de riqueza');

subplot(2,2,2)
plot(eta, gmed_inf,'-o', eta, gmed_for,'-s', eta, gmed_tot,'-^','LineWidth',1.4);
legend('$Med[a|inf]$','$Med[a|for]$','$Med[a]$','Interpreter','latex','Location','best');
xlabel('$\eta$','Interpreter','latex'); ylabel('Mediana de riqueza'); grid on; title('Medianas de riqueza');

subplot(2,2,3)
plot(eta, gini_inf,'-o', eta, gini_for,'-s', eta, gini_tot,'-^','LineWidth',1.4);
legend('$Gini_{inf}$','$Gini_{for}$','$Gini_{tot}$','Interpreter','latex','Location','best');
xlabel('$\eta$','Interpreter','latex'); ylabel('Gini'); grid on; title('Gini de riqueza');

subplot(2,2,4)
plot(eta, Cmean_inf,'-o', eta, Cmean_for,'-s', eta, Cmean_tot,'-^','LineWidth',1.4);
legend('$E[c|inf]$','$E[c|for]$','$E[c]$','Interpreter','latex','Location','best');
xlabel('$\eta$','Interpreter','latex'); ylabel('Media de consumo'); grid on; title('Consumo medio');
set(gcf,'Color','w');
