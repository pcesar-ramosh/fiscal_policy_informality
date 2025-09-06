%% Run Model Solution
%% [2] Many values for "RRA" (informal)
% this version reproduces graphs in the main tex
% Informality size: eta=64%
% I use this m-file in the main tex (May 2023, Feb 2024)

clear; close all;

%% Parameters
% eta_vector  = 0.64*ones(1,19);% = 0.64*ones(1,9); % For Informal Economy  % &&&&&
% ds = (0.3 - 0.15)/18; %= (0.3 - 0.15)/8; % espaciado uniforme (ds) entre valores de RRA (9 puntos ⇒ 8 intervalos.)
% sI_vector1  = 0.15:ds:0.3;% 0.15:ds:0.3; % RRA = s    vector de aversión al riesgo relativa (RRA) para los trabajadores informales.
% sF_vector1  = 0.15:ds:0.3; % 0.15:ds:0.3;
% % sF_vector1  = 0.30*ones(1,19); %= 0.15*ones(1,9); % For Formal Economy Define un vector de RRA constante para los trabajadores formales.
% % &&&&&

%% Parameters ajustables
n_agents = 10; % Inicial estaba con 9. also we add this in for structure
% Definición del rango para la aversión al riesgo relativa (RRA)
s_min = 0.15; %0.15 0.3;  % 0.15;  With this values the code works  + lambda of function
s_max = 0.30; %0.30 0.9; % 0.30;
% Vectores automáticos para trabajadores informales y formales (heterogéneos)
eta_vector   = 0.75 * ones(1, n_agents);   % 0.75 Stable (alta 85% - baja 20%)             % Mismo tamaño que los tipos de agentes
sI_vector1   = linspace(s_min, s_max, n_agents);        % RRA heterogénea informal
% sF_vector1   = linspace(s_min, s_max, n_agents);        % RRA heterogénea formal
% Alternativa: RRA constante para formales (descomenta si se desea usar)
 sF_vector1 = 0.30 * ones(1, n_agents);       % 0.30 Initial 0.15        % RRA constante formal
% Definir un set de parametros iniciales, 1 de esos en teroi definir
% calibrar los parametros, valores iniciales, del modelo base comparen los
% resultados de wealth, cons, policy function, cuando el RRA es mayor a 1
% (enteros) comparens es en la RRA, con cual de los dos nos quedamos en
% escenario base? y con eso poodemoms agregar el escenerario

%% Simulation (Solves the General Equilibrium Model for each input parameter combination)
%[r, ir, pop1_vector, statsMatrix1, statsCMatrix1, GDistribution1, a, Distribution1] = ...
%    huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1);

[r_opt, ir, pop1_vector, a, g_opt, c_opt] =...
    huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1);




%%=================================
% ANALYSIS
%%=================================

%% Preliminary
z1 = 0.33; %informal % ingreso informal
I= 700; %1340; %1000; % cuántos puntos va a tener tu grilla de activos.
amin = -0.3*z1; % 30% of the min income   límite inferior de la grilla de activos
amax = 5; % límite superior de la riqueza
da = (amax-amin)/(I-1); % es el tamaño del paso entre cada punto de la grilla.



%% Getting stats for each RRA value of the informal agent

for j=1:length(sI_vector1)

%% Total wealth distribution (for every level of "a")
g = g_opt{j};

g1 = g_opt{1}(:,1);
g2 = g_opt{1}(:,2);

G = g1 + g2; % # of TOTAL people per "a" level
                     % it is not expressed in "fraction" of total people
check3 = sum(G)*da; % = 1

%% People at constraint

pp_fraction = g/sum(sum(g)); % pp/TotalPop % OBS: is the number of persons in each asset state and work state (formal/informal)
pp_fraction_const = pp_fraction(1,:); % percentage of total people at the constraint

%% Distribution Statistics (wealth)
% wealth mean for: Inf | For
stats.gmean = sum(a.*g)./sum(g); 
%stats.gsd = [sqrt(var(a,g(:,1))) sqrt(var(a,g(:,2)))]; %SD
% Note: g1 and g2 must be +
g1_old = g1;
g2_old = g2;

g1(g1 < 0) = 0;
g2(g2 < 0) = 0;

stats.gsd = [sqrt(var(a,g1)) sqrt(var(a,g2))]; %SD
    D = [a a];%a.*g;
    W = g./sum(g);

%WW{jj} = W; % weights   

    stats.gmedian = [weightedMedian(D(:,1),W(:,1))...
                     weightedMedian(D(:,2),W(:,2))];

%skewness is a measure of the asymmetry of the probability distribution 
%negative skew commonly indicates that the tail is on the left side of the distribution
stats.gskew = skewness(g);        % skew Normal Distr = 0 (Symmetry)

%A positive value tells you that you have heavy-tails (i.e. a lot of data in your tails).
stats.gkurt = kurtosis(g);        % kurt Standard normal Dis = 3

% G stats:
G_old = G;
G(G < 0) = 0;

stats.Gmean = sum(a.*G)./sum(G); 
stats.Gsd = sqrt(var(a,G)); %SD
    D = a;%a.*G;
    W1G = G./sum(G);
stats.Gmedian = [weightedMedian(D,W1G)];

% Gini
%gini = ginicoeff(G/sum(G), D/sum(D));

% Stats (Wealth)
%statsMatrix(1,:) = [stats.gmean stats.gsd stats.gmedian stats.gskew stats.gkurt... % informal | formal
%                     stats.Gmean stats.Gsd stats.Gmedian... % Total
%                     p11 r pp_fraction_const]; % more info
r = r_opt(j);
statsMatrix(j,:) = [stats.gmean stats.gsd stats.gmedian stats.gskew stats.gkurt... % informal | formal
                     stats.Gmean stats.Gsd stats.Gmedian... % Total
                     r pp_fraction_const]; % more info

disp(stats)

%% Distribution Statistics (consumption)
% WRONG! Cons_distribution = c.*g; % the total consumption per "a" level
                          % for informal and formal
                          % "g" is the # of people per "a" level
                          
% Total consumption distribution (per level of "c")
    % GC = Cons_distribution(:,1) + Cons_distribution(:,2); % No correct!
    % GC = g(:,1) + g(:,2); % informal ppl + formal ppl
        % GC is "G" previously defined. So, I use "G"

% consumption mean for: Inf | For
c = c_opt{j};

statsC.gmean = sum(c.*g)./sum(g); 
statsC.gsd = [sqrt(var(c(:,1),g1)) sqrt(var(c(:,2),g2))]; %SD
    DC = c; %c.*g;
    W = g./sum(g);
statsC.gmedian = [weightedMedian(DC(:,1),W(:,1))...
                  weightedMedian(DC(:,2),W(:,2))];


%skewness is a measure of the asymmetry of the probability distribution 
%negative skew commonly indicates that the tail is on the left side of the distribution
statsC.gskew = skewness(c.*g);        % skew Normal Distr = 0 (Symmetry)

%A positive value tells you that you have heavy-tails (i.e. a lot of data in your tails).
statsC.gkurt = kurtosis(c.*g);        % kurt Standard normal Dis = 3

% GC stats:
%statsC.Gmean = sum(a.*GC)./sum(GC); % this is WRONG
ct = [c(:,1);c(:,2)];
gt = [g1;g2];

statsC.Gmean = sum(gt.*ct)/sum(gt);%sum(G.*GC)./sum(G); 
statsC.Gsd = sqrt(var(ct,gt)); %sqrt(var(GC,G)); %SD
    D1 = ct; %GC;%GC.*G;
    W1 = gt./sum(gt); %G./sum(G);
statsC.Gmedian = [weightedMedian(D1,W1)];

% Stats
statsCMatrix(j,:) = [statsC.gmean statsC.gsd statsC.gmedian...
                      statsC.gskew statsC.gkurt...
                      statsC.Gmean statsC.Gsd statsC.Gmedian];

disp(statsC)
end

%% Export Statistics Wealth and Consumption 

% Para statsMatrix (Estadísticas de riqueza)
% Nombres de columnas para statsMatrix
varnames_stats = {...
    'gmean_inf', 'gmean_form', ...
    'gsd_inf', 'gsd_form', ...
    'gmedian_inf', 'gmedian_form', ...
    'gskew_inf', 'gskew_form', ...
    'gkurt_inf', 'gkurt_form', ...
    'Gmean', 'Gsd', 'Gmedian', ...
    'InterestRate', 'ConstraintInf', 'ConstraintForm'};

% Convertir en tabla
T_stats = array2table(statsMatrix, 'VariableNames', varnames_stats);

% Para statsCMatrix (Estadísticas de consumo)
% Nombres de columnas para statsCMatrix
varnames_statsC = {...
    'cmean_inf', 'cmean_form', ...
    'csd_inf', 'csd_form', ...
    'cmedian_inf', 'cmedian_form', ...
    'cskew_inf', 'cskew_form', ...
    'ckurt_inf', 'ckurt_form', ...
    'Cmean', 'Csd', 'Cmedian'};

% Convertir en tabla
T_statsC = array2table(statsCMatrix, 'VariableNames', varnames_statsC);

%Exportar ambas tablas al archivo Excel
filename = 'summary_stats.xlsx';

% Escribir ambas tablas en diferentes hojas
writetable(T_stats, filename, 'Sheet', 'WealthStats');
writetable(T_statsC, filename, 'Sheet', 'ConsumptionStats');


%% Export results in Excel File

% Supongamos que ya ejecutaste el modelo y tienes:
% c_opt y g_opt como celdas de 1x20, cada una con matriz 700x2

% Define el número de agentes
%n_agents = length(c_opt);  % normalmente 20
%a = linspace(-0.3*0.33, 5, 700)';  % si no está definido en tu entorno actual

% Crear archivos Excel
filename_c = 'consumption_data.xlsx';
filename_g = 'wealth_data.xlsx';

% Loop para exportar cada combinación
for j = 1:n_agents
    % Extrae matrices
    c = c_opt{j};  % consumo
    g = g_opt{j};  % distribución de población

    % Crea una tabla para Excel
    T_c = table(a, c(:,1), c(:,2), 'VariableNames', {'Assets', 'c_informal', 'c_formal'});
    T_g = table(a, g(:,1), g(:,2), 'VariableNames', {'Assets', 'g_informal', 'g_formal'});

    % Nombre de la hoja: RRA_j (puedes cambiar el formato)
    sheet_name = ['RRA_', num2str(j)];

    % Escribir en archivos Excel   
    writetable(T_c, filename_c, 'Sheet', sheet_name);
    writetable(T_g, filename_g, 'Sheet', sheet_name);
end

%% Some Plots

for j = 1:n_agents
    figure;
    
    subplot(1,2,1)
    plot(a, c_opt{j}(:,1), 'r--', 'LineWidth', 1.5); hold on;
    plot(a, c_opt{j}(:,2), 'b--', 'LineWidth', 1.5);
    xlabel('Activos (a)'); ylabel('Consumo');
    legend('Informal', 'Formal');
    title(['Consumo óptimo - RRA informal #' num2str(j)]);
    xlim([1 5]);  % opcional
    grid on;
    
    subplot(1,2,2)
    plot(a, g_opt{j}(:,1), 'r--', 'LineWidth', 1.5); hold on;
    plot(a, g_opt{j}(:,2), 'b--', 'LineWidth', 1.5);
    xlabel('Activos (a)'); ylabel('Distribución de riqueza');
    legend('Informal', 'Formal');
    title(['Distribución g(a) - RRA informal #' num2str(j)]);
    xlim([0.01 0.6]);
    %ylim([0 5]);  % ajusta según tus curvas
    grid on;

end



figure;

% --- Subplot 1: Consumo óptimo (todos los agentes) ---
subplot(1,2,1);
hold on;  % importante: permite superposición
for j = 1:n_agents
    plot(a, c_opt{j}(:,1), 'r--', 'LineWidth', 1.2);  % informal
    plot(a, c_opt{j}(:,2), 'b--', 'LineWidth', 1.2); % formal
end
xlabel('Activos (a)');
ylabel('Consumo');
title('Consumo óptimo - todos los RRA informales');
legend({'Informal', 'Formal'}, 'Location', 'SouthEast');
xlim([1 5]);
grid on;

% --- Subplot 2: Distribución de riqueza (todos los agentes) ---
subplot(1,2,2);
hold on;  % permite superposición
for j = 1:n_agents
    g1 = g_opt{j}(:,1);
    g2 = g_opt{j}(:,2);

    % Visualmente eliminamos valores muy pequeños (ruido numérico)
    g1(g1 < 1e-5) = NaN;
    g2(g2 < 1e-5) = NaN;

    plot(a, g1, 'r--', 'LineWidth', 1.2);  % informal
    plot(a, g2, 'b--', 'LineWidth', 1.2); % formal
end
xlabel('Activos (a)');
ylabel('Distribución de riqueza');
title('Distribución g(a) - todos los RRA informales');
legend({'Informal', 'Formal'}, 'Location', 'Northeast');
xlim([0.01 0.6]);
grid on;





figure;

% --- Subplot 1: Optimal consumption (all agents) ---
subplot(1,2,1);
hold on;  % important: allows overlay
for j = 1:n_agents
    plot(a, c_opt{j}(:,1), 'r-', 'LineWidth', 1.2);  % informal
    plot(a, c_opt{j}(:,2), 'b-', 'LineWidth', 1.2); % formal
end
xlabel('Assets (a)');
ylabel('Consumption');
title('Optimal Consumption - RRA 3.15 - 5.30 ');
legend({'Informal', 'Formal'}, 'Location', 'SouthEast');
xlim([1 5]);
grid on;

% --- Subplot 2: Wealth distribution (all agents) ---
subplot(1,2,2);
hold on;  % allows overlay
for j = 1:n_agents
    g1 = g_opt{j}(:,1);
    g2 = g_opt{j}(:,2);

    % Visually remove very small values (numerical noise)
    g1(g1 < 1e-5) = NaN;
    g2(g2 < 1e-5) = NaN;

    plot(a, g1, 'r-', 'LineWidth', 1.2);  % informal
    plot(a, g2, 'b-', 'LineWidth', 1.2); % formal
end
xlabel('Assets (a)');
ylabel('Wealth distribution');
title('Walth Distribution - RRA 3.15 - 5.30');
legend({'Informal', 'Formal'}, 'Location', 'Northeast');
xlim([0.01 0.55]);
%xlim([0.01 4.99]);
grid on;



%% We do not run this

%{
%% Store 
C{jj} = c; %consumption per agent for the value of the state variable "a"
Saving{jj} = adot(:,:,ir); %saving per agent

% (1) Informal | Formal  (DISTRIBUTION)
    % (a) Wealth
    Distribution{jj} = g_r(:,:,ir);  % = g (the optimal): used for "a"
                        % column1: informal population per wealth level
                        % column2: formal population per wealth level

    % (b) Consumption
        % CDistribution{jj} = Cons_distribution; WRONG!
    CDistribution{jj} = g_r(:,:,ir); % = g (the optimal): used for "c"
    
% (2) Total: informal + formal (DISTRIBUTION)
    var_int = g_r(:,:,ir); % getting "g"

        % I do not need this :(. I have them in "Distribution"
        informal{jj} = var_int(:,1); % informal population per wealth level
        formal{jj}   = var_int(:,2); % formal population per wealth level
        
    % (a) TOTAL Wealth
    GDistribution{jj} = var_int(:,1)+ var_int(:,2); %Total ppl per "a" level
    
    % (b) TOTAL Consumption
        % GCdistribution{jj} = GC; %Total consumption distribution WRONG!
    GCdistribution{jj} = var_int(:,1)+ var_int(:,2); %Total consumption distribution

CasesSim(jj,:) = [jj pop1 taxI theta sI sF];
%}