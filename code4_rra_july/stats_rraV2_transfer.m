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
n_agents = 20; % Inicial estaba con 9. also we add this in for structure
% Definición del rango para la aversión al riesgo relativa (RRA)
s_min = 0.15;  % 0.15;  With this values the code works  + lambda of function
s_max = 0.30; % 0.30;
% Vectores automáticos para trabajadores informales y formales (heterogéneos)
eta_vector   = 0.75 * ones(1, n_agents);   % 0.65             % Mismo tamaño que los tipos de agentes
sI_vector1   = linspace(s_min, s_max, n_agents);        % RRA heterogénea informal
% sF_vector1   = linspace(s_min, s_max, n_agents);        % RRA heterogénea formal
% Alternativa: RRA constante para formales (descomenta si se desea usar)
sF_vector1 = 0.30 * ones(1, n_agents);       % Initial 0.15        % RRA constante formal
% Definir un set de parametros iniciales, 1 de esos en teroi definir
% calibrar los parametros, valores iniciales, del modelo base comparen los
% resultados de wealth, cons, policy function, cuando el RRA es mayor a 1
% (enteros) comparens es en la RRA, con cual de los dos nos quedamos en
% escenario base? y con eso poodemoms agregar el escenerario

%% Simulation (Solves the General Equilibrium Model for each input parameter combination)
%[r, ir, pop1_vector, statsMatrix1, statsCMatrix1, GDistribution1, a, Distribution1] = ...
%   huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1);

[r_opt, ir, pop1_vector, a, g_opt, c_opt] =...
    huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1);


% [r_opt, ir, pop1_vector, a, g_opt, c_opt, statsCMatrix1] = ...
%    huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1);

    % r	Vector con tasas de interés de equilibrio para cada RRA.
    % ir	Iteración interna del ajuste de r.
    % pop1_vector	Tamaño relativo de la población informal para cada simulación.
    % statsMatrix1	Estadísticas de distribución de riqueza (medias, medianas, Gini, etc.).
    % statsCMatrix1	Estadísticas de distribución del consumo.
    % GDistribution1	Distribución total de población por nivel de activos (a).
    % a	Vector de niveles de activos (grid).
    % Distribution1	Matriz de densidad poblacional informal/formal por a (g1 y g2).

%% a. Extract what we need
    % Extract: vector of parameter values
    eta = sI_vector1'; % the X axis is "Informal RRA"   Extrae vector columna de aversión informal.
    % Extract: moments of "wealth distribution" 
    eta_matrix = statsMatrix1; % contiene estadísticas  Extrae la matriz de estadísticas de riqueza.

%% b. Borrowing / Lending
% Cálculo de prestatarios y clasificación por activos
% Para cada nivel de RRA (9 simulaciones):
% total people per level of "a"
GDD = GDistribution1; % every column if for \tax1 value: 0 ---> 0.18

% informal / formal people per level of "a"
GDW = Distribution1; % every column if for \eta value: 0.2 ---> 0.9
                                                   % g1 and g2
da = 0.0000001; % 0.005104104104104;     % tamaño del paso de activos ??? Como se obtuvo?
% vector de la irqueza miniam a una a, un paso constante, grid, deberiamos
% de calcularlo, endogenamente 0.03 max 5, con linspace
N_borrowers = [];
frac_N = [];
borrowing = [];
lending =[];
g1T = [];
g2T = [];

% calcutating fraction of borrowers
    for i=1:n_agents % for i=1:19 % 
%     for i=1:19 % for i=1:9 % &&&&&%%-------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % - What is the fraction of the population that are net borrowers?
    % answer (for pop1_vector{1} = eta(i)):
        GA = [GDD{i} a]; %step1 distribución total (informal + formal) para RRA i
        GA(:,3) = GA(:,2) < 0;     %Restriccion? Podemos Modificarla?
        % % Marcamos con 1 si el nivel de activo es negativo (a<0)
        N_borrowers(i) = sum(GA(:,1).*GA(:,3));
        frac_N(i) = N_borrowers(i)/sum(GDD{1}); % this is the fraction of total
        % population that are net borrowers (a < 0) 
        % Same total pop across eta: 
            % 195.92 = sum(GDD{1})=sum(GDD{2})=...sum(GDD{9})
    
    % - What is the fraction of the informal/formal pop..
    %   that are net borrowers?
        GAIF = [GDW{i} a]; %step1: 
                           % column1(informal pop)
                           % column2(formal pop)
                           % column3(asset level)
        GAIF(:,4) = GAIF(:,3) < 0;     %indicator function: 1 (a<0) (prestatarios)
        GAIF(:,5) = GAIF(:,3) > 0;     %indicator function: 1 (a>0) (ahorradores)
        
        N_borrowersInf(i) = sum(GAIF(:,1).*GAIF(:,4)); % Informales con a<0
        N_borrowersFor(i) = sum(GAIF(:,2).*GAIF(:,4)); % Formales con a<0

%-------------------
% ppl: Informal    
    N_abarInf(i) = GAIF(1,1); % at a_bar  Informales (activo mínimo)
    N_anInf(i) = N_borrowersInf(i) - N_abarInf(i); %  a_bar < a <0 % Informales con a entre a_bar y 0
    N_apInf(i) = sum(GAIF(:,1).*GAIF(:,5)); % a>0  Informales con a > 0 (ahorradores)
    Total_info(i) = N_abarInf(i) + N_anInf(i) + N_apInf(i); % Total de informales
        %check
        Total_info_GAIF(i) = sum(GAIF(:,1)); % it should be = to "Total_info(i)" % Chequeo: debe coincidir

% ppl: Formal
    N_abarFor(i) = GAIF(1,2); % at a_bar
    N_anFor(i) = N_borrowersFor(i) - N_abarFor(i); %  a_bar < a <0
    N_apFor(i) = sum(GAIF(:,2).*GAIF(:,5)); % a>0
    Total_for(i) = N_abarFor(i) + N_anFor(i) + N_apFor(i);
        %check
        Total_for_GAIF(i) = sum(GAIF(:,2)); % it should be = to "Total_for(i)"        
% People: #  Guardar resumen por grupo
    people{i} = [N_abarInf(i), N_anInf(i), N_apInf(i), Total_info(i), Total_info_GAIF(i);...
        N_abarFor(i), N_anFor(i), N_apFor(i), Total_for(i), Total_for_GAIF(i)];
%-------------------
        %### Fracción de prestatarios dentro de su grupo y del total
        % this is the fraction of total of type
        frac_NI(i) = N_borrowersInf(i)/sum(GAIF(:,1)); % dentro de informales
        frac_NF(i) = N_borrowersFor(i)/sum(GAIF(:,2)); % dentro de formales
        % population that are net borrowers (a < 0)                                        
    
        % this is the fraction of Total pop
        frac_NIT(i) = N_borrowersInf(i)/sum(sum(GAIF(:,1:2))); % sobre el total
        frac_NFT(i) = N_borrowersFor(i)/sum(sum(GAIF(:,1:2))); % sobre el total
        % population that are net borrowers (a < 0)                                        

    % Asset Supply/Demand
        % Informal: asset supply (borrower):a<0
        IaB(i) = sum(GAIF(:,1).*GAIF(:,4).*a*da); % Informales prestatarios (a<0)
        % Informal: asset demand (lender):a>0
        IaL(i) = sum(GAIF(:,1).*GAIF(:,5).*a*da); % Informales prestamistas (a>0)
    
        % Formal: asset supply (borrower):a<0
        FaB(i) = sum(GAIF(:,2).*GAIF(:,4).*a*da);  % Formales prestatarios
        % Formal: asset demand (lender):a>0
        FaL(i) = sum(GAIF(:,2).*GAIF(:,5).*a*da); % Formales prestamistas
    
    % - check: total borrowing = -total lending
    % For for pop1_vector{1} = eta(i):
        GA(:,4) = GA(:,1).*GA(:,2); % total pop per "a" level % población total por cada "a"
        borrowing(i) = sum(GA(:,4).*GA(:,3)); % Total prestado (a<0)
        GA(:,5) = GA(:,2) > 0;     %indicator function: 1 (a>0) % Indicador de a > 0
        lending(i)   = sum(GA(:,4).*GA(:,5));  % Total ahorrado (a > 0)
        % - borrowing =~ lending (In equilibrium)
    
    % - inspection of g1 and g2
        g1 = GDW{i}(:,1); % distribución informal
        g2 = GDW{i}(:,2); % distribución formal
        g1T(i) = g1'*a*da; % Media de activos para informales
        g2T(i) = g2'*a*da; % Media de activos para formales
    end

% Storing in one matrix  ( guarda los principales indicadores de cada simulación.)
BL = [eta, N_borrowers', frac_N', frac_NI', frac_NF', frac_NIT',...
      frac_NFT',borrowing', lending', 100*eta_matrix(:,15)];

    % eta	Niveles de aversión al riesgo informal (sI_vector1').
    % N_borrowers'	Número total de prestatarios netos (a<0).
    % frac_N'	Fracción del total de población que son prestatarios.
    % frac_NI'	Fracción de informales que son prestatarios.
    % frac_NF'	Fracción de formales que son prestatarios.
    % frac_NIT'	Fracción de prestatarios informales sobre el total poblacional.
    % frac_NFT'	Fracción de prestatarios formales sobre el total poblacional.
    % borrowing'	Total de activos negativos (endeudamiento agregado).
    % lending'	Total de activos positivos (ahorro agregado).
    % 100*eta_matrix(:,15)	Índice de Gini multiplicado por 100 (columna 15 contiene el Gini).

%% Assets: demand / supply (to explain "r")
totalSupply = IaB' + FaB'; % Oferta agregada de activos (deuda)
totalDemand = IaL' + FaL'; % Demanda agregada de activos (ahorro)
    % Suma el total de activos de los agentes prestatarios (a < 0) → oferta
    % Suma el total de activos de los agentes ahorradores (a > 0) → demanda
AssetDS = [IaB' FaB' IaL' FaL' totalSupply totalDemand]; 

    % IaB'	Activos negativos de informales (deuda).
    % FaB'	Activos negativos de formales (deuda).
    % IaL'	Activos positivos de informales (ahorro).
    % FaL'	Activos positivos de formales (ahorro).
    % totalSupply	Oferta agregada de activos (endeudamiento).
    % totalDemand	Demanda agregada de activos (ahorro).

%% Graphs    
    run g2_asset.m      % assets(supply / demand) + r
    run g3_stat.m       % stats: median/Std: total/informal/formal
    run g4_borrower.m   % Population (borrowers): frac_N,...  

    % g2_asset.m	Gráfica de oferta vs demanda de activos y su relación con r.
    % g3_stat.m	Estadísticas de distribución (media, mediana, std) por grupo.
    % g4_borrower.m	Gráficas de fracciones de prestatarios por RRA o tipo.
    
