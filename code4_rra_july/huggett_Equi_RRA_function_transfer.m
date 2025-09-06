%function [r, ir, pop1_vector, statsMatrix, statsCMatrix, GDistribution, a, Distribution] =...
%         huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1)

function [r_opt, ir, pop1_vector, a, g_opt, c_opt] =...
         huggett_Equi_RRA_function_transfer(eta_vector,sI_vector1,sF_vector1)

    % eta_vector: vector con proporciones del sector informal
    % sI_vector1: aversión al riesgo relativa (RRA) para informales.
    % sF_vector1: RRA para formales.
    
    % r: tasa de interés de equilibrio.
    % ir: número de iteración donde se alcanza el equilibrio.
    % pop1_vector: copia de eta_vector.
    % statsMatrix, statsCMatrix: estadísticas de riqueza y consumo.
    % GDistribution: distribución total de activos.
    % a: grilla de activos.
    % Distribution: densidad de población por grupo (informal y formal).

%s = RRA; % same for formal and informal (converge: 0.25, 0.3)
 %% Model
%clear all; clc; close all;

tic;

pop1_vector = eta_vector;   %Informal population

for jj=1:size(pop1_vector,2)
    pop1 = pop1_vector(jj);

% Lambda & Jump Probability
p22 = 0.75; %0.75; %0.20 oK % 0.8155; % 0.8155;   0.62=0.83 %Prob[y2|y2] = exp(-lambda_2) | From lognormal.m (Luis Y.) Sensible
la2 = -log(p22); %intensity Formal==>Informal
    %From calibration, but we need that p11 changes with the "informality size"
    %p11 = 0.8944;    %Prob[y1|y1] = exp(-lambda_1) | From lognormal.m (Luis Y.)
    %la1 = -log(p11); %intensity Formal==>Informal

% Tasa de salto para informales (de informal a formal).
la1 =(1-pop1)*la2/pop1; %Prob[y1|y1] = exp(-lambda_1)
p11 = exp(-la1); % with pop1 = 0.64, we got p11 (0.8916) close to the calibrated value
% Probabilidad de que un agente informal permanezca en el estado informal.
la = [la1,la2]; %%% Sensibilad REVISAR POR PARTES

% Tax rate
taxF = 0.10; %0.17; %0.18; 0.15;
taxI = 0.00;
%%%phi = 0.05; %%% Transfers

% Premium
%theta_vector = [0 0.02 0.02]; %informal pays 1.5 times the "r" paid by formal
                              % theta = 20% (calibration)

% RRA (s) and discount rate (rho)                              
sI_vector    = sI_vector1; %RRA informal
sF_vector    = sF_vector1; %RRA formal

% Prima de riesgo: adicional que pagan los informales al endeudarse.
    %theta = theta_vector(jj);
    theta = 0.02; % risk premium  % prima de riesgo por endeudamiento informal
    sI = sI_vector(jj);
    sF = sF_vector(jj);

rho = 0.05; % tasa de descuento

% Income
z1 = 0.33; %informal % ingreso informal
z2 = 1;    %formal % ingreso formal
z = [z1,z2];



%% Discretization

% Los agentes eligen cuánto ahorrar o endeudarse en cada periodo. 
% Como no se puede resolver eso analíticamente se discretiza el espacio de decisiones en una grilla de activos, 
% que es simplemente un vector que contiene distintos posibles niveles de riqueza (o deuda).

r0 = 0.03;  % Guess Interest Rate
rmin = 0.01; % -0.01; % 
rmax = 0.04;

I= 700; %1340; %1000; % cuántos puntos va a tener tu grilla de activos.
amin = -0.3*z1; % 30% of the min income   límite inferior de la grilla de activos
amax = 5; % límite superior de la riqueza
a = linspace(amin,amax,I)'; % grilla de activos a
da = (amax-amin)/(I-1); % es el tamaño del paso entre cada punto de la grilla.

aa = [a,a];
zz = ones(I,1)*z;


maxit= 100; %0; %100;  Número máximo de iteraciones permitidas para resolver 
                    % la ecuación de Bellman (función de valor) por el método iterativo.
crit = 10^(-6); % Criterio de convergencia: cuando el cambio entre iteraciones de la 
                    % función de valor (V) es menor que 1e-6, se considera que la solución ha convergido.
Delta = 1000; % 10000; %1000; discretización temporal implícita del sistema de ecuaciones de HJB.

%  Inicialización de variables para derivadas y consumo:
dVf = zeros(I,2);  % matrices de derivadas hacia adelante (forward)
dVb = zeros(I,2);  % matrices de derivadas hacia atrás (backward)
c = zeros(I,2); % Inicializa la matriz de consumo c(a) para ambos agentes.

% Matriz de transición entre sectores (informal/formal):
Aswitch = [-speye(I)*la(1),speye(I)*la(1);
                  speye(I)*la(2),-speye(I)*la(2)];
% Esta es una matriz 2I × 2I que representa las transiciones entre estados ocupacionales:
% la(1) es la intensidad de cambio de informal a formal.
% la(2) es la intensidad de cambio de formal a informal.

Ir = 100; %40; %Number of simulations (to find "r")
% Número máximo de iteraciones para encontrar la tasa de interés de equilibrio r, 
% mercado de activos esté en equilibrio (oferta = demanda agregada de ahorro).

crit_S = 10^(-5); % Criterio de convergencia del exceso de ahorro o demanda neta de activos.


%%%%%%%%%%%%%%%%%%%%%%%%%%
% HJB EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL GUESS
r = r0;
%Gov = taxI*z(1)*(la2/(la1+la2)) + taxF*z(2)*(la1/(la1+la2));
phi = 0.13; %0.12; %0.13; %0.06; %0.06; %0.06; %0.04;  %0.05;                % proporción del ingreso informal que se transfiere
Transfer = phi * z1;         % transferencia fija para cada informal

popI = la2 / (la1 + la2);    % proporción informal en estado estacionario
popF = la1 / (la1 + la2);    % proporción formal en estado estacionario

% Recaudación total esperada
TaxRevenue = taxI*z1*popI + taxF*z2*popF;

% Gasto público común para todos (residual)
Gov_pub = TaxRevenue - Transfer * popI;

% Gasto total disponible por agente (consumido como bien público)
Gov = Gov_pub + Transfer;  % Correcto: bien público + transfer focalizada


    % Interest rate is different for formal and informal agent (in the borrowing side: a<0)
    % rr: vector interest rate for informal over "a"
    for i=1:size(a,1)
        if a(i)>=0
            rr(i,1) = r;         %informal: lending (=formal)
        else
            rr(i,1) = r + theta; %informal: borrowing
        end
    end
    check = [a rr]; % I expect "r+theta" for a<0

%v0(:,1) = (((1-taxI)*zz(:,1) + rr.*a).^(1-sI)/(1-sI) + Gov)/rho;
%v0(:,1) = (((1-taxI)*zz(:,1) + rr.*a + phi * z1).^(1-sI)/(1-sI)) / rho;
%v0(:,1) = (((1-taxI)*zz(:,1) + rr.*a + Transfer).^(1-sI)/(1-sI) + Gov)/rho;
v0(:,1) = (((1-taxI)*zz(:,1) + rr.*a + Transfer).^(1-sI)/(1-sI) + Gov)/rho;


%v0(:,2) = (((1-taxF)*zz(:,2) + r.*a).^(1-sF)/(1-sF)  + Gov)/rho;
%v0(:,2) = (((1-taxF)*zz(:,2) + r.*a).^(1-sF)/(1-sF) )/rho;
%v0(:,2) = (((1-taxF)*zz(:,2) + r.*a).^(1-sF)/(1-sF) + Gov)/rho;
v0(:,2) = (((1-taxF)*zz(:,2) + r.*a).^(1-sF)/(1-sF) + Gov)/rho;

for ir=1:Ir

r_r(ir)=r;
rmin_r(ir)=rmin;
rmax_r(ir)=rmax;
    
if ir>1
v0 = V_r(:,:,ir-1);
end

v = v0;

for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    
    %% (I) ----Forward/Backward Diff Approx----
    % forward difference    
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,1) = ((1-taxI)*z(1) + r.*amax).^(-sI); %will never be used, but impose state constraint a<=amax just in case
    dVf(I,2) = ((1-taxF)*z(2) + r.*amax).^(-sF); %will never be used, but impose state constraint a<=amax just in case
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,1) = ((1-taxI)*z(1) + (r+theta).*amin).^(-sI); %state constraint boundary condition
    dVb(1,2) = ((1-taxF)*z(2) + r.*amin).^(-sF); %state constraint boundary condition
        
    I_concave = dVb > dVf; %indicator whether value function is concave (problems arise if this is not the case)

        % rr: vector interest rate for informal over "a"
        for i=1:size(a,1)
            if a(i)>=0
                rr(i,1) = r;         %informal: lending (=formal)
            else
                rr(i,1) = r + theta; %informal: borrowing
            end
        end
    %consumption and savings with forward difference
    cf = [max(dVf(:,1),10^(-10)).^(-1/sI) max(dVf(:,2),10^(-10)).^(-1/sF)];
    ssf = [(1-taxI)*zz(:,1) (1-taxF)*zz(:,2)] + [rr r*ones(size(aa,1),1)].*aa - cf;
    %consumption and savings with backward difference
    cb = [max(dVb(:,1),10^(-10)).^(-1/sI) max(dVb(:,2),10^(-10)).^(-1/sF)];
    ssb = [(1-taxI)*zz(:,1) (1-taxF)*zz(:,2)] + [rr r*ones(size(aa,1),1)].*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = [(1-taxI)*zz(:,1) (1-taxF)*zz(:,2)] + [rr r*ones(size(aa,1),1)].*aa;
    
    % upwind method makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    
    c = cf.*If + cb.*Ib + c0.*I0;
    %u = [c(:,1).^(1-sI)/(1-sI) c(:,2).^(1-sF)/(1-sF)] + Gov.*ones(I,2);
    %u = [((c(:,1) + phi * z1).^(1 - sI)) / (1 - sI), (c(:,2).^(1 - sF)) / (1 - sF)];
    %u = [((c(:,1) + Transfer).^(1 - sI)) / (1 - sI), (c(:,2).^(1 - sF)) / (1 - sF)] + Gov.*ones(I,2);
    u = [((c(:,1) + Transfer).^(1 - sI)) / (1 - sI), (c(:,2).^(1 - sF)) / (1 - sF)] + Gov.*ones(I,2);


    %CONSTRUCT MATRIX
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    B = (1/Delta + rho)*speye(2*I) - A;

    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOKKER-PLANCK EQUATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
AT = A';
b = zeros(2*I,1);

%need to fix one value, otherwise matrix is singular
i_fix = 1;  %row 1
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%Solve linear system
gg = AT\b;

g_sum = gg'*ones(2*I,1)*da; %total sum = 1: informal + formal
gg = gg./g_sum;             %fraction wrt total

g = [gg(1:I),gg(I+1:2*I)]; %column1 (informal), column2 (formal)
                           %sum Row: total distribution per level of "a" 
                           %this is "wealth distribution" per level of "a"

check1 = g(:,1)'*ones(I,1)*da; % = pop1_vector (informal pop)  
check2 = g(:,2)'*ones(I,1)*da; % = 1 - pop1_vector (formal pop)  

g_r(:,:,ir) = g;
adot(:,:,ir) = [(1-taxI)*zz(:,1) (1-taxF)*zz(:,2)] + [rr r*ones(size(aa,1),1)].*aa - c;
V_r(:,:,ir) = V;

S(ir) = g(:,1)'*a*da + g(:,2)'*a*da;

%UPDATE INTEREST RATE
if S(ir)>crit_S
    disp('Excess Supply')
    rmax = r;
    r = 0.5*(r+rmin);
        disp(r)
elseif S(ir)<-crit_S;
    disp('Excess Demand')
    rmin = r;
    r = 0.5*(r+rmax);
elseif abs(S(ir))<crit_S;
    display('Equilibrium Found, Interest rate =')
    disp(r)
    break
end

%% Storage
r_opt(jj) = r;
g_opt{jj} = g;
c_opt{jj} = c;
end




%{

%% Total wealth distribution (for every level of "a")
G = g(:,1) + g(:,2); % # of TOTAL people per "a" level
                     % it is not expressed in "fraction" of total people
check3 = sum(G)*da; % = 1

%% People at constraint

pp_fraction = g/sum(sum(g)); % pp/TotalPop % OBS: is the number of persons in each asset state and work state (formal/informal)
pp_fraction_const = pp_fraction(1,:); % percentage of total people at the constraint

%% Distribution Statistics (wealth)
% wealth mean for: Inf | For
stats.gmean = sum(a.*g)./sum(g); 
stats.gsd = [sqrt(var(a,g(:,1))) sqrt(var(a,g(:,2)))]; %SD
    D = [a a];%a.*g;
    W = g./sum(g);

WW{jj} = W; % weights   

    stats.gmedian = [weightedMedian(D(:,1),W(:,1))...
                     weightedMedian(D(:,2),W(:,2))];

%skewness is a measure of the asymmetry of the probability distribution 
%negative skew commonly indicates that the tail is on the left side of the distribution
stats.gskew = skewness(g);        % skew Normal Distr = 0 (Symmetry)

%A positive value tells you that you have heavy-tails (i.e. a lot of data in your tails).
stats.gkurt = kurtosis(g);        % kurt Standard normal Dis = 3

% G stats:
stats.Gmean = sum(a.*G)./sum(G); 
stats.Gsd = sqrt(var(a,G)); %SD
    D = a;%a.*G;
    W1G = G./sum(G);
stats.Gmedian = [weightedMedian(D,W1G)];

% Gini
%gini = ginicoeff(G/sum(G), D/sum(D));

% Stats (Wealth)
statsMatrix(jj,:) = [stats.gmean stats.gsd stats.gmedian stats.gskew stats.gkurt... % informal | formal
                     stats.Gmean stats.Gsd stats.Gmedian... % Total
                     p11 r pp_fraction_const]; % more info
 
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
statsC.gmean = sum(c.*g)./sum(g); 
statsC.gsd = [sqrt(var(c(:,1),g(:,1))) sqrt(var(c(:,2),g(:,2)))]; %SD
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
gt = [g(:,1);g(:,2)];

statsC.Gmean = sum(gt.*ct)/sum(gt);%sum(G.*GC)./sum(G); 
statsC.Gsd = sqrt(var(ct,gt)); %sqrt(var(GC,G)); %SD
    D1 = ct; %GC;%GC.*G;
    W1 = gt./sum(gt); %G./sum(G);
statsC.Gmedian = [weightedMedian(D1,W1)];

% Stats
statsCMatrix(jj,:) = [statsC.gmean statsC.gsd statsC.gmedian...
                      statsC.gskew statsC.gkurt...
                      statsC.Gmean statsC.Gsd statsC.Gmedian];

disp(statsC)

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
end

