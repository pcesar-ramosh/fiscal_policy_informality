%% stats_Eta_v2 (ADAPTADO AL MODELO BASE)
clear; clc;

RRA = 4.20;
eta_blocks = {[0.2 0.3 0.4], [0.5 0.6 0.64], [0.7 0.8 0.9]};

[r1, ~, pop1_vector1, statsMatrix1, statsCMatrix1A, GDistribution1, a, Distribution1] = ...
    huggett_Equi_INFO_function(RRA, eta_blocks{1});

[r2, ~, pop1_vector2, statsMatrix2, statsCMatrix2,  GDistribution2, a2, Distribution2] = ...
    huggett_Equi_INFO_function(RRA, eta_blocks{2});

[r3, ~, pop1_vector3, statsMatrix3, statsCMatrix3,  GDistribution3, a3, Distribution3] = ...
    huggett_Equi_INFO_function(RRA, eta_blocks{3});

if any(a~=a2) || any(a~=a3), error('Los grids de activos no coinciden.'); end

eta       = [pop1_vector1(:); pop1_vector2(:); pop1_vector3(:)];
r_all     = [r1(:); r2(:); r3(:)];
eta_matrix   = [statsMatrix1; statsMatrix2; statsMatrix3];
statsCMatrix = [statsCMatrix1A; statsCMatrix2; statsCMatrix3];

GDD = [GDistribution1, GDistribution2, GDistribution3];
GDW = [Distribution1,   Distribution2,   Distribution3];

I  = numel(a);
da = (a(end)-a(1))/(I-1);
nE = numel(GDW);

N_borrowers   = zeros(nE,1);
frac_N        = zeros(nE,1);
N_borrowersInf= zeros(nE,1);
N_borrowersFor= zeros(nE,1);
frac_NI       = zeros(nE,1);
frac_NF       = zeros(nE,1);
frac_NIT      = zeros(nE,1);
frac_NFT      = zeros(nE,1);
IaB           = zeros(nE,1);
IaL           = zeros(nE,1);
FaB           = zeros(nE,1);
FaL           = zeros(nE,1);
borrowing     = zeros(nE,1);
lending       = zeros(nE,1);
g1T           = zeros(nE,1);
g2T           = zeros(nE,1);

for i = 1:nE
    Gi = GDD{i};
    isBorrow = (a<0); isLend = (a>0);
    N_borrowers(i) = sum(Gi(isBorrow));
    frac_N(i)      = N_borrowers(i)/sum(Gi);
    borrowing(i)   = sum(Gi(isBorrow).*a(isBorrow))*da;
    lending(i)     = sum(Gi(isLend)  .*a(isLend))  *da;

    gPair = GDW{i};
    N_borrowersInf(i) = sum(gPair(isBorrow,1));
    N_borrowersFor(i) = sum(gPair(isBorrow,2));
    massInf = sum(gPair(:,1)); massFor = sum(gPair(:,2)); massTot = massInf+massFor;
    frac_NI(i)  = N_borrowersInf(i)/max(massInf,eps);
    frac_NF(i)  = N_borrowersFor(i)/max(massFor,eps);
    frac_NIT(i) = N_borrowersInf(i)/max(massTot,eps);
    frac_NFT(i) = N_borrowersFor(i)/max(massTot,eps);

    IaB(i) = sum(gPair(isBorrow,1).*a(isBorrow))*da;
    IaL(i) = sum(gPair(isLend,1)  .*a(isLend))  *da;
    FaB(i) = sum(gPair(isBorrow,2).*a(isBorrow))*da;
    FaL(i) = sum(gPair(isLend,2)  .*a(isLend))  *da;

    g1T(i) = sum(gPair(:,1).*a)*da; 
    g2T(i) = sum(gPair(:,2).*a)*da;
end

totalSupply = IaB + FaB; 
totalDemand = IaL + FaL;
AssetDS     = [IaB, FaB, IaL, FaL, totalSupply, totalDemand];

% Gini total en columna 9 del statsMatrix adaptado
gini_tot = eta_matrix(:,9);

BL = [ eta, N_borrowers, frac_N, frac_NI, frac_NF, ...
       frac_NIT, frac_NFT, borrowing, lending, 100*gini_tot ];

disp('stats_Eta_v2 listo: BL, AssetDS, eta_matrix, statsCMatrix, r_all en workspace.');

% --- Gráficos auxiliares
run g2_asset.m
run g3_stat.m
run g4_borrower.m
