%% =============================================================
%  huggett_Equi_INFO_Jan24.m  (ADAPTADO AL MODELO BASE)
%  Corre el solver para eta = [0.20, 0.64] y deja variables
%  en el workspace con los nombres que usa tu main de gráficos.
% =============================================================
clearvars -except -regexp ^(runningFromMain)$
rng(123);

RRA = 4.20;             % RRA base (puedes cambiar)
eta_vector = [0.20, 0.64];

[r, ir, pop1_vector, statsMatrix, statsCMatrix, ...
 GDistribution, a, Distribution, informal, formal, ...
 C, Saving, CDistribution, GCdistribution, I, amin] = ...
    huggett_Equi_INFO_function(RRA, eta_vector);

% Nota: las variables quedan en el workspace y tu "main" de gráficos
% las utiliza directamente (C, Saving, Distribution, GDistribution, etc.).
