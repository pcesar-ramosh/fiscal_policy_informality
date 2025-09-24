function out = hugguet_covid_mastransf(params)
% Wrapper para el escenario COVID con distintas transferencias.
% Reutiliza el solver general del modelo base con impuestos/transferencias.
%
% ENTRADA:
%   params : struct con los mismos campos que huggett_base_covid_function
%
% SALIDA:
%   out    : struct con .r, .a, .g, .c, .s, .popI, .popF, .Ctot, .Y
%            .fiscal (Tl, Tc, Tr, G, B, rB, PB, BB)
%            .stats  (momentos incl. gini)
%            .borrowers (fracciones y volúmenes a<0/a>0)

out = huggett_base_covid_function(params);

end
