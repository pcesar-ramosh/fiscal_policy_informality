function q = wquantile(x, w, p, da)
% Cuantil ponderado continuo (aprox Riemann).
% x: malla ordenada, w: densidad marginal (suma w * da = 1), p en [0,1].
    if nargin < 4, da = 1; end
    w = w / (sum(w)*da);                    % normaliza
    F = cumsum(w) * da;                     % CDF
    q = interp1(F, x, p, 'linear', 'extrap');
end