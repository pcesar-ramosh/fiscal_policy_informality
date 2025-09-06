function G = gini_weighted_discrete(x, w)
% Gini ponderado para par (x,w) discreto (x>=0).
    x = x(:); w = w(:);
    w = w / sum(w);
    x = max(x, 0);
    mu = sum(x .* w);
    if mu <= 0
        G = 0; return
    end
    % Ordenar
    [xs, idx] = sort(x);
    ws = w(idx);
    % Curva de Lorenz
    F = cumsum(ws);
    S = cumsum(xs .* ws);
    L = S / mu;
    areaL = trapz(F, L);
    G = 1 - 2*areaL;
    G = min(max(G,0),1);
end