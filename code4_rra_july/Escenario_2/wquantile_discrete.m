function q = wquantile_discrete(x, w, p)
% Cuantil ponderado para datos discretos (x no necesariamente ordenado).
    [xs, idx] = sort(x(:));
    ws = w(:); ws = ws(idx);
    ws = ws / sum(ws);
    F  = cumsum(ws);
    q  = interp1(F, xs, p, 'linear', 'extrap');
end