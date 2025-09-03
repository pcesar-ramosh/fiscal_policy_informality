function q = wquantile_discrete(x, w, p)
    x = x(:); w = w(:);
    m = isfinite(x) & isfinite(w) & (w > 0);
    x = x(m); w = w(m);
    if isempty(x), q = NaN; return; end

    [x, k] = sort(x);
    w = w(k);
    F = cumsum(w) / sum(w);

    [Fu, idx] = unique(F, 'stable');
    xu = x(idx);
    if Fu(1) > 0, Fu = [0; Fu]; xu = [x(1); xu]; end
    if Fu(end) < 1, Fu = [Fu; 1]; xu = [xu; x(end)]; end
    q = interp1(Fu, xu, p, 'linear', 'extrap');
end
