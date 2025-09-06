function q = wquantile(x, w, p, da)
    if nargin < 4, da = 1; end
    x = x(:); w = w(:);
    w(~isfinite(w)) = 0; w = max(w,0);

    F = cumsum(w) * da;
    tot = F(end);
    if tot <= 0
        q = NaN; return
    end
    F = F / tot;

    [Fu, idx] = unique(F, 'stable');
    xu = x(idx);
    if Fu(1) > 0, Fu = [0; Fu]; xu = [x(1); xu]; end
    if Fu(end) < 1, Fu = [Fu; 1]; xu = [xu; x(end)]; end
    q = interp1(Fu, xu, p, 'linear', 'extrap');
end