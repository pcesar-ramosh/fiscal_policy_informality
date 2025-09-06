function G = gini_weighted_discrete(x, w)
    x = x(:); w = w(:);
    m = isfinite(x) & isfinite(w) & w > 0;
    x = x(m); w = w(m);
    if isempty(x), G = NaN; return; end
    [x, k] = sort(x);
    w = w(k);

    p = w / sum(w);
    mu = sum(p .* x); if mu<=0, G = NaN; return; end
    cum_p = cumsum(p);
    cum_y = cumsum(p .* x) / mu;
    cum_y_lag = [0; cum_y(1:end-1)];
    G = 1 - sum(p .* (cum_y + cum_y_lag));
end
