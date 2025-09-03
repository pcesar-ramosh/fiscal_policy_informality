function G = gini_weighted(x, w, da)
    if nargin < 3, da = 1; end
    x = x(:); w = w(:);
    p = max(w,0) * da;
    s = sum(p); if s<=0, G = NaN; return; end
    p = p / s;

    [x, k] = sort(x);
    p = p(k);

    mu = sum(p .* x); if mu<=0, G = NaN; return; end
    cum_p = cumsum(p);
    cum_y = cumsum(p .* x) / mu;
    cum_y_lag = [0; cum_y(1:end-1)];
    G = 1 - sum(p .* (cum_y + cum_y_lag));
end
