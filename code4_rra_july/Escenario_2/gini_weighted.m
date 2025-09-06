function G = gini_weighted(x, w, da)
% Gini ponderado para variable continua x (malla ordenada) y densidad w.
% Fórmula basada en área de la Curva de Lorenz.
    if nargin < 3, da = 1; end
    % Ordenar por x, asegurar no-negatividad
    [x, idx] = sort(x(:));
    w = w(:); w = w(idx);
    w = w / (sum(w)*da);

    % Ingresos (o riqueza) no negativos
    x = max(x, 0);

    mu = sum(x .* w) * da;
    if mu <= 0
        G = 0; return
    end

    % Curva de Lorenz (trapezoidal)
    F = cumsum(w) * da;                 % proporción población
    S = cumsum(x .* w) * da;            % masa de recursos
    L = S / (mu);                       % Lorenz

    % Área bajo Lorenz por trapecios en F
    areaL = trapz(F, L);
    G = 1 - 2*areaL;
    % Acotar por robustez numérica
    G = min(max(G,0),1);
end