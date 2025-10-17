function out = solve_two_type_huggett_fiscal_tauC(cfg, tau_c)
% Wrapper para analizar el escenario con un τ_c específico.
% Reutiliza tu solver base: solve_two_type_huggett_fiscal_Bfixed(cfg)
%
% Uso:
%   out = solve_two_type_huggett_fiscal_tauC(cfg0, 0.22);

    if nargin<2 || isempty(tau_c)
        error('Debes proporcionar el valor de tau_c (IVA).');
    end
    alt = cfg;
    alt.tau_c = tau_c;
    out = solve_two_type_huggett_fiscal_Bfixed(alt);
end
