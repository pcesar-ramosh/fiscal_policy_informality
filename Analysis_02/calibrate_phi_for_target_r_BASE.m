function [phi_star, S_at_target, dbg] = calibrate_phi_for_target_r_BASE( ...
    share_I, sI_vec, sF_vec, trans, cfg, phi_lo, phi_hi, tol)

% Calibra phi tal que S_total(r_target; phi)=0
% Estrategia:
%  (1) Barrido en malla [phi_lo, phi_max] para localizar cambio de signo
%  (2) Bisección en el sub-intervalo que cambia de signo
%
% Salidas:
%  - phi_star: solución
%  - S_at_target: S_total en la solución (≈0)
%  - dbg: struct con rejilla y valores S(phi) para diagnóstico

if nargin < 7 || isempty(tol), tol = 1e-5; end
if nargin < 6 || isempty(phi_hi), phi_hi = 0.08; end
if nargin < 5, error('cfg requerido'); end

rT = cfg.r_target;
cfg_try = cfg;

% ---------- (1) BARRIDO EN MALLA AMPLIO ----------
phi_max = max(phi_hi, 0.6);     % rango amplio por defecto
Ngrid   = 41;                   % ~cada 0.015 si phi_max=0.6
phi_grid = linspace(max(0,phi_lo), phi_max, Ngrid);
S_grid   = nan(Ngrid,1);

for i=1:Ngrid
    cfg_try.phi = phi_grid(i);
    [~, S_grid(i)] = huggett_S_given_r_BASE(share_I, sI_vec, sF_vec, trans, cfg_try, rT);
end

dbg = struct('phi_grid',phi_grid(:), 'S_grid',S_grid(:));

% Buscar primer par contiguo con cambio de signo
idx = find(S_grid(1:end-1).*S_grid(2:end) <= 0, 1, 'first');

if isempty(idx)
    % No hay cambio de signo en todo el rango. Explicar y abortar.
    msg = sprintf(['No se encontró cambio de signo de S_total(r=%.4f) en phi ∈ [%.3f, %.3f].\n' ...
                   'Valores extremos: S(%.3f)=%.3e, S(%.3f)=%.3e.\n' ...
                   'Interpretación: con su calibración actual (tau_l=%.3f, tau_c=%.3f, theta=%.3f, a_min=%.2f),\n' ...
                   'a r=%.4f la economía quiere %s ahorro neto para todo el rango de phi probado.\n\n' ...
                   'Sugerencias para habilitar el bracketing SIN cambiar la estructura del modelo:\n' ...
                   '  • Ampliar el rango de phi (p.ej. hasta 0.8 o 1.0)\n' ...
                   '  • Hacer a_min más laxo (más negativo), p.ej. -2.0\n' ...
                   '  • Reducir theta (spread), p.ej. 0.005 → 0.002\n' ...
                   '  • Reducir tau_c (IVA) si es alto, p.ej. 0.08 → 0.04 o 0.00\n' ...
                   '  • Aumentar tau_l formal si se desea reducir recursos formales (0.15 → 0.20)\n'], ...
                   rT, phi_grid(1), phi_grid(end), phi_grid(1), S_grid(1), phi_grid(end), S_grid(end), ...
                   cfg.tau_l, cfg.tau_c, cfg.theta, cfg.amin_abs, rT, ternary_str(all(S_grid>0),'EXCESO POSITIVO (oferta)','EXCESO NEGATIVO (demanda)'));
    error(msg);
end

% ---------- (2) BISECCIÓN EN EL SUB-INTERVALO ----------
a = phi_grid(idx); b = phi_grid(idx+1);
Sa = S_grid(idx);  Sb = S_grid(idx+1);

for it=1:60
    m  = 0.5*(a+b);
    cfg_try.phi = m;
    [~, Sm] = huggett_S_given_r_BASE(share_I, sI_vec, sF_vec, trans, cfg_try, rT);
    if abs(Sm) < tol
        phi_star   = m;
        S_at_target= Sm;
        return
    end
    if sign(Sm) == sign(Sa)
        a = m; Sa = Sm;
    else
        b = m; Sb = Sm;
    end
end

phi_star    = 0.5*(a+b);
cfg_try.phi = phi_star;
[~, S_at_target] = huggett_S_given_r_BASE(share_I, sI_vec, sF_vec, trans, cfg_try, rT);

end

% --------- util chiquita ---------
function s = ternary_str(tf, s1, s2)
if tf, s = s1; else, s = s2; end
end
