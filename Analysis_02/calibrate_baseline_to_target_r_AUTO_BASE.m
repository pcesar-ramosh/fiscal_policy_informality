function [cfg_out, phi_star, S_at_target] = calibrate_baseline_to_target_r_AUTO_BASE( ...
    share_I, sI_vec, sF_vec, trans, cfg_in, auto_opts, tol_phi)
% Busca la calibración "mínima" que permita S_total(r_target)=0 ajustando:
%   1) a_min (más laxo), 2) theta (más bajo), 3) tau_c (más bajo),
%   4) tau_l (más alto), 5) ampliar phi_max; y luego calibra phi por bisección.
%
% Devuelve cfg_out (con parámetros elegidos + phi*), phi_star y S_at_target.

if nargin < 7 || isempty(tol_phi), tol_phi = 1e-5; end
cfg_out = cfg_in;

% Listas por defecto si no se pasan
if ~isfield(auto_opts,'amin_list'),  auto_opts.amin_list  = [cfg_in.amin_abs, -2.0, -2.5]; end
if ~isfield(auto_opts,'theta_list'), auto_opts.theta_list = [cfg_in.theta, 0.005, 0.002]; end
if ~isfield(auto_opts,'tauc_list'),  auto_opts.tauc_list  = [cfg_in.tau_c, 0.04, 0.00]; end
if ~isfield(auto_opts,'taul_list'),  auto_opts.taul_list  = [cfg_in.tau_l, 0.20]; end
if ~isfield(auto_opts,'phimax_list'),auto_opts.phimax_list= [0.60, 0.80, 1.00]; end

% Intento en orden lexicográfico de “menor intervención”
for a_min = auto_opts.amin_list
    for th = auto_opts.theta_list
        for tc = auto_opts.tauc_list
            for tl = auto_opts.taul_list
                for phimax = auto_opts.phimax_list

                    cfg_try = cfg_in;
                    cfg_try.amin_abs = a_min;
                    cfg_try.theta    = th;
                    cfg_try.tau_c    = tc;
                    cfg_try.tau_l    = tl;

                    try
                        [phi_star, S_at_target] = calibrate_phi_for_target_r_BASE( ...
                            share_I, sI_vec, sF_vec, trans, cfg_try, 0.0, phimax, tol_phi);

                        % ¡Éxito!
                        cfg_try.phi = phi_star;
                        cfg_out = cfg_try;
                        return

                    catch
                        % No bracketing aún: sigue probando
                        continue
                    end
                end
            end
        end
    end
end

% Si llegamos aquí, ningún combo produjo sign-change para phi ∈ [0, max(phimax_list)]
% Elegimos el “mejor” combo (mínimo |S| en el extremo superior) y avisamos.
best = struct('absS',inf,'cfg',cfg_in,'phi',NaN,'S',NaN);
for a_min = auto_opts.amin_list
    for th = auto_opts.theta_list
        for tc = auto_opts.tauc_list
            for tl = auto_opts.taul_list
                phimax = auto_opts.phimax_list(end);
                cfg_try = cfg_in;
                cfg_try.amin_abs = a_min; cfg_try.theta = th;
                cfg_try.tau_c = tc;       cfg_try.tau_l = tl;
                cfg_try.phi   = phimax;

                [~, S_end] = huggett_S_given_r_BASE(share_I, sI_vec, sF_vec, trans, cfg_try, cfg_try.r_target);
                if abs(S_end) < best.absS
                    best.absS = abs(S_end);
                    best.cfg  = cfg_try;
                    best.phi  = phimax;
                    best.S    = S_end;
                end
            end
        end
    end
end

error(['AUTO-calibración no logró bracketing ni con combos ampliados.\n' ...
       'Mejor intento: amin=%.2f, theta=%.3f, tau_c=%.3f, tau_l=%.3f, phi=%.2f -> S=%.3e\n' ...
       'Sugerencias: permitir aún más deuda (a_min más negativo), reducir theta más,\n' ...
       'o bajar tau_c; también puedes mover (z1,z2) o ρ si fuera necesario.'], ...
       best.cfg.amin_abs, best.cfg.theta, best.cfg.tau_c, best.cfg.tau_l, best.phi, best.S);
end
