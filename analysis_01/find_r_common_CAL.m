function [r_star_common, history] = find_r_common_CAL(eta_vector, sI_vector1, sF_vector1, cfg, rL, rU)
% Búsqueda por bisección de r común que vacía el mercado agregado:
% S_total(r*) = 0

if nargin < 5, rL = -0.01; end
if nargin < 6, rU =  0.06; end

history = struct('r',[],'S',[]);
for it=1:50
    rM = 0.5*(rL + rU);
    [~, S_tot] = huggett_S_given_r_CAL(eta_vector, sI_vector1, sF_vector1, cfg, rM);
    history.r(end+1,1) = rM;  %#ok<AGROW>
    history.S(end+1,1) = S_tot;

    if abs(S_tot) < 1e-5, break; end
    if S_tot > 0
        rU = rM;
    else
        rL = rM;
    end
end
r_star_common = rM;
end
