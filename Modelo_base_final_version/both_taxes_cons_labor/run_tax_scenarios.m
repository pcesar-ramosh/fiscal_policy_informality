function R = run_tax_scenarios(kind, cfg0, grid_vals)
% RUN_TAX_SCENARIOS  Ejecuta escenarios para IVA ('VAT') o laboral ('LAB')
% y arma estructuras de resultados + tablas (stats, fiscal, assets, summary).
%
% Requiere en el path:
%   - solve_two_type_huggett_fiscal_VATdifferential.m   (para 'VAT')
%   - solve_two_type_huggett_fiscal_LABtax.m            (para 'LAB')
%
% Uso:
%   R = run_tax_scenarios('VAT', cfg0, [0.18 0.20 0.22]);
%   R = run_tax_scenarios('LAB', cfg0, [0.12 0.15 0.18 0.22]);

    sc = struct('name',{},'cfg',{},'sol',{},'param_value',{});
    labels = cell(numel(grid_vals),1);

    switch upper(kind)
        case 'VAT'
            solver = @(cfg,val,extra) solve_two_type_huggett_fiscal_VATdifferential(setfield(cfg,'tau_c',val), val);
            for k=1:numel(grid_vals)
                sc(k).name = sprintf('VAT_%02d', round(100*grid_vals(k)));
                sc(k).cfg  = cfg0; sc(k).cfg.tau_c = grid_vals(k);
                sc(k).param_value = grid_vals(k);
                labels{k} = sprintf('\\tau_c=%.2f', grid_vals(k));
                sc(k).sol = solver(sc(k).cfg, grid_vals(k), []);
                print_banner(sc(k).name, sc(k).sol);
            end

        case 'LAB'
            solver = @(cfg,val,extra) solve_two_type_huggett_fiscal_LABtax(setfield(cfg,'tau_l',val), val);
            for k=1:numel(grid_vals)
                sc(k).name = sprintf('LAB_%02d', round(100*grid_vals(k)));
                sc(k).cfg  = cfg0; sc(k).cfg.tau_l = grid_vals(k);
                sc(k).param_value = grid_vals(k);
                labels{k} = sprintf('\\tau_\\ell=%.2f', grid_vals(k));
                sc(k).sol = solver(sc(k).cfg, grid_vals(k), []);
                print_banner(sc(k).name, sc(k).sol);
            end

        otherwise
            error('kind debe ser ''VAT'' o ''LAB''.');
    end

    % ----- Tablas (I/F y total) -----
    T_stats  = scenarios_stats_table(sc);
    T_fiscal = scenarios_fiscal_table(sc);
    T_assets = scenarios_assets_table(sc);
    T_sum    = scenarios_summary_table(sc, labels, grid_vals);

    % ----- Salida -----
    R.kind = upper(kind);
    R.sc = sc;
    R.labels = labels;
    R.solver_handle = solver;   % usado por funciones de plotting (MPC y mercado)
    R.extra = [];               % placeholder
    R.T_stats   = T_stats;
    R.T_fiscal  = T_fiscal;
    R.T_assets  = T_assets;
    R.T_summary = T_sum;
end

% ========================== Helpers de tablas =============================

function T = scenarios_stats_table(sc)
    % 19 columnas: scenario, r, popI, popF, Y, C,
    %              wmean_I/F/T, giniW_I/F/T, cmean_I/F/T, giniC_I/F/T, p11_rep
    n = numel(sc);
    rows = cell(n, 19);
    for k=1:n
        s = sc(k).sol; S = s.stats;
        rows(k,:) = { sc(k).name, s.r, s.popI, s.popF, s.Y, s.Ctot, ...
                      S.wealth_mean(1), S.wealth_mean(2), S.wealth_mean(3), ...
                      S.giniW(1),       S.giniW(2),       S.giniW(3), ...
                      S.cons_mean(1),   S.cons_mean(2),   S.cons_mean(3), ...
                      S.giniC(1),       S.giniC(2),       S.giniC(3), ...
                      S.p11 };
    end
    T = cell2table(rows, 'VariableNames', ...
        {'scenario','r','popI','popF','Y','C', ...
         'wmean_I','wmean_F','wmean_T','giniW_I','giniW_F','giniW_T', ...
         'cmean_I','cmean_F','cmean_T','giniC_I','giniC_F','giniC_T','p11_rep'});
end

function T = scenarios_fiscal_table(sc)
    % Detecta si hay desagregación TcI/TcF (propio de VAT no-neutral)
    f0 = sc(1).sol.fiscal;
    hasTcSplit = isfield(f0,'TcI') && isfield(f0,'TcF');

    n = numel(sc);
    if hasTcSplit
        rows = cell(n, 11);
        for k=1:n
            f = sc(k).sol.fiscal;
            rows(k,:) = { sc(k).name, f.Tl, f.Tc, f.TcI, f.TcF, f.Tr, f.G, f.rB, f.PB, f.B, f.BB };
        end
        T = cell2table(rows, 'VariableNames', ...
            {'scenario','Tl','Tc','TcI','TcF','Tr','G','rB','PB','B','BB'});
    else
        rows = cell(n, 9);
        for k=1:n
            f = sc(k).sol.fiscal;
            rows(k,:) = { sc(k).name, f.Tl, f.Tc, f.Tr, f.G, f.rB, f.PB, f.B, f.BB };
        end
        T = cell2table(rows, 'VariableNames', ...
            {'scenario','Tl','Tc','Tr','G','rB','PB','B','BB'});
    end
end

function T = scenarios_assets_table(sc)
    % 5 columnas: scenario, A_I, A_F, A_priv, B_public
    n = numel(sc);
    rows = cell(n, 5);
    for k=1:n
        s = sc(k).sol; da = s.a(2)-s.a(1);
        A_I   = sum(s.g(:,1).*s.a)*da;
        A_F   = sum(s.g(:,2).*s.a)*da;
        A_prv = A_I + A_F;
        rows(k,:) = { sc(k).name, A_I, A_F, A_prv, s.fiscal.B };
    end
    T = cell2table(rows, 'VariableNames', {'scenario','A_I','A_F','A_priv','B_public'});
end

function T = scenarios_summary_table(sc, labels, grid_vals)
    % 11 columnas: scenario(label bonito), policy_val, r, Y, C, G, PB, BB, giniW_T, giniC_T, eta_real
    n = numel(sc);
    rows = cell(n, 11);
    for k=1:n
        s = sc(k).sol;
        rows(k,:) = { labels{k}, grid_vals(k), s.r, s.Y, s.Ctot, s.fiscal.G, ...
                      s.fiscal.PB, s.fiscal.BB, s.stats.giniW(3), s.stats.giniC(3), ...
                      s.popI/(s.popI+s.popF) };
    end
    T = cell2table(rows, 'VariableNames', ...
        {'scenario','policy_val','r','Y','C','G','PB','BB','giniW_T','giniC_T','eta_real'});
end

% ============================ Helper banner ===============================

function print_banner(nm, s)
    fprintf('[%s] r*=%.4f | eta=%.3f | B=%.3f | G=%.3f | Y=%.3f | C=%.3f\n', ...
        nm, s.r, s.popI/(s.popI+s.popF), s.fiscal.B, s.fiscal.G, s.Y, s.Ctot);
end
