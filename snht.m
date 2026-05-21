function [T_max, breakpoint_time, p_value, PrctDiff] = snht(input_data, var_name, alpha)
% SNHT - Standard Normal Homogeneity Test (Alexandersson, 1986)
%
% Accepte un vecteur numérique OU une timetable MATLAB. Quand une timetable
% est fournie, la rupture est exprimée en temps (datetime/duration) plutôt
% qu'en indice. Si le test n'est pas significatif, aucune rupture n'est
% renvoyée (breakpoint_time = NaT) et aucun graphique n'est tracé.
%
% SYNTAXE
%   [T_max, bp_time, p_value, T_series] = snht(tt)
%   [T_max, bp_time, p_value, T_series] = snht(tt, var_name)
%   [T_max, bp_time, p_value, T_series] = snht(tt, var_name, alpha)
%   [T_max, bp_time, p_value, T_series] = snht(data_vector)
%   [T_max, bp_time, p_value, T_series] = snht(data_vector, [], alpha)
%
% ENTREES
%   input_data - Timetable (recommande) OU vecteur numerique (nx1 ou 1xn)
%   var_name   - Nom de la variable dans la timetable a analyser
%                (optionnel : si absent, prend la 1ere variable numerique)
%   alpha      - Niveau de signification (defaut : 0.05)
%
% SORTIES
%   T_max           - Statistique maximale du test
%   breakpoint_time - Temps de la rupture (datetime si timetable, indice si
%                     vecteur). Vaut NaT / NaN si non significatif.
%   p_value         - p-valeur estimee par bootstrap parametrique (1000 tirages)
%   T_series        - Vecteur des statistiques T(k), k = 1...n-1
%
% EXEMPLES
%   % --- Avec une timetable ---
%   t    = datetime(2000,1,1) + caldays(0:99)';
%   vals = [randn(50,1); randn(50,1) + 1.5];
%   tt   = timetable(t, vals, 'VariableNames', {'Temperature'});
%   [T, bp, pv] = snht(tt, 'Temperature');
%
%   % --- Avec un vecteur ---
%   data = [randn(50,1); randn(50,1) + 1.5];
%   [T, bp, pv] = snht(data);
%
% REFERENCE
%   Alexandersson, H. (1986). A homogeneity test applied to precipitation data.
%   Journal of Climatology, 6(6), 661-675.

    %% --- 1. Analyse des arguments d'entree ---
    is_timetable = isa(input_data, 'timetable');

    if nargin < 2 || isempty(var_name)
        var_name = [];
    end
    if nargin < 3
        alpha = 0.05;
    end

    %% --- 2. Extraction des donnees et de l'axe temporel ---
    if is_timetable
        tt = input_data;

        % Selection de la variable a analyser
        if isempty(var_name)
            vars  = tt.Properties.VariableNames;
            found = false;
            for v = 1:numel(vars)
                if isnumeric(tt.(vars{v}))
                    var_name = vars{v};
                    found    = true;
                    break;
                end
            end
            if ~found
                error('snht: aucune variable numerique trouvee dans la timetable.');
            end
        else
            if ~ismember(var_name, tt.Properties.VariableNames)
                error('snht: variable "%s" introuvable dans la timetable.', var_name);
            end
        end
%% remove NAN
ind_nan=isnan(tt.(var_name));
tt=tt(~ind_nan,:) ;

        data          = double(tt.(var_name));
        time_vec      = tt.Properties.RowTimes;   % datetime ou duration
        var_name_plot = var_name;

    else
        data          = double(input_data(:));
        time_vec      = [];
        var_name_plot = 'Valeur';
    end

    data = data(:);
    n    = numel(data);

    if n < 10
        error('snht: la serie doit contenir au moins 10 observations.');
    end

    %% --- 3. Standardisation ---
    mu    = mean(data, 'omitnan');
    sigma = std(data,  'omitnan');

    if sigma == 0
        error('snht: la serie est constante, le test ne peut pas etre applique.');
    end

    z = (data - mu) / sigma;

    %% --- 4. Calcul vectorise de T(k) via sommes cumulees ---
    cumz     = cumsum(z);
    T_series = zeros(n - 1, 1);

    for k = 1:(n - 1)
        z1_bar   = cumz(k) / k;
        z2_bar   = (cumz(n) - cumz(k)) / (n - k);
        T_series(k) = k * z1_bar^2 + (n - k) * z2_bar^2;
    end

    %% --- 5. Statistique maximale et indice de rupture ---
    [T_max, bp_idx] = max(T_series);

    %% --- 6. Bootstrap parametrique pour la p-valeur ---
    n_boot = 1000;
    T_boot = zeros(n_boot, 1);

    for b = 1:n_boot
        z_sim = randn(n, 1);
        cs    = cumsum(z_sim);
        T_sim = zeros(n - 1, 1);
        for k = 1:(n - 1)
            m1       = cs(k) / k;
            m2       = (cs(n) - cs(k)) / (n - k);
            T_sim(k) = k * m1^2 + (n - k) * m2^2;
        end
        T_boot(b) = max(T_sim);
    end

    p_value   = mean(T_boot >= T_max);
    is_signif = (p_value <= alpha);

    %% --- 7. Point de rupture exprime en temps ---
    if is_timetable
        if is_signif
            breakpoint_time = time_vec(bp_idx);
        else
            if isa(time_vec, 'datetime')
                breakpoint_time = NaT;
            else
                breakpoint_time = NaN;
            end
        end
    else
        if is_signif
            breakpoint_time = bp_idx;
        else
            breakpoint_time = NaN;
        end
    end
    %% median diff before and after break point
PrctDiff=(prctile(data(bp_idx:end),[10 50 90])-prctile(data(1:bp_idx),[10 50 90]).*100./prctile(data,[10 50 90]));

    %% --- 8. Affichage texte ---
    fprintf('\n=== Resultats du SNHT ===\n');
    fprintf('Variable analysee     : %s\n', var_name_plot);
    fprintf('Taille de la serie    : %d\n', n);
    fprintf('Statistique T_max     : %.4f\n', T_max);
    fprintf('p-valeur (bootstrap)  : %.4f\n', p_value);

    if is_signif
        if is_timetable
            fprintf('Point de rupture      : %s  (indice %d)\n', ...
                string(breakpoint_time), bp_idx);
        else
            fprintf('Point de rupture      : indice %d\n', bp_idx);
        end
        fprintf('Decision (alpha=%.2f)  : RUPTURE SIGNIFICATIVE detectee\n', alpha);
    else
        fprintf('Decision (alpha=%.2f)  : Pas de rupture significative\n', alpha);
    end
    fprintf('=========================\n\n');

    %% --- 9. Visualisation uniquement si significatif ---
    if ~is_signif
        return;
    end

    %figure('Name', 'SNHT - Standard Normal Homogeneity Test', 'NumberTitle', 'off');

    % Axes des abscisses
    if is_timetable
        x_data   = time_vec;
        x_break  = breakpoint_time;
        x_T      = time_vec(1:end-1);
        x_label  = 'Temps';
        bp_label = sprintf('Rupture : %s', string(breakpoint_time));
        bp_title = string(breakpoint_time);
    else
        x_data   = (1:n)';
        x_break  = bp_idx;
        x_T      = (1:n-1)';
        x_label  = 'Indice temporel';
        bp_label = sprintf('Rupture (k=%d)', bp_idx);
        bp_title = sprintf('k = %d', bp_idx);
    end

    % Panneau 1 : serie temporelle avec rupture
    % % subplot(2, 1, 1);
    % % plot(x_data, data, 'b-o', 'MarkerSize', 3, 'LineWidth', 1);
    % % hold on;
    % % xline(x_break, 'r--', 'LineWidth', 1.5, ...
    % %       'Label', bp_label, 'LabelVerticalAlignment', 'bottom');
    % % xlabel(x_label);
    % % ylabel(var_name_plot);
    % % title(sprintf('Serie temporelle : %s', var_name_plot));
    % % grid on;
    % % 
    % % % Panneau 2 : statistique T(k)
    % % subplot(2, 1, 2);
    % % plot(x_T, T_series, 'k-', 'LineWidth', 1.2);
    % % hold on;
    % % plot(x_break, T_max, 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    % % xlabel(x_label);
    % % ylabel('T(k)');
    % % title(sprintf('Statistique T(k)  —  T_{max} = %.2f  a  %s', T_max, bp_title));
    % % grid on;
    % % 
    % % sgtitle(sprintf('SNHT  |  p = %.4f  |  n = %d  |  \\alpha = %.2f', ...
    % %         p_value, n, alpha));
end
