function TT = read_spo_actris(filepath)
% READ_SPO_ACTRIS  Lit un fichier de données aérosols ACTRIS au format SPO
%   et retourne une timetable MATLAB avec des NaN pour les valeurs manquantes.
%
%   TT = read_spo_actris(filepath)
%
%   Structure du fichier (4 lignes d'en-tête) :
%     Ligne 1 : longueurs d'onde (ex. ,467,528,652,...)
%     Ligne 2 : descriptions longues des variables
%     Ligne 3 : noms courts des colonnes  <-- utilisée pour VariableNames
%               Première colonne = "DateTimeUTC"
%     Ligne 4+ : données  "YYYY-MM-DD hh:mm:ss", val1, val2, ...
%
%   Remplacement par NaN :
%     * Valeurs 99999.99  et  99999.999  (flags de données manquantes)
%     * Cellules vides  (champs ,, consécutifs)
%     * Toute valeur non convertible en double
%
%   Arguments :
%     filepath  - Chemin complet vers le fichier (string ou char)
%                 Ex : 'M:\pay-prod\Aerosol_actris_trend\data\spo\spo_1979_2025'
%
%   Retour :
%     TT  - timetable indicée par DateTimeUTC (datetime, format UTC)
%
%   Exemple d'utilisation :
%     TT = read_spo_actris('M:\pay-prod\Aerosol_actris_trend\data\spo\spo_1979_2025');
%     head(TT)
%     summary(TT)

    % -----------------------------------------------------------------------
    % 0.  Constantes
    % -----------------------------------------------------------------------
    N_HEADER_LINES = 3;               % lignes d'en-tête à sauter
    VARNAMES_LINE  = 3;               % ligne (1-based) portant les noms courts
    MISSING_FLAGS  = [99999.99, 99999.999];  % valeurs sentinelles → NaN
    MISSING_TOL    = 1e-4;            % tolérance pour la comparaison flottante

    % -----------------------------------------------------------------------
    % 1.  Vérification de l'existence du fichier
    % -----------------------------------------------------------------------
    if ~isfile(filepath)
        error('read_spo_actris:fileNotFound', ...
              'Fichier introuvable :\n  %s', filepath);
    end

    % -----------------------------------------------------------------------
    % 2.  Lecture des en-têtes pour extraire les noms de variables
    % -----------------------------------------------------------------------
    fid = fopen(filepath, 'r', 'n', 'UTF-8');
    if fid == -1
        error('read_spo_actris:cannotOpen', ...
              'Impossible d''ouvrir le fichier :\n  %s', filepath);
    end

    headerLines = cell(N_HEADER_LINES, 1);
    for k = 1:N_HEADER_LINES
        headerLines{k} = fgetl(fid);
    end
    fclose(fid);

    % -----------------------------------------------------------------------
    % 3.  Extraction des noms de colonnes (ligne VARNAMES_LINE)
    %     Format attendu : "DateTimeUTC,BaB_A11,BaG_A11,..."
    % -----------------------------------------------------------------------
    varLine  = strtrim(headerLines{VARNAMES_LINE});
    allNames = strsplit(varLine, ',');
    allNames = strtrim(allNames);

    if ~strcmpi(allNames{1}, 'DateTimeUTC')
        warning('read_spo_actris:unexpectedHeader', ...
            'Première colonne attendue "DateTimeUTC", trouvée "%s".', allNames{1});
    end

    dataVarNames = allNames(2:end);           % tous les noms sauf la date
    nVars        = numel(dataVarNames);

    % Noms valides et uniques pour MATLAB
    dataVarNames = matlab.lang.makeValidName(dataVarNames);
    dataVarNames = matlab.lang.makeUniqueStrings(dataVarNames);

    % -----------------------------------------------------------------------
    % 4.  Lecture des données avec textscan
    %     - délimiteur virgule
    %     - champs vides → NaN  (EmptyValue + TreatAsMissing)
    %     - 1 colonne string (datetime) + nVars colonnes numériques
    % -----------------------------------------------------------------------
    fprintf('Lecture du fichier en cours...\n');

    fid = fopen(filepath, 'r', 'n', 'UTF-8');

    % Sauter les lignes d'en-tête
    for k = 1:N_HEADER_LINES
        fgetl(fid);
    end

    % Format : %s pour la date, %f pour chaque variable
    fmt = ['%s', repmat(' %f', 1, nVars)];

    C = textscan(fid, fmt, ...
        'Delimiter',      ',', ...
        'EmptyValue',     NaN, ...     % champ vide → NaN immédiatement
        'CollectOutput',  false);
    fclose(fid);

    % -----------------------------------------------------------------------
    % 5.  Conversion de la colonne DateTimeUTC
    % -----------------------------------------------------------------------
    dateStrings = C{1};
    timeVec     = datetime(dateStrings, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % Pour affecter explicitement le fuseau UTC, décommenter la ligne suivante :
    % timeVec.TimeZone = 'UTC';

    nRows = numel(timeVec);

    % -----------------------------------------------------------------------
    % 6.  Construction de la matrice de données + remplacement des sentinelles
    % -----------------------------------------------------------------------
    dataMatrix = NaN(nRows, nVars);

    for j = 1:nVars
        col = C{j + 1};            % vecteur double pour la j-ième variable

        % Remplacer chaque valeur sentinelle par NaN
        for flag = MISSING_FLAGS
            col(abs(col - flag) < MISSING_TOL) = NaN;
        end

        dataMatrix(:, j) = col;
    end

    % -----------------------------------------------------------------------
    % 7.  Assemblage de la timetable
    % -----------------------------------------------------------------------
    T  = array2table(dataMatrix, 'VariableNames', dataVarNames);
    TT = table2timetable(T, 'RowTimes', timeVec);

    % -----------------------------------------------------------------------
    % 8.  Rapport de lecture
    % -----------------------------------------------------------------------
    nNaN   = sum(sum(ismissing(TT)));
    pctNaN = 100 * nNaN / (nRows * nVars);

    fprintf('Lecture terminée : %d lignes x %d variables.\n', nRows, nVars);
    fprintf('Valeurs NaN      : %d (%.1f %% des cellules).\n', nNaN, pctNaN);
    fprintf('Période couverte : %s  -->  %s\n', ...
            string(TT.Time(1)), string(TT.Time(end)));
end