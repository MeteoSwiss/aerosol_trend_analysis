function TTout = timetable2num(TTin)
%TIMETABLE2NUM Convertit les variables textuelles d'une timetable en données numériques.
%   TTout = TIMETABLE2NUM(TTin)
%   - TTin : timetable d'entrée
%   - TTout : timetable avec colonnes converties en double si possible

    TTout = TTin;  % copie de travail
    varNames = TTin.Properties.VariableNames;

    % Expression régulière pour extraire un nombre
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*)|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*))(?<suffix>.*)';

    for v = 1:numel(varNames)

        col = TTin.(varNames{v});

        % On ne traite que les colonnes de texte ou de cellules
        if ~(iscell(col) || isstring(col) || ischar(col))
            continue
        end

        % Convertir en cellstr pour traitement uniforme
        col = cellstr(col);
        n = numel(col);
        numericCol = NaN(n,1);

        for i = 1:n
            value = col{i};

            % Extraction via regexp
            try
                result = regexp(value, regexstr, 'names');
                numbers = result.numbers;

                % Vérification des virgules
                if contains(numbers, ',')
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'))
                        continue
                    end
                    numbers = strrep(numbers, ',', '');
                end

                % Conversion
                tmp = textscan(numbers, '%f');
                numericCol(i) = tmp{1};

            catch
                % Valeur non convertible → NaN
                continue
            end
        end

        % Remplacement dans la timetable
        TTout.(varNames{v}) = numericCol;
    end
end
