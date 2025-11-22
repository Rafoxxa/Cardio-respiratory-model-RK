%% FUNCIÓN 1: reduce_by_param_names
%% Elimina parámetros del tensor, matriz de sensibilidad y lista de nombres
%% según una lista de nombres a mantener (o a excluir)
%% =========================================================================
function [STensor_red, sens_matrix_red, pars_to_sens_red] = ...
    reduce_by_param_names(STensor, sens_matrix, pars_to_sens, varargin)
% REDUCE_BY_PARAM_NAMES  Filtra parámetros por nombre en las tres estructuras
%
% Uso:
%   [S_red, M_red, P_red] = reduce_by_param_names(STensor, sens_matrix, pars_to_sens, ...
%                                       'keep', {'k1','alpha','R0'}, ...
%                                       'exclude', {'deadspace','fuzzy'});
%
% Opciones:
%   'keep'     : cell con nombres de parámetros que SÍ quieres mantener
%   'exclude'  : cell con nombres de parámetros que quieres eliminar
%   'regex'    : expresión regular (ej: '^k[0-9]+$' para todos los k1,k2,...)
%
% Si no se da ninguna opción → devuelve todo (útil para pruebas)

    p = inputParser;
    addParameter(p, 'keep',    {}, @(x) iscellstr(x) || isstring(x));
    addParameter(p, 'exclude', {}, @(x) iscellstr(x) || isstring(x));
    addParameter(p, 'regex',   '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    
    keep_names    = p.Results.keep;
    exclude_names = p.Results.exclude;
    regex_pat     = p.Results.regex;
    
    % Convertir todo a cellstr para uniformidad
    pars_cell = cellstr(pars_to_sens);
    
    % Construir máscara lógica de parámetros a mantener
    mask = true(size(pars_cell));
    
    if ~isempty(exclude_names)
        exclude_names = cellstr(exclude_names);
        [~, idx_ex] = ismember(exclude_names, pars_cell);
        idx_ex = idx_ex(idx_ex > 0);
        mask(idx_ex) = false;
    end
    
    if ~isempty(keep_names)
        keep_names = cellstr(keep_names);
        [~, idx_keep] = ismember(keep_names, pars_cell);
        if any(idx_keep == 0)
            missing = keep_names(idx_keep == 0);
            warning('Los siguientes parámetros no existen y serán ignorados: %s', ...
                strjoin(missing, ', '));
        end
        mask = false(size(pars_cell));
        mask(idx_keep(idx_keep > 0)) = true;
    end
    
    if ~isempty(regex_pat)
        mask = mask & ~cellfun('isempty', regexp(pars_cell, regex_pat, 'once'));
    end
    
    if ~any(mask)
        error('La selección dejó 0 parámetros. Revisa los nombres.');
    end
    
    % Aplicar máscara (dimensión 1 = parámetros)
    STensor_red       = STensor(mask, :, :);
    sens_matrix_red   = sens_matrix(mask, :);
    pars_to_sens_red  = pars_cell(mask);
    
    fprintf('Reducción completada: %d → %d parámetros\n', length(pars_cell), sum(mask));
end