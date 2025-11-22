%% =========================================================================
function [STensor_red, pars_red, idx_kept] = reduce_STensor_by_names(STensor, pars_original, pars_to_keep)
% REDUCE_STENSOR_BY_NAMES  Elimina parámetros del STensor según nombres
%
% Entradas:
%   STensor       : tensor original [nPar × nVar × nTime] o [nPar × nVar]
%   pars_original : cell con todos los nombres originales (mismo orden que STensor)
%   pars_to_keep  : cell con los nombres de parámetros que QUIERES CONSERVAR
%
% Salidas:
%   STensor_red   : tensor reducido (solo parámetros seleccionados)
%   pars_red      : nombres de los parámetros conservados (subconjunto)
%   idx_kept      : índices originales de los parámetros mantenidos (útil para debug)
%
% Ejemplo:
%   [S_red, p_red] = reduce_STensor_by_names(STensor, pars_full, pars_filtered);

    % Convertir a cellstr por seguridad
    pars_original = cellstr(pars_original);
    pars_to_keep  = cellstr(pars_to_keep);
    
    % Encontrar coincidencias exactas
    [found, idx_kept] = ismember(pars_to_keep, pars_original);
    
    if ~all(found)
        missing = pars_to_keep(~found);
        error(['Los siguientes parámetros no se encontraron en pars_original:\n   ' ...
               strjoin(missing, ', ') '\nRevisa ortografía o mayúsculas/minúsculas.']);
    end
    
    % Reordenar idx_kept para mantener el orden de pars_to_keep
    idx_kept = idx_kept(found);  % ya está en orden correcto por ismember
    
    % Reducir el tensor (soporta 2D o 3D)
    if ndims(STensor) == 3
        STensor_red = STensor(idx_kept, :, :);
    elseif ndims(STensor) == 2
        STensor_red = STensor(idx_kept, :);
    else
        error('STensor debe ser 2D (matriz) o 3D (tensor)');
    end
    
    pars_red = pars_to_keep;  % ya está filtrado y en el orden deseado
    
    % Reporte bonito
    fprintf('Reducción exitosa: %d → %d parámetros\n', length(pars_original), length(pars_red));
    if length(pars_red) <= 20
        disp('Parámetros conservados:');
        disp(pars_red);
    end
end