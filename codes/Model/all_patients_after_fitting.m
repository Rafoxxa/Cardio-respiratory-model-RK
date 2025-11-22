function bestpars_table = all_patients_after_fitting()
    set(0, 'DefaultFigureVisible', 'off');

    % Sujetos de interés (orden tal como usaste antes)
    sujetos = [1, 2, 3,4];
    num_sujetos = numel(sujetos);

    % Primera iteración: obtener estructura base (sujeto 1 en tu esquema)
    setup0 = ForwardFittingModel("hard:only-loading", 1);
    variable_names = setup0.xnames_fitting;
    keys_pars = setup0.pars.keys;
    names = keys_pars(setup0.idx_optpars);
    num_parametros = numel(names);

    % Preallocaciones para recoger datos por sujeto
    data_matrix = nan(num_parametros, num_sujetos); % valores optimizados (por sujeto)
    lb_matrix = nan(num_parametros, num_sujetos);   % lower bounds por sujeto
    ub_matrix = nan(num_parametros, num_sujetos);   % upper bounds por sujeto
    optpars0_all = cell(1, num_sujetos);           % valores iniciales por sujeto (para comprobar igualdad)
    old_pars_cell = cell(1, num_sujetos);          % mapas de pars antiguos por sujeto
    errores_H = []; errores_N = [];

    % Recolectar datos para cada sujeto
    for sidx = 1:num_sujetos
        subj = sujetos(sidx);
        setup = ForwardFittingModel("hard:only-loading", subj);

        % recoger best x (asegurar orientación columna)
        bestx = setup.out_solver.x;
        if isrow(bestx)
            bestx = bestx';
        end
        if numel(bestx) ~= num_parametros
            error('Número de parámetros optimizados (bestx) no coincide con names para sujeto %d.', subj);
        end
        data_matrix(:, sidx) = bestx;

        % bounds específicos por sujeto (se asume vector columna del mismo largo)
        lb_tiny = setup.lb_old_tiny; % debería ser vector de length num_parametros
        ub_tiny = setup.ub_old_tiny;
        if numel(lb_tiny) ~= num_parametros || numel(ub_tiny) ~= num_parametros
            error('lb_old_tiny/ub_old_tiny no coinciden en dimensiones para sujeto %d.', subj);
        end
        lb_matrix(:, sidx) = lb_tiny(:);
        ub_matrix(:, sidx) = ub_tiny(:);

        optpars0_all{sidx} = setup.optpars_0_old{end};       

        % pars antiguos (map)
        updated_pars_old = setup.updated_pars_old;
        old_pars_cell{sidx} = mergeUpdatedPars(updated_pars_old);

        % errores (JH, JN) con reorden que usabas
        JH = setup.JH; JN = setup.JN;
        JJH = [JH(3) JH(4) JH(5) JH(6) JH(1) JH(2) JH(7) JH(10) JH(8) JH(9)];
        JJN = [JN(3) JN(4) JN(5) JN(6) JN(1) JN(2) JN(7) JN(10) JN(8) JN(9)];
        errores_H(:, end+1) = JJH(:);
        errores_N(:, end+1) = JJN(:);
    end

    % Comprobar que valores iniciales son consistentes (si todos iguales, usar primero)
    % Si no iguales conserva la primera columna y alerta por consola (no detiene ejecución)
    num_parametros = numel(optpars0_all{1});  % cantidad de parámetros
    valores_iniciales = NaN(num_parametros,1); % inicializar como NaN
    
    for pidx = 1:num_parametros
        % Extraer el valor de cada sujeto para este parámetro
        vals = cellfun(@(c)c(pidx), optpars0_all);
        vals = cell2mat(vals);
        if all(vals == vals(1))   % todos iguales
            valores_iniciales(pidx) = vals(1);
        else
            valores_iniciales(pidx) = NaN; % distinto entre sujetos
        end
    end

    
    num_parametros = size(data_matrix,1);
    num_sujetos = size(data_matrix,2);
    
    % Construir matriz de valores iniciales por sujeto
    optpars0_matrix = nan(num_parametros, num_sujetos);
    for sidx = 1:num_sujetos
        optpars0_matrix(:, sidx) = cell2mat(optpars0_all{sidx});
    end
    
    % Calcular distancia euclidiana columna a columna
    % Nota: aquí asumimos que optpars0_matrix y data_matrix no contienen NaN; 
    % si aparece NaN, indica que hay un error previo
    dist_per_col = sqrt(sum((data_matrix - optpars0_matrix).^2, 1));


    % Calcular promedio y desviación estándar sobre los valores optimizados (solo sujetos)
    Promedio = mean(data_matrix, 2);
    Desviacion_Estandar = std(data_matrix, 0, 2);

    
    % Construcción de nombres de columnas en el orden requerido:
    % Parametro | LB_S1 | S1 | UB_S1 | ... (cada sujeto) ... | Promedio | Desviacion_Estandar | Valores_iniciales_comunes
    col_names = cell(1, 1 + num_sujetos*3 + 3);
    col_names{1} = 'Parametros';
    idx = 2;
    sujeto_colnames = cell(1, num_sujetos);
    for sidx = 1:num_sujetos
        subj = sujetos(sidx);
        col_names{idx} = sprintf('LB_S%g', subj); idx = idx + 1;
        col_names{idx} = sprintf('S%g', subj);   idx = idx + 1;
        col_names{idx} = sprintf('UB_S%g', subj); idx = idx + 1;
        sujeto_colnames{sidx} = sprintf('S%g', subj);
    end
    col_names{idx} = 'Promedio'; idx = idx + 1;
    col_names{idx} = 'Desviacion_Estandar'; idx = idx + 1;
    col_names{idx} = 'Valores_iniciales_comunes';

    % Crear tabla base con parámetros optimizados
    tabla_parametros = table();
    tabla_parametros.Parametros = names(:);

    % Añadir columnas LB/S/UB por sujeto
    for sidx = 1:num_sujetos
        subj = sujetos(sidx);
        lb_col = lb_matrix(:, sidx);
        val_col = data_matrix(:, sidx);
        ub_col = ub_matrix(:, sidx);

        tabla_parametros.(sprintf('LB_S%g', subj)) = lb_col;
        tabla_parametros.(sprintf('S%g', subj)) = val_col;
        tabla_parametros.(sprintf('UB_S%g', subj)) = ub_col;
    end

    % Añadir Promedio, Desviación y Valores iniciales
    tabla_parametros.Promedio = Promedio;
    tabla_parametros.Desviacion_Estandar = Desviacion_Estandar;
    tabla_parametros.Valores_iniciales_comunes = valores_iniciales;

    % Añadir fila adicional 'Distancia al p0' (valores por columna: sujetos + inicial)
    % Crear fila con NaN para campos LB/S/UB/Promedio/Std/Inicial y rellenar S* y Valores_iniciales_comunes específicas
    % Construir fila como celda para concatenar fácilmente
    % Primero determinar la cantidad de columnas en la tabla actual
    current_vars = tabla_parametros.Properties.VariableNames;
    nvars = numel(current_vars);

    % Convertir tabla_parametros a celdas antes de agregar fila heterogénea
    tabla_cells = table2cell(tabla_parametros);
    
    % Crear fila de distancia
    dist_row = cell(1, nvars);
    dist_row{1} = 'Distancia al p0';
    for sidx = 1:num_sujetos
        colname_val = sprintf('S%g', sujetos(sidx));
        vidx = find(strcmp(current_vars, colname_val));
        if ~isempty(vidx)
            dist_row{vidx} = dist_per_col(sidx);
        end
    end
    vidx_init = find(strcmp(current_vars,'Valores_iniciales_comunes'));
    dist_row{vidx_init} = dist_per_col(end);
    
    % Concatenar fila
    tabla_cells = [tabla_cells; dist_row];
    
    % Reconstruir tabla
    tabla_parametros = cell2table(tabla_cells, 'VariableNames', current_vars);

    % Ajustar Parametros columna (ya la pusimos como string antes)
    % Asegurarse de que Parametros es cellstr, no cell of numeric
    tabla_parametros.Parametros = cellstr(tabla_parametros.Parametros);

    % --- Construcción de filas de parámetros "old" ---
    old_rows = table();
    % Construir plantilla de columnas con NaN compatibles tipos (usar cell temp y luego convertir)
    for subj_idx = 1:num_sujetos
        subj = sujetos(subj_idx);
        % nothing to do: columns created already in tabla_parametros, we'll create temp tables matching same varnames
    end

    % Recolectar todas las claves antiguas para crear filas individuales
    all_old_keys = {};
    for subj_idx = 1:num_sujetos
        old_map = old_pars_cell{subj_idx};
        if isempty(old_map)
            continue;
        end
        keys_i = keys(old_map);
        all_old_keys = [all_old_keys; keys_i(:)];
    end
    all_old_keys = unique(all_old_keys);

    % Crear tabla vacía con mismas columnas (inicialmente as cell)
    empty_row = cell(1, nvars);
    empty_row(:) = {NaN};
    empty_row{1} = ''; % Parametro
    temp_cells_old = cell(0, nvars); % acumular filas como celdas
    for k = 1:numel(all_old_keys)
        key = all_old_keys{k};
        % fila base
        row = empty_row;
        row{1} = [key '_old_value'];
        % rellenar solo las columnas S<subj> si tienen valor
        for subj_idx = 1:num_sujetos
            old_map = old_pars_cell{subj_idx};
            subj = sujetos(subj_idx);
            if isempty(old_map)
                continue;
            end
            if isKey(old_map, key)
                val = old_map(key);
                % ubicar columna S<subj>
                colname_val = sprintf('S%g', subj);
                vidx = find(strcmp(current_vars, colname_val));
                if ~isempty(vidx)
                    row{vidx} = val;
                end
                % IMPORTANTE: según tu petición, no rellenar LB_* ni UB_* para filas old (quedarán NaN)
            end
        end
        temp_cells_old = [temp_cells_old; row];
    end

    % Guardar columna de parámetros como estaba originalmente
    param_col = tabla_parametros.Parametros;
    
    % Elimina temporalmente la columna para no convertirla erróneamente
    tabla_parametros.Parametros = [];
    
    % Normalizar columnas numéricas: convertir celdas que contienen scalars a double
    vars = tabla_parametros.Properties.VariableNames;
    nrows = height(tabla_parametros);
    
    for i = 1:numel(vars)
        col = tabla_parametros.(vars{i});
        
        % Si ya es numeric (double), nada que hacer
        if isnumeric(col)
            continue;
        end
        
        % Si la columna es cell, intentamos convertirla a vector double
        if iscell(col)
            newcol = nan(nrows,1); % por defecto NaN
            for r = 1:nrows
                entry = col{r};
                if isempty(entry)
                    % 0x0 double u otro vacío -> NaN
                    newcol(r) = NaN;
                elseif isnumeric(entry)
                    % numeric: aceptar solo escalares
                    if isscalar(entry)
                        newcol(r) = double(entry);
                    else
                        % vector/array numérico: no sabemos qué tomar -> NaN
                        newcol(r) = NaN;
                    end
                else
                    % no numérico (string, cell, map, etc.) -> NaN
                    newcol(r) = NaN;
                end
            end
            tabla_parametros.(vars{i}) = newcol;
            continue;
        end
        
        % Si la columna es string / cellstr, dejar Parametros u otros textos
        if isstring(col) || ischar(col) || iscellstr(col)
            % Mantener, pero para seguridad convertir Parametros a string
            if strcmp(vars{i}, 'Parametros')
                tabla_parametros.Parametros = string(col);
            end
            continue;
        end
        
        % Otros tipos (timetable, categorical, etc.), dejarlos tal cual
    end

    
    % Reincorporar la columna original de parámetros al inicio de la tabla
    tabla_parametros = addvars(tabla_parametros, param_col, 'Before', 1, ...
                               'NewVariableNames', 'Parametros');


    % Convertir temp_cells_old a tabla (si hay filas)
    if ~isempty(temp_cells_old)
        old_rows = cell2table(temp_cells_old, 'VariableNames', current_vars);
        % Asegurar Parametros como cellstr
        old_rows.Parametros = cellstr(old_rows.Parametros);
        % Concatenar con tabla_parametros
        tabla_parametros = [tabla_parametros; old_rows];
    end

    % Mostrar la tabla resultante
    disp('Tabla de Parámetros por Sujeto:');
    disp('=====================================');
    disp(tabla_parametros);

    set(0, 'DefaultFigureVisible', 'on');

    % Mostrar en uitable y guardar imagen
    f = figure('Visible','on');
    uit = uitable(f, ...
        'Data', table2cell(tabla_parametros), ...
        'ColumnName', tabla_parametros.Properties.VariableNames, ...
        'RowName', [], ...
        'FontSize', 12, ...
        'Units','normalized', ...
        'Position', [0 0 1 1]);
    drawnow;
    %frame = getframe(f);
    %imwrite(frame.cdata,'tabla_parametros.png');
    %close(f);

    % --- Tablas de errores (igual que antes) ---
    Tabla_N = array2table(errores_N', 'VariableNames', variable_names);
    Tabla_N.Sujeto = sujetos(:);
    Tabla_N = movevars(Tabla_N, 'Sujeto', 'before', 1);

    Tabla_H = array2table(errores_H', 'VariableNames', variable_names);
    Tabla_H.Sujeto = sujetos(:);
    Tabla_H = movevars(Tabla_H, 'Sujeto', 'before', 1);

    Tabla_summary = table();
    Tabla_summary.Variable = variable_names(:);
    Tabla_summary.N_mean = mean(errores_N, 2);
    Tabla_summary.N_std = std(errores_N, 0, 2);
    Tabla_summary.H_mean = mean(errores_H, 2);
    Tabla_summary.H_std = std(errores_H, 0, 2);

    % Mostrar y guardar tablas de errores
    disp('Tabla de Errores N:');
    disp(Tabla_N);
    disp('Tabla de Errores H:');
    disp(Tabla_H);
    disp('Tabla resumen de errores:');
    disp(Tabla_summary);

    % Guardar tablas de errores como figuras (opcional)
  % --- Tabla de errores N ---
    fN = figure('Name','Errores N','NumberTitle','off','Visible','on');
    uitable(fN, ...
        'Data', table2cell(Tabla_N), ...
        'ColumnName', Tabla_N.Properties.VariableNames, ...
        'RowName', [], ...
        'Units','normalized', ...
        'Position',[0 0 1 1]);
    drawnow;
    
    % --- Tabla de errores H ---
    fH = figure('Name','Errores H','NumberTitle','off','Visible','on');
    uitable(fH, ...
        'Data', table2cell(Tabla_H), ...
        'ColumnName', Tabla_H.Properties.VariableNames, ...
        'RowName', [], ...
        'Units','normalized', ...
        'Position',[0 0 1 1]);
    drawnow;
    
    % --- Tabla resumen ---
    fSummary = figure('Name','Resumen de Errores','NumberTitle','off','Visible','on');
    uitable(fSummary, ...
        'Data', table2cell(Tabla_summary), ...
        'ColumnName', Tabla_summary.Properties.VariableNames, ...
        'RowName', [], ...
        'Units','normalized', ...
        'Position',[0 0 1 1]);
    drawnow;


    % Retornar la tabla resultante
    bestpars_table = tabla_parametros;

    exportSingleTableToExcel('Tabla_Parametros', tabla_parametros, 'Parametros');
    exportSingleTableToExcel('Tabla_J_N', Tabla_N, 'Errores_N');
    exportSingleTableToExcel('Tabala_J_H', Tabla_H, 'Errores_H');
    exportSingleTableToExcel('Tabala_Summary', Tabla_summary, 'Resumen');

end

%% Funciones auxiliares (mantuve tus originales, con leves defensas)
function all_params = mergeUpdatedPars(updated_pars_old)
    % Combina una celda de containers.Map devolviendo un único containers.Map
    all_params = containers.Map();
    if isempty(updated_pars_old)
        return;
    end
    for i = 1:numel(updated_pars_old)
        entry = updated_pars_old{i};
        if isempty(entry)
            continue;
        end
        keys_i = entry.keys;
        values_i = entry.values;
        if numel(keys_i) ~= numel(values_i)
            error('Número de keys y values no coincide en índice %d.', i);
        end
        for k = 1:numel(keys_i)
            all_params(keys_i{k}) = values_i{k};
        end
    end


   

end


 function exportSingleTableToExcel(file_base, T, sheet_name)
    % EXPORTSINGLETABLETOEXCELPDF Exporta una tabla a Excel y luego a PDF.
    %
    % file_base: nombre base del archivo (sin extensión)
    % T: tabla a exportar
    % sheet_name: nombre de la hoja en Excel

    % Nombre del archivo Excel y PDF
    excel_file = [file_base '.xlsx'];
    

    % 1) Guardar la tabla en Excel
    writetable(T, excel_file, 'Sheet', sheet_name, 'WriteRowNames', false);
    fprintf('Tabla "%s" exportada a Excel: %s\n', sheet_name, excel_file);
   
end
