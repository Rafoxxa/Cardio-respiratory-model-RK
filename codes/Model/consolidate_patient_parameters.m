function consolidate_patient_parameters(patient_ids, output_filename)
    % Consolida los parámetros de múltiples pacientes desde archivos Excel
    % 
    % Inputs:
    %   patient_ids: vector con los IDs de los pacientes (ej: [1, 2, 3, 4])
    %   output_filename: nombre del archivo Excel de salida (opcional)
    %
    % Ejemplo de uso:
    %   consolidate_patient_parameters([1, 2, 3, 4])
    %   consolidate_patient_parameters([1, 2, 3, 4], 'resumen_pacientes.xlsx')
    
    if nargin < 2
        output_filename = 'parametros_consolidados.xlsx';
    end
    
    n_patients = length(patient_ids);
    
    % Read data from each patient file
    all_data = cell(n_patients, 1);
    param_names = {};
    
    fprintf('Leyendo archivos de pacientes...\n');
    
    for i = 1:n_patients
        patient_id = patient_ids(i);
        filename = sprintf('estimate_ohmnewton_parametros_paciente_%d.xlsx', patient_id);
        
        try
            % Read the Excel file
            [~, ~, raw_data] = xlsread(filename, 'Comparación');
            
            % Skip header row and extract data
            if i == 1
                % Get parameter names from first file
                param_names = raw_data(2:end, 1);
                % Remove empty cells
                param_names = param_names(~cellfun(@isempty, param_names));
            end
            
            % Store all data for this patient
            all_data{i} = raw_data;
            fprintf('  ✓ Paciente %d cargado\n', patient_id);
            
        catch ME
            error('Error al leer archivo del paciente %d: %s', patient_id, ME.message);
        end
    end
    
    fprintf('\nProcesando datos...\n');
    
    % Initialize matrices for values
    n_params = length(param_names);
    old_values = zeros(n_params, 1);
    new_values = zeros(n_params, n_patients);
    
    % Extract values from each patient
    for i = 1:n_patients
        raw_data = all_data{i};
        
        for j = 1:n_params
            param_name = param_names{j};
            
            % Find the row for this parameter
            row_idx = find(strcmp(raw_data(:, 1), param_name));
            
            if ~isempty(row_idx)
                row_idx = row_idx(1); % Take first match
                
                % Get old value (only from first patient, as it's common)
                if i == 1
                    old_val = raw_data{row_idx, 2};
                    if isnumeric(old_val)
                        old_values(j) = old_val;
                    else
                        old_values(j) = NaN;
                    end
                end
                
                % Get new value for this patient
                new_val = raw_data{row_idx, 3};
                if isnumeric(new_val)
                    new_values(j, i) = new_val;
                else
                    new_values(j, i) = NaN;
                end
            else
                if i == 1
                    old_values(j) = NaN;
                end
                new_values(j, i) = NaN;
            end
        end
    end
    
    % Calculate statistics
    mean_values = mean(new_values, 2, 'omitnan');
    std_values = std(new_values, 0, 2, 'omitnan');
    
    % Create consolidated table
    n_cols = 3 + n_patients + 2; % Param name + Old value + patients + mean + std
    consolidated_data = cell(n_params + 1, n_cols);
    
    % Headers
    consolidated_data{1, 1} = 'Parámetro';
    consolidated_data{1, 2} = 'Valor Inicial';
    for i = 1:n_patients
        consolidated_data{1, 2 + i} = sprintf('Paciente %d', patient_ids(i));
    end
    consolidated_data{1, 3 + n_patients} = 'Promedio';
    consolidated_data{1, 4 + n_patients} = 'Desv. Estándar';
    
    % Fill data
    for i = 1:n_params
        consolidated_data{i + 1, 1} = param_names{i};
        
        % Old value
        if isnan(old_values(i))
            consolidated_data{i + 1, 2} = 'N/A';
        else
            consolidated_data{i + 1, 2} = old_values(i);
        end
        
        % New values for each patient
        for j = 1:n_patients
            consolidated_data{i + 1, 2 + j} = new_values(i, j);
        end
        
        % Statistics
        consolidated_data{i + 1, 3 + n_patients} = mean_values(i);
        consolidated_data{i + 1, 4 + n_patients} = std_values(i);
    end
    
    % Write to Excel
    fprintf('\nGuardando archivo consolidado...\n');
    try
        writecell(consolidated_data, output_filename, 'Sheet', 'Datos Consolidados');
        
        % Add summary sheet
        summary_data = cell(6, 2);
        summary_data{1, 1} = 'RESUMEN';
        summary_data{2, 1} = 'Número de pacientes';
        summary_data{2, 2} = n_patients;
        summary_data{3, 1} = 'IDs de pacientes';
        summary_data{3, 2} = strjoin(arrayfun(@num2str, patient_ids, 'UniformOutput', false), ', ');
        summary_data{4, 1} = 'Número de parámetros';
        summary_data{4, 2} = n_params;
        summary_data{5, 1} = 'Fecha de generación';
        summary_data{5, 2} = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        
        writecell(summary_data, output_filename, 'Sheet', 'Resumen');
        
        fprintf('✓ Archivo guardado exitosamente: %s\n', output_filename);
        fprintf('  - Hoja 1: Datos Consolidados (parámetros de todos los pacientes)\n');
        fprintf('  - Hoja 2: Resumen (información general)\n');
        
    catch ME
        error('Error al guardar el archivo Excel: %s', ME.message);
    end
    
    % Display summary statistics
    fprintf('\n=== ESTADÍSTICAS GENERALES ===\n');
    fprintf('Total de pacientes procesados: %d\n', n_patients);
    fprintf('Total de parámetros: %d\n', n_params);
    fprintf('Archivo generado: %s\n\n', output_filename);
    
end