function generate_percentages_table(percentages_file, output_filename)
    % Genera una tabla que mapea los índices de x con los porcentajes usados
    % y los parámetros correspondientes en estimate_newton_ohm
    %
    % Inputs:
    %   percentages_file: archivo .mat con la estructura 'percentages'
    %   output_filename: nombre del archivo Excel de salida (opcional)
    %
    % Ejemplo de uso:
    %   generate_percentages_table('Ohmpercentages.mat')
    %   generate_percentages_table('Ohmpercentages.mat', 'tabla_porcentajes.xlsx')
    
    if nargin < 2
        output_filename = 'tabla_porcentajes.xlsx';
    end
    
    % Load percentages
    fprintf('Cargando archivo de porcentajes...\n');
    try
        data = load(percentages_file);
        if isfield(data, 'percentages')
            perc_struct = data.percentages;
            if isstruct(perc_struct) && isfield(perc_struct, 'percentages')
                perc = perc_struct.percentages;
            else
                perc = perc_struct;
            end
        else
            % Try to find the first numeric array in the file
            fields = fieldnames(data);
            perc = data.(fields{1});
        end
        fprintf('✓ Archivo cargado exitosamente\n');
    catch ME
        error('Error al cargar el archivo: %s', ME.message);
    end
    
    % Define the mapping based on estimate_newton_ohm.m
    % x(index) -> Variable name
    mapping = {
        % Index, Variable
        1, 'P_m_p';
        2, 'P_h_p';
        3, 'P_p_p';
        4, 'P_s_p';
        5, 'P_e_p';
        6, 'P_b_p';
        7, 'V_m_Total';
        8, 'V_h_Total';
        9, 'V_p_Total';
        10, 'V_s_Total';
        11, 'V_e_Total';
        12, 'V_b_Total';
        13, 'V_vc_Total';
        14, 'V_m_v';
        15, 'V_h_v';
        16, 'V_p_v';
        17, 'V_s_v';
        18, 'V_e_v';
        19, 'V_b_v';
        20, 'V_m_v_un';
        21, 'V_h_v_un';
        22, 'V_p_v_un';
        23, 'V_s_v_un';
        24, 'V_e_v_un';
        25, 'V_b_v_un';
        26, 'V_m_p_un';
        27, 'V_h_p_un';
        28, 'V_p_p_un';
        29, 'V_s_p_un';
        30, 'V_e_p_un';
        31, 'V_b_p_un';
        32, 'Q_m';
        33, 'Q_h';
        34, 'Q_p';
        35, 'Q_s';
        36, 'Q_e';
        37, 'Q_b';
        38, 'Q_vc';
    };
    
    % Create table
    n_params = size(mapping, 1);
    table_data = cell(n_params + 1, 3);
    
    % Headers
    table_data{1, 1} = 'Índice x';
    table_data{1, 2} = 'Variable';
    table_data{1, 3} = 'Porcentaje';
    
    % Fill data
    for i = 1:n_params
        table_data{i + 1, 1} = mapping{i, 1};
        table_data{i + 1, 2} = mapping{i, 2};
        
        % Get percentage value
        idx = mapping{i, 1};
        if idx <= length(perc)
            table_data{i + 1, 3} = perc(idx);
        else
            table_data{i + 1, 3} = NaN;
        end
    end
    
    % Save to Excel
    fprintf('\nGuardando tabla de porcentajes...\n');
    try
        writecell(table_data, output_filename, 'Sheet', 'Mapeo de Porcentajes');
        
        fprintf('✓ Archivo guardado exitosamente: %s\n', output_filename);
        fprintf('  - Tabla de mapeo: Índice x -> Variable -> Porcentaje\n');
        
    catch ME
        error('Error al guardar el archivo Excel: %s', ME.message);
    end
    
    % Display summary
    fprintf('\n=== RESUMEN ===\n');
    fprintf('Total de parámetros: %d\n', n_params);
    fprintf('Longitud del vector de porcentajes: %d\n', length(perc));
    fprintf('\n');
    
end