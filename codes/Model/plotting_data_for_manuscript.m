%plotting_data_for_manuscript
%DATA
% % 
% [t1_N, y1_N, poly1VO2_N, poly1VCO2_N, poly1fiO2_N, ~] = data_preprocessing(1, "normoxia", ".", 0, "ori");
% [t4_N, y4_N, poly4VO2_N, poly4VCO2_N, poly4fiO2_N, ~] = data_preprocessing(4, "normoxia", ".", 0, "ori");
% [t5_N, y5_N, poly5VO2_N, poly5VCO2_N, poly5fiO2_N, ~] = data_preprocessing(5, "normoxia", ".", 0, "ori");
% [t6_N, y6_N, poly6VO2_N, poly6VCO2_N, poly6fiO2_N, ~] = data_preprocessing(6, "normoxia", ".", 0, "ori");
% 
% [t1_H, y1_H, poly1VO2_H, poly1VCO2_H, poly1fiO2_H, ~] = data_preprocessing(1, "hipoxia", "mix", 0, "ori");
% [t4_H, y4_H, poly4VO2_H, poly4VCO2_H, poly4fiO2_H, ~] = data_preprocessing(4, "hipoxia", "mix", 0, "ori");
% [t5_H, y5_H, poly5VO2_H, poly5VCO2_H, poly5fiO2_H, ~] = data_preprocessing(5, "hipoxia", "mix", 0, "ori");
% [t6_H, y6_H, poly6VO2_H, poly6VCO2_H, poly6fiO2_H, ~] = data_preprocessing(6, "hipoxia", "mix", 0, "ori");


% hipoxia = 1;
% normoxia = 1 - hipoxia;
% 
% 
% 
% 
% variables_to_show = [11, 12];
% width = 1;
% height = 2;
% with_mask = 0;
% 
% y_cells_H = {y1_H', y4_H', y5_H', y6_H'};
% t_cells_H = {t1_H', t4_H', t5_H', t6_H'};
% 
% y_cells_N = {y1_N', y4_N', y5_N', y6_N'};
% t_cells_N = {t1_N, t4_N, t5_N, t6_N};
% 
% VO2_cells_H = {poly1VO2_H, poly4VO2_H, poly5VO2_H, poly6VO2_H};
% VCO2_cells_H = {poly1VCO2_H, poly4VCO2_H, poly5VCO2_H, poly6VCO2_H};
% 
% VO2_cells_N = {poly1VO2_N, poly4VO2_N, poly5VO2_N, poly6VO2_N};
% VCO2_cells_N = {poly1VCO2_N, poly4VCO2_N, poly5VCO2_N, poly6VCO2_N};
% 
% VT1_N = {791.00, 705.00, 668.00, 766.00}; %t@VT1
% VT1_H = {2864.00, 2721.00, 2799.00, 3184.00}; %t@VT1
% 
% tmax_N = 2300;
% tmax_H = 3400;
% 
% variable_names = ["dVE", "VT", "TI", "Tresp","PAO2", "PACO2", "HR", "PS", "PD", "PM", "VO2", "VCO2", "Pow"];
% subject_ids = [1,4,5,6];
% units_table = readtable("variables_units.xlsx");

% if hipoxia
%    plot_data(y_cells_H, variables_to_show ,t_cells_H, variable_names , {}, units_table,  VO2_cells_H, VCO2_cells_H, VT1_H, subject_ids, tmax_H, with_mask);
% else
%    plot_data(y_cells_N, variables_to_show ,t_cells_N, variable_names , {}, units_table,  VO2_cells_N, VCO2_cells_N, VT1_N, subject_ids, tmax_N, with_mask);
% end

%average = 0;

%plot_masked_comparison(y_cells_N, t_cells_N, y_cells_H, t_cells_H, variables_to_show, variable_names, units_table, width, height);

%plot_poly(t_cells_H, VCO2_cells_H, subject_ids);

%SIMULATIONS
subject_ids = [1,4,5,6];

% autonomic_vars = ["fh_s", "fp_s", "fv_s", "fv"];
% local_vars = ["xO2_b", "xO2_h", "xO2_rm", "xO2_am"];
% resp_vars = ["HR"];
% gas_vars = ["PACO2"];
% 
% autonomic_colors = [
%             0.6, 0.8, 0.6;   % Verde claro
%             0.3, 0.6, 0.3;   % Verde medio
%             0.1, 0.4, 0.1;   % Verde oscuro
%             0.5, 0.25, 0     % Marrón
%         ];
% local_colors = [
%     0.1216, 0.4667, 0.7059;  % Azul (similar al default de MATLAB)
%     1.0000, 0.4980, 0.0549;  % Naranja
%     0.1725, 0.6275, 0.1725;  % Verde
%     0.8392, 0.1529, 0.1569;  % Rojo
% ];
% common_vars = ["aCO2", "vCO2"];
% title_name = "Gaseous variables and gas concentrations";
% plot_simulations(subject_ids, gas_vars, [], title_name, common_vars);

thresh_arr = [0.6, 0.8, 0.6, 0.5];
plot_LSA(subject_ids, "all", thresh_arr);
%The Sall was computed in the shared_params code


function plot_LSA(subject_ids, plot_type, thresh_arr)
    figure;
    for idx = 1:length(subject_ids)
        subplot(2,2,idx)
        parameter_analysis_fun_plot(subject_ids(idx), plot_type, thresh_arr(idx));
        hold on;
    end

end

function plot_simulations(subjects, vars, colors, title_name, common)

    units_table = readtable("variables_units.xlsx");    
    for idx = 1:length(subjects)

        path_n = sprintf("../Simulations/only_simulation/%d/1200_sec_normoxia-28-07-2025.mat", subjects(idx));
        path_h = sprintf("../Simulations/only_simulation/%d/3300_sec_hipoxia-28-07-2025.mat", subjects(idx));
        
        sim_n = load(path_n);
        sim_h = load(path_h);

        struct_vars_n = sim_n.struct_vars;
        time_n = sim_n.t;

        struct_vars_h = sim_h.struct_vars;
        time_h = sim_h.t;
        
        
        custom_plot("multiple_to_show", {vars, struct_vars_n, time_n, units_table, colors, common}); 
        str_title_n = sprintf("Subject %d: %s during normoxia", subjects(idx), title_name); 
        %title("holi");
        title(char(str_title_n));


        custom_plot("multiple_to_show", {vars, struct_vars_h, time_h, units_table, colors, common}); 
        str_title_h = sprintf("Subject %d: %s during hypoxia", subjects(idx), title_name);
        title(char(str_title_h));
    end
end

function plot_poly(time_cell, vo_cell,subject_ids)
    for index = 1:length(vo_cell)
        subplot(2,2,index);
        bestPolynomialFit(time_cell{index}, vo_cell{index}, 8, "VCO2(l/min)", 2);
        name_title = sprintf("subject %d", subject_ids(index));
        title(name_title);
    end

end









%set(findall(gca, 'Type', 'Line'), 'LineWidth', 0.1);
function plot_data(y_cells, V, t_cells, nombres, together, units_table, vo2_cells, vco2_cells, vt1_times, subject_ids, tmax, with_mask)
% y_cells: celda de matrices VxT, una por sujeto
% V: vector de índices de variables a graficar
% t_cells: celda de vectores de tiempo, uno por sujeto
% nombres: nombres de las variables (cell array de strings)
% together: cell array agrupando variables por subplot
% units_table: tabla con columnas 'Variable' y 'MeasureUnit'
% vo2_cells, vco2_cells: no se usan, mantenidos por compatibilidad
% vt1_times: celda con tiempos VT1 por sujeto (en segundos)
% subject_ids: vector con los identificadores reales de cada sujeto (ej: [1 4 5 6])

    if ~iscell(y_cells)
        error("y_cells debe ser una celda de matrices VxT.");
    end
    if ~iscell(t_cells)
        error("t_cells debe ser una celda de vectores de tiempo.");
    end
    if length(y_cells) ~= length(t_cells)
        error("y_cells y t_cells deben tener el mismo largo.");
    end
    if isempty(together)
        together = num2cell(V(:)');
    end
    if nargin < 10 || isempty(subject_ids)
        subject_ids = 1:length(y_cells);  % valores por defecto si no se entrega
    end
    if length(subject_ids) ~= length(y_cells)
        error("subject_ids debe tener la misma longitud que y_cells.");
    end

    T_units = units_table;
    num_subjects = length(y_cells);
    num_plots = numel(together);

    base_colores = lines(num_subjects);
    line_styles = {'-', '--', ':', '-.'};
    y1 = y_cells{1};
    t1 = t_cells{1};

    power_series = y1(13, :);
    time_for_power = t1;
    figure;
    for i = 1:num_plots
        subplot(num_plots, 1, i);
        hold on;

        vars = together{i};
        num_vars = length(vars);

        if num_vars > numel(line_styles)
            warning('Más variables que estilos de línea disponibles; algunos estilos se repetirán.');
        end

        for s = 1:num_subjects
            y = y_cells{s};
            t = t_cells{s};
            subj_id = subject_ids(s);

            % Obtener máscara (última columna)
            mask = y(end - 1, :);  % vector lógico por variable (una fila por var)
            mask = logical(mask);  % asegúrate de que sea lógico

            for v = 1:num_vars
                var_idx = vars(v);
                style_idx = mod(v-1, length(line_styles)) + 1;
                style = line_styles{style_idx};

                legend_name = sprintf('Sujeto %d', subj_id);

                % Datos originales
                data = y(var_idx, 1:end);
                time = t;

                % Aplicar máscara (pon NaN donde mask = 0)
                %valid_mask = logical(mask(var_idx));
                if with_mask
                    data(~mask) = NaN;   
                end

                plot(time, data', ...
                     'Color', base_colores(s,:), ...
                     'LineStyle', style, ...
                     'LineWidth', 1.2, ...
                     'DisplayName', legend_name);
            end

            

            % Línea vertical en tiempo de VT1 (sin leyenda)
            if ~isempty(vt1_times) && ~isempty(vt1_times{s}) && ~isnan(vt1_times{s})
                t_vt1 = vt1_times{s};

                xline(t_vt1, '--', num2str(subj_id), ...
                      'Color', base_colores(s,:), ...
                      'LabelOrientation', 'horizontal', ...
                      'Alpha', 0.8, ...
                      'LineWidth', 1.5, ...
                      'HandleVisibility', 'off');
            end

            legend('show', 'Location', 'northwest');
        end

        yyaxis right
        ax = gca;
        ax.YColor = 'k';  % Cambia el color del eje derecho a negro
        plot(time_for_power, power_series, 'Color', [0.3 0.3 0.3], 'DisplayName', "power-ladder");
        ylabel('W')  % o la unidad que corresponda
        yyaxis left  % volver al eje original para futuros plots


        hold off;
        legend('show');
        grid on;

        if isfinite(tmax)
            xlim([0 tmax]);
        end

        nombre_vars = nombres(vars);
        unidades = strings(size(vars));
        for k = 1:length(vars)
            try
                unidades(k) = find_unit(T_units, nombre_vars{k});
            catch
                unidades(k) = "";
            end
        end

        if all(unidades == unidades(1)) && unidades(1) ~= ""
            ylabel(unidades(1));
        else
            ylabel('Valor');
        end

        title(strjoin(nombre_vars, ', '));
        xlabel('Tiempo [s]');
    end
    
end


function plot_masked_comparison(y_cells_N, t_cells_N, y_cells_H, t_cells_H, variables_to_show, variable_names, units_table, width, height)
    % Colores por condición
    color_N = [0.00, 0.45, 0.74]; % azul normoxia
    color_H = [0.85, 0.33, 0.10]; % rojo hipoxia

    num_subjects = length(y_cells_N);
    if length(y_cells_H) ~= num_subjects
        error("Cantidad de sujetos en normoxia y hipoxia no coincide.");
    end

    % === Interpolar y alinear ===
    % Tiempo de referencia: primer t de normoxia sujeto 1
    

    % Crear malla común de tiempo (hasta 1s después del máximo compartido)
    t_grid = 0:1:3000; % puedes ajustar este rango si quieres
    

    % === Normoxia ===
   
    num_vars = length(variables_to_show);
    figure; 
    for v = 1:num_vars
        var_idx = variables_to_show(v);
        subplot(width, height, v);
        hold on;
        Y_interp_N = nan(num_subjects, length(t_grid));
        for i = 1:num_subjects
            y = y_cells_N{i};
            t = t_cells_N{i};
            mask = logical(y(end, :));
            mask = mask & logical(y(end - 1,:));
            t_masked = t(mask);
            t_masked = t_masked - t_masked(1);
            y_masked = y(var_idx, mask);
            
    
            [t_unique, ia] = unique(t_masked);       % <-- Evita duplicados
            y_unique = y_masked(ia);
    
            Y_interp_N(i, :) = interp1(t_unique, y_unique, t_grid, 'linear', NaN);
        end
        mean_N = mean(Y_interp_N, 1, 'omitnan');
        std_N  = std(Y_interp_N, 0, 1, 'omitnan');
    
        mean_N = mean_N(~isnan(mean_N));
        std_N = std_N(~isnan(mean_N));
        t_grid_N = t_grid(~isnan(mean_N));
    
        
    
    
    
        % === Hipoxia ===
        Y_interp_H = nan(num_subjects, length(t_grid));
        for i = 1:num_subjects
            y = y_cells_H{i};
            t = t_cells_H{i};
            mask = logical(y(end, :));
            power = y(end - 2, :);
            t_masked = t(mask);
            t_masked = t_masked - t_masked(1);
            y_masked = y(var_idx, mask);
            power_masked = power(mask);
           
    
            [t_unique, ia] = unique(t_masked);       % <-- Evita duplicados
            y_unique = y_masked(ia);
            power_unique = power_masked(ia);
        
            Y_interp_H(i, :) = interp1(t_unique, y_unique, t_grid, 'linear');
            power_interp = interp1(t_unique, power_unique, t_grid, 'linear');
        end
        mean_H = mean(Y_interp_H, 1, 'omitnan');
        std_H  = std(Y_interp_H, 0, 1, 'omitnan');
    
        mean_H = mean_H(~isnan(mean_H));
        std_H = std_H(~isnan(mean_H));
        t_grid_H = t_grid(~isnan(mean_H));
        power_H = power_interp(~isnan(mean_H));
        
    
    
        % === Plot ===
        %figure;
        %hold on;
    
        % Normoxia: sombreado ±1 std
        yyaxis left

        fill([t_grid_N, fliplr(t_grid_N)], ...
             [mean_N + std_N, fliplr(mean_N - std_N)], ...
            color_N, 'FaceAlpha', 0.3, 'EdgeColor', 'none',  'DisplayName', 'std Normoxia');
        
        plot(t_grid_N, mean_N, 'Color', color_N, 'LineWidth', 2, ...
             'DisplayName', 'Normoxia', 'LineStyle','-');
    
        % Hipoxia: sombreado ±1 std
        fill([t_grid_H, fliplr(t_grid_H)], ...
             [mean_H + std_H, fliplr(mean_H - std_H)], ...
             color_H, 'FaceAlpha', 0.2, 'EdgeColor', 'none',  'DisplayName', 'std Hipoxia');
        plot(t_grid_H, mean_H, 'Color', color_H, 'LineWidth', 2, ...
             'DisplayName', 'Hipoxia', 'LineStyle','-');
    
        yyaxis right        
        plot(t_grid_H, power_H, 'DisplayName', 'Potencia', 'Color', 'black');        
        ax = gca;
        ax.YColor = [0 0 0];
        
    
        xlim([0, 1000]);
    
        % === Formato gráfico ===
        xlabel('Tiempo alineado [s]');
        ylabel(sprintf('Variable %d', var_idx));
        title(sprintf('Comparación promedio ± std (var %d)', var_idx));
        legend('Location', 'best');
        grid on;
    
        var_name = variable_names(var_idx);
        try
            unidad = find_unit(units_table, var_name);
        catch
            unidad = "";
        end
        yyaxis left
        ylabel(unidad);
        title(var_name);
        legend('Location', 'best');
        ax = gca;
        ax.YColor = [0 0 0];
        grid on;
        yyaxis right
        ylabel('W')

   end
end





function [t_common, avg_series] = average_series(y_cells, t_cells, var_idx, with_mask, tmax)
    % Encuentra un tiempo común interpolado (ej: 1000 puntos)
    t_common = linspace(0, tmax, 1000);
    all_series = zeros(length(y_cells), length(t_common));

    for i = 1:length(y_cells)
        y = y_cells{i}(var_idx,:);
        t = t_cells{i};
        if with_mask
            mask = logical(y_cells{i}(end-1,:));
            y(~mask) = NaN;
        end
        y_interp = interp1(t, y, t_common, 'linear', NaN);
        all_series(i,:) = y_interp;
    end
    avg_series = nanmean(all_series, 1);
end








function unit = find_unit(table, var)
    row_index = strcmp(table.Variable, var);
    unit = table.MeasureUnit{row_index};    
end