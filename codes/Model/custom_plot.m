function custom_plot(type_of_plot, custom_args)       

 if strcmp(type_of_plot, 'vars_to_show')
    
    t = custom_args{3};
    struct_vars = custom_args{2};
    vars_to_show = custom_args{1};
    units_table = custom_args{4};

    %  vars_to_show = list_of_names;  %With 2 options, both variables will be plotted with their own axis, with more, there will be just plot's superposition.
    plot_vars(t, struct_vars, vars_to_show, units_table); 
    %figure; 
    %IF YOU WANT TO CONTINUE PLOTTING
    % vars_to_show = {" ", " "};
    % plot_vars(t, struct_vars, vars_to_show, units_table);
    % figure;
 elseif strcmp(type_of_plot, 'optim-boundries')
     lower_bounds = custom_args{1};
     initial_values = custom_args{2};
     upper_bounds = custom_args{3};
     final_values = custom_args{4};
     names = custom_args{5};
     boundary_saturation_table(lower_bounds, initial_values, upper_bounds, final_values, names);

     
 
 elseif strcmp(type_of_plot, 'sim_vs_exp')
     disp('aloha');
    %figure;  
    t_sim = custom_args{1};
    t_exp = custom_args{2};
    struct_vars = custom_args{3};
    X_exp = custom_args{4};
    Power = X_exp(:,11);
    notnan_mask = X_exp(:, 12);
    X_exp = X_exp(:, 1:10);    
    list_of_names = custom_args{5};
    units_table = custom_args{6};
    width = custom_args{7};
    height = custom_args{8};
    simulation_filename = custom_args{9};
    old_mode = custom_args{10};
    patient_number = custom_args{11};
    fitting_status = custom_args{12};
    J = custom_args{13};
    state = custom_args{14};
    
    if numel(custom_args) > 14
        reverse = custom_args{15};
    else
        reverse = 0;
    end

    

    if strcmp(old_mode, 'only-exp')
        plot_simulated_vs_experimental(t_sim, t_exp, [0], X_exp, list_of_names, units_table, width, height, state);
    else

    
        if ~strcmp(simulation_filename, '')
           loaded_simulation = load(simulation_filename);
           cell_l_s = fieldnames(loaded_simulation);
           cell_1_name = cell_l_s{1};
           cell_2_name = cell_l_s{2};
           try
               cell_3_name = cell_l_s{3};
           catch
               cell_3_name = 0;
           end

           if strcmp(class(loaded_simulation.(cell_1_name)), 'struct')
               struct_vars_name = cell_1_name;
               t_name = cell_2_name;
           elseif numel(loaded_simulation.(cell_1_name)) > 10
               t_name = cell_1_name;
               if strcmp(class(loaded_simulation.(cell_2_name)), 'struct')
                   struct_vars_name = cell_2_name;
                   J_name = cell_3_name;
               else
                   struct_vars_name = cell_3_name;
                   J_name = cell_2_name;
               end
           
           else
               J_name = cell_1_name;
               if strcmp(class(loaded_simulation.(cell_2_name)), 'struct')
                   struct_vars_name = cell_2_name;
                   t_name = cell_3_name;
               else
                   struct_vars_name = cell_3_name;
                   t_name = cell_2_name;
               end

               
           end

           struct_vars = loaded_simulation.(struct_vars_name);
           t_sim = loaded_simulation.(t_name);
        end
        
        if strcmp(old_mode, 'on')
            X_sim = [struct_vars.dVE', struct_vars.V', struct_vars.TI', struct_vars.Tresp', struct_vars.PAO2', struct_vars.PACO2', 60 * struct_vars.HR', struct_vars.P_sa', struct_vars.P_sa', struct_vars.mean_P_sa'];    
            %time = linspace(0, 1260, size(struct_vars.dVE,2));
            time = t_sim;
            X_sim(:,2) = data_processing("volume", X_sim(:, 2), time'); %volume
            [pm, ps, pd] = data_processing("pressure", X_sim(:, 8), time');
            X_sim(:, 8) = ps;
            X_sim(:, 9) = pd;
            
            X_exp(:, 7) = X_exp(:, 7)*60;
        else
            X_sim = [struct_vars.dVE', struct_vars.VT', struct_vars.TI', struct_vars.Tresp', struct_vars.PAO2', struct_vars.PACO2', struct_vars.HR' * 60, struct_vars.ps', struct_vars.pd', struct_vars.pm'];    
            X_exp(:, 7) = X_exp(:, 7)*60;
            
        end    

        t_exp = t_exp(t_exp <= t_sim(end));
        X_exp = X_exp(t_exp <= t_sim(end),:);
        Power = Power(t_exp <= t_sim(end),:);
        notnan_mask = notnan_mask(t_exp <= t_sim(end),:);

        MRtCO2 = struct_vars.MRtCO2;
        MRtO2 = struct_vars.MRtO2;

        
        plot_simulated_vs_experimental(t_sim, t_exp, X_sim, X_exp, list_of_names, units_table, width, height, Power, patient_number, fitting_status,J, MRtO2, MRtCO2, reverse, state, notnan_mask);
    end
elseif strcmp(type_of_plot,'LSA-plot')
    
    sens_matrix = custom_args{1};
    pars_to_sens = custom_args{2};   
    variables_of_interest = custom_args{3}; 
    idx_variable_of_interest = custom_args{4};
    if length(custom_args) < 5
        patient_idx = 'all';
        title_text = 'All: Local Sensitivity Analysis integrated over time (sum sqr)';
    else 
        patient_idx = custom_args{5};
        title_text = sprintf('Subject %d: Local Sensitivity Analysis integrated over time (sum sqr)', patient_idx);
    end
    
    
   
    %figure;
    % Convert cell array to categorical if needed
    if iscell(pars_to_sens)
        pars_categorical = categorical(pars_to_sens);
    else
        pars_categorical = categorical(cellstr(pars_to_sens));
    end
    
    b = bar(pars_categorical, sens_matrix, 'stacked'); 
    
    num_variables = size(sens_matrix, 2);      
    colors = jet(num_variables);  % Change 'parula' to 'jet', 'hsv', etc., if preferred    
    
    for i = 1:num_variables
        set(b(i), 'FaceColor', colors(i, :));
    end
    
    
    % Title, labels, and grid
    
    title(title_text);
    xlabel('Parameters');
    ylabel('Sensitivities');
    grid on;
    
    % Add legend to identify variables
    if iscell(variables_of_interest)
        legend_vars = variables_of_interest(idx_variable_of_interest);
    else
        legend_vars = cellstr(variables_of_interest(idx_variable_of_interest));
    end
    legend(legend_vars, 'Location', 'best');

elseif strcmp(type_of_plot,'ident-plot')  % Fixed typo: strmp -> strcmp

    correlation_matrix = custom_args{1};
    pars_to_sens_red = custom_args{2};    
    if length(custom_args) < 4
        patient_idx = 'all';
        title_text = 'All: Correlation between parameters';
    else
        patient_idx = custom_args{4};
        title_text = sprintf('Subject %d: Correlation between parameters', patient_idx);
    end
    

    %figure;
    imagesc(abs(correlation_matrix));

    colorbar; % Add a color scale
    %caxis([-1, 1]); % Normalize between -1 and 1
    colormap jet; % Use a color map for better visualization

    % Set axis labels
    set(gca, 'XTick', 1:length(pars_to_sens_red));
    set(gca, 'YTick', 1:length(pars_to_sens_red));
    if iscell(pars_to_sens_red)
        set(gca, 'XTickLabel', pars_to_sens_red);
        set(gca, 'YTickLabel', pars_to_sens_red);
    else
        set(gca, 'XTickLabel', cellstr(pars_to_sens_red));
        set(gca, 'YTickLabel', cellstr(pars_to_sens_red));
    end
    % Rotate x-axis labels (manual approach for 2017)
    hx = get(gca, 'XLabel');
    set(hx, 'Units', 'data');
    h_labels = get(gca, 'XTickLabel');
    delete(get(gca, 'XLabel'));
    for i = 1:length(h_labels)
        text(i, 0.5, h_labels{i}, 'HorizontalAlignment', 'right', 'Rotation', 45, ...
             'Units', 'data', 'VerticalAlignment', 'top');
    end
    
    title(title_text);
    

elseif strcmp(type_of_plot,'multiple_to_show')
    t = custom_args{3};
    struct_vars = custom_args{2};
    vars_to_show = custom_args{1};
    units_table = custom_args{4};
    colors = custom_args{5};
    common = custom_args{6};
    multiple_to_show = vars_to_show;
    %common = {"MRtO2", "MRtCO2"};
    plot_mutiple_vs_common(t, struct_vars,units_table, multiple_to_show, common, colors);  
 
% elseif strcmp(type_of_plot, "same_units")
%     show_same_units_vars(t, struct_vars, units_table, 'blood_pressure');
%  elseif strcmp(type_of_plot, "interest_variables")
%     show_interest_variables(t, struct_vars, units_table);
 end
 


%% FUNCTIONS

%Plotting

    function plot_simulated_vs_experimental(t_sim, t_exp, X_sim, X_exp, Xnames, Xunits, width, height, Power, patient_number, fitting_status, J, MRtO2, MRtCO2, reverse, state, notnan_mask)
    t_exp = t_exp;
    x_exp = X_exp;
    x_sim = X_sim;
    variable_names = Xnames;
    units_table = Xunits;

    
    % Create the figure and 5x2 grid of subplots
    %figure;
    JJ = [J(3) J(4) J(5) J(6) J(1) J(2) J(7) J(10) J(8) J(9)]; %small adjustment because order discrepancies
    if reverse
        ww = width;
        width = height;
        height = ww;
    end
    
    % Loop through each variable
    for i = 1:size(x_exp, 2)
        subplot(width, height, i); % Define the subplot position
        hold on;
        %yyaxis left
        if i >= 8
            var = x_exp(logical(notnan_mask),i);
            time_ = t_exp(logical(notnan_mask));
        else
            if patient_number == 1 && strcmp(state, 'Normoxia') && i >= 7
                var = x_exp(logical(notnan_mask),i);
                time_ = t_exp(logical(notnan_mask));
            else
                var = x_exp(:,i);
                time_ = t_exp;                
            end

        end
        yyaxis right
        plot(t_exp, Power, 'Color', 'black', 'LineWidth', 1.5); % mostaza
        ylabel('Power (W)');
    
        yyaxis left
        
        plot(time_, var, 'Color', [0.4, 0.6, 0.8], 'LineWidth', 1, ...
             'LineStyle', '-'); % l�nea s�lida azul suave
        
        % --- Configuraci�n de ejes uniformes ---
        ax = gca;
        ax.YColor = 'k'; % eje izquierdo negro
        ax.YAxis(2).Color = 'k'; % eje derecho negro
        ax.XColor = 'k';
        ax.Box = 'on'; % mantiene bordes visibles
        
        xlim([0 Inf]);
        ylim([0 Inf]);
        if iscell(variable_names)             
            %title_ = sprintf('%s, J: %.3f', variable_names{i}, JJ(i));% Set the title to the variable name
            title_ = sprintf('%s', variable_names{i});% Set the title to the variable name
            title(title_); 
        else
            %title_ = sprintf('%s:%d', variable_names(i), JJ(i));
            title_ = sprintf('%s', variable_names(i));
            title(title_); 
        end
        
        xlabel('Time (s)'); % Label the x-axis
        unit = find_unit(units_table, variable_names{i});
        if iscell(variable_names)
            yl = strcat(variable_names{i},'(', unit, ')');
        else
            yl = strcat(variable_names(i),'(', unit, ')');
        end
        ylabel(yl); % Label the y-axis with the variable name
        grid on; % Add grid for better visualization
    end

    for i = 1:size(x_sim, 2)
        subplot(width, height, i); % Define the subplot position
        hold on;
        
        if strcmp(fitting_status, 'actual')
           plot(t_sim, x_sim(:, i),  'Color', [1.0, 0.5, 0.0], 'LineWidth', 1, 'LineStyle', '-', 'DisplayName', 'actual'); % Plot the variable against time
            
        elseif strcmp(fitting_status, 'old')
            plot(t_sim, x_sim(:, i), 'Color', [0.5, 0.0, 0.6], 'LineWidth', 1, 'LineStyle', '-', 'DisplayName', 'old'); % Plot the variable against time
        else
            plot(t_sim, x_sim(:, i), 'DisplayName', 'data'); 

        end
        grid on; % Add grid for better visualization
    end
    for i = 1:size(x_sim, 2)
        subplot(width, height, i);
        %lgd = legend('show');
        %lgd.Position = [0.8 0.9 0.15 0.1]; 
    end
    
    
    %for i = 1:size(x_sim, 2)
    %    subplot(width, height, i);
    %    hold on;
        % Manual yyaxis implementation for MATLAB 2017
        %ax1 = gca;
        %ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', ...
                   %'YAxisLocation', 'right', 'XTick', [], 'Box', 'off');
        %axes(ax2);
    %    yyaxis right;
    %    plot(t_sim, MRtO2, 'Color', [0.9, 0.5, 0.5]);
    %    set(gca, 'YColor', 'red');
        %hold on;
        %plot(t_exp, MRtCO2);

        %set(ax2, 'YColor', [0.5 0.5 0.5]);
        %axes(ax1); % Return to original axis
        %yyaxis left;
    %end

    
    % Adjust the layout to prevent overlap
    % sgtitle not available in 2017, use manual title
    title_name = sprintf('Sim vs Exp: S%d - %s', patient_number, state);
    annotation('textbox', [0 0.9 1 0.1], 'String', title_name, ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14);
end

function plot_vars(t,rt, var_names,units_table)
    
    if length(var_names) == 2
    
       % Manual yyaxis implementation for MATLAB 2017
       ax1 = gca;
       plot(t, real(rt.(var_names{1})));
       unit1 = find_unit(units_table, var_names{1});
       ylabel(unit1);
       set(ax1, 'YColor', 'b');

       ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', ...
                  'YAxisLocation', 'right', 'XTick', [], 'Box', 'off');
       axes(ax2);
       plot(t, real(rt.(var_names{2})), 'r');
       unit2 = find_unit(units_table, var_names{2});
       ylabel(unit2);
       set(ax2, 'YColor', 'r');

    else
        for i = 1: length(var_names)
            if iscell(var_names)
                plot(t, real(rt.(var_names{i})))
            else
                plot(t, real(rt.(var_names(i))))
            end
            hold on 
        end
       
    end 
    title('Plot');
    if iscell(var_names)
        legend_labels = var_names;
    else
        legend_labels = cellstr(var_names);
    end
    legend(legend_labels);
    % Replace underscores in legend
    legend_h = legend;
    legend_str = get(legend_h, 'String');
    for i = 1:length(legend_str)
        legend_str{i} = strrep(legend_str{i}, '_', ' ');
    end
    set(legend_h, 'String', legend_str);
    xlabel('time');
 
    
end

function plot_mutiple_vs_common(t, rt,units_table, multiple_to_show, common, colors_)
     
    figure;
    if isempty(colors_)
        colors_ = lines(length(multiple_to_show));
    end
    
    % Manual yyaxis left
    ax1 = gca;
    for var_idx = 1:length(multiple_to_show)        
        if iscell(multiple_to_show)
            plot(t, real(rt.(multiple_to_show{var_idx})), '-', 'Color', colors_(var_idx, :));
            unit = find_unit(units_table, multiple_to_show{var_idx});
        else
            plot(t, real(rt.(multiple_to_show(var_idx))), '-', 'Color', colors_(var_idx, :));
            unit = find_unit(units_table, multiple_to_show(var_idx));
        end
        ylabel(unit);
        ylim([0, 70]);
        hold on;
    end
    xlabel('t(s)');
    set(ax1, 'YColor', 'black');

    % Manual yyaxis right
    ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', ...
               'YAxisLocation', 'right', 'XTick', [], 'Box', 'off');
    axes(ax2);
    for common_idx = 1:length(common)
        if iscell(common)
            plot(t, real(rt.(common{common_idx})), 'Color', 'black');
        else
            plot(t, real(rt.(common(common_idx))), 'Color', 'black');
        end
        hold on;
    end
    
    % Create combined legend
    if iscell(multiple_to_show) && iscell(common)
        legend_vars = [multiple_to_show, common];
    elseif iscell(multiple_to_show) && ~iscell(common)
        legend_vars = [multiple_to_show, cellstr(common)];
    elseif ~iscell(multiple_to_show) && iscell(common)
        legend_vars = [cellstr(multiple_to_show), common];
    else
        legend_vars = [cellstr(multiple_to_show), cellstr(common)];
    end
    legend(legend_vars);
    
    if iscell(common)
        unit = find_unit(units_table, common{1});
    else
        unit = find_unit(units_table, common(1));
    end
    ylabel(unit);
    ylim([0,2]);
    set(ax2, 'YColor', 'black');
end

function show_same_units_vars(t, struct_vars, units_table, label)

    if strcmp(label, 'blood_pressure')
        vars_to_show = {'P_sa', 'P_sp'};
        plot_vars(t, struct_vars, vars_to_show, units_table);
        figure;
    elseif strcmp(label, 'O2')
        vars_to_show = {'P_1O2', 'P_2O2', 'P_3O2', 'P_4O2', 'P_5O2', 'PAO2', 'PaO2'};
        plot_vars(t, struct_vars, vars_to_show, units_table);
        figure;
    elseif strcmp(label, 'CO2')
        vars_to_show = {'P_1CO2', 'P_2CO2', 'P_3CO2', 'P_4CO2', 'P_5CO2', 'PACO2', 'PaCO2', 'mean_PbCO2', 'PvbCO2'};
        plot_vars(t, struct_vars, vars_to_show, units_table); 
        figure;  
    elseif strcmp(label, 'volume')
        vars_to_show = {'V_total_e_v', 'V_total_s_v', 'V_total_b_v', 'V_total_h_v', 'V_total_rm_v', 'V_total_am_v', 'V_total_vc', 'V_total_lv', 'V_total_la', 'V_total_rv', 'V_total_ra', 'V_total_pa', 'V_total_pp', 'V_total_pv'};
        plot_vars(t, struct_vars, vars_to_show, units_table);
        figure;
    elseif strcmp(label, 'flow')
        vars_to_show = {'Q_sa', 'Q_pa'};
        plot_vars(t, struct_vars, vars_to_show, units_table);
        figure;
    end   


end

function show_interest_variables(t, struct_vars, units_table)
    arteries = {'Q_sa', 'P_sa'};
    atriums = {'V_total_ra', 'V_total_la'};
    ventricles = {'V_total_lv', 'V_total_rv'};
    lungs = {'V_total_pa', 'V_total_pv', 'V_total_pp'};
    systemic_venous = {'V_total_e_v','V_total_b_v', 'V_total_rm_v', 'V_total_am_v', 'V_total_h_v', 'V_total_s_v'};
    %CHANGE THE ONES YOU WANT

    
    
    plot_vars(t, struct_vars, arteries, units_table);
    title('Arteries');
    figure;      
    
    plot_vars(t, struct_vars, atriums, units_table);
    title('Atriums');
    figure;
    
    plot_vars(t, struct_vars, ventricles, units_table);
    title('Ventricles');
    figure;
    
    plot_vars(t, struct_vars, lungs, units_table);
    title('Lungs');
    figure;    

    plot_vars(t, struct_vars, systemic_venous, units_table);
    title('Systemic Venous');
    
    

end

    function boundary_saturation_table(lower_bounds, initial_values, upper_bounds, final_values, param_names)
% BOUNDARY_SATURATION_TABLE - Genera tabla con análisis de saturación
% Entradas:
%   lower_bounds - vector con límites inferiores
%   initial_values - vector con valores iniciales  
%   upper_bounds - vector con límites superiores
%   final_values - vector con valores finales
%   param_names - (opcional) cell array o string array con nombres de parámetros

% Datos de ejemplo si no se proporcionan entradas
if nargin < 4
    lower_bounds = [8.7000 22.5000 0.5687 0.1998 0.0515 0.2000 23.8900 0.0650 3.5300 10.0000];
    initial_values = [17.4000 45.0000 1.1375 0.3997 0.1030 0.4000 47.7800 0.1300 7.0600 20.0000];
    upper_bounds = [174.0000 450.0000 11.3750 3.9970 1.0300 4.0000 477.8000 1.3000 70.6000 200.0000];
    final_values = [15.6313 450.0000 3.4914 0.1998 0.0515 3.9075 477.8000 0.3145 66.8010 199.9710];
    param_names = {'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'kappa', 'lambda'};
end

n_params = length(lower_bounds);

% Crear nombres de parámetros si no se proporcionan
if nargin < 5 || isempty(param_names)
    param_names = arrayfun(@(x) sprintf('P%d', x), 1:n_params, 'UniformOutput', false);
else
    % Convertir a cell array si es necesario
    if ischar(param_names) || isstring(param_names)
        param_names = cellstr(param_names);
    end
    % Verificar que el tamaño coincida
    if length(param_names) ~= n_params
        warning('El número de nombres de parámetros no coincide con el número de parámetros. Usando nombres por defecto.');
        param_names = arrayfun(@(x) sprintf('P%d', x), 1:n_params, 'UniformOutput', false);
    end
end

% Detectar saturación (tolerancia ajustable)
tolerance = 1e-6;
saturated_lower = abs(final_values - lower_bounds) <= tolerance;
saturated_upper = abs(final_values - upper_bounds) <= tolerance;

% Crear estado de saturación
saturation_status = cell(n_params, 1);
saturation_status(~saturated_lower & ~saturated_upper) = {'OK'};
saturation_status(saturated_lower) = {'SATURADO_INFERIOR'};
saturation_status(saturated_upper) = {'SATURADO_SUPERIOR'};

% Crear la tabla
T = table(param_names', lower_bounds', initial_values', upper_bounds', final_values', saturation_status, ...
          'VariableNames', {'Parametro', 'Limite_Inferior', 'Valor_Inicial', 'Limite_Superior', 'Valor_Final', 'Estado'});

% Mostrar la tabla
disp(' ');
disp('=============== AN�?LISIS DE SATURACIÓN DE PAR�?METROS ===============');
disp(T);

% Mostrar resumen estadístico
fprintf('\n==================== RESUMEN ====================\n');
fprintf('Total de parámetros: %d\n', n_params);
fprintf('Parámetros OK: %d\n', sum(~saturated_lower & ~saturated_upper));
fprintf('Saturados en límite inferior: %d\n', sum(saturated_lower));
fprintf('Saturados en límite superior: %d\n', sum(saturated_upper));

if sum(saturated_lower) > 0
    fprintf('\nParámetros saturados en límite INFERIOR:\n');
    sat_lower_indices = find(saturated_lower);
    for i = 1:length(sat_lower_indices)
        idx = sat_lower_indices(i);
        fprintf('  %s: %.6f (límite: %.6f)\n', param_names{idx}, final_values(idx), lower_bounds(idx));
    end
end

if sum(saturated_upper) > 0
    fprintf('\nParámetros saturados en límite SUPERIOR:\n');
    sat_upper_indices = find(saturated_upper);
    for i = 1:length(sat_upper_indices)
        idx = sat_upper_indices(i);
        fprintf('  %s: %.6f (límite: %.6f)\n', param_names{idx}, final_values(idx), upper_bounds(idx));
    end
end

fprintf('=====================================================\n');

end


function unit = find_unit(table, var)
    if iscell(var)
        var = var{1};
    end
    row_index = strcmp(table.Variable, var);
    unit = table.MeasureUnit{row_index};    
end

end