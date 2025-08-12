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
 elseif strcmp(type_of_plot, 'sim_vs_exp')
    %figure;  
    t_sim = custom_args{1};
    t_exp = custom_args{2};
    struct_vars = custom_args{3};
    X_exp = custom_args{4};
    Power = X_exp(:,11);
    X_exp = X_exp(:, 1:10);    
    list_of_names = custom_args{5};
    units_table = custom_args{6};
    width = custom_args{7};
    height = custom_args{8};
    simulation_filename = custom_args{9};
    old_mode = custom_args{10};

    

    if strcmp(old_mode, 'only-exp')
        plot_simulated_vs_experimental(t_sim, t_exp, [0], X_exp, list_of_names, units_table, width, height);
    else

    
        if ~strcmp(simulation_filename, '')
           loaded_simulation = load(simulation_filename);
           cell_l_s = fieldnames(loaded_simulation);
           struct_vars_name = cell_l_s{1};
           t_name = cell_l_s{2};

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
        
        plot_simulated_vs_experimental(t_sim, t_exp, X_sim, X_exp, list_of_names, units_table, width, height, Power);
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

function plot_simulated_vs_experimental(t_sim, t_exp, X_sim, X_exp, Xnames, Xunits, width, height, Power)
    t_exp = t_exp;
    x_exp = X_exp;
    x_sim = X_sim;
    variable_names = Xnames;
    units_table = Xunits;
    % Create the figure and 5x2 grid of subplots
    %figure;
    
    % Loop through each variable
    for i = 1:size(x_exp, 2)
        subplot(width, height, i); % Define the subplot position
        hold on;
        %yyaxis left
        plot(t_exp, x_exp(:, i)); % Plot the variable against time
        %yyaxis right
        %plot(t_exp, Power);
        %yyaxis left
        xlim([0 Inf]);
        ylim([0 Inf]);
        if iscell(variable_names)
            title(variable_names{i}); % Set the title to the variable name
        else
            title(variable_names(i)); 
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
        
        plot(t_sim, x_sim(:, i)); % Plot the variable against time
        grid on; % Add grid for better visualization
    end
    
    for i = 1:size(x_sim, 2)
        subplot(width, height, i);
        hold on;
        % Manual yyaxis implementation for MATLAB 2017
        ax1 = gca;
        ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', ...
                   'YAxisLocation', 'right', 'XTick', [], 'Box', 'off');
        axes(ax2);
        plot(t_exp, Power, 'Color', [0.5 0.5 0.5]);
        set(ax2, 'YColor', [0.5 0.5 0.5]);
        axes(ax1); % Return to original axis
    end

    
    % Adjust the layout to prevent overlap
    % sgtitle not available in 2017, use manual title
    annotation('textbox', [0 0.9 1 0.1], 'String', 'Sim vs Exp', ...
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

function unit = find_unit(table, var)
    if iscell(var)
        var = var{1};
    end
    row_index = strcmp(table.Variable, var);
    unit = table.MeasureUnit{row_index};    
end

end