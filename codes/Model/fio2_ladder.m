% 
% 
% vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");
% 
% patient_idx = 5;
% 
% hipoxic_iterations = 1:20;
% 
% fiO2_init = 0.21;
% fiO2_input = fiO2_init - hipoxic_iterations * 0.005;
% all_vars = cell(1, length(hipoxic_iterations));
% vars_matrix = zeros(length(hipoxic_iterations), 10); 
% 
% for hi = hipoxic_iterations    
%     [setup] = set_up("fiO2_ladder", patient_idx, "hipoxia", "mix", "dt", 0.1,"simulation_time", 200, "settling_time", 40);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%     s = setup;
%     s.pars("MRO2") = s.VO2_ladder_points(2)/1000; %from 2 - 6 points
%     s.pars("MRCO2") = s.VCO2_ladder_points(2)/1000; %from 2 - 6 points
%     %s.pars("tauMR") = 100;
%     %s.pars("tauMRv") = 100;
%     s.pars("fO2") = fiO2_input(hi) * 100;      
% 
%     [t, x_dot, x_vars, x_keys, index] = s.run_ode_fun(s.model, s.pars, s.init, s.simulation_time, s.dt);
%     struct_vars = arrange_results(x_dot, x_vars, x_keys, t);
%     all_vars{hi} = struct_vars;
% 
%     pout = data_processing("pressure", struct_vars.P_sa, t);
%     vout = data_processing("volume", struct_vars.V, t);
% 
%     vars_of_interest = [struct_vars.dVE', vout', struct_vars.TI', struct_vars.Tresp', struct_vars.PAO2', struct_vars.PACO2', 60 * struct_vars.HR', pout{2}', pout{3}', struct_vars.mean_P_sa'];    
%     final_values = mean(vars_of_interest(end-500:end,:),1);
%     %final_max_values = maxk(vars_of_interest(end-500:end,:), 5);
%     %final_min_values = mink(vars_of_interest(end-500:end,:), 5);
%     %final_values = mean([final_max_values; final_min_values], 1);  % Average of max and min for stability
%     vars_matrix(hi,:) = final_values';
% end
% 
% % Assume: vars_matrix exists in workspace from current run
% [N, M] = size(vars_matrix);
% 
% % Create subplot layout (only once)
% rows = ceil(sqrt(M));
% cols = ceil(M / rows);
% 
% % If figure doesn't exist or axes mismatch, create one and store handles
% if ~exist('subplot_axes', 'var') || length(subplot_axes) ~= M || any(~isgraphics(subplot_axes))
%     figure(100); clf;  % Fixed figure number to preserve across runs
%     subplot_axes = gobjects(1, M);  % Preallocate handles
% 
%     for i = 1:M
%         subplot_axes(i) = subplot(rows, cols, i);
%         hold(subplot_axes(i), 'on');
%         title(s.xnames_fitting{i}, 'Interpreter', 'none');
%         xlabel('FiO_2');
%         ylabel(s.xnames_fitting{i});
%         grid on;
%     end
% end
% 
% % === Use MRO2 and MRCO2 values for labeling ===
% % === Use MRO2 and MRCO2 values for labeling ===
% MRO2_val = round(s.pars("MRO2"), 2);
% MRCO2_val = round(s.pars("MRCO2"), 2);
% run_label = sprintf('[%.2f, %.2f]', MRO2_val, MRCO2_val);
% 
% 
% 
% % === Manage label → color map in base workspace ===
% color_list = lines(20);  % or use any other colormap
% if evalin('base', 'exist(''run_colors_map'', ''var'')')
%     run_colors_map = evalin('base', 'run_colors_map');
%     run_labels = evalin('base', 'run_labels');
% else
%     run_colors_map = containers.Map();
%     run_labels = {};
% end
% 
% % Assign a color if this label is new
% if ~isKey(run_colors_map, run_label)
%     color_idx = mod(length(run_labels), size(color_list, 1)) + 1;
%     run_colors_map(run_label) = color_list(color_idx, :);
%     run_labels{end+1} = run_label;
% end
% 
% % Save back to base
% assignin('base', 'run_colors_map', run_colors_map);
% assignin('base', 'run_labels', run_labels);
% 
% this_color = run_colors_map(run_label);
% 
% % === Plotting ===
% for i = 1:M
%     if ~isgraphics(subplot_axes(i))
%         subplot_axes(i) = subplot(rows, cols, i);
%         hold(subplot_axes(i), 'on');
%         title(s.xnames_fitting{i}, 'Interpreter', 'none');
%         xlabel('FiO_2');
%         ylabel(s.xnames_fitting{i});
%         grid on;
%     end
% 
%     axes(subplot_axes(i));
%     plot(fiO2_input, vars_matrix(:, i), 'o-', ...
%         'Color', this_color, ...
%         'DisplayName', run_label);
%     % If variable is (almost) constant, fix y-axis limits
%     % if max(vars_matrix(:, i)) - min(vars_matrix(:, i)) < 1e-6  % tolerance for floating point
%     %     ymid = mean(vars_matrix(:, i));
%     %     ylim([ymid - 0.01, ymid + 0.01]);  % flat ±0.01 range
%     %     yticks(ymid);                      % show just one tick
%     % end
% 
% end
% 
% for i = 1:M
%     axes(subplot_axes(i));
%     legend('show');
% end
% 
% 
% 
% %plot(fiO2_input, vars_matrix);

vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");

patient_idx = 5;
hipoxic_iterations = 1:20;

fiO2_init = 0.21;
fiO2_input = fiO2_init - hipoxic_iterations * 0.005;

num_points = 6;   % VO2/VCO2 ladder has 6 points
M = 10;           % Number of variables of interest

% Initialize subplot layout
rows = ceil(sqrt(M));
cols = ceil(M / rows);
figure(100); clf;
subplot_axes = gobjects(1, M);

for i = 1:M
    subplot_axes(i) = subplot(rows, cols, i);
    hold(subplot_axes(i), 'on');
    xlabel('FiO_2');
    ylabel(sprintf('Var %d', i));  % Placeholder, overwritten later
    grid on;
end

color_list = lines(num_points);

% Get VO2 and VCO2 ladder values only once
s0 = set_up("fiO2_ladder", patient_idx, "hipoxia", "mix", "dt", 0.1, ...
            "simulation_time", 200, "settling_time", 40);
VO2_ladder = s0.VO2_ladder_points / 1000;
VCO2_ladder = s0.VCO2_ladder_points / 1000;

for point_idx = 1:num_points
    all_vars = cell(1, length(hipoxic_iterations));
    vars_matrix = zeros(length(hipoxic_iterations), M); 

    for hi = hipoxic_iterations
        s = set_up("fiO2_ladder", patient_idx, "hipoxia", "mix", "dt", 0.1, ...
                   "simulation_time", 200, "settling_time", 40);
        
        s.pars("MRO2") = VO2_ladder(point_idx);
        s.pars("MRCO2") = VCO2_ladder(point_idx);
        s.pars("fO2") = fiO2_input(hi) * 100;

        [t, x_dot, x_vars, x_keys, index] = s.run_ode_fun(s.model, s.pars, s.init, s.simulation_time, s.dt);
        struct_vars = arrange_results(x_dot, x_vars, x_keys, t);
        all_vars{hi} = struct_vars;

        pout = data_processing("pressure", struct_vars.P_sa, t);
        vout = data_processing("volume", struct_vars.V, t);

        vars_of_interest = [struct_vars.dVE', vout', struct_vars.TI', struct_vars.Tresp', ...
                            struct_vars.PAO2', struct_vars.PACO2', 60 * struct_vars.HR', ...
                            pout{2}', pout{3}', struct_vars.mean_P_sa'];    
        final_values = mean(vars_of_interest(end-500:end, :), 1);
        vars_matrix(hi, :) = final_values;
    end

    % Generate label and color
    label_str = sprintf('[%.3f, %.3f]', VO2_ladder(point_idx), VCO2_ladder(point_idx));
    this_color = color_list(point_idx, :);

    % Plot results
    for i = 1:M
        axes(subplot_axes(i));
        if point_idx == 1
            title(s.xnames_fitting{i}, 'Interpreter', 'none');
            ylabel(s.xnames_fitting{i});
        end
        plot(fiO2_input, vars_matrix(:, i), 'o-', ...
             'Color', this_color, ...
             'DisplayName', label_str);
    end
end

% Show legends
for i = 1:M
    axes(subplot_axes(i));
    legend('show');
end

