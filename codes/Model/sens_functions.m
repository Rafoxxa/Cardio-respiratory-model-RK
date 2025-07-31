function out = sens_functions(mode, fun, setup)
    s = setup;

    if mode == "single"
        if fun == "run_base"
            out = run_base(s);       
        elseif fun == "run_LSA" 
            [sensitivities, error, perturbed] = run_LSA(s);    
            out = {sensitivities, error, perturbed}; 
        elseif fun == "build_LSA_matrix"
            out = build_LSA_matrix(s, s.variables_of_interest, s.sensitivities, s.perturbed);
        elseif fun == "filter_params_by_class"
            [reduced_sens, reduced_pars] = filter_params_by_class(s);
            out = {reduced_sens, reduced_pars};
        elseif fun == "sens_threshold"
            [sens_reduced, pars_to_sens_red] = sens_threshold(s);
            out = {sens_reduced, pars_to_sens_red};
        elseif fun == "show_p_space"
            show_p_space(s);
            out = "";       
        end
    elseif mode == "saving"
        if s.params_sample_size == 0
            
            results_base = run_base(s);
            s.results_base = results_base;
            [sens, error, perturbed] = run_LSA(s);
            
            sens_matrix = build_LSA_matrix(s, s.variables_of_interest, sens, perturbed); 
            s.sens_final_time_matrix = sens_matrix;        
            save_sens(s);
            out = {sens_matrix, sensitivities, error, perturbed, pars_to_sens};
        else
            for i = 1:s.params_sample_size
                % Create a copy of the parameter set
                params_base = s.pars;
                pars_free2move = s.pars_free2move;
                noise = 0.5; 
                params_base_noised = add_noise(params_base, pars_free2move, noise); 
                format long;
                disp(sqrt(sum(((params_base_noised.values - params_base.values)./(eps + params_base.values)).^2)));
                results_base = run_base(s);
                s.results_base = results_base;
                [sens, error, perturbed] = run_LSA(s, params_base);
                
                sens_matrix = build_LSA_matrix(s, s.variables_of_interest, sens, perturbed); 
                s.sens_final_time_matrix = sens_matrix;  
                s.new_pars =   params_base;      
                save_sens(s);
            end
        end   

    elseif mode == "loading"
        [sens_matrix, pars_to_sens] = load_sens(s);
        out = {sens_matrix, pars_to_sens};
    

    elseif mode == "p-space"
         out = show_p_space(s);     
    end    

    function save_sens(s)
        write_path_all = s.sensitivity_write_all_filename;
        write_path = s.sensitivity_write_filename;

        sensitivities = s.sensitivities;
        sens_final_time_matrix = s.sens_final_time_matrix;
        pars_to_sens = s.pars_to_sens;
        new_pars = s.new_pars;
        
        if s.params_sample_size == 0
        save(write_path_all, "sensitivities");  
        save(write_path, "sens_final_time_matrix", "pars_to_sens"); 
        else
            % Save the parameters and sensitivities for each sample
            for i = 1:s.params_sample_size
                filename = sprintf('%s_pspace_sampling_%d.mat', write_path, i);
                save(filename, "sens_final_time_matrix", "pars_to_sens", "new_pars"); 
            end
        end
    end

    function [sens_matrix, pars_to_sens] = load_sens(s)
        load_path = s.sensitivity_load_filename;
        sens_out_struct = load(load_path);
        sens_matrix = sens_out_struct.sens_final_time_matrix;
        pars_to_sens = sens_out_struct.pars_to_sens;

    end


    function results_base_ =  run_base(s)
        params_base = s.pars;
        [t, ~, x_base, ~, ~] = s.run_ode_fun(s.model, params_base, s.init, s.simulation_time, s.dt);
        results_base = x_base;
        results_base_ = data_processing("add-desired", {results_base, s.init}, t);
    end

    function [sens, error, perturbed] = run_LSA(s, params_base)
        pars_to_sens = s.pars_to_sens;        
        perturbed = cell(1, s.n_params_sens);
        error = cell(1, s.n_params_sens);
        sensitivities = s.sensitivities;
        epsilon = s.epsilon;
        percentages = s.percentages;
        model = s.model;
        run_ode_fun = s.run_ode_fun;
        init = s.init;
        dt = s.dt;
        simulation_time = s.simulation_time;
        results_base_ = s.results_base;

        for i = 1:s.n_params_sens
            
            try
                % Create a copy of the parameter set
                %all_global = zeros(15, round(10 * s.simulation_time/dt) + 1) + 0;
                params_perturbed = params_base;               
                
                % Perturb the i-th selected parameter by epsilon            
                param_key = pars_to_sens(i);  % Get the actual index of the 
                          
                params_perturbed(param_key) = params_base(param_key)*(1+ epsilon);        
                params_perturbed = estimate_newton_ohm(percentages, params_perturbed);
                % Run the ODE solver with the perturbed parameter set
                [t, ~, x_perturbed, ~, ~] = run_ode_fun(model, params_perturbed, init, simulation_time, dt);            
                
                
                 
                x_perturbed_ = data_processing("add-desired", {x_perturbed, init}, t); 
                matchSize_output = data_processing("match-size", {x_perturbed_, results_base_}, "-");
                x_perturbed__ = matchSize_output{1};
                results_base__ = matchSize_output{2};                   
                sensitivities{i} = (x_perturbed__ - results_base__)/(epsilon) ./ results_base__;
                error{i} = 0;
                if size(x_perturbed__, 1) <= 1
                    disp("error-size");
                    error{i} = 2;
                end
                perturbed{i} = x_perturbed__;
            catch
                disp("error-simulation");
                sensitivities{i} = 0;
                perturbed{i} = 0;
                error{i} = 1;
            end
        end

        sens = sensitivities;





    end

    function sens_final_time_matrix = build_LSA_matrix(s, variables_of_interest, sensitivities, perturbed)
        init = s.init;

        init("pm") = 0;
        init("ps") = 0;
        init("pd") = 0;
        init("VTidal") = 0;
        init("BF") = 0;
       
       num_variables = numel(variables_of_interest);       
       sens_final_time_matrix = zeros(s.n_params_sens, num_variables);       

       % Loop over each variable of interest
       for v = 1:num_variables
            key = variables_of_interest{v};
            target_variable_idx = find(init.keys == key);
            
            if numel(target_variable_idx) ~= 1
                error(['Key "' key '" should correspond to a single variable index.']);
            end
            
            % Calculate sensitivities for the current variable
            for i = 1:s.n_params_sens 
                if size(perturbed{i},1) > 1
                    sensitivity_data = sensitivities{i}(target_variable_idx, :);
                    N = numel(sensitivity_data);
                    sens_final_time_matrix(i, v) = 1/N * sqrt(sum(sensitivity_data.^2));
                else
                     sens_final_time_matrix(i, v) = 0;
                     disp(s.pars_to_sens(i));
                end


            end
        end

    end

    function [reduced_sens, reduced_pars] = filter_params_by_class(s)

        excelData = readtable('parameters_classified_ori.xlsx');  
        isArbitrary = strcmp(excelData.Class, 'Arbitrary');
        arbitraryParams = excelData.Parameter(isArbitrary);   
        
        [~, idxInSens] = ismember(arbitraryParams, s.pars_to_sens);    
        
        validIdx = idxInSens > 0;
        idxInSens = idxInSens(validIdx);  
        
        reduced_sens = s.sens_final_time_matrix(idxInSens, :);
        reduced_pars = s.pars_to_sens(idxInSens);
    end

    function [sens_reduced, pars_to_sens_red] = sens_threshold(s)
        sens_final_time_matrix = s.sens_matrix;
        pars_to_sens = s.pars_to_sens;
        sens_threshold = s.sens_threshold;
        idx_variable_of_interest = s.idx_variable_of_interest;

        if numel(idx_variable_of_interest) > 1
            mask = find(sum(sens_final_time_matrix(:, idx_variable_of_interest),2) > sens_threshold);
        else
            mask = find(sens_final_time_matrix(:, idx_variable_of_interest) > sens_threshold);
        end   
        
        %mask1D = find(pars_to_sens == ["GVdead"    "I_0_h_s"    "Wp_p_s"    "Wp_v_s"    "dPmax"    "f_ab_max"    "f_ab_min"    "k_isc_v_s"    "phi_min"    "tau_V_u_s_v"    "x_h_s"    "x_v_s"]);

        
        pars_to_sens_red = pars_to_sens(mask);      
        sens_reduced = sens_final_time_matrix(mask, idx_variable_of_interest); 
    end

    function pars = add_noise(pars, pars_free2move, noise)
        rng("shuffle");
        for i = 1:length(pars_free2move)
            key = pars_free2move(i);
            if sign(pars(key)) >= 0
                pars(key) = pars(key) * (1 + noise * randn(1));
            else
                pars(key) = pars(key) * (1 - noise * randn(1));    
            end    
        end
    end

    function out = show_p_space(s)
        files = s.p_space_files;
        len_files = length(files);
        load(files(1), "sens_final_time_matrix");
        
        big_sens_matrix  = zeros(size(sens_final_time_matrix, 1), len_files);
     
        for idx_file = 1:len_files
            load(files(idx_file), "sens_final_time_matrix");
            sens_matrix_idx = sens_final_time_matrix;
            sens_matrix_1d = sum(sens_matrix_idx, 2);            
            big_sens_matrix(:,idx_file) = sens_matrix_1d;
        end
        %for par_idx = 1:size(big_sens_matrix, 1)
        %    plot_distribution(big_sens_matrix(par_idx, :));
        %end
        %plot_distribution(big_sens_matrix, new_pars); 
        miniplot_distribution(big_sens_matrix);
        out = big_sens_matrix;


    end

    function plot_distribution(data, list)
        % Example matrix: each row is a 1D dataset
    % (Replace this with your actual data)
    M = data;
    
    list_of_pars =  {"GVdead"    "I_0_h_s"    "Wp_p_s"    "Wp_v_s"    "dPmax"    "f_ab_max"    "f_ab_min"    "k_isc_v_s"    "phi_min"    "tau_V_u_s_v"    "x_h_s"    "x_v_s"};

    [num_datasets, N] = size(M);

    % Determine global x-axis range, excluding zeros
    M_nonzero = M;
    M_nonzero(M_nonzero == 0) = NaN;  % Temporarily replace 0 with NaN
    x_vals = linspace(min(M_nonzero,[],'all','omitnan'), ...
                      max(M_nonzero,[],'all','omitnan'), 200);
    
    % Bandwidth using Silverman's rule (exclude zeros)
    all_nonzero_vals = M(M ~= 0);
    h = 1.06 * std(all_nonzero_vals) * numel(all_nonzero_vals)^(-1/5);
    
    % Plot
    figure; hold on;
    colors = lines(num_datasets);
    
    for i = 1:num_datasets
        % Remove zero entries from this row
        data = M(i, M(i,:) ~= 0);
        if isempty(data), continue; end
    
        mu = mean(data);
    
        % Manual Gaussian KDE
        pdf_vals = zeros(size(x_vals));
        for j = 1:numel(data)
            pdf_vals = pdf_vals + exp(-(x_vals - data(j)).^2 / (2*h^2));
        end
        pdf_vals = pdf_vals / (numel(data) * h * sqrt(2*pi));
    
        % Plot KDE curve
        plot(x_vals, pdf_vals, 'Color', colors(i,:), 'LineWidth', 2, ...
            'DisplayName', list_of_pars{i});
    
        % Plot mean line (no label)
        xline(mu, '--', 'Color', colors(i,:));
    end
    
    % Show only KDE curves in legend
    legend('show');
    xlabel('Value');
    ylabel('Estimated Density');
    title('Gaussian KDEs (Ignoring Zero Values, No Mean Labels)');
    grid on;


    end

    function miniplot_distribution(data)
    M = data;
    
    % Names for each row
    list_of_pars =  {"GVdead"    "I_0_h_s"    "Wp_p_s"    "Wp_v_s"    "dPmax"    "f_ab_max"    "f_ab_min"    "k_isc_v_s"    "phi_min"    "tau_V_u_s_v"    "x_h_s"    "x_v_s"};
    



    % Step 1: Clean values (0 < x ≤ 0.013)
    M_clean = NaN(size(M));
    for i = 1:size(M,1)
        valid = M(i,:) > 0 & M(i,:) <= 0.013;
        M_clean(i,1:sum(valid)) = M(i,valid);
    end

    % Step 2: Plot KDEs
    all_clean_vals = M_clean(~isnan(M_clean));
    x_vals = linspace(min(all_clean_vals), max(all_clean_vals), 200);
    h = 1.06 * std(all_clean_vals) * numel(all_clean_vals)^(-1/5);

    figure; hold on;
    colors = lines(size(M_clean,1));

    for i = 1:size(M_clean,1)
        data = M_clean(i, ~isnan(M_clean(i,:)));
        if isempty(data), continue; end

        mu = mean(data);

        % Gaussian KDE
        pdf_vals = zeros(size(x_vals));
        for j = 1:numel(data)
            pdf_vals = pdf_vals + exp(-(x_vals - data(j)).^2 / (2*h^2));
        end
        pdf_vals = pdf_vals / (numel(data) * h * sqrt(2*pi));

        plot(x_vals, pdf_vals, 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', list_of_pars{i});

        xline(mu, '--', 'Color', colors(i,:));  % Mean line without label
    end

    legend('show');
    xlabel('Value');
    ylabel('Estimated Density');
    title('KDEs (0 < x ≤ 0.013, Gaussian Kernel, No Toolbox)');
    grid on;

    % Step 3: Manual Kruskal-Wallis test

    % Flatten data
    data_all = [];
    group_labels = [];
    for i = 1:size(M_clean,1)
        x = M_clean(i, ~isnan(M_clean(i,:)));
        data_all = [data_all; x(:)];
        group_labels = [group_labels; repmat(i, numel(x), 1)];
    end

    % Manual tied ranking
    ranks = manual_tied_rank(data_all);

    % Compute Kruskal-Wallis statistic
    k = size(M_clean,1);
    N = numel(data_all);
    group_sizes = accumarray(group_labels, 1);
    group_ranks = accumarray(group_labels, ranks);

    H = (12 / (N*(N+1))) * sum((group_ranks.^2) ./ group_sizes) - 3*(N+1);

    % p-value from chi-squared distribution with (k-1) df
    p = 1 - chi2cdf(H, k - 1);

    fprintf('\nManual Kruskal-Wallis test (no toolbox):\n');
    fprintf('H statistic = %.4f\n', H);
    fprintf('p-value     = %.4f\n', p);
    if p < 0.05
        disp('→ Statistically significant differences found between distributions.');
    else
        disp('→ No significant differences found between distributions.');
    end


function ranks = manual_tied_rank(x)
    % Returns tied ranks of a vector x (1-based)
    [sorted, idx] = sort(x);
    ranks = zeros(size(x));
    i = 1;
    while i <= length(x)
        j = i;
        while j < length(x) && sorted(j) == sorted(j+1)
            j = j + 1;
        end
        avg_rank = (i + j) / 2;
        ranks(idx(i:j)) = avg_rank;
        i = j + 1;
    end
end


    end

    








end