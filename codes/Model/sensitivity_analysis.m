    %Only plotting from sensitivities from loaded
    %Only plotting from sens matrix from loaded
    %Run sensitivity analysis and save it with the actual date
    %Plot after running sensitivities
    %Selecting parameters

    
    
    currentDate = datetime('today');  
    formattedDate = datestr(currentDate, 'dd-mm-yyyy');
    vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");

    
    %% Driver constants for the code
    only_compute = 0;
    only_matrix = 0;
    only_plot = 1;
    load_sens = 1;
    save_sens = 1 - load_sens;

    state = "hipoxia-exercise";
    if only_matrix
        matrix = "with_matrix";
    else
        matrix = "";
    end
    file_load_name = "SensMatrix_07-04-2025.mat"; %"Sens_matrix.mat"; %"Sens_R_percentual.mat";
    file_save_name =  sprintf('Sens_R_percentual_%s_%s_%s.mat',state, matrix, formattedDate);

    %Analysis tools
    variables_of_interest = {'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal', 'dVE'};
    sens_threshold = 0.5;
    idx_variable_of_interest = [1:numel(variables_of_interest)];

    

    %% Computing things

    %Loading

    if load_sens
        folder_name = "../Sens_analysis/";
        sens_name = strcat(folder_name, file_load_name); 
        sens_out_struct = load(sens_name);
        %if only_matrix
            %sensitivities = sens_out_struct.sensitivities;
            sens_final_time_matrix = sens_out_struct.sens_final_time_matrix;
            pars_to_sens = sens_out_struct.pars_to_sens;
        %elseif only_plot            
        %    sens_final_time_matrix = sens_out_struct.sens_final_time_matrix;
        %end
    end

    %if only_compute || only_matrix

        % Specific parameters for sensitivity analysis
        [pars, init, taus] = load_global_easy();
        percentages = load("OhmNewton_percentages.mat").percentages;
        %Estimation of cardiovascular circuit components
        pars = estimate_newton_ohm(percentages, pars);
        
        %params_base = pars;
        %pars_to_sens = pars.keys;
        %Taking parameters from Dynamic Fitting paper
        %pars_to_sens = ["GTsym", "GTvagal", "G_R_e_p", "T0", "I0_met", "kmet", "PaCO2_n", "P_n", "phi_max", "K_E_lv", "K_E_rv", "KR_lv", "KR_rv", "R_sa", "A2", "alpha2", "C1", "C2", "K2", "MRbCO2", "Vtissue_CO2", "Kbg", "KcCO2", "KpCO2", "KpO2", "lambda1", "lambda2", "n", "dPmax", "Vdead", "El", "Rrs", "Ecw"];
        %pars_to_sens = ["Wb_h_s", "Wb_v_s", "Wc_h_s", "Wc_p_s", "Wc_v_s", "Wc_v", "Wp_h_s", "Wp_p_s", "Wp_v_s", "Wp_v", "Wt_h_s", "Wt_p_s", "Wt_v_s", "Wt_v", "gcc_h_s", "gcc_p_s", "gcc_v_s", "G_Emax_lv", "G_Emax_rv", "G_R_am_p", "G_R_e_p", "G_R_rm_p", "G_R_s_p", "G_V_u_am_v", "G_V_u_e_v", "G_V_u_rm_v", "G_V_u_s_v", "GTsym", "GTvagal", "KpCO2", "KcCO2", "KpO2"]; %autonomic
        %pars_to_sens = ["GTsym", "GTvagal"];
        n_params_sens = length(pars_to_sens);
    %end    

    if only_compute

        %Set up

        type_of_input = 7;
        control_on = 1;
        dt = 0.1;
        simulation_time = 1100;
        epsilon = 1e-3;

        VO2_external = 1;
        patient_idx = 1;
        hipoxia_state = "hipoxia";
        ascend_state = "exercise";

        if VO2_external
        [~, ~, VO2_poly, VCO2_poly, fO2_poly] = data_preprocessing(patient_idx, hipoxia_state, ascend_state,0);
            %VO2
            pars("MRO2_poly_0") = VO2_poly(1);
            pars("MRO2_poly_1") = VO2_poly(2);
            pars("MRO2_poly_2") = VO2_poly(3);
            pars("MRO2_poly_3") = VO2_poly(4);
            pars("MRO2_poly_4") = VO2_poly(5);
            pars("MRO2_poly_5") = VO2_poly(6);
            pars("MRO2_poly_6") = VO2_poly(7);
            pars("MRO2_poly_7") = VO2_poly(8);
            pars("MRO2_poly_8") = VO2_poly(9);
            
            %fiO2
        if hipoxia_state == "hipoxia"
        
            if ascend_state == "ascend"
                pars("fiO2_poly_0") = fO2_poly(1);
                pars("fiO2_poly_1") = fO2_poly(2);  %This can change depending on the best polynomial fit 
                pars("fiO2_poly_2") = fO2_poly(3);
                pars("fiO2_poly_3") = fO2_poly(4);
                pars("fiO2_poly_4") = fO2_poly(5);
            elseif ascend_state == "exercise" 
                pars("fiO2_poly_0") = fO2_poly(1);
                pars("fiO2_poly_1") = fO2_poly(2);
            end

        
        end
        
            %VCO2
        
            pars("MRCO2_poly_0") = VCO2_poly(1);
            pars("MRCO2_poly_1") = VCO2_poly(2);
            pars("MRCO2_poly_2") = VCO2_poly(3);
            pars("MRCO2_poly_3") = VCO2_poly(4);
            pars("MRCO2_poly_4") = VCO2_poly(5);
            pars("MRCO2_poly_5") = VCO2_poly(6);
            pars("MRCO2_poly_6") = VCO2_poly(7);
            pars("MRCO2_poly_7") = VCO2_poly(8);
            pars("MRCO2_poly_8") = VCO2_poly(9);
        end

    
        model = @(varargin) model_vec_hipoxia(varargin{:});
        run_ode_fun = @(varargin) run_ode_vec_hipoxia(varargin{:});
    
        
        pars('tau') = taus('tau_gases');
        pars('dt') = dt; 
        pars('type_of_input') = type_of_input;

        params_base = pars;
    
        %load('../simulations_saves/90sec_simulation.mat', 'x_vars');
        preloaded_vars = load('../simulations_saves/90sec_simulation.mat');
        init_values_loaded = preloaded_vars.x_vars(:, end);
        init_keys_loaded = fieldnames(preloaded_vars.struct_vars);
    
        init_keys = init.keys;
        for i = 1:length(init_values_loaded)
             if init_keys_loaded(i) ~= "vO2" && init_keys_loaded(i) ~= "PAO2" && init_keys_loaded(i) ~= "P_1O2" && init_keys_loaded(i) ~= "P_2O2" && init_keys_loaded(i) ~= "P_3O2" && init_keys_loaded(i) ~= "P_4O2" && init_keys_loaded(i) ~= "P_5O2" && init_keys_loaded(i) ~= "MRtO2"
                 key_i = init_keys(init_keys == init_keys_loaded(i));
                 init(key_i)= init_values_loaded(i);
             end
         end
    
        global all_global
        all_global = zeros(15, round(10 * simulation_time/dt) + 1) + 0;

        %To compute
        
        tic;
        sensitivities = cell(1, n_params_sens);
    
        [~, ~, x_base, ~, ~] = run_ode_fun(model, params_base, init, taus, simulation_time, dt, control_on);
        results_base = x_base;
        results_base_ = add_desired_variables(results_base, init);      
       
        

        for i = 1:n_params_sens
            
            %try
                % Create a copy of the parameter set
                all_global = zeros(15, round(10 * simulation_time/dt) + 1) + 0;
                params_perturbed = params_base;
                
                
                % Perturb the i-th selected parameter by epsilon            
                param_key = pars_to_sens(i);  % Get the actual index of the 
                          
                params_perturbed(param_key) = params_base(param_key)*(1+ epsilon);        
                params_perturbed = estimate_newton_ohm(percentages, params_perturbed);
                % Run the ODE solver with the perturbed parameter set
                [~, ~, x_perturbed, ~, ~] = run_ode_fun(model, params_perturbed, init, taus, simulation_time, dt, control_on);            
                perturbed_original{i} = x_perturbed;
                %sens = (x_perturbed - results_base)/(epsilon * params_perturbed(param_key)) ./ x_perturbed;
                x_perturbed_ = add_desired_variables(x_perturbed, init);        
                [x_perturbed__, results_base__] = matchSize(x_perturbed_, results_base_);
                sensitivities{i} = (x_perturbed__ - results_base__)/(epsilon) ./ results_base__;
                if size(x_perturbed__, 1) <= 1
                    disp("error");
                end
                perturbed{i} = x_perturbed__;
            %catch
                %sensitivities{i} = 0;
                %perturbed{i} = 0;
            %end
        end

        time_f = toc;
        disp(time_f);

        
        
        
    
        %Identificability analysis -- work in progress
    
        % sens = sensitivities
        % C= (sens'*sens)^-1;
        % 
        % c=[];
        % for i= 1:length(C)
        %     for j=1:length(C)
        %     c(i,j)= C(i,j)/sqrt(C(i,i)*C(j,j));
        %     end
        % end
    end


  

   %% Building sens_matrix  

   if only_matrix
       %adding post computed variables
       %pm; ps; pd; VT; BF
       
       init("pm") = 0;
       init("ps") = 0;
       init("pd") = 0;
       init("VTidal") = 0;
       init("BF") = 0;

       
       num_variables = numel(variables_of_interest);
       
       sens_final_time_matrix = zeros(n_params_sens, num_variables);
       

       % Loop over each variable of interest
       for v = 1:num_variables
            key = variables_of_interest{v};
            target_variable_idx = find(init.keys == key);
            
            if numel(target_variable_idx) ~= 1
                error(['Key "' key '" should correspond to a single variable index.']);
            end
            
            % Calculate sensitivities for the current variable
            for i = 1:n_params_sens 
                if sum(size(sensitivities{i})) > 2
                    sensitivity_data = sensitivities{i}(target_variable_idx, :);
                    N = numel(sensitivity_data);
                    sens_final_time_matrix(i, v) = 1/N * sqrt(sum(sensitivity_data.^2));
                else
                     sens_final_time_matrix(i, v) = 0;
                end


        end
        end
    end
        
    if save_sens 
        name = file_save_name;
        if only_matrix
            save(name, "sensitivities", "sens_final_time_matrix");  
          
        elseif only_compute
            save(name, "sensitivities");
        end        
        

    end    


    % Plotting
    only_plot = 1;
    sens_final_time_matrix = reduced_sens;
    pars_to_sens = pars_to_sens(idxInSens);
    if only_plot
        figure;
        %sens_threshold = 0.09;
        %idx_variable_of_interest = [3];
        if numel(idx_variable_of_interest) > 1
            mask = find(sum(sens_final_time_matrix(:, idx_variable_of_interest),2) > sens_threshold);
        else
            mask = find(sens_final_time_matrix(:, idx_variable_of_interest) > sens_threshold);
        end

        pars_to_sens_red = pars_to_sens(mask);            
        red_pars = numel(pars_to_sens_red) > 1;
        if red_pars
            
            b = bar(categorical(pars_to_sens_red), sens_final_time_matrix(mask, idx_variable_of_interest), 'stacked');
        else
            
            b = bar(categorical(pars_to_sens), sens_final_time_matrix(:, idx_variable_of_interest), 'stacked');
        end
        
        
        
        % Number of variables (stacked bars)
        num_variables = size(sens_final_time_matrix(:, idx_variable_of_interest), 2);  
        
        % Generate unique colors using a colormap
        colors = jet(num_variables);  % Change 'parula' to 'jet', 'hsv', etc., if preferred
        
        % Apply colors to each stacked segment
        for i = 1:num_variables
            b(i).FaceColor = 'flat';    
            b(i).CData = repmat(colors(i, :), size(b(i).XData, 1), 1); % Apply color to each bar segment
        end
        
        % Title, labels, and grid
        title('Sensibilidad para Diferentes Variables global (sum sqr)');
        xlabel('Parámetros');
        ylabel('Sensibilidad');
        grid on;
        
        % Add legend to identify variables
        legend(variables_of_interest(idx_variable_of_interest), 'Location', 'best');
    end  

    %To fitting

    % max_fiting_parameters = 35;
    % threshold_step = max(sum(sens_final_time_matrix(:, idx_variable_of_interest),2))/30;
    % sens_threshold = max(sum(sens_final_time_matrix(:, idx_variable_of_interest),2))/2;
    % 
    % pars_not_to_fit_table = readtable('pars_not_to_fit.xlsx', 'VariableNamingRule', 'preserve');
    % pars_not_to_fit = string(pars_not_to_fit_table{:, 1}); % Replace with actual column name
    % pars_needed_to_fit_table = readtable('pars_needed_to_fit.xlsx', 'VariableNamingRule', 'preserve');
    % pars_needed_to_fit = string(pars_needed_to_fit_table{:, 1}); % Replace with actual column name


    
    % pars_to_fit = pars_to_sens; %init condition
    % 
    % while size(pars_to_fit,1) > max_fiting_parameters
    % 
    %     if numel(idx_variable_of_interest) > 1
    %             mask = find(sum(sens_final_time_matrix(:, idx_variable_of_interest),2) > sens_threshold);
    %         else
    %             mask = find(sens_final_time_matrix(:, idx_variable_of_interest) > sens_threshold);
    %     end
    % 
    %     pars_to_sens_red = pars_to_sens(mask);  
    % 
    %     pars_to_fit = setdiff(pars_to_sens_red, pars_not_to_fit, 'stable');
    %     pars_to_fit= union(pars_to_fit, pars_needed_to_fit, 'stable');
    % 
    % 
    %     sens_threshold = sens_threshold + threshold_step; 
    % 
    % end
    % 
    % 
    % file_save_name_pars_to_fit =  sprintf('../Fitting/pars2fit/parameters_to_fit_%s_%s.mat',state, formattedDate);
    % save(file_save_name_pars_to_fit, "pars_to_fit");


    % Load data from Excel
    excelData = readtable('parameters_classified_ori.xlsx');  

    isArbitrary = strcmp(excelData.Class, 'Arbitrary');
    arbitraryParams = excelData.Parameter(isArbitrary);   
    
    [~, idxInSens] = ismember(arbitraryParams, pars_to_sens);    
    
    validIdx = idxInSens > 0;
    idxInSens = idxInSens(validIdx);  
    
    reduced_sens = sens_final_time_matrix(idxInSens, :);



    

    function [pm, ps, pd] = compute_presion(presion)
    % Crear un vector de tiempo basado en el tamaño de 'presion'
    n = length(presion);  % Número de puntos en la curva
    tiempo = linspace(0, n-1, n);  % Suponiendo un paso de tiempo uniforme
    
    % Derivadas de la curva de presión
    dp = diff(presion)./diff(tiempo);  % Primera derivada
    d2p = diff(dp)./diff(tiempo(1:end-1));  % Segunda derivada
    
    % Encontrar picos sistólicos (máximos locales) y diastólicos (mínimos locales)
    sistole_indices = find(dp(1:end-1) > 0 & dp(2:end) <= 0 & d2p < 0) + 1; % Máximos locales
    diastole_indices = find(dp(1:end-1) < 0 & dp(2:end) >= 0 & d2p > 0) + 1; % Mínimos locales
    
    % Extraer valores de presión sistólica y diastólica
    presion_sistolica = presion(sistole_indices);
    presion_diastolica = presion(diastole_indices);
    
    % Interpolación de los valores encontrados para tener un vector continuo
    ps = interp1(tiempo(sistole_indices), presion_sistolica, tiempo, 'linear', 'extrap');
    pd = interp1(tiempo(diastole_indices), presion_diastolica, tiempo, 'linear', 'extrap');
    
    % Calcular la presión media
    pm = trapz(tiempo, presion) / n;  % Integral dividida por la longitud del vector
    pm = ones(size(pd)) * pm;
    %[pm, ps, pd] = compute_presion(presion);
    
    end

    function [VT] = compute_VT(V)
        % Crear un vector de tiempo basado en el tamaño de 'presion'
        n = length(V);  % Número de puntos en la curva
        tiempo = linspace(0, n-1, n);  % Suponiendo un paso de tiempo uniforme
        
        % Derivadas de la curva de presión
        dp = diff(V)./diff(tiempo);  % Primera derivada
        d2p = diff(dp)./diff(tiempo(1:end-1));  % Segunda derivada
        
        % Encontrar picos sistólicos (máximos locales) y diastólicos (mínimos locales)
        VT_indices = find(dp(1:end-1) > 0 & dp(2:end) <= 0 & d2p < 0) + 1; % Máximos locales
        volumen_tidal = V(VT_indices);

        VT = interp1(tiempo(VT_indices), volumen_tidal, tiempo, 'linear', 'extrap');
    
    end

    function [BF] = compute_BF(TI, TE)
        BF = 1./(TI + TE);
    end

    function [Xout] = add_desired_variables(Xin, init)
        
        key_press = "P_sa";
        key_V = "V";
        key_TI = "TI";
        key_TE = "TE";
        
        press_idx = find(init.keys == key_press);
        V_idx = find(init.keys == key_V);
        TI_idx = find(init.keys == key_TI);
        TE_idx = find(init.keys == key_TE);

        presion = Xin(press_idx, :);
        V = Xin(V_idx, :);
        TI = Xin(TI_idx, :);
        TE = Xin(TE_idx, :);        

        [pm, ps, pd] = compute_presion(presion);
        [VT] = compute_VT(V);
        [BF] = compute_BF(TI, TE);

        Xout = [Xin; pm; ps; pd; VT; BF];
    end

    function [array1, array2] = matchSize(array1, array2)
        sz1 = size(array1, 2);
        sz2 = size(array2, 2);

        if sz1 == sz2
            return;
        end

        targetLength = min(sz1, sz2);
        array1 = array1(:, targetLength);
        array2 = array2(:, targetLength);

        %if sz1 > sz2            
            %array1 = downsampleArray(array1, targetLength);
        %else            
            %array2 = downsampleArray(array2, targetLength);
        %end
    end

    function outArray = downsampleArray(inArray, newSize)
        indices = round(linspace(1, length(inArray), newSize));
        outArray = inArray(indices);
    end



   %%plotting in time
   %plotting_in_time = 1;
   %plot(results_base_(find(init.keys == "dVE"), :));
   %hold on;
   %if plotting_in_time
   %
   %    for i=1:size(perturbed, 2)
   %        try
   %        param_effect = perturbed_original{i}; 
   %        plot(param_effect(find(init.keys == "dVE"), :));
   %        hold on 
   %        catch
   %            disp(i)
   %            pars_to_sens(i)
   %        end
   %    end
   %end