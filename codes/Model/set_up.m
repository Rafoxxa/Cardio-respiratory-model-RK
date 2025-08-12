function [setup_out] = ...
    set_up(case_of_use, patient_idx, hipoxia_state, ascend_state, varargin)

    model             = 'model';
    run_ode_fun       = 'run_ode_fun';
    pars              = 'pars';
    init              = 'init';
    type_of_input     = 7;
    control_on        = 1;
    only_plot         = 0;
    dt                = 0.01;
    simulation_time   = 2200;
    settling_time     = 20;
    pars_from_fitting = 0;
    fitting_mat_file  = '';
    VO2_external      = 1;    
    solver            = 'pattern_search';
    lb                = 'lb';
    ub                = 'ub';
    simulation_filename          = 'simulation_filename';
    fitting_filename   = 'fitting_filename';
    xnames_fitting    = 'xnames_fitting';
    percentages       = 'percentages';
    texp              = 'texp';
    yexp              = 'yexp';
    epsilon           = 1e-3; 
    requestedDate    = '';   
    fitting_filename = 'fitting_filename';
    fitting_folder = 'fitting_folder';
    best_fitting_filename = 'best_fitting_filename';    
    params_sample_size = 0;
    

    if strcmp(case_of_use, 'simulation')
        dt = 0.01;
        settling_time = 10;
        simulation_folder = 'only_simulation';

    elseif strcmp(case_of_use, 'fitting')
        dt = 0.1;
        settling_time = 11;
        simulation_folder = 'simulation_after_fitting';

    elseif strcmp(case_of_use, 'sens')
        dt = 0.1;
        settling_time = 11;        
        simulation_folder = 'sens_simulation';
        epsilon           = 1e-2; 

    elseif strcmp(case_of_use, 'fiO2_ladder')
        dt = 0.1;
        settling_time = 11;        
        simulation_folder = 'fiO2_ladder';
        epsilon           = 1e-2; 
        VO2_external = 0;
        %hipoxia_state = 'niether';
        type_of_input = 0;
        simulation_time = 500;

        
    end

    if strcmp(hipoxia_state, 'normoxia') && ~strcmp(case_of_use, 'fiO2_ladder')
        type_of_input = 6;
        simulation_time = 1200;%1200;
        
    elseif strcmp(hipoxia_state, 'hipoxia') && ~strcmp(case_of_use, 'fiO2_ladder')
        type_of_input = 7;
        if strcmp(ascend_state, 'ascend')
            simulation_time = 2200;
        elseif strcmp(ascend_state, 'exercise')
            simulation_time = 1200;
        elseif strcmp(ascend_state, 'mix')
            simulation_time = 2200 + 1100;
        end
        

    end
    
    defaults.model             = model;
    defaults.run_ode_fun       = run_ode_fun;
    defaults.pars              = pars;
    defaults.init              = init;
    defaults.type_of_input     = type_of_input;
    defaults.control_on        = 1;
    defaults.only_plot         = 0;
    defaults.dt                = dt;
    defaults.simulation_time   = simulation_time;
    defaults.settling_time     = settling_time;
    defaults.pars_from_fitting = pars_from_fitting;
    defaults.fitting_mat_file  = '../Fittings/Fitting-21-02-2025.mat';
    defaults.VO2_external      = VO2_external;
    defaults.patient_idx       = patient_idx;
    defaults.solver            = 'pattern_search';
    defaults.lb               = 'lb';
    defaults.ub               = 'ub';
    defaults.idx_optpars      = 'idx_optpars';
    defaults.optpars_0        = 'optpars_0';
    defaults.pars_values      = 'pars_values';
    defaults.init_keys        = 'init_keys';
    defaults.epsilon          = epsilon;
    defaults.requestedDate    = requestedDate;
    defaults.params_sample_size = params_sample_size;

    

    % Create input parser
    p = inputParser;
    addParameter(p, 'model', defaults.model);
    addParameter(p, 'run_ode_fun', defaults.run_ode_fun);
    addParameter(p, 'pars', defaults.pars);
    addParameter(p, 'init', defaults.init);
    addParameter(p, 'type_of_input', defaults.type_of_input);
    addParameter(p, 'control_on', defaults.control_on);
    addParameter(p, 'only_plot', defaults.only_plot);
    addParameter(p, 'dt', defaults.dt);
    addParameter(p, 'simulation_time', defaults.simulation_time);
    addParameter(p, 'settling_time', defaults.settling_time);
    addParameter(p, 'pars_from_fitting', defaults.pars_from_fitting);
    addParameter(p, 'fitting_mat_file', defaults.fitting_mat_file);
    addParameter(p, 'VO2_external', defaults.VO2_external);
    addParameter(p, 'patient_idx', defaults.patient_idx);
    addParameter(p, 'solver', defaults.solver);
    addParameter(p, 'lb', defaults.lb);
    addParameter(p, 'ub', defaults.ub);
    addParameter(p, 'idx_optpars', defaults.idx_optpars);
    addParameter(p, 'optpars_0', defaults.optpars_0);
    addParameter(p, 'pars_values', defaults.pars_values);
    addParameter(p, 'init_keys', defaults.init_keys);
    addParameter(p, 'epsilon', defaults.epsilon);
    addParameter(p, 'requestedDate', defaults.requestedDate);
    addParameter(p, 'params_sample_size', defaults.params_sample_size);



    parse(p, varargin{:});
    opts = p.Results;

    % Extract values
    model             = opts.model;
    run_ode_fun       = opts.run_ode_fun;
    pars              = opts.pars;
    init              = opts.init;
    type_of_input     = opts.type_of_input;
    control_on        = opts.control_on;
    only_plot         = opts.only_plot;
    dt                = opts.dt;
    simulation_time   = opts.simulation_time;
    settling_time     = opts.settling_time;
    pars_from_fitting = opts.pars_from_fitting;
    fitting_mat_file  = opts.fitting_mat_file;
    VO2_external      = opts.VO2_external;
    patient_idx       = opts.patient_idx;
    solver            = opts.solver;
    lb                = opts.lb;
    ub                = opts.ub;
    idx_optpars       = opts.idx_optpars;
    optpars_0         = opts.optpars_0;
    pars_values       = opts.pars_values;
    init_keys         = opts.init_keys;
    epsilon           = opts.epsilon;
    requestedDate     = opts.requestedDate;
    params_sample_size = opts.params_sample_size;


    


    %Select model
    model = @(varargin) model_vec_hipoxia(varargin{:});
    run_ode_fun = @(varargin) run_ode_vec_hipoxia(varargin{:});


    %Loadings
    [pars, init, taus] = load_global_easy();
    data = struct('p1', [70.8, 177, 2], 'p4', [85.5, 173, 2], 'p5', [77.5, 185, 2], 'p6', [68, 175, 2]);
    pars('BW') =  data.(sprintf('p%d', patient_idx))(1);
    pars('Hgt') =  data.(sprintf('p%d', patient_idx))(2);
    pars('Gender') =  data.(sprintf('p%d', patient_idx))(3);
    percentages = load('OhmNewton_percentages.mat');
    percentages = percentages.percentages;
    
    %Estimation of cardiovascular circuit components
    pars = estimate_newton_ohm(percentages, pars);
    pars_keys = keys(pars);

    %Sensitivity analysis
    pars_not_to_sens = load_pars_not_to_sens();
    pars_free2move = setdiff(pars_keys, pars_not_to_sens);     
    pars_to_sens = pars_free2move;
    %pars_to_sens =  {'GVdead'    'I_0_h_s'    'Wp_p_s'    'Wp_v_s'    'dPmax'    'f_ab_max'    'f_ab_min'    'k_isc_v_s'    'phi_min'    'tau_V_u_s_v'    'x_h_s'    'x_v_s'};
    pars_to_sens = {'GTsym'}; %, 'G_R_e_p', 'T0', 'I0_met', 'kmet'];
    %pars_to_sens = {'GTsym', 'GTvagal', 'G_R_e_p', 'T0', 'I0_met', 'kmet', 'PaCO2_n', 'P_n', 'phi_max', 'K_E_lv', 'K_E_rv', 'KR_lv', 'KR_rv', 'R_sa', 'A2', 'alpha2', 'C1', 'C2', 'K2', 'MRbCO2', 'Vtissue_CO2', 'Kbg', 'KcCO2', 'KpCO2', 'KpO2', 'lambda1', 'lambda2', 'n', 'dPmax', 'Vdead', 'El', 'Rrs', 'Ecw'};
    n_params_sens = length(pars_to_sens);
    variables_of_interest = {'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal', 'dVE'};
    idx_variable_of_interest = [1:numel(variables_of_interest)];

    %Load fitting parameters
    if pars_from_fitting
        %disp(patient_idx)
        if strcmp(fitting_mat_file, 'last')
            basePath = sprintf('../Fitting/parsFitted/%d', patient_idx);
            formattedDate = getLatestFittingDateStr(basePath);
            disp(formattedDate);
            fitting_mat_file = sprintf('Fitting-%s/best.mat', formattedDate);
            
        end
        fitting_mat_path = sprintf('../Fitting/parsFitted/%d/%s', patient_idx, fitting_mat_file);

        pars_struct = load(fitting_mat_path);
        updated_pars = pars_struct.updated_pars;
        
        pars_keys_updated = keys(updated_pars);
        for i = 1:length(pars_keys_updated)
            key_ = pars_keys_updated{i};
            pars(key_) = updated_pars(key_);
        end

        simulation_folder = 'simulation_after_fitting';
    end

    %read fast data from each pacient
    %if strcmp(hipoxia_state, 'normoxia') || strcmp(hipoxia_state, 'hipoxia')
    fast_data_filename = sprintf('../fast_data/%d/%s_data_preprocessed.mat', patient_idx, hipoxia_state);
    load(fast_data_filename, 'texp', 'yexp', 'VO2_poly', 'VCO2_poly', 'fO2_poly', 'basal', 'VO2_ladder_points', 'VCO2_ladder_points', 'AT');
    disp(basal(1))
    disp(basal(2))
    disp(AT);
    %end
    pars('AT') = AT;
    
    pars('MRtO2_basal') = basal(1); %this is making trouble
    pars('MRtCO2_basal') = basal(2); %this is making trouble

    VO2_ladder_points_ = VO2_ladder_points;
    VCO2_ladder_points_ = VCO2_ladder_points;

    %disp('MRtO2_basal:');
    %disp(basal(1));

    %save input coefficients in pars
    if VO2_external
        %[~, ~, VO2_poly, VCO2_poly, fO2_poly] = data_preprocessing(patient_idx, hipoxia_state, ascend_state,0);
        %VO2
        pars('MRO2_poly_0') = VO2_poly(1);
        pars('MRO2_poly_1') = VO2_poly(2);
        pars('MRO2_poly_2') = VO2_poly(3);
        pars('MRO2_poly_3') = VO2_poly(4);
        pars('MRO2_poly_4') = VO2_poly(5);
        pars('MRO2_poly_5') = VO2_poly(6);
        pars('MRO2_poly_6') = VO2_poly(7);
        pars('MRO2_poly_7') = VO2_poly(8);
        pars('MRO2_poly_8') = VO2_poly(9);
        
        %fiO2
        if strcmp(hipoxia_state, 'hipoxia')
            
            if strcmp(ascend_state, 'ascend') || strcmp(ascend_state, 'mix')
                pars('fiO2_poly_0') = fO2_poly(1);
                pars('fiO2_poly_1') = fO2_poly(2);  %This can change depending on the best polynomial fit 
                pars('fiO2_poly_2') = fO2_poly(3);
                pars('fiO2_poly_3') = fO2_poly(4);
                pars('fiO2_poly_4') = fO2_poly(5);
            elseif strcmp(ascend_state, 'exercise') 
                pars('fiO2_poly_0') = fO2_poly(1);
                pars('fiO2_poly_1') = fO2_poly(2);
            end

            
        end
        %VCO2

        pars('MRCO2_poly_0') = VCO2_poly(1);
        pars('MRCO2_poly_1') = VCO2_poly(2);
        pars('MRCO2_poly_2') = VCO2_poly(3);
        pars('MRCO2_poly_3') = VCO2_poly(4);
        pars('MRCO2_poly_4') = VCO2_poly(5);
        pars('MRCO2_poly_5') = VCO2_poly(6);
        pars('MRCO2_poly_6') = VCO2_poly(7);
        pars('MRCO2_poly_7') = VCO2_poly(8);
        pars('MRCO2_poly_8') = VCO2_poly(9);
    end

    %

    %save hyperparams in pars
    pars('tau') = taus('tau_gases');
    pars('dt') = dt; 
    pars('type_of_input') = type_of_input;
    pars('settling_time') = settling_time;
    units_table = readtable('variables_units.xlsx');

    %Load initial conditions in init
    if ~only_plot
        preloaded_vars = load('../Simulations/only_simulation/90sec_simulation.mat');
        init_values_loaded = preloaded_vars.x_vars(:, end);
        init_keys_loaded = fieldnames(preloaded_vars.struct_vars);

        init_keys = keys(init);


        % Redefine the values in the dictionary using the new vector
        for i = 1:length(init_values_loaded)
            if ~strcmp(init_keys_loaded{i}, 'vO2') && ~strcmp(init_keys_loaded{i}, 'PAO2') && ~strcmp(init_keys_loaded{i}, 'P_1O2') && ~strcmp(init_keys_loaded{i}, 'P_2O2') && ~strcmp(init_keys_loaded{i}, 'P_3O2') && ~strcmp(init_keys_loaded{i}, 'P_4O2') && ~strcmp(init_keys_loaded{i}, 'P_5O2') && ~strcmp(init_keys_loaded{i}, 'MRtO2') && ~strcmp(init_keys_loaded{i}, 'aO2') && ~strcmp(init_keys_loaded{i}, 'vO2')  && ~strcmp(init_keys_loaded{i}, 'PvbCO2') && ~strcmp(init_keys_loaded{i}, 'PbCO2') && ~strcmp(init_keys_loaded{i}, 'PCSFCO2')  && ~strcmp(init_keys_loaded{i}, 'mean_PbCO2') && ~strcmp(init_keys_loaded{i}, 'TI') && ~strcmp(init_keys_loaded{i}, 'TE') && ~strcmp(init_keys_loaded{i}, 'a1') && ~strcmp(init_keys_loaded{i}, 'a2') && ~strcmp(init_keys_loaded{i}, 'a0') 
                key_i_idx = find(strcmp(init_keys, init_keys_loaded{i}));
                if ~isempty(key_i_idx)
                    key_i = init_keys{key_i_idx(1)};
                    init(key_i) = init_values_loaded(i);
                end
            end
        end
        
    end

    %Simulation savings
    currentDate = datetime('today');  
    formattedDate = datestr(currentDate, 'dd-mm-yyyy');
    simulation_filename = sprintf('../Simulations/%s/%d/%d_sec_%s-%s.mat',simulation_folder, patient_idx, simulation_time, hipoxia_state,formattedDate);
    %Sensitivity    
    sensitivity_write_all_filename = sprintf('../Sens_analysis/%d/sensitivities_%s.mat', patient_idx,formattedDate); 
    sensitivity_write_filename = sprintf('../Sens_analysis/%d/SensMatrix_%s.mat', patient_idx,formattedDate);
    sensitivity_load_filename = sprintf('../Sens_analysis/%d/SensMatrix_%s.mat', patient_idx,requestedDate);


    %optimization hyperparameters
    %[texp, yexp, ~, ~, ~] = data_preprocessing(patient_idx, hipoxia_state, ascend_state); %its better to run it always
    xnames_fitting = {'dVE', 'V', 'TI', 'Tresp', 'PAO2', 'PACO2', 'HR', 'PS', 'PD', 'PM'};  
    
    if strcmp(case_of_use, 'fitting')
        currentDate = datetime('today');  
        formattedDate = datestr(currentDate, 'dd-mm-yyyy');
        if ~strcmp(requestedDate, '')
            if strcmp(requestedDate, 'last')
                basePath = sprintf('../Fitting/parsFitted/%d', patient_idx);
                formattedDate = getLatestFittingDateStr(basePath);
            else
                formattedDate = requestedDate;    
            end
                
        end
        fitting_folder = sprintf('../Fitting/parsFitted/%d/Fitting-%s/', patient_idx, formattedDate);              
        if ~exist(fitting_folder, 'dir')
            mkdir(fitting_folder);
        end 
        
        timestamp = datestr(now, 'yyyymmdd_HHMMSSFFF');
        fitting_filename = sprintf('%s%s.mat', fitting_folder, timestamp);
        best_fitting_filename = sprintf('%sbest.mat', fitting_folder);
        
        [lb, ub] = load_optim_boundries(pars, patient_idx);
        %Small size pars domain for optim solver
        idx_optpars = find(ub ~= lb);
        pars_values = values(pars);
        optpars_0 = pars_values(idx_optpars);  

    end

    sensitivities = cell(1, n_params_sens);
    

    

    
    
    
    setup_out.model = model;
    setup_out.run_ode_fun = run_ode_fun;
    setup_out.pars = pars;
    setup_out.init = init;
    setup_out.simulation_time = simulation_time;
    setup_out.settling_time = settling_time;
    setup_out.dt = dt;
    setup_out.control_on = control_on;
    setup_out.only_plot = only_plot;
    setup_out.units_table = units_table;
    setup_out.solver = solver;
    setup_out.ub = ub;
    setup_out.lb = lb;
    setup_out.idx_optpars = idx_optpars;
    setup_out.optpars_0 = optpars_0;
    setup_out.pars_values = pars_values;
    setup_out.init_keys = init_keys;

    setup_out.fitting_filename = fitting_filename;
    setup_out.fitting_folder = fitting_folder;
    setup_out.best_fitting_filename = best_fitting_filename;
    setup_out.simulation_filename = simulation_filename;

    setup_out.xnames_fitting = xnames_fitting;
    setup_out.percentages = percentages;
    setup_out.texp = texp;
    setup_out.yexp = yexp;
    setup_out.pars_to_sens = pars_to_sens;
    setup_out.pars_free2move = pars_free2move;
    setup_out.params_sample_size = params_sample_size;
    setup_out.n_params_sens = n_params_sens;
    setup_out.sensitivity_write_all_filename = sensitivity_write_all_filename;
    setup_out.sensitivity_write_filename = sensitivity_write_filename;
    setup_out.sensitivity_load_filename = sensitivity_load_filename;
    setup_out.sensitivities = sensitivities;
    setup_out.variables_of_interest = variables_of_interest;
    setup_out.epsilon = epsilon;
    setup_out.idx_variable_of_interest = idx_variable_of_interest;

    setup_out.VO2_ladder_points = VO2_ladder_points_;
    setup_out.VCO2_ladder_points = VCO2_ladder_points_;

    

    function [lb, ub] = load_optim_boundries(pars, patient_idx)
        lb = pars;
        ub = pars;
        %list_of_pars = {'C2', 'G_R_e_p', 'KpCO2', 'T0', 'MRbCO2', 'lambda1'};    
        try    
           cell_of_pars = load_pars_to_fit(patient_idx);      %instead of this code we should use pars2fit files.
        catch
           error('parameters to fit not correctly saved');
        end
        
        
        for i = 1:length(cell_of_pars)  %10 y 0.1
            key = cell_of_pars{i}; 
            if sign(pars(key)) >= 0
                lb(key) = pars(key) * 0.5;
                ub(key) = pars(key) * 10;
            else
                lb(key) = pars(key) * 10;
                ub(key) = pars(key) * 0.5;    
            end    
        end

        lb = values(lb);
        ub = values(ub);
    end

function latestDateStr = getLatestFittingDateStr(basePath)
%GETLATESTFITTINGDATESTR Returns the latest dd-MM-yyyy string from Fitting-dd-MM-yyyy folders,
% excluding today's date. If the latest folder is empty, checks the next latest, and so on.

    % Get list of folders
    folders = dir(basePath);
    folders = folders([folders.isdir]);  % Keep only directories
    folders = folders(~ismember({folders.name}, {'.', '..'}));

    % Initialize
    validDates = datetime.empty;
    dateStrings = {};

    for k = 1:length(folders)
        name = folders(k).name;

        % Look for pattern: Fitting-dd-MM-yyyy
        if strncmp(name, 'Fitting-', length('Fitting-'))
            if length(name) > length('Fitting-')
                dateStr = name(length('Fitting-')+1:end);
            else
                dateStr = '';
            end
            try
                Dt = datetime(dateStr, 'InputFormat', 'dd-MM-yyyy');
                %if Dt ~= datetime('today')  % Exclude today's date
                    validDates(end+1) = Dt;
                    dateStrings{end+1} = dateStr;
                %end
            catch
                % Skip invalid formats
            end
        end
    end

    % Sort dates in descending order
    if ~isempty(validDates)
        [sortedDates, sortIdx] = sort(validDates, 'descend');
        sortedDateStrings = dateStrings(sortIdx);

        % Loop over sorted dates and check folder contents
        for i = 1:length(sortedDates)
            folderPath = fullfile(basePath, ['Fitting-' sortedDateStrings{i}]);
            contents = dir(folderPath);
            % Exclude . and .. entries
            contents = contents(~ismember({contents.name}, {'.', '..'}));

            if ~isempty(contents)
                latestDateStr = sortedDateStrings{i};
                return;
            end
        end
    end

    % No valid folder with files found
    latestDateStr = '';
    warning('No valid Fitting-dd-MM-yyyy folders with files found in %s (excluding today)', basePath);
end


        
        
        



end