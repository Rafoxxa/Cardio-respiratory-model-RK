% set_up.m
% Prepara y devuelve una estructura 'setup_out' con todos los parámetros y
% configuraciones necesarias para ejecutar el modelo y la simulación.
% Entradas:
%  - case_of_use: cadena que indica el propósito ('simulation', 'fitting', etc.)
%  - patient_idx: índice del paciente (entero)
%  - hipoxia_state: 'normoxia'|'hipoxia'
%  - ascend_state: describe el tipo de subida ('ascend','exercise','mix', ...)
%  - varargin: parámetros opcionales pasados mediante name-value pairs
% Salida:
%  - setup_out: estructura con funciones, parámetros, condiciones iniciales,
%    nombres de archivos y otros metadatos necesarios.

function [setup_out] = ...
    set_up(case_of_use, patient_idx, hipoxia_state, ascend_state, varargin)

      % ---------------------------------------------------------------------
    % 1) Valores por defecto: variables internas y parámetros "hard-coded"
    %    Se definen valores por defecto que pueden ser sobrescritos por
    %    argumentos en varargin (mediante inputParser más abajo).
    % ---------------------------------------------------------------------

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
    simulation_filename  = 'initial-conditions';
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
    type_of_optim = 'pattern_search';
    initial_condition_filename = 'initial-conditions';%'old-classic';
    pars2loadfolder = 'last';
    time_from_data = 0;
    load_averaged_fitted_values = 0;
    param_index_sens = 1;
    load_rb = 0;
    estimated_newton = 1;

     
    % ---------------------------------------------------------------------
    % 2) Ajustes específicos según el use-case (case_of_use)
    %    Aquí se cambian dt, settling_time, simulation_folder y otros
    %    parámetros en función de la finalidad (simulación, fitting, sens, ...).
    % ---------------------------------------------------------------------
    
    

    if strcmp(case_of_use, 'compute-initials')
        dt = 0.1;
        settling_time = 0;
        simulation_folder = 'only_simulation';
        
    
    elseif strcmp(case_of_use, 'simulation')
        dt = 0.1;
        settling_time = 10;
        simulation_folder = 'only_simulation';
        load_averaged_fitted_values = 0;

    elseif strcmp(case_of_use, 'fitting')
        dt = 0.1;
        settling_time = 11;
        simulation_folder = 'simulation_after_fitting';

    elseif strcmp(case_of_use, 'sens')
        dt = 0.1;
        %settling_time = 11;  
        %simulation_time = 10;
        settling_time = 10;
        simulation_folder = 'sens_simulation';
        epsilon           = 1e-2; 
        load_averaged_fitted_values = 0;
        pars_from_fitting = 1;


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

     % ---------------------------------------------------------------------
    % 3) Ajustes según hipoxia_state y ascend_state
    %    Determinan el tipo de entrada y duración de la simulación.
    % ---------------------------------------------------------------------

    if strcmp(hipoxia_state, 'normoxia') && ~strcmp(case_of_use, 'fiO2_ladder') %&& ~strcmp(case_of_use, 'sens')
        type_of_input = 6;
        simulation_time = 1200;%1200;
        
    elseif strcmp(hipoxia_state, 'hipoxia') && ~strcmp(case_of_use, 'fiO2_ladder') %&& ~strcmp(case_of_use, 'sens')
        type_of_input = 7;
        if strcmp(ascend_state, 'ascend')
            simulation_time = 2200;
        elseif strcmp(ascend_state, 'exercise')
            simulation_time = 1200;
        elseif strcmp(ascend_state, 'mix')
            simulation_time = 2200 + 1100;
        end
        

    end

    % ---------------------------------------------------------------------
    % 4) Construcción de la estructura 'defaults' usada por inputParser
    %    Reúne los valores por defecto que se pueden sobrescribir.
    % ---------------------------------------------------------------------
    
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
    defaults.type_of_optim = type_of_optim;
    defaults.initial_noise = 0.5;
    defaults.iterations = 'None';
    defaults.initial_condition_filename = initial_condition_filename; 
    defaults.simulation_filename = simulation_filename; 
    defaults.pars2loadfolder = pars2loadfolder;
    defaults.time_from_data = time_from_data;
    defaults.param_index_sens = param_index_sens;
    defaults.load_rb = load_rb;
    defaults.estimated_newton = estimated_newton;

     % ---------------------------------------------------------------------
    % 5) inputParser: permite pasar parámetros opcionales vía name-value pairs
    % ---------------------------------------------------------------------  

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
    addParameter(p, 'type_of_optim', defaults.type_of_optim);
    addParameter(p, 'initial_noise', defaults.initial_noise);
    addParameter(p, 'iterations', defaults.iterations);
    addParameter(p, 'initial_condition_filename', defaults.initial_condition_filename);
    addParameter(p, 'simulation_filename', defaults.simulation_filename);
    addParameter(p, 'pars2loadfolder', defaults.pars2loadfolder);
    addParameter(p, 'time_from_data', defaults.time_from_data);
    addParameter(p, 'param_index_sens', defaults.param_index_sens);
    addParameter(p,'load_rb', defaults.load_rb);
    addParameter(p,'estimated_newton', defaults.estimated_newton);
    




    parse(p, varargin{:});
    opts = p.Results;

        % ---------------------------------------------------------------------
    % 6) Extraer valores finales (posibles sobrescrituras) desde opts
    % ---------------------------------------------------------------------

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
    initial_noise = opts.initial_noise;
    iterations = opts.iterations;
    initial_condition_filename = opts.initial_condition_filename;
    simulation_filename = opts.simulation_filename;
    pars2loadfolder = opts.pars2loadfolder;
    time_from_data = opts.time_from_data;
    param_index_sens = opts.param_index_sens;
    load_rb = opts.load_rb;
    estimated_newton = opts.estimated_newton;

    


    % ---------------------------------------------------------------------
    % 7) Selección de funciones del modelo (se descargan las versiones vectorizadas)
    %    Se asignan funciones anónimas que llaman a model_vec_hipoxia y run_ode_vec_hipoxia.
    % ---------------------------------------------------------------------

    %Select model
    model = @(varargin) model_vec_hipoxia(varargin{:});
    run_ode_fun = @(varargin) run_ode_vec_hipoxia(varargin{:});



    % ---------------------------------------------------------------------
    % 8) Carga inicial de parámetros, condiciones y constantes (load_global_easy)
    %    load_global_easy devuelve containers.Map para pars y init, y taus.
    % ---------------------------------------------------------------------

    %Loadings
    [pars, init, taus] = load_global_easy();
    data = struct('p1', [70.8, 177, 2], 'p4', [85.5, 173, 2], 'p5', [77.5, 185, 2], 'p6', [68, 175, 2]);
    pars('BW') =  data.(sprintf('p%d', patient_idx))(1);
    pars('Hgt') =  data.(sprintf('p%d', patient_idx))(2);
    pars('Gender') =  data.(sprintf('p%d', patient_idx))(3);
    percentages = load('OhmNewton_percentages.mat');
    percentages = percentages.percentages;
    
    
    % ---------------------------------------------------------------------
    % 9) Estimación por método Newton/OHM (opcional)
    %    estimate_newton_ohm puede modificar algunos parámetros basales.
    % ---------------------------------------------------------------------
    %Estimation of cardiovascular circuit components
    if estimated_newton
        pars = estimate_newton_ohm(percentages, pars, patient_idx);
    end

    % ---------------------------------------------------------------------
    % 10) Preparación para análisis de sensibilidad
    %    Se define qué parámetros están disponibles para sensibilidad, etc.
    % ---------------------------------------------------------------------

    pars_keys = keys(pars);
    %Sensitivity analysis
    pars_not_to_sens = load_pars_not_to_sens();
    pars_free2move = setdiff(pars_keys, pars_not_to_sens);     
    pars_to_sens = pars_free2move;
    disp(size(pars_to_sens));
    pars_to_sens = {pars_to_sens{param_index_sens}};
    %pars_to_sens =  {'GVdead'    'I_0_h_s'    'Wp_p_s'    'Wp_v_s'};%    'dPmax'    'f_ab_max'    'f_ab_min'    'k_isc_v_s'    'phi_min'    'tau_V_u_s_v'    'x_h_s'    'x_v_s'};
    %pars_to_sens = {'A','AT','C_ra','Kbg','LCTV','MO2_bp','P_n','PaO2_ac_n','V_unstressed_la','Wp_p_s','beta1','f_ac_max','gcc_p_s','kcc_v_s','tau_V_u_s_v','tau_w','ub_TI','ub_a1','ub_tau','vO2_am_n','vO2_b_n','vO2_e_n','vO2_h_n','vO2_rm_n','vO2_s_n','x_h_s','x_p_s','x_v_s'};

    %pars_to_sens = {'GTsym', 'G_R_e_p'}; %, 'G_R_e_p', 'T0', 'I0_met', 'kmet'];
    %pars_to_sens = {'GTsym', 'GTvagal', 'G_R_e_p', 'T0', 'I0_met', 'kmet', 'PaCO2_n', 'P_n', 'phi_max', 'K_E_lv', 'K_E_rv', 'KR_lv', 'KR_rv', 'R_sa', 'A2', 'alpha2', 'C1', 'C2', 'K2', 'MRbCO2', 'Vtissue_CO2', 'Kbg', 'KcCO2', 'KpCO2', 'KpO2', 'lambda1', 'lambda2', 'n', 'dPmax', 'Vdead', 'El', 'Rrs', 'Ecw'};
    n_params_sens = length(pars_to_sens);
    variables_of_interest = {'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal', 'dVE'};
    idx_variable_of_interest = [1:numel(variables_of_interest)];

    
    % ---------------------------------------------------------------------
    % 11) Carga de parámetros ajustados (fitting) si corresponde
    %    Si pars_from_fitting se activa, carga updated_pars desde archivos .mat
    %    y sobrescribe valores en 'pars'.
    % ---------------------------------------------------------------------
    %Load fitting parameters
    updated_pars_old = {};
    optpars_0_old = {};
    try
    if pars_from_fitting
        %disp(patient_idx)
        if strcmp(fitting_mat_file, 'last')
            % Busca la última carpeta con fitting válido y carga el best.mat
            basePath = sprintf('../Fitting/parsFitted/%d', patient_idx);
            formattedDate = getLatestFittingDateStr(basePath);
            disp(formattedDate);
            fitting_mat_file = sprintf('Fitting-%s/best.mat', formattedDate);

            fitting_mat_path = sprintf('../Fitting/parsFitted/%d/%s', patient_idx, fitting_mat_file);

            pars_struct = load(fitting_mat_path);
            updated_pars = pars_struct.updated_pars;
    
            pars_keys_updated = keys(updated_pars);
            
            for in = 1:length(pars_keys_updated)
                keyy = pars_keys_updated{in};
                pars(keyy) = updatedpars(keyy);
            end
                
        
        %fitting_mat_path = sprintf('../Fitting/parsFitted/%d/%s', patient_idx, fitting_mat_file);
        
        else
            % Si se provee una lista de fechas/archivos, se itera sobre ellos
            fitting_mat_files = fitting_mat_file;
            for file_idx = 1:length(fitting_mat_files)
                disp(file_idx);
                fitting_mat_file = fitting_mat_files{file_idx};
                fitting_mat_path = sprintf('../Fitting/parsFitted/%d/Fitting-%s/best.mat', patient_idx, fitting_mat_file);
                
                pars_struct = load(fitting_mat_path);
                updated_pars = pars_struct.updated_pars;
                
                pars_keys_updated = keys(updated_pars);
                disp(pars_keys_updated);
                disp(updated_pars.values);
                for in = 1:length(pars_keys_updated)
                    key_ = pars_keys_updated{in};
                    pars(key_) = updated_pars(key_);  % Las claves repetidas se sobrescriben
                end
                updated_pars_old{file_idx} = updated_pars;

            if file_idx == length(fitting_mat_files) - 1 %tomo los olds
                
                if strcmp(pars2loadfolder, 'last')
                    p2lf = '';
                else
                    p2lf = pars2loadfolder;
                end
                [lb_old, ub_old] = load_optim_boundries(pars, patient_idx, p2lf);
                %Small size pars domain for optim solver
                % Preparar optpars_0 con la subselección de parámetros a optimizar
                pars_values = pars.values;
                idx_optpars_old = find(~cellfun(@isequal, ub_old, lb_old));                
                optpars_0_old{file_idx} = pars_values(idx_optpars_old);  
                
                lb_old_mat = cell2mat(lb_old);
                ub_old_mat = cell2mat(ub_old);
                
                lb_old_tiny = lb_old_mat(idx_optpars_old);
                ub_old_tiny = ub_old_mat(idx_optpars_old);
                
                
               

            end
            end


            disp('pars loaded')
        end
    
    
        

        %simulation_folder = 'simulation_after_fitting';
    end
    catch
        disp('continue');
    end

    disp('check');

    %read fast data from each pacient
    %if strcmp(hipoxia_state, 'normoxia') || strcmp(hipoxia_state, 'hipoxia')


    % ---------------------------------------------------------------------
    % 12) Carga de datos del paciente y cálculo de condiciones iniciales
    %    Se lee el archivo preprocesado y se toman valores iniciales promedios.
    % ---------------------------------------------------------------------
    
    fast_data_filename = sprintf('../fast_data/%d/%s_data_preprocessed.mat', patient_idx, hipoxia_state);
    load(fast_data_filename, 'texp', 'yexp', 'VO2_poly', 'VCO2_poly', 'fO2_poly', 'TI_poly', 'Tresp_poly', 'basal', 'VO2_ladder_points', 'VCO2_ladder_points', 'AT', 'simulation_time_from_data');
    initials = take_initials(yexp);
    disp(simulation_time_from_data);

    if time_from_data
        simulation_time = simulation_time_from_data;
    end


    

    

    disp(basal(1))
    disp(basal(2))
    disp(AT);
    %end
    pars('AT') = AT;

    % Guardar metabolismos basales en pars (algunos sobrescritos luego)    
    pars('MRtO2_basal') = basal(1);%0.4;%basal(1); %this is making trouble
    pars('MRtCO2_basal') = basal(2);%0.33;%basal(2); %this is making trouble
    pars('MRO2') = pars('MRtO2_basal');
    pars('MRCO2') = pars('MRtCO2_basal'); 

    VO2_ladder_points_ = VO2_ladder_points;
    VCO2_ladder_points_ = VCO2_ladder_points;

    %disp('MRtO2_basal:');
    %disp(basal(1));


     % ---------------------------------------------------------------------
    % 13) Guardar polinomios externos en pars si VO2_external habilitado
    %    Estos coeficientes son usados por la entrada externa del modelo.
    % ---------------------------------------------------------------------
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

        pars('TI_poly_0') = TI_poly(1);
        pars('TI_poly_1') = TI_poly(2);
        pars('TI_poly_2') = TI_poly(3);
        pars('TI_poly_3') = TI_poly(4);
        pars('TI_poly_4') = TI_poly(5);
        pars('TI_poly_5') = TI_poly(6);
        pars('TI_poly_6') = TI_poly(7);
        pars('TI_poly_7') = TI_poly(8);
        pars('TI_poly_8') = TI_poly(9);
        pars('TI_poly_9') = TI_poly(10);
        pars('TI_poly_10') = TI_poly(11);
        pars('TI_poly_11') = TI_poly(12);
        pars('TI_poly_12') = TI_poly(13);
        pars('TI_poly_13') = TI_poly(14);
        pars('TI_poly_14') = TI_poly(15);
        pars('TI_poly_15') = TI_poly(16);

        pars('Tresp_poly_0') = Tresp_poly(1);
        pars('Tresp_poly_1') = Tresp_poly(2);
        pars('Tresp_poly_2') = Tresp_poly(3);
        pars('Tresp_poly_3') = Tresp_poly(4);
        pars('Tresp_poly_4') = Tresp_poly(5);
        pars('Tresp_poly_5') = Tresp_poly(6);
        pars('Tresp_poly_6') = Tresp_poly(7);
        pars('Tresp_poly_7') = Tresp_poly(8);
        pars('Tresp_poly_8') = Tresp_poly(9);
        pars('Tresp_poly_9') = Tresp_poly(10);
        pars('Tresp_poly_10') = Tresp_poly(11);
        pars('Tresp_poly_11') = Tresp_poly(12);
        pars('Tresp_poly_12') = Tresp_poly(13);
        pars('Tresp_poly_13') = Tresp_poly(14);
        pars('Tresp_poly_14') = Tresp_poly(15);
        pars('Tresp_poly_15') = Tresp_poly(16);

    end

    %
      % ---------------------------------------------------------------------
    % 14) Guardar hiperparámetros en pars (tau, dt, type_of_input, etc.)
    % ---------------------------------------------------------------------

    %save hyperparams in pars
    pars('tau') = taus('tau_gases');
    pars('dt') = dt; 
    pars('type_of_input') = type_of_input;
    pars('settling_time') = settling_time;
    units_table = readtable('variables_units.xlsx');
    disp('stop');

     % ---------------------------------------------------------------------
    % 15) Cargar condiciones iniciales pre-calc. y reasignarlas en init
    %    Lee un fichero .mat con x_vars y struct_vars y actualiza init (containers.Map)
    % ---------------------------------------------------------------------

    %Load initial conditions in init

    if ~only_plot
        if strcmp(initial_condition_filename, 'old-classic')
            preloaded_vars = load('../Simulations/only_simulation/90sec_simulation.mat');
        else
            path_file = sprintf('../Simulations/only_simulation/%d/%s.mat', patient_idx, initial_condition_filename);
            preloaded_vars = load(path_file);
        end
        init_values_loaded = preloaded_vars.x_vars(:, end);
        init_keys_loaded = fieldnames(preloaded_vars.struct_vars);

        init_keys = keys(init);


        % Redefine the values in the dictionary using the new vector
        if strcmp(initial_condition_filename, 'old-classic')
            for i = 1:length(init_values_loaded)
                if ~strcmp(init_keys_loaded{i}, 'vO2') && ~strcmp(init_keys_loaded{i}, 'PAO2') && ~strcmp(init_keys_loaded{i}, 'P_1O2') && ~strcmp(init_keys_loaded{i}, 'P_2O2') && ~strcmp(init_keys_loaded{i}, 'P_3O2') && ~strcmp(init_keys_loaded{i}, 'P_4O2') && ~strcmp(init_keys_loaded{i}, 'P_5O2') && ~strcmp(init_keys_loaded{i}, 'MRtO2') && ~strcmp(init_keys_loaded{i}, 'aO2') && ~strcmp(init_keys_loaded{i}, 'vO2')  && ~strcmp(init_keys_loaded{i}, 'PvbCO2') && ~strcmp(init_keys_loaded{i}, 'PbCO2') && ~strcmp(init_keys_loaded{i}, 'PCSFCO2')  && ~strcmp(init_keys_loaded{i}, 'mean_PbCO2') && ~strcmp(init_keys_loaded{i}, 'TI') && ~strcmp(init_keys_loaded{i}, 'TE') && ~strcmp(init_keys_loaded{i}, 'a1') && ~strcmp(init_keys_loaded{i}, 'a2') && ~strcmp(init_keys_loaded{i}, 'a0') 
                    key_i_idx = find(strcmp(init_keys, init_keys_loaded{i}));
                    if ~isempty(key_i_idx)
                        key_i = init_keys{key_i_idx(1)};
                        init(key_i) = init_values_loaded(i);
                    end
                end
            end
            init('MRtO2') = 0.4;
            init('MRtCO2') = 0.33;
            %{"dVE", "VT", "TI", "Tresp", "PAO2", "PACO2", "HR", "PS", "PD", "PM"}
            init('dVE') = initials(1);  
            init('TI') = initials(3);
            init('Tresp') = initials(4);
            init('PAO2') = initials(5);
            init('PACO2') = initials(6);
            if initials(7) > 0
                init('Theart') = 1/initials(7);
            end
            if initials(9) > 0
                init('P_sa') = initials(9);
            end
        
        else
            for i = 1:length(init_values_loaded)                
                key_i_idx = find(strcmp(init_keys, init_keys_loaded{i}));
                if ~isempty(key_i_idx)
                    key_i = init_keys{key_i_idx(1)};
                    init(key_i) = init_values_loaded(i);
                end
                init('vO2') = 0.185;
            end                

        end
        
        
    end

    % ---------------------------------------------------------------------
    % 16) Ajustes adicionales luego de estimated_newton
    % ---------------------------------------------------------------------

    %init('MRtO2') = pars('MRtO2_basal');
    %init('MRtCO2') = pars('MRtCO2_basal');
    %init('V_total_e_v') = 1000;
    if estimated_newton
        correction_I0 = @(I0) (I0 * (1 - 0.3) + 0.3 - pars('MRtCO2_basal'))/(pars('AT')- pars('MRtCO2_basal'));
        pars('I_0_h_s') = correction_I0(pars('I_0_h_s'));
        pars('I_0_v_s') = correction_I0(pars('I_0_v_s'));
        pars('I_0_p_s') = correction_I0(pars('I_0_p_s'));
        pars('I0_met') = correction_I0(pars('I0_met'));
    end

    if load_averaged_fitted_values
        pars = load_averaged_fitted_values_fun(pars);
    end
    
    %pars('PaO2_ac_n') = 45;
    %pars('f_ab_max') = 100;%47.78;
    %pars('kev') = 7.06;
    %pars('Wt_v_s') = 0.4;
    pars('R_p_p_n') = 0.0894;
    if estimated_newton
        pars('phi_max') = 13;
    end

    % ---------------------------------------------------------------------
    % 17) Construcción de nombres de ficheros de salida (simulaciones, fittings)
    % ---------------------------------------------------------------------

    %Simulation savings
    currentDate = datetime('today');  
    formattedDate = datestr(currentDate, 'dd-mm-yyyy');
    if strcmp(case_of_use, 'compute-initials')
        simulation_filename = sprintf('../Simulations/%s/%d/%s.mat',simulation_folder, patient_idx, simulation_filename);
    else
        simulation_filename = sprintf('../Simulations/%s/%d/%d_sec_%s-%s.mat',simulation_folder, patient_idx, simulation_time, hipoxia_state,formattedDate);
    end

    % ---------------------------------------------------------------------
    % 18) Rutas y nombres para sensibilidad y optimización
    % ---------------------------------------------------------------------
    
    %Sensitivity    
    sensitivity_write_all_filename = sprintf('../Sens_analysis/%d/sensitivities_%s_%s.mat', patient_idx, hipoxia_state, formattedDate); 
    %sensitivity_write_filename = sprintf('../Sens_analysis/%d/SensMatrix_%s_%s.mat', patient_idx, hipoxia_state ,formattedDate);
    sensitivity_write_filename = sprintf('../Sens_analysis/%d/SensTensor_%s_%s.mat', patient_idx, pars_to_sens{1},formattedDate);
    sensitivity_load_filename = sprintf('../Sens_analysis/%d/SensTensor_%s_%s.mat', patient_idx, pars_to_sens{1}, requestedDate);
    %sensitivity_load_filename = sprintf('../Sens_analysis/%d/SensMatrix_%s_%s.mat', patient_idx, hipoxia_state, requestedDate);


    %optimization hyperparameters
    %[texp, yexp, ~, ~, ~] = data_preprocessing(patient_idx, hipoxia_state, ascend_state); %its better to run it always
    xnames_fitting = {'dVE', 'VT', 'TI', 'Tresp', 'PAO2', 'PACO2', 'HR', 'PS', 'PD', 'PM'};  
    
    if strcmp(case_of_use, 'fitting')
        currentDate = datetime('today');  
        formattedDate = datestr(currentDate, 'dd-mm-yyyy');
        formattedDate = sprintf('Fitting-%s',formattedDate);
        if ~strcmp(requestedDate, '')
            if strcmp(requestedDate, 'last')
                basePath = sprintf('../Fitting/parsFitted/%d', patient_idx);
                formattedDate = getLatestFittingDateStr(basePath);
                formattedDate = sprintf('Fitting-%s',formattedDate);
            else   %comentar si se quiere entregar un fecha específica para optimizar
                formattedDate = sprintf('Fitting-%s', requestedDate);  %Esto solo sirve para local, para calcular el best.mat de una fecha en particular  
            end
                
        end
        fitting_folder = sprintf('../Fitting/parsFitted/%d/%s/', patient_idx, formattedDate);              
        if ~exist(fitting_folder, 'dir')
            mkdir(fitting_folder);
        end 
        
        timestamp = datestr(now, 'yyyymmdd_HHMMSSFFF');
        fitting_filename = sprintf('%s%s.mat', fitting_folder, timestamp);
        best_fitting_filename = sprintf('%sbest.mat', fitting_folder);
        
        [lb, ub] = load_optim_boundries(pars, patient_idx, pars2loadfolder);
        %Small size pars domain for optim solver
        idx_optpars = find(~cellfun(@isequal, ub, lb));
        pars_values = values(pars);
        optpars_0 = pars_values(idx_optpars);  
        disp('optpars_0');
        disp(optpars_0);
        

    end

    sensitivities = cell(1, n_params_sens);

    % ---------------------------------------------------------------------
    % 19) Identificador único para nodos de cálculo (e.g., cluster)
    % ---------------------------------------------------------------------
    
    %id for particle-swarm
    timestamp = round(posixtime(datetime('now')) * 1000);    
    % Número aleatorio (6 dígitos)
    randomNum = randi([100000, 999999]);    
    % Combinar como números
    node_id = timestamp * 1000000 + randomNum; % Multiplico por 1M para hacer espacio
       
    
    % ---------------------------------------------------------------------
    % 20) Montar la estructura de salida 'setup_out' con todo lo necesario
    % ---------------------------------------------------------------------
    
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

    setup_out.node_id = node_id;
    setup_out.patient_idx = patient_idx;
    
    setup_out.iterations = iterations;
    setup_out.initial_noise = initial_noise;
    setup_out.hipoxia_state = hipoxia_state;
    setup_out.load_rb = load_rb;

    try 
        setup_out.updated_pars_old = updated_pars_old;
        setup_out.optpars_0_old = optpars_0_old;
        setup_out.lb_old_tiny = lb_old_tiny;
        setup_out.ub_old_tiny = ub_old_tiny;
    catch
        setup_out.updated_pars_old = 'error';
    end
    

    % ---------------------------------------------------------------------
    % 21) Funciones locales anidadas: load_optim_boundries
    %    Devuelve lb y ub como celdas (values) para uso por optimizadores.
    % ---------------------------------------------------------------------

    function [lb, ub] = load_optim_boundries(pars, patient_idx, pars2loadfolder)
        lb = containers.Map(pars.keys, pars.values);
        ub = containers.Map(pars.keys, pars.values);
        %list_of_pars = {'C2', 'G_R_e_p', 'KpCO2', 'T0', 'MRbCO2', 'lambda1'};    
        %try    
           cell_of_pars = load_pars_to_fit(patient_idx, pars2loadfolder);      %instead of this code we should use pars2fit files.
        %catch
        %   error('parameters to fit not correctly saved');
        %end
        
        
        for i = 1:length(cell_of_pars)  %10 y 0.1
            key = cell_of_pars{i}; 
            [big, little] = boundries_factor(key, pars);
            if sign(pars(key)) >= 0
                lb(key) = pars(key) * little;
                ub(key) = pars(key) * big;
                
            else
                lb(key) = pars(key) * big;
                ub(key) = pars(key) * little;    
                
            end  
            disp('-----')
            disp(key);
            disp(pars(key));
            disp(lb(key));
            disp(ub(key));
        end

        lb = values(lb);
        ub = values(ub);
    end

    % ---------------------------------------------------------------------
    % 22) Función auxiliar boundries_factor: devuelve factores para límites
    % ---------------------------------------------------------------------
function  [b,l, pars] = boundries_factor(key, pars)
    %b = 10;
    b = 2.5;
    l = 0.05;
    
    switch key
        case 'PaO2_ac_n'            
            b = 60/pars(key); %mmHg
            l = 40/pars(key);
        case 'f_ab_max'
            b = 100/pars(key); %Hz(o spikes/s)
        case 'aO2_n'
            b = 1.2;
            l = 0.8;
        case 'R_p_p_n'
            b = 2;
            l = 0.5;
        case 'P_n'
            b = 100/pars(key);
            l = 85/pars(key);
        case 'C_sa'
            b = 4;
            l = 0.1;
        case 'R_sa'
            b = 10;
            l = 0.1;
        case 'L_sa'
            b = 1;
            l = 0.01;
        case 'GTsym'
            b = 5;
            l = 0.1;
        case 'phi_max'
            b = 1.1;
            l = 0.9;
       
        
        
    end
    
end
    % ---------------------------------------------------------------------
    % 23) getLatestFittingDateStr: devuelve la última carpeta Fitting-... con archivos
    % ---------------------------------------------------------------------
function latestDateStr = getLatestFittingDateStr(basePath)

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

function initials = take_initials(yexp)
    initials = mean(yexp(1:30, :), 1);
    
end




        
        
        



end