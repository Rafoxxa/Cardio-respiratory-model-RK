function [out_solver] = exec_solver(type_of_solver, setup)
    rng('shuffle');
    out_solver = struct();
    % Execute the solver
    %
    % :param type_of_solver: Type of solver to use
    % :param filename: Name of the file to save the results
    % :param time_exp: Experimental time vector
    % :param y_exp: Experimental data
    % :param init: Initial conditions
    % :param dt: Time step
    % :param settling_time: Settling time
    % :param pars: Parameters
    % :param idx_optpars: Indexes of the parameters to optimize
    % :param optpars_0: Initial guess for the parameters
    % :param lb: Lower boundaries for the parameters
    % :param ub: Upper boundaries for the parameters
    % :param percentages: Percentages of the parameters to optimize
    %
    % :return results: Results of the optimization
    % :return updated_pars: Updated parameters
    %
%for pattern search

    if strcmp(type_of_solver, 'find-best-solution')
        % Folder where the .mat files are located
    folderPath = setup.fitting_folder;  % <-- replace with actual path
    files = dir(fullfile(folderPath, '*.mat'));

    minFval = Inf;
    bestX = [];
    otherJvsiter = {};
    param_space = [];
    J_space = [];

    cnt = 1;
    for i = 1:length(files)
        data = load(fullfile(folderPath, files(i).name));
        
        % Check that both x and fval exist
        if isfield(data, 'x') && isfield(data, 'fval') && isfield(data, 'JvsIter')
        param_space = [param_space, data.JvsIter.bestx];
        J_space = [J_space,  data.JvsIter.bestfval]; 
        otherJvsiter{cnt} = data.JvsIter.bestfval;
            if data.fval < minFval
                minFval = data.fval;
                bestX = data.x;
                bestJvsiter = data.JvsIter.bestfval;                
            end
        cnt = cnt + 1;
        else
            %warning('File %s is missing "x" or "fval"', files(i).name);
        end
        
    end

    % Result
    fprintf('Minimum fval found: %g\n', minFval);
    disp('Corresponding x:');
    disp(bestX);

    pars = setup.pars_list{1};
    pars_keys = keys(pars);
    pars_keys_updated = pars_keys(setup.idx_optpars);
    updated_pars = containers.Map(pars_keys_updated, num2cell(bestX'));
    x = bestX;
    fval = minFval;


    out_solver = struct('x', bestX, 'fval', minFval);    
    save(setup.best_fitting_filename, 'updated_pars', 'fval', 'x');  

    %Show Iterations over each nuclei
    figure;
    hold on;
    plotHandles = gobjects(length(otherJvsiter) + 1, 1); % uno mÃ¡s para bestJvsiter
    legendLabels = cell(length(otherJvsiter) + 1, 1);
    
    % Plots de otherJvsiter
    for index = 1:length(otherJvsiter)  
        
        h = plot(otherJvsiter{index}, '-o');
        plotHandles(index) = h(1); % en caso de que sean varios, solo toma el primero
        legendLabels{index} = num2str(index + 3);
        
    end    
    
    % Plot adicional: bestJvsiter
    plotHandles(end) = plot(bestJvsiter, '-o', 'LineWidth', 2); % por ejemplo, lÃ­nea negra
    legendLabels{end} = 'best';
    
    
    legend(plotHandles, legendLabels);
    xlabel('iterations');
    ylabel('J'); 
    disp('end');
    uistack(plotHandles(end), 'top');

    %Checking boundries
    lb = cell2mat(setup.lb);
    lb_red = lb(setup.idx_optpars);
    ub = cell2mat(setup.ub);
    ub_red = ub(setup.idx_optpars);

    upper = abs(cell2mat(setup.optpars_0)) * 0.5;
    initial = abs(cell2mat(setup.optpars_0));
    lower = abs(cell2mat(setup.optpars_0)) * 10;
    keys_pars = setup.pars.keys;
    
    names = keys_pars(setup.idx_optpars);
    custom_plot('optim-boundries', {upper, initial, lower, abs(bestX), names});


    % ========== PCA VISUALIZATION OF PARAMETER SPACE ==========
    % PCA analysis of parameter space with J value coloring
    if size(param_space, 2) > 1 && size(param_space, 1) > 1
    fprintf('\nPerforming PCA analysis of parameter space...\n');
    
    % Transpose param_space so each row is an observation (solution)
    
    param_space = param_space(:, J_space < 0);
    J_space = J_space(J_space < 0);
    param_data = param_space';
    
    
    % Standardize the data (center and scale)
    param_centered = param_data - mean(param_data, 1);
    param_std = param_centered ./ std(param_centered, 0, 1);
    
    % Perform PCA
    [coeff, score, latent, ~, explained] = pca(param_std);
    
    % Find minimum J value and its index
    [min_J, min_idx] = min(J_space);
    
    % Create PCA visualization
    figure('Name', 'PCA Parameter Space Analysis', 'Position', [100, 100, 1200, 400]);
    
    % Subplot 1: PCA scatter plot (2D projection)
    subplot(1, 3, 1);
    scatter(score(:, 1), score(:, 2), 50, J_space, 'filled');
    hold on;
    % Highlight the minimum
    scatter(score(min_idx, 1), score(min_idx, 2), 200, 'red', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 2);
    colorbar;
    colormap(jet);
    xlabel(sprintf('PC1 (%.1f%% variance)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%% variance)', explained(2)));
    title('PCA Parameter Space (colored by J value)');
    grid on;
    legend('Solutions', 'Best Solution', 'Location', 'best');
    
    % Subplot 2: Explained variance
    subplot(1, 3, 2);
    bar(explained(1:min(10, length(explained))));
    xlabel('Principal Component');
    ylabel('Explained Variance (%)');
    title('PCA Explained Variance');
    grid on;
    
    % Subplot 3: Parameter contributions to PC1 and PC2
    subplot(1, 3, 3);
    n_params = min(size(coeff, 1), 10); % Show max 10 parameters
    x_pos = 1:n_params;
    bar_width = 0.35;
    
    bar(x_pos - bar_width/2, coeff(1:n_params, 1), bar_width, 'DisplayName', 'PC1');
    hold on;
    bar(x_pos + bar_width/2, coeff(1:n_params, 2), bar_width, 'DisplayName', 'PC2');
    
    xlabel('Parameter Index');
    ylabel('Loading Coefficient');
    title('Parameter Loadings on PC1 & PC2');
    legend('Location', 'best');
    grid on;
    xticks(x_pos);
    
    % Print PCA summary
    fprintf('PCA Summary:\n');
    fprintf('- Total variance explained by PC1-PC2: %.1f%%\n', sum(explained(1:2)));
    fprintf('- Minimum J value: %.6g at solution %d\n', min_J, min_idx);
    fprintf('- PC1 coordinates of best solution: %.3f\n', score(min_idx, 1));
    fprintf('- PC2 coordinates of best solution: %.3f\n', score(min_idx, 2));
    end
    





    
elseif strcmp(type_of_solver, 'pattern_search')
        % Define the objective function

        logStruct.iteration = [];
        logStruct.bestfval = [];
        logStruct.meshsize = [];
        logStruct.bestx = [];

        
        lower_boundries = setup.lb;
        upper_boundries = setup.ub;
        
        objFunParams.texp_list =  setup.texp_list;
        objFunParams.yexp_list = setup.yexp_list;
        objFunParams.dt = setup.dt;
        objFunParams.settling_time =  setup.settling_time;
        objFunParams.init =  setup.init;
        objFunParams.pars_list = setup.pars_list;
        objFunParams.idx_optpars = setup.idx_optpars;
        objFunParams.percentages = setup.percentages;
        objFunParams.simulation_time_list =  setup.simulation_time_list;       
        objFunParams.lb = lower_boundries;
        objFunParams.ub = upper_boundries;   
        objFunParams.xnames_fitting = setup.xnames_fitting;

        % Bounds for the optimization
        lb = cell2mat(lower_boundries(setup.idx_optpars)); % Lower boundaries
        ub = cell2mat(upper_boundries(setup.idx_optpars)); % Upper boundaries
        % Initial guesses for the optimization
        %num_start_points = 4; % Number of starting points
        initial_point = cell2mat(setup.optpars_0) .* (1 + rand(1, length(setup.optpars_0)) * 0.5); % Add small random perturbations
        initial_point = max(lb, min(ub, initial_point));
        disp(sprintf('SEED: %s', mat2str(initial_point, 2)));
        %Folder to store results in case of a crush
        objFunParams.initial_point = initial_point;
        objFunParams.lb = cell2mat(lower_boundries);
        objFunParams.ub = cell2mat(upper_boundries);

        objFun = @(args) obj_fun(args, objFunParams);

        % p = gcp('nocreate'); 
        % if isempty(p)
        % pool = parpool; % Start parallel pool
        % addAttachedFiles(pool, {'obj_fun.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m', 'estimate_newton_ohm.m', 'data_processing.m'}); % Attach necessary files
        % end

        % Options for patternsearch
        options = optimoptions('patternsearch', ...
        'OutputFcn', @saving_optim_info, ...
        'UseParallel', false, ...
        'Display', 'iter', ...
        'PollMethod', 'GPSPositiveBasis2N', ...
        'MaxFunctionEvaluations', 5);%, ...
        %'MaxIterations', 5, ... % Example maximum number of iterations
       % 'OutputFcn', @(optimValues, optold, flag) saveResultsOutputFcn(optimValues, optold, flag, folderName)); % Example max time (1 hour)

        % Initialize results storage
        %results = struct('x', [], 'fval', [], 'exitflag', [], 'output', []);       


        % Run patternsearch for each starting point
        %initial = tic;
        %for i = 1:num_start_points %parfor i = 1:num_start_points
        %fprintf('Running patternsearch from starting point %d\n', i);
        [x, fval, exitflag, output] = patternsearch(objFun, initial_point, [], [], [], [], lb, ub, options);
        %x;
        %fval;
        %exitflag = exitflag;
        %output = output;
       
        JvsIter = logStruct;
        save(setup.fitting_filename, 'x', 'fval', 'JvsIter'); 
        %out_solver.JvsIter = logStruct;

        %figure;
        %plot(logStruct.iteration, logStruct.bestfval, '-o');
        %xlabel('Iteration');
        %ylabel('Best Function Value');
        %title('Optimization Progress');

        % Display the best result
        % [Jmin, bestIdx] = min([results.fval]);
        % fprintf('Best solution found:\n');
        % disp(results(bestIdx).x);
        % disp(['Objective value: ', num2str(results(bestIdx).fval)]);

        % pars = setup.pars_list{1};
        % pars_keys = keys(pars);
        % pars_keys_updated = pars_keys(setup.idx_optpars);
        % updated_pars = containers.Map(pars_keys_updated, x');
        % save(setup.fitting_filename, 'updated_pars', 'fval');    

 
    
    
    elseif strcmp(type_of_solver, 'particle_swarm')
        
        % Para PSO: UseParallel = true es MÁS efectivo que múltiples puntos
        lower_boundries = setup.lb;
        upper_boundries = setup.ub;
        
        objFunParams.texp_list =  setup.texp_list;
        objFunParams.yexp_list = setup.yexp_list;
        objFunParams.dt = setup.dt;
        objFunParams.settling_time =  setup.settling_time;
        objFunParams.init =  setup.init;
        objFunParams.pars_list = setup.pars_list;
        objFunParams.idx_optpars = setup.idx_optpars;
        objFunParams.percentages = setup.percentages;
        objFunParams.simulation_time_list =  setup.simulation_time_list;       
        objFunParams.lb = lower_boundries;
        objFunParams.ub = upper_boundries;   
        objFunParams.xnames_fitting = setup.xnames_fitting;

        % Bounds for the optimization
        lb = cell2mat(lower_boundries(setup.idx_optpars)); 
        ub = cell2mat(upper_boundries(setup.idx_optpars)); 
        initial_point = cell2mat(setup.optpars_0); % Add small random perturbations
        objFunParams.initial_point = initial_point;
        objFunParams.lb = cell2mat(lower_boundries);
        objFunParams.ub = cell2mat(upper_boundries);
        objFun = @(args) obj_fun(args, objFunParams);

        
       
            % Setup parallel pool ANTES de configurar PSO
            %pool_success = setupParallelPool();
            
            % Configurar tamaño de swarm basado en dimensiones y workers
            nvars = length(lb);
            %if pool_success
            %    num_workers = gcp().NumWorkers;
                % Hacer swarm múltiplo del número de workers para eficiencia
             %   swarm_size = 20;%max(20, ceil(4*nvars/num_workers)*num_workers);
               % use_parallel = true;
               % fprintf('Using parallel PSO with %d workers, swarm size: %d\n', num_workers, swarm_size);
            %else
                swarm_size = 4;%max(10, 4*nvars);
                use_parallel = false;
                fprintf('Using serial PSO, swarm size: %d\n', swarm_size);
            %end
            
            % Initialize logging
            logStruct.iteration = [];
            logStruct.bestfval = [];
            logStruct.bestx = [];
            
            options = optimoptions('particleswarm', ...
           'MaxIterations', 50, ...               
           'SwarmSize', swarm_size, ...            
           'FunctionTolerance', 1e-6, ...          
           'MaxStallIterations', 20, ...           
           'UseParallel', use_parallel, ...        % Basado en disponibilidad del pool
           'OutputFcn', @saving_pso_info, ...      
           'Display', 'iter');                     

            fprintf('Starting PSO with SEED bounds: [%s] to [%s]\n', ...
                mat2str(lb, 2), mat2str(ub, 2));
                
            [x, fval, exitflag, output] = particleswarm(objFun, nvars, lb, ub, options);


        % Save results with logging structure
        JvsIter = logStruct;
        save(setup.fitting_filename, 'x', 'fval', 'JvsIter', 'exitflag', 'output'); 
        
        % Display results
        fprintf('Best parameter set found: %s\n', mat2str(x, 4));
        fprintf('Objective function value: %.6f\n', fval);
        
        out_solver = struct('x', x, 'fval', fval, 'JvsIter', logStruct);

    elseif strcmp(type_of_solver, 'particle_swarm_island')
        
        
        
        is_island_mode = logical(1);
        if is_island_mode
            node_id = setup.node_id;
            total_nodes = 10; % Por defecto 10 nodos
            if isfield(setup, 'total_nodes')
                total_nodes = setup.total_nodes;
            end
            
            % Crear directorio para migración
            migration_base_dir = sprintf('migration_data/p%d', setup.patient_idx);
            if ~exist(migration_base_dir, 'dir'), mkdir(migration_base_dir); end
            
            fprintf('=== ISLAND PSO MODE ===\n');
            fprintf('Node %d/%d starting...\n', node_id, total_nodes);
        else
            fprintf('=== SINGLE NODE PSO MODE ===\n');
        end
        
        % Configuración base del PSO
        lower_boundries = setup.lb;
        upper_boundries = setup.ub;
        
        objFunParams.texp_list =  setup.texp_list;
        objFunParams.yexp_list = setup.yexp_list;
        objFunParams.dt = setup.dt;
        objFunParams.settling_time =  setup.settling_time;
        objFunParams.init =  setup.init;
        objFunParams.pars_list = setup.pars_list;
        objFunParams.idx_optpars = setup.idx_optpars;
        objFunParams.percentages = setup.percentages;
        objFunParams.simulation_time_list =  setup.simulation_time_list;       
        objFunParams.lb = lower_boundries;
        objFunParams.ub = upper_boundries;   
        objFunParams.xnames_fitting = setup.xnames_fitting;

        lb = cell2mat(lower_boundries(setup.idx_optpars)); 
        ub = cell2mat(upper_boundries(setup.idx_optpars)); 
        nvars = length(lb);
        
        % Configuración del punto inicial
        initial_point = cell2mat(setup.optpars_0);
        
        if is_island_mode
            % Variación por nodo para diversidad inicial
            node_variation = 0.2 * (rand(1, length(initial_point)) - 0.5) * mod(node_id, total_nodes)/total_nodes;
            initial_point = initial_point .* (1 + node_variation);
            initial_point = max(lb, min(ub, initial_point));
            
            % Semilla específica por nodo
            rng(1000 + mod(node_id,100) * 137);
            
            % Parámetros optimizados para island
            swarm_size = 4; % Menos partículas por nodo
            max_iterations = 50; % Más iteraciones
            migration_interval = 6;
            
            % Modificar filename para incluir node_id
            [pathstr, name, ext] = fileparts(setup.fitting_filename);
            setup.fitting_filename = fullfile(pathstr, sprintf('%s_node_%d%s', name, node_id, ext));
            
        else
            % Configuración normal para nodo único
            swarm_size = max(10, 4*nvars);
            max_iterations = 100;
        end
        
        objFunParams.initial_point = initial_point;
        objFunParams.lb = cell2mat(lower_boundries);
        objFunParams.ub = cell2mat(upper_boundries);
        objFun = @(args) obj_fun(args, objFunParams);
        
        % Logging structure
        logStruct.iteration = [];
        logStruct.bestfval = [];
        logStruct.bestx = [];
        if is_island_mode
            logStruct.node_id = node_id;
            logStruct.migrations = [];
        end
        
        
        
        % Configurar opciones PSO
        if is_island_mode
            % Configuración para island model
            options = optimoptions('particleswarm', ...
                'MaxIterations', max_iterations, ...
                'SwarmSize', swarm_size, ...
                'FunctionTolerance', 1e-6, ...
                'MaxStallIterations', 15, ...
                'UseParallel', false, ...
                'OutputFcn', @adaptive_pso_output, ...
                'Display', 'iter', ...
                'SelfAdjustmentWeight', 1.49 + 0.05*randn(), ...
                'SocialAdjustmentWeight', 1.49 + 0.05*randn(), ...
                'InertiaRange', [0.1 + 0.02*randn(), 1.1 + 0.05*randn()]);
        else
            % Tu configuración original
            swarm_size = max(10, 4*nvars);
            use_parallel = false;
            fprintf('Using serial PSO, swarm size: %d\n', swarm_size);
            
            options = optimoptions('particleswarm', ...
               'MaxIterations', 100, ...               
               'SwarmSize', swarm_size, ...            
               'FunctionTolerance', 1e-6, ...          
               'MaxStallIterations', 20, ...           
               'UseParallel', use_parallel, ...        
               'OutputFcn', @adaptive_pso_output, ...      
               'Display', 'iter');
        end

        fprintf('Starting PSO with bounds: [%s] to [%s]\n', ...
            mat2str(lb, 2), mat2str(ub, 2));
        if is_island_mode
            fprintf('Swarm size: %d, Migration every %d iterations\n', swarm_size, migration_interval);
        end
            
        [x, fval, exitflag, output] = particleswarm(objFun, nvars, lb, ub, options);

        % Guardar resultados
        JvsIter = logStruct;
        save(setup.fitting_filename, 'x', 'fval', 'JvsIter', 'exitflag', 'output'); 
        
        fprintf('Best parameter set found: %s\n', mat2str(x, 4));
        fprintf('Objective function value: %.6f\n', fval);
        
        out_solver = struct('x', x, 'fval', fval, 'JvsIter', logStruct);
    
    elseif strcmp(type_of_solver, 'local_solver')
        % %Testing fitting with local solvers
        options = optimset('fmincon');
        options.Algorithm = 'sqp';
        options.Display = 'off';
        options.Diagnostics = 'off';
        options.FinDiffRelStep = 0.02; %This is a key element to make the solver run
        options.UseParallel = 1;
        options.RelLineSrchBnd = 1;
        %options.InitTrustRegionRadius = 10; %none effect

        options.MaxIter = 100;
        %options.OutputFcn = @outfun;
        lower_boundries = lb;
        upper_boundries = ub;

        %global history;
        %history.fval = []; % Store objective function values
        %history.x = [];    % Store solution at each iteration
        initial_time = tic;
        


        [optpars, J] = fmincon(@(args) obj_fun(args, time_exp, y_exp, init.values, dt, settling_time, init.keys, pars.values, idx_optpars), optpars_0, [], [], [], [], lower_boundries(idx_optpars), upper_boundries(idx_optpars), [], options);
        elapsed_time = toc;
        
        %pars(idx_optpars) = optpars;
        %save(filename, 'pars');
    elseif strcmp(type_of_solver, 'surrogate')
    % Initialize logging structure for surrogate optimization
    logStruct.iteration = [];
    logStruct.bestfval = [];
    logStruct.bestx = [];
    
    lower_boundries = setup.lb;
    upper_boundries = setup.ub;
    
    objFunParams.texp_list = setup.texp_list;
    objFunParams.yexp_list = setup.yexp_list;
    objFunParams.dt = setup.dt;
    objFunParams.settling_time = setup.settling_time;
    objFunParams.init = setup.init;
    objFunParams.pars_list = setup.pars_list;
    objFunParams.idx_optpars = setup.idx_optpars;
    objFunParams.percentages = setup.percentages;
    objFunParams.simulation_time_list = setup.simulation_time_list;       
    objFunParams.lb = lower_boundries;
    objFunParams.ub = upper_boundries;   
    objFunParams.xnames_fitting = setup.xnames_fitting;
    
    % Bounds for the optimization
    lb = cell2mat(lower_boundries(setup.idx_optpars)); % Lower boundaries
    ub = cell2mat(upper_boundries(setup.idx_optpars)); % Upper boundaries
    
    % Initial point for reference (opcional para puntos iniciales)
    initial_point = cell2mat(setup.optpars_0) + rand(1, length(setup.optpars_0)) * 0.01;
    initial_point = max(lb, min(ub, initial_point));
    disp(sprintf('SEED: %s', mat2str(initial_point, 2)));
    
    objFunParams.initial_point = initial_point;
    objFunParams.lb = cell2mat(lower_boundries);
    objFunParams.ub = cell2mat(upper_boundries);

    objFun = @(args) obj_fun(args, objFunParams);

    % Configurar opciones para custom_surrogateopt
    options_custom = struct();    
    options_custom.MinSampleDistance = 1e-6;
    options_custom.MinSurrogatePoints = max(10, length(lb) + 1);
    options_custom.MaxFunctionEvaluations = 20 * options_custom.MinSurrogatePoints;
    options_custom.Display = 'iter';
    options_custom.FunctionTolerance = 1e-6;
    options_custom.MaxStallIterations = 20;
    
    % Función de output personalizada para logging
    options_custom.OutputFcn = @saving_custom_surrogate_info;
    
    % Usar puntos iniciales si se desean (opcional)
    if isfield(setup, 'initial_points') && ~isempty(setup.initial_points)
        options_custom.InitialPoints = setup.initial_points;
    else
        % Generar algunos puntos iniciales incluyendo el punto de referencia
        n_initial = options_custom.MinSurrogatePoints;
        initial_points = zeros(n_initial, length(lb));
        
        % Primer punto: el punto de referencia
        initial_points(1, :) = initial_point;
        
        % Resto de puntos: aleatorios dentro de los bounds
        for i = 2:n_initial
            initial_points(i, :) = lb + rand(1, length(lb)) .* (ub - lb);
        end
        
        options_custom.InitialPoints = initial_points;
    end

    % Ejecutar custom_surrogateopt
    fprintf('Running custom surrogate optimization...\n');
    [x, fval, exitflag, output, log_custom] = custom_surrogateopt(objFun, lb, ub, options_custom);

    % Combinar logging structures (tu logging + el del algoritmo)
    JvsIter = logStruct;
    if ~isempty(log_custom) && isfield(log_custom, 'bestfval')
        % Usar datos del algoritmo si están disponibles
        JvsIter.iteration = log_custom.iteration;
        JvsIter.bestfval = log_custom.bestfval;
        JvsIter.bestx = log_custom.bestx;
        JvsIter.fval_history = log_custom.fval_history;
    end
    
    % Guardar resultados
    save(setup.fitting_filename, 'x', 'fval', 'JvsIter', 'exitflag', 'output', 'log_custom'); 
    
    % Mostrar resultados
    fprintf('\n=== CUSTOM SURROGATE OPTIMIZATION RESULTS ===\n');
    fprintf('Best parameter set found: %s\n', mat2str(x, 4));
    fprintf('Objective function value: %.6f\n', fval);
    fprintf('Function evaluations: %d\n', output.funcCount);
    fprintf('Iterations: %d\n', output.iterations);
    fprintf('Exit flag: %d\n', exitflag);
    fprintf('Algorithm: %s\n', output.algorithm);
    fprintf('==============================================\n');
    
    out_solver = struct('x', x, 'fval', fval, 'JvsIter', JvsIter, 'exitflag', exitflag, 'output', output);


    elseif strcmp(type_of_solver, 'multistart')
        % Define options for fmincon
        options = optimset('fmincon');
        options.Algorithm = 'sqp';
        options.Display = 'iter';
        options.Diagnostics = 'off';
        options.FinDiffRelStep = 0.02; % Key element to make the solver run
        %options.UseParallel = true;    % Enable parallel computation
        options.RelLineSrchBnd = 1;
        options.MaxIter = 10;
        options.MaxTime = 500;
        options.OutputFcn =  @(x, optimValues, state) myOutputFcn(x, optimValues, state, 3);
        
        
        % Define lower and upper boundaries
        lower_boundries = lb;
        upper_boundries = ub;
        
        % Create the problem structure for fmincon
        problem = createOptimProblem('fmincon', ...
            'objective', @(args) obj_fun(args, time_exp, y_exp, init.values, dt, settling_time, init.keys, pars.values, idx_optpars), ...
            'x0', optpars_0, ...
            'lb', lower_boundries(idx_optpars), ...
            'ub', upper_boundries(idx_optpars), ...
            'options', options);
        
        % Initialize the MultiStart solver
        ms = MultiStart('UseParallel', true, 'Display', 'iter', 'MaxTime',1200);
        
        % Run the MultiStart solver
        num_start_points = 2; % Set the number of starting points
        tic;
        %parpool;
        %results = run(ms, problem, num_start_points);
        [optpars, J] = run(ms, problem, num_start_points);
        elapsed_time = toc;
        disp(elapsed_time);


    
    
    elseif strcmp(type_of_solver, 'global_search')

        % Configure fmincon options
        options = optimoptions('fmincon', ...
            'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'Diagnostics', 'off', ...
            'FiniteDifferenceStepSize', 0.02, ... % Key element to make the solver run
            'UseParallel', true, ...
            'RelLineSrchBnd', 1, ...
            'MaxIter', 50);
        
        % Define bounds
        lower_boundries = lb;
        upper_boundries = ub;
        
        % Create the problem structure for fmincon
        problem = createOptimProblem('fmincon', ...
            'objective', @(args) obj_fun(args, time_exp, y_exp, init.values, dt, settling_time, init.keys, pars.values, idx_optpars), ...
            'x0', optpars_0, ...
            'lb', lower_boundries(idx_optpars), ...
            'ub', upper_boundries(idx_optpars), ...
            'options', options);
        
        % Initialize GlobalSearch object
        gs = GlobalSearch('Display', 'iter', 'NumTrialPoints', 100, 'NumStageOnePoints', 20);
        
        % Run GlobalSearch
        tic;
        [optpars, J, exitflag, output] = run(gs, problem);
        elapsed_time = toc;
        disp(elapsed_time);
        % Display results
        fprintf('Optimization completed in %.2f seconds\n', elapsed_time);
        fprintf('Objective Function Value: %.4f\n', J);
        fprintf('Optimal Parameters: %s\n', mat2str(optpars));
        fprintf('Exit Flag: %d\n', exitflag);
        disp('GlobalSearch Output:');
        disp(output);

    
    
    end

    %function stop = outfun(x, optimValues, state)
    %    stop = false; % Do not stop the optimization prematurely
    %    switch state
    %        case 'iter'
    %            % Append current solution and objective value to history
    %            history.x = [history.x; x];
    %            history.fval = [history.fval; optimValues.fval];
    %    end
    %end

    function [stop, options, optchanged] = saving_optim_info(optimValues, options, state)
        
        stop = false;
        optchanged = false;

        switch state
            case 'init'
                logStruct.iteration = [];
                logStruct.bestfval = [];
                logStruct.meshsize = [];
                logStruct.bestx = [];

            case 'iter'
                logStruct.iteration(end+1) = optimValues.iteration;
                logStruct.bestfval(end+1) = optimValues.fval;
                logStruct.meshsize(end+1) = optimValues.meshsize;
                logStruct.bestx(:,end+1) = optimValues.x;

                fprintf('Iter %d | fval = %.4f\n', ...
                    optimValues.iteration, optimValues.fval);

            case 'done'
                assignin('base', 'ps_log', logStruct);  % store in base workspace
                fprintf('Optimization complete. Log saved to variable `ps_log`.\n');
        end
    end


    function stop = myOutputFcn(~, optimValues, state, threshold)           

        stop = false; % Default: do not stop

        if strcmp(state, 'iter') % Check at each iteration
            if optimValues.fval < threshold
                stop = true; % Stop if the objective is less than threshold
                fprintf('Stopping: Objective function value %.4f is below threshold %.4f\n', optimValues.fval, threshold);
                
            end
        end

        
    end



    function folderName = createResultsFolder()
        % Define the parent folder outside the current directory
        parentFolder = fullfile('..', 'optimization_savings'); % One level up

        % Ensure the parent folder exists
        if ~exist(parentFolder, 'dir')
            mkdir(parentFolder);
        end

        % Create a timestamped subfolder inside optimization_savings
        timestamp = datestr(now, 'yyyy_mm_dd_HHMMSS');
        folderName = fullfile(parentFolder, ['optimization_data_', timestamp]);

        % Ensure the folder is created
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
    end

% Crear función de output personalizada para migración
    % Función de output adaptativa
    function [stop, options, optchanged] = saving_custom_surrogate_info(optimValues, options, state)
    stop = false;
    optchanged = false;
    
    switch state
        case 'init'
            logStruct.iteration = [];
            logStruct.bestfval = [];
            logStruct.bestx = [];
            fprintf('Custom surrogate optimization initialized\n');
            
        case 'iter'
            logStruct.iteration(end+1) = optimValues.iteration;
            logStruct.bestfval(end+1) = optimValues.fval;
            logStruct.bestx(:,end+1) = optimValues.x(:);
            
            fprintf('Custom Surrogate - Iter %d | fval = %.6f | FuncCount = %d\n', ...
                optimValues.iteration, optimValues.fval, optimValues.funccount);
            
        case 'done'
            fprintf('Custom surrogate optimization completed\n');
    end
    end
        function stop = adaptive_pso_output(optimValues, state)
            stop = false;
            
            switch state
                case 'init'
                    if is_island_mode
                        fprintf('Node %d: Initializing Island PSO\n', node_id);
                    end
                    
                case 'iter'
                    % Log básico
                    logStruct.iteration(end+1) = optimValues.iteration;
                    logStruct.bestfval(end+1) = optimValues.bestfval;
                    logStruct.bestx(:,end+1) = optimValues.bestx;
                    
                    if is_island_mode
                        fprintf('Node %d | Iter %d | Best fval = %.6f\n', ...
                            node_id, optimValues.iteration, optimValues.bestfval);
                        
                        % Checkpoint cada 5 iteraciones
                        if mod(optimValues.iteration, 5) == 0
                            checkpoint_file = sprintf('%s/checkpoint_node_%d_iter_%d.mat', ...
                            migration_base_dir, node_id, optimValues.iteration);                             
                            save(checkpoint_file, 'optimValues', 'logStruct');
                        end
                        
                        % Migración periódica
                        if mod(optimValues.iteration, migration_interval) == 0 && optimValues.iteration > 0
                            perform_migration(optimValues, node_id, total_nodes);
                        end
                    else
                        fprintf('Iter %d | fval = %.6f\n', ...
                            optimValues.iteration, optimValues.bestfval);
                    end
                    
                case 'done'
                    if is_island_mode
                        fprintf('Node %d: Optimization completed\n', node_id);
                        % Guardar resultado final del nodo
                        final_file = sprintf('%s/final_node_%d.mat', migration_base_dir, node_id);
                        save(final_file, 'optimValues', 'logStruct');
                    else
                        fprintf('PSO optimization complete.\n');
                    end
            end
        end
    function perform_migration(optimValues, node_id, total_nodes)
        try
            % Guardar estado actual de este nodo
            migration_file = sprintf('%s/node_%d_migration.mat', migration_base_dir, node_id);
            migration_data = struct();
            migration_data.bestx = optimValues.bestx;
            migration_data.bestfval = optimValues.bestfval;
            migration_data.iteration = optimValues.iteration;
            migration_data.timestamp = now;
            migration_data.node_id = node_id;
            
            save(migration_file, 'migration_data');
            
            % Leer datos de otros nodos
            other_nodes_data = [];
            for other_node = 1:total_nodes
                if other_node ~= node_id
                    other_file = sprintf('migration_data/node_%d_migration.mat', other_node);
                    if exist(other_file, 'file')
                        try
                            other_data = load(other_file);
                            % Solo datos recientes (últimas 3 horas)
                            time_diff = abs(now - other_data.migration_data.timestamp);
                            if time_diff < 0.125
                                other_nodes_data = [other_nodes_data; other_data.migration_data];
                            end
                        catch
                            % Ignorar archivos corruptos
                        end
                    end
                end
            end
            
            % Log de migración
            if ~isempty(other_nodes_data)
                logStruct.migrations(end+1) = optimValues.iteration;
                [best_external_fval, ~] = min([other_nodes_data.bestfval]);
                
                fprintf('Node %d | Migration: Found %d other nodes | Best external: %.6f\n', ...
                    node_id, length(other_nodes_data), best_external_fval);
            else
                fprintf('Node %d | Migration: No other nodes found yet\n', node_id);
            end
            
        catch ME
            fprintf('Node %d | Migration error: %s\n', node_id, ME.message);
        end
    end


    
    function [stop, optnew, changed] = saveResultsOutputFcn(optimValues, optold, flag, folderName)
        stop = false;  % Ensure optimization continues
        optnew = optimValues; % Typically unchanged
        changed = false; % We are not modifying optimization values
    
        % Ensure the folder exists
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        
        % Get worker ID (for parallel execution)
        workerID = getCurrentWorkerID();
        
        % Define the filename inside the results folder
        filename = fullfile(folderName, sprintf('optimization_results_worker_%d.mat', workerID));
    
        % Only save every 10 iterations
        if mod(optimValues.iteration, 10) == 0
            % Load previous results if file exists
            if exist(filename, 'file')
                try
                    load(filename, 'results');
                catch
                    results = [];
                end

            else
                results = [];  % Initialize if file doesn't exist
            end
    
            % Store current iteration results
            current_result = struct( ...
                'iteration', optimValues.iteration, ...
                'x', optimValues.x, ...
                'fval', optimValues.fval ...
            );
    
            % Append to results array
            results = [results; current_result];
    
            % Save results to a worker-specific file inside the folder
            save(filename, 'results');
    
            % Display a message
            fprintf('Worker %d saved iteration %d in %s\n', workerID, optimValues.iteration, filename);
        end
    end     
    
    
    function workerID = getCurrentWorkerID()
        % Check if running in parallel
        pool = gcp('nocreate');
        if isempty(pool)
            workerID = 0; % Running in serial mode
        else
            try
                t = getCurrentTask();
                workerID = t.ID; % Get parallel worker ID
            catch
                workerID = 1; % Fallback for older MATLAB versions
            end
        end
    end




end