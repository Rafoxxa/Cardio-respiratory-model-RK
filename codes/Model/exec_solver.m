function [out_solver] = exec_solver(type_of_solver, setup)
    rng("shuffle");
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

    if type_of_solver == "find-best-solution"
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
    pars_keys = pars.keys;
    pars_keys_updated = pars_keys(setup.idx_optpars);
    updated_pars = dictionary(pars_keys_updated, bestX');
    x = bestX;
    fval = minFval;


    out_solver = struct('x', bestX, 'fval', minFval);    
    save(setup.best_fitting_filename, "updated_pars", "fval", "x");  

    %Show Iterations over each nuclei
    figure;
    hold on;
    plotHandles = gobjects(length(otherJvsiter) + 1, 1); % uno más para bestJvsiter
    legendLabels = cell(length(otherJvsiter) + 1, 1);
    
    % Plots de otherJvsiter
    for index = 1:length(otherJvsiter)  
        
        h = plot(otherJvsiter{index}, '-o');
        plotHandles(index) = h(1); % en caso de que sean varios, solo toma el primero
        legendLabels{index} = num2str(index + 3);
        
    end    
    
    % Plot adicional: bestJvsiter
    plotHandles(end) = plot(bestJvsiter, '-o', 'LineWidth', 2); % por ejemplo, línea negra
    legendLabels{end} = 'best';
    
    
    legend(plotHandles, legendLabels);
    xlabel("iterations");
    ylabel("J"); 
    disp("end");
    uistack(plotHandles(end), 'top');





    
elseif type_of_solver == "pattern_search"
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
        lb = lower_boundries(setup.idx_optpars); % Lower boundaries
        ub = upper_boundries(setup.idx_optpars); % Upper boundaries
        % Initial guesses for the optimization
        %num_start_points = 4; % Number of starting points
        initial_point = setup.optpars_0' + rand(1, length(setup.optpars_0)) * 0.01; % Add small random perturbations
        initial_point = max(lb', min(ub', initial_point));
        disp(sprintf("SEED: %s", mat2str(initial_point, 2)));
        %Folder to store results in case of a crush
        objFunParams.initial_point = initial_point;

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
        save(setup.fitting_filename, "x", "fval", "JvsIter"); 
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
        % pars_keys = pars.keys;
        % pars_keys_updated = pars_keys(setup.idx_optpars);
        % updated_pars = dictionary(pars_keys_updated, x');
        % save(setup.fitting_filename, "updated_pars", "fval");    

 
    elseif type_of_solver == "surrogate"
        objFunParams.texp_list =  setup.texp_list;
        objFunParams.yexp_list = setup.yexp_list;
        objFunParams.dt = setup.dt;
        objFunParams.settling_time =  setup.settling_time;
        objFunParams.init =  setup.init;
        objFunParams.pars_list = setup.pars_list;
        objFunParams.idx_optpars = setup.idx_optpars;
        objFunParams.percentages = setup.percentages;
        objFunParams.simulation_time_list =  setup.simulation_time_list;       
        
        
        objFun = @(args) obj_fun(args, objFunParams);

        lower_boundries = setup.lb;
        upper_boundries = setup.ub;

        % Bounds for the optimization
        lb = lower_boundries(setup.idx_optpars); % Lower boundaries
        ub = upper_boundries(setup.idx_optpars); % Upper boundaries

        % Initial guesses for the optimization
        %num_start_points = 4; % Number of starting points
        %initial_point = setup.optpars_0' + rand(1, length(setup.optpars_0)) * 0.01; % Add small random perturbations

        options = optimoptions('surrogateopt', ...
       'MaxFunctionEvaluations', 40, ...      % total expensive evals allowed
       'UseParallel', false, ...               % use parpool if available     
       'Display', 'iter');                    % show progress in Command Window

        [x, fval, exitflag, output] = surrogateopt(objFun, lb, ub, options);

        % Display results
        disp('Best parameter set found:');
        disp(x);
        disp('Objective function value:');
        disp(fval);
        
        save(setup.fitting_filename, "x", "fval"); 

    elseif type_of_solver == "local_solver"
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
        %save(filename, "pars");


    elseif type_of_solver == "multistart"
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


    
    
    elseif type_of_solver == "global_search"

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
            workerID = getCurrentTask().ID; % Get parallel worker ID
        end
    end




end