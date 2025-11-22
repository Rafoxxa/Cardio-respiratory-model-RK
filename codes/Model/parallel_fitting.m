function setup = parallel_fitting(idx, solver_option, date, initial_noise, iterations, pars2loadfolder)

    if ischar(idx)
        idx = str2double(idx);
    end
    if nargin > 2
        if strcmp(solver_option, 'find-best-solution') || strcmp(solver_option, 'find-best-solution-local') 
            pars_from_fitting = 1;
            requestedDate = date{end}; % If the solver option is to find the best solution, we use the date provided
            fitting_mat_file = date;
            iterations = 'None';
            initial_noise = 0.5;
        elseif strcmp(solver_option, 'pattern_search')
            pars_from_fitting = 1;
            requestedDate = ''; % Otherwise, we don't need a specific date
            if strcmp(date, 'last')
                fitting_mat_file = 'last'; % If the date is 'last', we set fitting_mat_file to 'last'
            else
                if isempty(date)
                    fitting_mat_file = '';
                    requestedDate = '';
                    pars_from_fitting = 0;
                else
                    fitting_mat_file = date;
                    requestedDate = date{end};
                end               
            end
               
        end
    else
        requestedDate = ''; %Requested date is only to have the name of the folder based on the date in the setup output. Just in case we want to use that date from setup output in other functions 
        fitting_mat_file = '';
        pars_from_fitting = 0; % Default value if not specified
        pars2loadfolder = 'last';

    end
    
    patient_array = [1,4,5,6];
    patient_idx = patient_array(idx);
    
 
    
    %

    %fitting_mat_file is actually the file from which the last fitted parameters come from.
    
    % Setting up
    %vectorize_dicts('run_ode.m', 'model_basic.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m');    
    %patient_idx = 5;
    setup_out_normoxia = set_up('fitting', patient_idx, 'normoxia', '-', 'requestedDate', requestedDate, 'fitting_mat_file', fitting_mat_file, 'pars_from_fitting', pars_from_fitting, 'type_of_optim', solver_option, 'initial_noise', initial_noise, 'iterations', iterations, 'pars2loadfolder', pars2loadfolder, 'time_from_data', 1);
    setup_out_hipoxia = set_up('fitting', patient_idx, 'hipoxia', 'mix', 'requestedDate', requestedDate, 'fitting_mat_file', fitting_mat_file, 'pars_from_fitting', pars_from_fitting, 'type_of_optim', solver_option, 'initial_noise', initial_noise, 'iterations', iterations, 'pars2loadfolder', pars2loadfolder,  'time_from_data', 1);          
    setup = setup_out_normoxia;

    setup.texp_list = {setup_out_normoxia.texp, setup_out_hipoxia.texp};
    setup.yexp_list = {setup_out_normoxia.yexp, setup_out_hipoxia.yexp};
    setup.pars_list = {setup_out_normoxia.pars, setup_out_hipoxia.pars}; 
    setup.simulation_time_list = {setup_out_normoxia.simulation_time, setup_out_hipoxia.simulation_time};
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Run fitting and internally save it in fitting folder.    
    out_solver = exec_solver(solver_option, setup);
    setup.out_solver = out_solver;
    

    
    %%%%%%%%%%%%%%%%%%%%%%%
end