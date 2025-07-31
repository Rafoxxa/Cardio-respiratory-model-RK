function setup = parallel_fitting(patient_idx, solver_option, date)
    if nargin > 2
        if solver_option == "find-best-solution"
            pars_from_fitting = 0;
            requestedDate = date; % If the solver option is to find the best solution, we use the date provided
            fitting_mat_file = "";
        elseif solver_option == "pattern_search"
            pars_from_fitting = 1;
            requestedDate = ""; % Otherwise, we don't need a specific date
            if date == "last"
                fitting_mat_file = "last"; % If the date is 'last', we set fitting_mat_file to 'last'
            else
                fitting_mat_file = sprintf("%s/best.mat", date);
            end
               
        end
    else
        requestedDate = ""; %Requested date is only to have the name of the folder based on the date in the setup output. Just in case we want to use that date from setup output in other functions 
        fitting_mat_file = "";
        pars_from_fitting = 0; % Default value if not specified
    end
 
    
    %

    %fitting_mat_file is actually the file from which the last fitted parameters come from.
    
    % Setting up
    vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");    
    %patient_idx = 5;
    setup_out_normoxia = set_up("fitting", patient_idx, "normoxia", "-", "requestedDate", requestedDate, "fitting_mat_file", fitting_mat_file, "pars_from_fitting", pars_from_fitting, "simulation_time", 500);
    setup_out_hipoxia = set_up("fitting", patient_idx, "hipoxia", "mix", "requestedDate", requestedDate, "fitting_mat_file", fitting_mat_file, "pars_from_fitting", pars_from_fitting, "simulation_time", 500);         
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