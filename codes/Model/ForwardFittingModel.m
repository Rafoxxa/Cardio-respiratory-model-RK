    
    %%%%%%%%%%%%%%%%%%%%%%%

    % Setting up and run
    %setup = parallel_fitting(5, "pattern_search");
    vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");    
    patient_idx = 5;
    % setup_out_normoxia = set_up("fitting", patient_idx, "normoxia", "-", "simulation_time", 500);
    % setup_out_hipoxia = set_up("fitting", patient_idx, "hipoxia", "mix", "simulation_time", 500);         
    % setup = setup_out_normoxia;
    % 
    % setup.texp_list = {setup_out_normoxia.texp, setup_out_hipoxia.texp};
    % setup.yexp_list = {setup_out_normoxia.yexp, setup_out_hipoxia.yexp};
    % setup.pars_list = {setup_out_normoxia.pars, setup_out_hipoxia.pars}; 
    % setup.simulation_time_list = {setup_out_normoxia.simulation_time, setup_out_hipoxia.simulation_time};


    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Take the best solution   
    %solver_output = exec_solver("find-best-solution", setup);
    setup = parallel_fitting(patient_idx, "find-best-solution", "last");
    out_solver = setup.out_solver;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%

    % Run simulation    
    parsfitted_filename = setup.best_fitting_filename;
    parts = split(parsfitted_filename, filesep);
    lastTwo = parts(end-1:end);
    path_parsfitted = fullfile(lastTwo{1}, lastTwo{2});
    setup_out_normoxia_sim = set_up("simulation", patient_idx, "normoxia", "-", "pars_from_fitting", 1, 'fitting_mat_file', path_parsfitted);
    sn = setup_out_normoxia_sim;
    setup_out_hipoxia_sim = set_up("simulation", patient_idx, "hipoxia", "mix", "pars_from_fitting", 1, 'fitting_mat_file', path_parsfitted);
    sh = setup_out_hipoxia_sim;

    global all_global;
    all_global = zeros(15, round(10 * sn.simulation_time/sn.dt) + 1) + 0;  %This array saves all the data used for delays and for integration 
    [t_normoxia, x_dot, x_vars, x_keys, ~] = sn.run_ode_fun(sn.model, sn.pars, sn.init, sn.simulation_time, sn.dt);
    struct_vars_normoxia = arrange_results(x_dot, x_vars, x_keys, t_normoxia);
    all_global = zeros(15, round(10 * sh.simulation_time/sh.dt) + 1) + 0;  %This array saves all the data used for delays and for integration 
    [t_hipoxia, x_dot, x_vars, x_keys, ~] = sh.run_ode_fun(sh.model, sh.pars, sh.init, sh.simulation_time, sh.dt);
    struct_vars_hipoxia = arrange_results(x_dot, x_vars, x_keys, t_hipoxia);

    %%%%%%%%%%%%%%%%%%%%%%%%

    % Save simulation    
    save(sn.simulation_filename, "struct_vars_normoxia", "t_normoxia");
    save(sh.simulation_filename, "struct_vars_hipoxia", "t_hipoxia");

    %%%%%%%%%%%%%%%

    % Plot
    figure;    
    old_mode = "off";  %that means, the struct_vars have the old version
    time_sim = t_normoxia;
    X_sim = struct_vars_normoxia;
    custom_plot("sim_vs_exp", {time_sim, sn.texp, X_sim, sn.yexp, setup.xnames_fitting, setup.units_table, 5, 2, sn.simulation_filename, old_mode});
    original_normoxia_simulation_filename = "../Simulations/only_simulation/5/1200_sec_normoxia-24-04-2025-p.mat";
    custom_plot("sim_vs_exp", {"-", sn.texp, "-", sn.yexp, setup.xnames_fitting, setup.units_table, 5, 2, original_normoxia_simulation_filename, old_mode});
    
    figure;
    time_sim = t_hipoxia; 
    X_sim = struct_vars_hipoxia;
    custom_plot("sim_vs_exp", {time_sim, sh.texp, X_sim, sh.yexp, setup.xnames_fitting, setup.units_table, 5, 2, sh.simulation_filename, old_mode});
    original_hipoxia_simulation_filename = "../Simulations/only_simulation/5/3300_sec_hipoxia-24-04-2025-p.mat";
    custom_plot("sim_vs_exp", {"-", sh.texp, "-", sh.yexp, setup.xnames_fitting, setup.units_table, 5, 2, original_hipoxia_simulation_filename, old_mode});


    
    sh_filename = sh.simulation_filename;
    sn_filename = sn.simulation_filename;
    name = sprintf("../plot_data/%d/logJ_is:%.4f.mat", patient_idx, out_solver.fval);
    if ~isfolder(sprintf("../plot_data/%d", patient_idx))
        mkdir(sprintf("../plot_data/%d", patient_idx));
    end
    save(name, "sh_filename", "sn_filename", "original_hipoxia_simulation_filename", "original_normoxia_simulation_filename");



    

    


        
    

