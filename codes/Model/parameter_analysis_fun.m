function parameter_analysis_fun(p)
    vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");
    patients = [1, 4, 5, 6];
    %% SENSITIVITY ANALYSIS
    patient_idx = patients(p);    
    setup_n = set_up("sens", patient_idx, "normoxia", ".");
    setup_h = set_up("sens", patient_idx, "hipoxia", "mix");
    setup = {setup_n, setup_h};
    sens_functions("saving", "-", setup);
end