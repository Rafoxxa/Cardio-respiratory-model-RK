%Compute initials 4 subjects

vectorize_dicts('run_ode.m', 'model_basic.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m');
patient_list = [1 4 5 6];
for j = 1:4
    patient_idx = patient_list(j);
    [setup] = set_up("compute-initials", patient_idx, "normoxia", "mix", "dt", 0.1,'pars_from_fitting', 1, 'fitting_mat_file', 'last', 'simulation_time', 200);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
    s = setup;
    %s.pars('f_ab_max') = 60;
    correction_I0 = @(I0) (I0 * (1 - 0.3) + 0.3 - s.pars('MRtCO2_basal'))/(s.pars('AT')- s.pars('MRtCO2_basal'));
    s.pars('I_0_h_s') = correction_I0(s.pars('I_0_h_s'));
    s.pars('I_0_v_s') = correction_I0(s.pars('I_0_v_s'));
    s.pars('I_0_p_s') = correction_I0(s.pars('I_0_p_s'));
    s.pars('I0_met') = correction_I0(s.pars('I0_met'));
    %s.pars('P_n') = 100;
    %s.pars('kes') = s.pars('kes') * 0.5; %reduje el efecto de fab sobre el HR
    %s.pars('Wb_v_s') = 0;
    %s.pars('fes_max') = 100;
    
    %s.simulation_time = 100;
    %s.pars('V_tot') = s.pars('V_tot')  + 2200;    
    
    % Ejecutar simulación ODE
    [t, x_dot, x_vars, x_keys, index] = s.run_ode_fun(s.model, s.pars, s.init, s.simulation_time, s.dt);
    
    % Organizar resultados
    struct_vars = arrange_results(x_dot, x_vars, x_keys, t);
    
    % Guardar resultados
    save(s.simulation_filename, 'struct_vars', 't', 'x_vars');
end