%INSTRUCTIONS: In README.md there is an example of the model basic
%simulation that is executed here.

%Clear and setup
rng(2);  
%clc
clear -global delays_global
clear -global all_global
clear -global externals_global


vectorize_dicts('run_ode.m', 'model_basic.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m');

% State definition and patient selection: 
state = 'normoxia';
patient_idx = 5;
[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1,'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'time_from_data', 1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%[setup] = set_up("simulation", patient_idx, state, 'mix', 'pars_from_fitting', 0, "dt", 0.1, 'time_from_data', 1, 'estimated_newton', 1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1, 'initial_condition_filename', 'initial-conditions', 'time_from_data', 0, 'simulation_time', 500);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%[setup] = set_up('simulation', patient_idx, state, 'mix', 'dt', 0.1);
s = setup;
s.simulation_time = 100;

% Run ODE
[t, x_dot, x_vars, x_keys, index] = s.run_ode_fun(s.model, s.pars, s.init, s.simulation_time, s.dt);

% Organize results
struct_vars = arrange_results(x_dot, x_vars, x_keys, t);

% Save results in 'Simulations' folder
%save(s.simulation_filename, 'struct_vars', 't');

figure;
%Crear cell array para parámetros de plotting
plot_params = {t, s.texp, struct_vars, s.yexp, s.xnames_fitting, s.units_table, 5, 2, '', 'off',  patient_idx, 'actual', [1 2 3 4 5 6 7 8 9 10],state, 0};
custom_plot('sim_vs_exp', plot_params);

