%INSTRUCTIONS: In README.md there is an example of the model basic
%simulation that is executed here.

%Clear and setup
rng(2);  
%clc
clear -global delays_global
clear -global all_global
clear -global externals_global


vectorize_dicts('run_ode.m', 'model_basic.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m');

state = 'normoxia';
patient_idx = 5;
[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1,'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'time_from_data', 1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%[setup] = set_up("simulation", patient_idx, state, 'mix', 'pars_from_fitting', 0, "dt", 0.1, 'time_from_data', 1, 'estimated_newton', 1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1, 'initial_condition_filename', 'initial-conditions', 'time_from_data', 0, 'simulation_time', 500);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
%[setup] = set_up('simulation', patient_idx, state, 'mix', 'dt', 0.1);
s = setup;
s.dt = 0.01;
%s.pars('gO2_p') = 10 * s.pars('gO2_p');
disp('alto');
%s.pars('f_ab_max') = 60;
% correction_I0 = @(I0) (I0 * (1 - 0.3) + 0.3 - s.pars('MRtCO2_basal'))/(s.pars('AT')- s.pars('MRtCO2_basal'));
% s.pars('I_0_h_s') = correction_I0(s.pars('I_0_h_s'));
% s.pars('I_0_v_s') = correction_I0(s.pars('I_0_v_s'));
% s.pars('I_0_p_s') = correction_I0(s.pars('I_0_p_s'));
% s.pars('I_0_v') = correction_I0(s.pars('I_0_v'));
% pars('I0_met') = correction_I0(s.pars('I0_met'));F
%s.pars('P_n') = 100;
%s.pars('kes') = s.pars('kes') * 0.5; %reduje el efecto de fab sobre el HR
%s.pars('Wb_v_s') = s.pars('Wb_v_s') + sqrt(eps); %* 1.001;   %-1.3422
%s.pars('fes_max') = 60;
%s.pars('phi_max') = 13;
%This are basals for normoxia
%s.pars('MRtO2_basal') = s.pars('MRtO2_basal') * 2.5;

%s.pars('MRtCO2_basal') = 0.28;
%xperturbed_Wb_v_s.mat
s.simulation_time = 100;
s.pars('settling_time') = 0;
%s.pars('vO2_b_n') = 0.14; %0.03; %0.4 ;%0.6 %1  %parece que valores un poco más altos estabilizan la PD, pero en ejercicio se vuelve loco.. muy alto
%s.pars('vO2_am_n') = 0.08;
%s.pars('vO2_s_n') = 0.06;
%s.pars('L_sa') = 0.2;
%s.pars('gO2_b') = 0.01;
%s.pars('gO2_am') = 10;
%s.pars('vO2_am_n') = 0.25;

%s.pars('V_tot') = s.pars('V_tot')  + 2200; 
%s.simulation_time = 500;


%s.pars('R_sa') = 0.1;
%s.pars('C_sa') = 0.8;
%s.pars('L_sa') = 0.022;
%s.pars('GTsym') = s.pars('GTsym');



% Ejecutar simulación ODE
[t, x_dot, x_vars, x_keys, index] = s.run_ode_fun(s.model, s.pars, s.init, s.simulation_time, s.dt);

% Organizar resultados
struct_vars = arrange_results(x_dot, x_vars, x_keys, t);

% Guardar resultados
%save(s.simulation_filename, 'struct_vars', 't');


%Plotting - usando cell arrays y char arrays para compatibilidad
%custom_plot('vars_to_show', {['dVE'], struct_vars, t, s.units_table}); 




%figure;
%Crear cell array para parámetros de plotting
%plot_params = {t, s.texp, struct_vars, s.yexp, s.xnames_fitting, s.units_table, 5, 2, '', 'off',  patient_idx, 'actual', [1 2 3 4 5 6 7 8 9 10],state, 0};
%custom_plot('sim_vs_exp', plot_params);

