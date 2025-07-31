
close all
clear all
rng(2);  
clc
clear -global delays_global
clear -global all_global
clear -global externals_global

type_of_input = 1;
control_on = 1;
only_plot = 0;
dt = 0.01;
simulation_time = 45;   
pars_from_fitting = 0;
fitting_mat_file = "Fitting-04-12-24.mat";
VO2_external = 1;
patient_idx = 1;
hipoxia_state = "hipoxia";
ascend_state = "ascend";

if type_of_input == 2
    %simulation_time = 1260;
    simulation_time = 32;
elseif type_of_input == 150
    simulation_time = 200; 
end

model = @(varargin) model_hipoxia_vec_test(varargin{:});
run_ode_fun = @(varargin) run_ode_hipoxia_vec(varargin{:});
structure_data_fun = @(varargin) structure_data(varargin{:});

%Loadings

[pars, init, taus] = load_global_easy();

if pars_from_fitting
    pars_struct = load(fitting_mat_file);
    pars = pars_struct.pars;
end

%Input coefficients ()
if VO2_external
    [~, ~, VO2_poly, ~, fO2_poly] = data_preprocessing(patient_idx, hipoxia_state, ascend_state);
    %VO2
    pars("MRO2_poly_0") = VO2_poly(1);
    pars("MRO2_poly_1") = VO2_poly(2);
    pars("MRO2_poly_2") = VO2_poly(3);
    pars("MRO2_poly_3") = VO2_poly(4);
    pars("MRO2_poly_4") = VO2_poly(5);
    pars("MRO2_poly_5") = VO2_poly(6);
    pars("MRO2_poly_6") = VO2_poly(7);
    pars("MRO2_poly_7") = VO2_poly(8);
    pars("MRO2_poly_8") = VO2_poly(9);
    
    %VCO2
    pars("fiO2_poly_0") = fO2_poly(1);
    pars("fiO2_poly_1") = fO2_poly(2);
    pars("fiO2_poly_2") = fO2_poly(3);
    pars("fiO2_poly_3") = fO2_poly(4);
    pars("fiO2_poly_4") = fO2_poly(5);
end


%
pars('tau') = taus('tau_gases');
pars('dt') = dt; 
pars('type_of_input') = type_of_input;
units_table = readtable("variables_units.xlsx");

if ~only_plot
    preloaded_vars = load('../simulations_saves/90sec_simulation.mat');
    init_values_loaded = preloaded_vars.x_vars(:, end);
    init_keys_loaded = fieldnames(preloaded_vars.struct_vars);

    init_keys = init.keys;


    % Redefine the values in the dictionary using the new vector
    for i = 1:length(init_values_loaded)
         if init_keys_loaded(i) ~= "vO2" && init_keys_loaded(i) ~= "PAO2" && init_keys_loaded(i) ~= "P_1O2" && init_keys_loaded(i) ~= "P_2O2" && init_keys_loaded(i) ~= "P_3O2" && init_keys_loaded(i) ~= "P_4O2" && init_keys_loaded(i) ~= "P_5O2" && init_keys_loaded(i) ~= "MRtO2"
             key_i = init_keys(init_keys == init_keys_loaded(i));
             init(key_i)= init_values_loaded(i);
         end
     end
    
end    



num_values = 20;
studied_param_name = "GTsym";
% Rsa_values = linspace(pars("R_sa") * 1, pars("R_sa") * 1.5, num_values);
% Rsa_Csa_ratio = 0.1;
%Pn_values = linspace(pars("P_n") * 0.8, pars("P_n") * 1.5, num_values);
studied_param_values = linspace(pars(studied_param_name) * 0.1, pars(studied_param_name) * 1.8, num_values);
%studied_param_name = "T0";
% K_E_lv_values = linspace(pars("K_E_lv") * 1, pars("K_E_lv") * 1.5, num_values);
% fab_0_values = linspace(pars("fab_0") * 1, pars("fab_0") * 1.5, num_values);
%Wb_h_s_values = linspace(pars("Wb_h_s") * 0.1, pars("Wb_h_s") * 1, num_values);
%T0_values = linspace(pars("T0") * 0.1, pars("T0") * 1, num_values);
%GTvagal_values = linspace(pars("GTvagal") * 0.1, pars("GTvagal") * 1, num_values);

%studied_param_values = T0_values;




% Initialize the results matrix
results_matrix = zeros(num_values, length(init.keys) + 3);
init_keys = init.keys;
init_keys_modified = [init_keys; "PS"; "PM"; "PD"];

% Run the simulations in parallel
for i = 1:num_values
    % Create a copy of the parameters for this iteration
    pars_i = pars;
    
    % Update the parameter value
    % pars_i("R_sa") = Rsa_values(i);
    % pars_i("C_sa") = Rsa_values(i)/Rsa_Csa_ratio;
    
    pars_i(studied_param_name) = studied_param_values(i);

     %pars_i("Wb_h_s") = Wb_h_s_values(i); 
     %pars_i("T0") = T0_values(i); 
     %pars_i("GTvagal") = GTvagal_values(i); 
    % pars_i("V_u_am_v_0") = V_u_am_v_0_values(i); 
    % pars_i("V_u_s_v_0") = V_u_s_v_0_values(i); 
    
    
    % Run the model
    
        
        try
        all_global = zeros(15, round(10 * simulation_time/dt) + 1) + 0;  %This array saves all the data used for delays and for integration 

        [t, x_dot, x_vars, ~, index] = run_ode_fun(model, pars_i, init, taus, simulation_time, dt, control_on);
        
        % Extract P_sa from x_vars
        P_sa_index = find(strcmp(init_keys, 'P_sa'));
        P_sa = x_vars(P_sa_index, :);
        
        % Compute PS, PM, and PD using compute_presion
        [PM, PS, PD] = compute_presion(P_sa, t);
        
        % Include PS, PM, and PD in x_vars
        x_vars = [x_vars; PS; PM; PD];
        avg_values = mean(x_vars(:, end-14:end), 2);
        catch
    
        avg_values = zeros(size(init_keys_modified));
        end
    
    
    % Store the results
    results_matrix(i, :) = avg_values';
end

mask = sum(results_matrix ~= 0, 2) > 0;
results_matrix_non_zero = results_matrix(mask, :);





% List of variables to plot
variables_to_plot = ["Theart"]; % Add more variables as needed

% Plot the results matrix for the specified variables with units
figure;
hold on;
for var = variables_to_plot
    variable_index = find(strcmp(init_keys_modified, var));
    %unit = find_unit(units_table, var);
    if ~isempty(variable_index)
        plot(studied_param_values(mask), results_matrix_non_zero(:, variable_index), 'DisplayName', var);
    else
        disp(['Variable ', var, ' not found in x_keys']);
    end
end
xlabel(studied_param_name);
ylabel('Values');
title(['Variation of variables with ', studied_param_name]);
legend show;
hold off;

variables_to_plot = ["PS", "PD", "mean_P_sa"]; % Add more variables as needed

% Plot the results matrix for the specified variables with units
figure;
hold on;
for var = variables_to_plot
    variable_index = find(strcmp(init_keys_modified, var));
    %unit = find_unit(units_table, var);
    if ~isempty(variable_index)
        plot(studied_param_values(mask), results_matrix_non_zero(:, variable_index), 'DisplayName', var);
    else
        disp(['Variable ', var, ' not found in x_keys']);
    end
end
xlabel(studied_param_name);
ylabel('Values');
title(['Variation of variables with ', studied_param_name]);
legend show;
hold off;



function unit = find_unit(table, var)
    row_index = strcmp(table.Variable, var);
    unit = table.MeasureUnit{row_index};    
end


function [pm, ps, pd] = compute_presion(presion, t)
    % Crear un vector de tiempo basado en el tamaño de 'presion'
    n = length(presion);  % Número de puntos en la curva
    tiempo = t;  % Suponiendo un paso de tiempo uniforme
    
    % Derivadas de la curva de presión
    dp = diff(presion)./diff(tiempo);  % Primera derivada
    d2p = diff(dp)./diff(tiempo(1:end-1));  % Segunda derivada
    
    % Encontrar picos sistólicos (máximos locales) y diastólicos (mínimos locales)
    sistole_indices = find(dp(1:end-1) > 0 & dp(2:end) <= 0 & d2p < 0) + 1; % Máximos locales
    diastole_indices = find(dp(1:end-1) < 0 & dp(2:end) >= 0 & d2p > 0) + 1; % Mínimos locales
    
    % Extraer valores de presión sistólica y diastólica
    presion_sistolica = presion(sistole_indices);
    presion_diastolica = presion(diastole_indices);
    
    % Interpolación de los valores encontrados para tener un vector continuo
    ps = interp1(tiempo(sistole_indices), presion_sistolica, tiempo, 'linear', 'extrap');
    pd = interp1(tiempo(diastole_indices), presion_diastolica, tiempo, 'linear', 'extrap');
    
    % Calcular la presión media
    pm = trapz(tiempo, presion) / n;  % Integral dividida por la longitud del vector
    pm = ones(size(pd)) * pm;
    %[pm, ps, pd] = compute_presion(presion);
    
    end
