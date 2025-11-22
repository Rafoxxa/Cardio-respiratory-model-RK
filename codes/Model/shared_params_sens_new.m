%% =========================================================================
%%  ANÁLISIS DE SENSIBILIDAD GLOBAL + IDENTIFICABILIDAD (Pipeline completo)
%%  Fecha: Noviembre 2025
%% =========================================================================

clear; clc; close all;

%% -------------------------------------------------------------------------
%% 1. CONFIGURACIÓN INICIAL
%% -------------------------------------------------------------------------
vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");

state = 'normoxia';
patient_idx = 5;
[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1,'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'time_from_data', 1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");

T_target       = 45053;                                      % Longitud temporal fija
vars_needed    = {'PAO2','PACO2','pd','ps','pm','Theart','VTidal','dVE'};
extra_vars     = {'pm', 'ps', 'pd', 'VTidal', 'BF'};

% Función auxiliar para extraer variables y reducir a tiempo fijo
extract_and_reduce = @(results_base_) ...
    extract_vars_fixed_time(results_base_, vars_needed, extra_vars, setup);

%% -------------------------------------------------------------------------
%% 2. CARGA Y PREPROCESADO DE TRAZAS BASE (Normoxia + Hipoxia)
%% -------------------------------------------------------------------------
loadedN = load('results_base_local_normoxia_eps-5_mrfChanges.mat');
rb_normoxia = extract_and_reduce(loadedN.results_base_);

loadedH = load('results_base_local_hipoxia_eps-5_mrfChanges.mat');
rb_hipoxia = extract_and_reduce(loadedH.results_base_);

% Concatenar en el tiempo: [variables × sujetos × tiempo]
rbt = cat(3, rb_normoxia, rb_hipoxia);
rbt_fixed = rbt;

% Interpolación a longitud fija T_target
% [NS, NV, T_orig] = size(rbt);
% rbt_fixed = zeros(NS, NV, T_target);
% 
% t_orig = linspace(0, 1, T_orig);
% t_new  = linspace(0, 1, T_target);
% 
% for s = 1:NS
%     for v = 1:NV
%         rbt_fixed(s,v,:) = interp1(t_orig, squeeze(rbt(s,v,:)), t_new, 'linear');
%     end
% end
% 
% fprintf('rbt_fixed creado → tamaño: %s\n', mat2str(size(rbt_fixed)));

%% -------------------------------------------------------------------------
%% 3. CARGA DEL TENSOR DE SENSIBILIDAD
%% -------------------------------------------------------------------------
loadedTensor = load('../Sens_analysis/STensorEps-5C_20-11-2025.mat');
STensor      = loadedTensor.STensor.tensor;
pars_to_sens = loadedTensor.STensor.parameters;

% Eliminar variables no deseadas (ajustar según necesidad)
STensor(:,[7,8],:) = [];        % Ejemplo: quitar VTidal y dVE del tensor

fprintf('\n--- Tensor de sensibilidad cargado ---\n');
fprintf('Tamaño: %d parámetros × %d variables × %d puntos temporales\n', size(STensor));

%% -------------------------------------------------------------------------
%% 4. CÁLCULO DE MATRIZ DE SENSIBILIDAD GLOBAL (RMS normalizado)
%% -------------------------------------------------------------------------
N = size(STensor, 3);
sens_matrix_r = (1/sqrt(N)) * sqrt( sum( (STensor ./ (rbt_fixed + 0.01)).^2 , 3) );  % [P × V]

fprintf('Matriz de sensibilidad calculada → tamaño: %s\n', mat2str(size(sens_matrix_r)));

[STensor_new, sens_matrix_r, pars_to_sens] = reduce_by_param_names(STensor, sens_matrix_r, pars_to_sens, 'exclude', {'tiny_dt', 'Tsys_0'});

%% -------------------------------------------------------------------------
%% 5. VISUALIZACIÓN INICIAL (Matriz completa)
%% -------------------------------------------------------------------------
setup.variables_of_interest    = {'PAO2','PACO2','pd','ps','pm','Theart','VTidal','dVE'};
setup.idx_variable_of_interest = 1:8;

figure('Name','Sensibilidad global completa','NumberTitle','off');
custom_plot("LSA-plot", {sens_matrix_r, pars_to_sens, ...
                        setup.variables_of_interest, setup.idx_variable_of_interest});
title('A) Matriz de sensibilidad completa');

%% -------------------------------------------------------------------------
%% 6. APLICAR UMBRAL DE SENSIBILIDAD
%% -------------------------------------------------------------------------
sens_threshold = 8500;   % Ajustar según distribución (antes era ~0.4 para cardio)

setup.sens_threshold = sens_threshold;
setup.sens_matrix    = sens_matrix_r;
setup.pars_to_sens   = pars_to_sens;

out                  = sens_functions("single", "sens_threshold", setup);
sens_reduced         = out{1};
pars_to_sens_reduced = out{2};

figure('Name','Después del umbral','NumberTitle','off');
custom_plot("LSA-plot", {sens_reduced, pars_to_sens_reduced, ...
                        setup.variables_of_interest, 1:numel(setup.variables_of_interest)});
title('B) Parámetros sensibles (|S| > threshold)');

%% -------------------------------------------------------------------------
%% 7. FILTRADO POR CLASE DE PARÁMETRO (opcional)
%% -------------------------------------------------------------------------
setup.sens_final_time_matrix = sens_reduced;
setup.pars_to_sens           = pars_to_sens_reduced;

filtered_output        = sens_functions("single", "filter_params_by_class", setup);
sens_matrix_filtered   = filtered_output{1};
pars_to_sens_filtered  = filtered_output{2};

%subplot(2,2,3);
custom_plot("LSA-plot", {sens_matrix_filtered, pars_to_sens_filtered, ...
                        setup.variables_of_interest, 1:numel(setup.variables_of_interest)});
title('C) Filtrado por clase de parámetro');

[STensor_red, pars_red, idx_kept] = reduce_STensor_by_names(STensor_new, pars_to_sens, pars_to_sens_filtered);

%% -------------------------------------------------------------------------
%% 8. ANÁLISIS DE IDENTIFICABILIDAD (Correlación entre parámetros sensibles)
%% -------------------------------------------------------------------------
ident_args.sens_matrix      = sens_reduced;         % o sens_matrix_filtered
ident_args.pars_list        = pars_to_sens_reduced; % o pars_to_sens_filtered
ident_args.corr_threshold  = 0.95;

IDENT_output = ident_functions("compute-corr", "-", ident_args);

subplot(2,2,4);
custom_plot("ident-plot", IDENT_output);
title('D) Identificabilidad (correlación > 0.95)');

%% -------------------------------------------------------------------------
%% 9. (Opcional) VISUALIZACIÓN DE TRAZAS PERTURBADAS - EJEMPLO
%% -------------------------------------------------------------------------
[pDim, vDim, tDim] = size(STensor);
perturbed = zeros(pDim, vDim, tDim);

figure('Name','Ejemplo de perturbación por parámetro','Position',[100 100 1200 800]);
for v = 1:vDim
    subplot(4,2,v);
    hold on;
    plot(squeeze(rbt_fixed(1,v,:)), 'k', 'LineWidth', 2);  % Línea base
    for p = 1:max(5,pDim)  % Solo los 5 parámetros más sensibles como ejemplo
        perturbed = squeeze(STensor(p,v,:))/100000 + squeeze(rbt_fixed(1,v,:));
        plot(perturbed, 'LineWidth', 1.5);
    end
    title(sprintf('Variable %d', v));
    xlabel('Tiempo'); ylabel('Valor');
    hold off;
end

disp('¡Análisis completado con éxito!');


%% =========================================================================
