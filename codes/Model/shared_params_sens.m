vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");


[pDim, vDim, tDim] = size(STensor);

perturbed = zeros(pDim, vDim, tDim);

%[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1,'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'time_from_data', 1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
c= 0;
figure;
for v = 1:vDim
    

    for p = 1:pDim
        
        perturbed(p,v,:) = STensor(p,v,:)/(100000) + rbt_fixed(1, v,:);
        

        if sum(sum(sum(perturbed))) == 0
            disp('alerta')
        end

        % Graficar el resultado de esta combinación (descomenta si quieres)
        hold on;
        subplot(4, 2, v);
        plot(squeeze(perturbed(p,v,:)));
        title(sprintf('Perturbed p=%d, v=%d', p, v));
        xlabel('t'); ylabel('value');
        
    end
    hold on;
    subplot(4, 2, v);
    plot(squeeze(rbt_fixed(1,v,:)), 'r', 'LineWidth', 2.5);
    hold off;
    disp('end of subjetc')
end


disp('ende');


%% ===============================
%       CONFIGURACIÓN GENERAL
% ================================
T_target = 45050;   % tamaño fijo deseado en el tiempo
vars_needed = {'PAO2','PACO2','pd','ps','pm','Theart','VTidal','dVE'};
extra_vars = {'pm', 'ps', 'pd', 'VTidal', 'BF'};

%% ===============================
%       FUNCION AUXILIAR
% ================================
extract_and_reduce = @(results_base_) ...
    extract_vars_fixed_time(results_base_, vars_needed, extra_vars, setup);

%% ===============================
%   1) CARGAR NORMOXIA E HIPÓXIA
% ================================
loadedN = load('results_base_local_normoxia.mat');
rb_normoxia = extract_and_reduce(loadedN.results_base_);

loadedH = load('results_base_local_hipoxia.mat');
rb_hipoxia = extract_and_reduce(loadedH.results_base_);

%% ===============================
%   2) CONCATENAR EN EL TIEMPO
% ================================
rbt = cat(3, rb_normoxia, rb_hipoxia);  
% tamaño: [271 × 8 × T_total]

%% ===============================
%   3) INTERPOLAR A T_target
% ================================
[NS, NV, T_orig] = size(rbt);
rbt_fixed = zeros(NS, NV, T_target);

t_orig = linspace(0,1,T_orig);
t_new  = linspace(0,1,T_target);

for s = 1:NS
    for v = 1:NV
        rbt_fixed(s,v,:) = interp1(t_orig, squeeze(rbt(s,v,:)), t_new, 'linear');
    end
end

%% ===============================
%   4) rbt_fixed es el resultado final
% ================================
size(rbt_fixed)   % --> [271 8 45050]



%% SENSITIVITY ANALYSIS
%Set up
% sensMatrices = {};
% cnt = 0;
% patients = [1,4,5,6];
% window_name = ["normoxia", "hipoxia"];
% 
% for i = 1:length(patients)
%     sens_matrix = 0;
%     for j = 1:length(window_name)   
%         weight = 0.33;
%         mode = "loading"; %"loading"        
%         setup = set_up("sens", patients(i), "hipoxia", "mix");
%         basePath = sprintf('../Sens_analysis/%d/%s/', patients(i), window_name(j));
%         latestFullPath = getLatestFittingDateStr(basePath);
%         setup.sensitivity_load_filename = latestFullPath;
%         %{'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal', 'dVE'};
%         setup.variables_of_interest = {'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal', 'dVE'}; %{'pd', 'ps', 'pm', 'Theart'};  %selecting only cardiovascular parameters
%         setup.idx_variable_of_interest = [1 2 3 4];
% 
%         LSA_output = sens_functions(mode, "-", setup);
%         sens_matrix_ = LSA_output{1};
%         pars_to_sens = LSA_output{2}; 
%         if window_name(j) == "normoxia"
%             sens_matrix_ = mapping_2017_to_2025_pars_indexes(pars_to_sens, sens_matrix_);
%             sens_matrix_ = weight * sens_matrix_;
%         else
%             sens_matrix_ = (1 - weight)/(length(window_name) - 1) * sens_matrix_;
%         end
%         sens_matrix = sens_matrix_  + sens_matrix;
% 
%     end
% 
% 
%     %sens_matrix(:,1:6) = 0;
%     sensMatrices{i} = sens_matrix;
% end
% 
% 
% sens_threshold = 0.23;
% nCases = numel(sensMatrices);
% 
% % 1) detectar orientación
% rows = zeros(1,nCases); cols = zeros(1,nCases);
% for i = 1:nCases
%     [r,c] = size(sensMatrices{i});
%     rows(i) = r; cols(i) = c;
% end
% 
% % heurística: si la dimensión media de filas es mayor que la de columnas,
% % asumimos que las matrices están en formato (parametros x variables) => hay que transponer
% transposeFlag = mean(rows) > mean(cols);
% 
% % 2) construir versión consistente (variables x parametros)
% Scell = cell(1,nCases);
% for i = 1:nCases
%     S = sensMatrices{i};
%     if transposeFlag
%         S = S';   % ahora filas = variables, columnas = parametros
%     end
%     Scell{i} = S;
% end
% SS = (Scell{1} + Scell{2} + Scell{3} + Scell{4})/4;
% SSS = sum(SS, 1);
% 
% [f,xi] = ksdensity(SSS);
% plot(xi, f, 'LineWidth', 2);
% xlabel('Value');
% ylabel('Density');
% title('Kernel Density Estimate (Smooth Histogram)');
% hold on
% p80 = prctile(SSS, 80);
% xline(p80, '--r', 'Percentile 80');
% 
% disp('halt');
% 
% 
% % sanity check: verificar que todas las matrices tengan mismas columnas (mismo nº parámetros)
% nParams = size(Scell{1},2);
% for i = 2:nCases
%     if size(Scell{i},2) ~= nParams
%         error('Número de parámetros inconsistente entre casos después de transponer.');
%     end
% end
% 
% % 3) máscara de parámetros activos por caso
% activeMask = false(nCases, nParams);
% for i = 1:nCases
%     S = Scell{i};
%     % un parámetro "sobrevive" en este caso si en alguna variable |S| > threshold
%     activeMask(i,:) = any(abs(S) > sens_threshold, 1);
% end
% 
% % 4) elegir el caso más restrictivo (el que deja menos parámetros)
% survivorsCount = sum(activeMask,2);
% [~, idxMin] = min(survivorsCount);
% finalMask = activeMask(idxMin,:);   % máscara lógica 1 x nParams
% 
% fprintf('Parametros originales: %d\n', nParams);
% fprintf('Parametros sobrevivientes por caso: %s\n', mat2str(survivorsCount'));
% fprintf('Caso mas restrictivo: %d -> quedan %d parametros\n', idxMin, sum(finalMask));
% 
% % 5) reducir todas las matrices y apilar (stack)
% reduced = cell(1,nCases);
% for i = 1:nCases
%     reduced{i} = Scell{i}(:, finalMask); % filas = variables, cols = parametros finales
% end
% S_stack = vertcat(reduced{:});  % filas = sum(n_i), columnas = n_final_params
% 
% common_averaged = 0;
% if common_averaged

%% Load results_base

%% 1) Cargar archivo
load('all_last_structs.mat');  % Debe traer "data"

%subject_names = {'subj1','subj4','subj5','subj6'};
subject_names = {'subj5'};
%cond_names    = {'normoxia','hipoxia'};
cond_names    = {'normoxia', 'hipoxia'};

%vars_needed = {'PAO2','PACO2','pd','ps','pm','Theart','VT','dVE'};
vars_needed = {'PAO2','PACO2','pd','ps','Theart','VT','dVE'};
%vars_needed = {'pd','ps','pm'};
n_target_points = 500;  % largo temporal deseado por bloque

all_blocks = [];  % Aquí concatenaremos

for s = 1:numel(subject_names)
    subj = subject_names{s};

    for c = 1:numel(cond_names)
        cond = cond_names{c};

        % Extraer struct_vars de este sujeto/condición
        S = data.(subj).(cond).struct_vars;

        % --- 2) Extraer solo las 8 variables en el orden deseado ---
        block = zeros(numel(vars_needed), n_target_points);

        for v = 1:numel(vars_needed)
            varName = vars_needed{v};

            if isfield(S, varName)
                vec = S.(varName)(:);  % asegurar columna
            else
                warning('Variable %s no encontrada en %s.%s, rellenando con ceros.', ...
                        varName, subj, cond);
                if isfield(S, varName)
                    vec = S.(varName)(:);  % asegurar columna
                else
                    warning('Variable %s no encontrada en %s.%s, rellenando con ceros.', ...
                            varName, subj, cond);
                    fn = fieldnames(S);               % obtener nombres de campos
                    vec = zeros(length(S.(fn{1})),1); % usar el primero como referencia de largo
                end

            end

            % --- 3) Submuestrear a 500 puntos (ventana centrada) ---
            n_timepoints = length(vec);
            idx = round(linspace(1, n_timepoints, n_target_points));
            win = max(1, round(n_timepoints / n_target_points / 2));

            vec_sub = zeros(1, n_target_points);
            for k = 1:n_target_points
                low  = max(1, idx(k) - win);
                high = min(n_timepoints, idx(k) + win);
                vec_sub(k) = mean(vec(low:high));
            end

            block(v, :) = vec_sub;
        end

        % --- 4) Concatenar en el eje del tiempo ---
        all_blocks = [all_blocks, block];

    end
end

%% Resultado
results_base = all_blocks;  % 8 x 4000
disp(size(results_base))  % debería mostrar [8 4000]

results_base = reshape(results_base, [1, size(results_base,1), size(results_base,2)]);

% Replicar 270 veces en la 1ra dimensión
results_base = repmat(results_base, [271, 1, 1]);







%%





loadedTensor = load('../Sens_analysis/STensor_16-11-2025.mat');
STensor = loadedTensor.STensor.tensor;
STensor(:,[7,8],:) = [];
%STensor(:,[1,2,6,7,8,9, 10],:) = [];

pars_to_sens = loadedTensor.STensor.parameters;

fprintf('\n--- Loaded Sensitivity Tensor ---\n');
fprintf('Tensor size: %d params x %d vars x %d timepoints\n', size(STensor));

%% === 2) Compute sensitivity matrix from tensor ===
N = size(STensor,3);
sens_matrix_r = (1/N) * sqrt(sum( (STensor./(rbt_fixed + 0.01)).^2, 3));  % P x V

fprintf('Sensitivity matrix computed: %d x %d\n', size(sens_matrix_r));

%% === 3) Wrap into cell array to emulate previous multi-subject pipeline ===
%sensMatrices = { sens_matrix };


    %sens_matrix = sensMatrices;%simple_merge_S(sensMatrices);
    %subplot(2,2,1);
    setup.variables_of_interest = {'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'VTidal', 'dVE'};
    setup.idx_variable_of_interest = [1 2 3 4 5 6 7 8];
    custom_plot("LSA-plot", {sens_matrix_r__, pars_to_sens, setup.variables_of_interest, setup.idx_variable_of_interest});
    setup.variables_of_interest = {'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'VTidal', 'dVE'};
    %setup.variables_of_interest = {'pd', 'ps', 'pm'};
    setup.idx_variable_of_interest = [1 2 3 4 5 6 7 8];
    %setup.idx_variable_of_interest = [1 2 3];
    custom_plot("LSA-plot", {sens_matrix_r, pars_to_sens, setup.variables_of_interest, setup.idx_variable_of_interest});
    title("A");
    hold on;
    
    % % Apply sensitivity threshold
    sens_threshold = 8.5;%0.375;%0.6931; %for cardiovascular %0.4; %for all %0.1;
    idx_variable_of_interest = [1:numel(setup.variables_of_interest)];
    setup.idx_variable_of_interest = idx_variable_of_interest;
    setup.sens_threshold = sens_threshold;
    setup.sens_matrix = sens_matrix_r__;
    setup.pars_to_sens = pars_to_sens;
    out = sens_functions("single", "sens_threshold", setup);
    sens_reduced = out{1};
    pars_to_sens_reduced = out{2};
    %subplot(2,2,2);
    figure;
    custom_plot("LSA-plot", {sens_reduced, pars_to_sens_reduced,setup.variables_of_interest, idx_variable_of_interest});
    title("B");
    hold on;
    
    
    % %Apply fitler by class
    setup.sens_final_time_matrix = sens_reduced; %sens_matrix;%
    setup.pars_to_sens = pars_to_sens_reduced;%pars_to_sens;%
    filtered_output = sens_functions("single", "filter_params_by_class", setup);
    sens_matrix_filtered = filtered_output{1};
    pars_to_sens_filtered = filtered_output{2};
    subplot(2,2,3);
    custom_plot("LSA-plot", {sens_matrix_filtered, pars_to_sens_filtered, setup.variables_of_interest, idx_variable_of_interest});
    title("C")
    hold on;
    
    % %% IDENTIFIABILITY ANALYSIS
    % 
    ident_args.sens_matrix = sens_reduced; %sens_matrix_filtered;
    ident_args.pars_list = pars_to_sens_reduced; %pars_to_sens_filtered;
    ident_args.corr_threshold = 0.95; % Set the correlation threshold for the analysis
    IDENT_output = ident_functions("compute-corr", "-", ident_args);
    subplot(2,2,4);
    custom_plot("ident-plot", IDENT_output);
    title("D")
%end



% [sharedParamsIndices, sensMatricesOut] = selectSharedSensitiveParamIndices(sensMatrices, 0.32);
% 
% pars_filtered_not_ident = {};
% pars_filtered_ident = {};
% sensMatricesFiltered = {};
% 
% for i = 1:length(patients)
%     setup = set_up("sens", patients(i), "hipoxia", "mix");
%     pars_keys = setup.pars.keys;
%     sharedParams = pars_keys(sharedParamsIndices);
%     setup.sens_final_time_matrix = sensMatricesOut{i};
%     setup.pars_to_sens = sharedParams;
%     filtered_output = sens_functions("single", "filter_params_by_class", setup);
%     sens_matrix_filtered = filtered_output{1};
%     pars_to_sens_filtered = filtered_output{2};
%     sensMatricesFiltered{i} = sens_matrix_filtered;
%     pars_filtered_not_ident{i} = pars_to_sens_filtered;
% 
%     ident_args.sens_matrix = sens_matrix_filtered; %sens_reduced; %sens_matrix_filtered;
%     ident_args.pars_list = pars_to_sens_filtered; %pars_to_sens_reduced; %pars_to_sens_filtered;
%     ident_args.corr_threshold = 0.98; % Set the correlation threshold for the analysis
%     IDENT_output = ident_functions("compute-corr", "-", ident_args);
%     pars_filtered_ident{i} = IDENT_output{2};
% 
% 
% 
% 
% 
% end


%custom_plot("ident-plot", IDENT_output);


function latestFullPath = getLatestFittingDateStr(basePath)
    %GETLATESTFITTINGDATESTR Returns the full path containing the latest dd-MM-yyyy date
    % found in folders or files (e.g., 'Fitting-dd-MM-yyyy' or 'Sens_matrix_dd-MM-yyyy'),
    % excluding today's date.

    % Get list of all items in basePath
    items = dir(basePath);
    items = items(~ismember({items.name}, {'.', '..'}));  % Exclude . and ..

    validDates = datetime.empty;
    fullPaths = {};

    % Regular expression to find dd-MM-yyyy
    datePattern = '\d{2}-\d{2}-\d{4}';

    for k = 1:length(items)
        name = items(k).name;
        tokens = regexp(name, datePattern, 'match');
        
        if ~isempty(tokens)
            dateStr = tokens{1};  % Take first match
            try
                Dt = datetime(dateStr, 'InputFormat', 'dd-MM-yyyy');
                if Dt ~= datetime('today')  % Exclude today's date
                    validDates(end+1) = Dt;
                    fullPaths{end+1} = fullfile(basePath, name);  % Save full path
                end
            catch
                % Skip invalid dates
            end
        end
    end

    if ~isempty(validDates)
        [~, idx] = max(validDates);
        latestFullPath = fullPaths{idx};  % Return full path including the date
    else
        latestFullPath = '';
        warning('No valid dated items found in %s (excluding today)', basePath);
    end
end

function [SensMatrix] = simple_merge_S(sensMatrices)
    SensMatrix = zeros(size(sensMatrices{1}));
    

    for idx = 1:numel(sensMatrices)
        matrix = sensMatrices{idx};
        SensMatrix = SensMatrix + matrix.^2;
                
    end
    SensMatrix = sqrt(SensMatrix)/numel(sensMatrices);
end




function [sharedParamIndices, sensMatricesOut] = selectSharedSensitiveParamIndices(sensMatrices, threshold)
% sensMatrices: cell array de N matrices de sensibilidad {N x 1}, cada una [numVars x numParams]
% threshold: valor mínimo de sensibilidad normalizada (por ej. 0.2)

    N = numel(sensMatrices);
    numParams = size(sensMatrices{1}, 1);
    sharedParamsMask = true(numParams, 1); % inicializar como todos verdaderos

    for i = 1:N
        S = sensMatrices{i}; % [numVars x numParams]
        
        % 1. Normalización min-max por parámetro
        minVal = min(S, [], 2);
        maxVal = max(S, [], 2);
        range = maxVal - minVal;
        range(range == 0) = 1; % evitar división por cero
        S_norm = (S - minVal) ./ range;
        

        % 2. Sumar sensibilidad por parámetro
        sumSens = sum(S_norm, 2)/size(S,2);

        % 3. Crear máscara de parámetros sensibles (sobre threshold)
        paramMask = sumSens > threshold;

        % 4. Intersección acumulativa
        sharedParamsMask = sharedParamsMask & paramMask;


        %5. Acá habría que hacer un análisis de identificabilidad

    end


    % Convertir máscara lógica a índices
    sharedParamIndices = find(sharedParamsMask);
    sensMatricesOut = {};
    for i = 1:N
        S = sensMatrices{i};
        sensMatricesOut{i} = S(sharedParamIndices, :);
    end

end

function [sens_matrix_to_2025] = mapping_2017_to_2025_pars_indexes(pars2sens2017, sens_matrix2017)
    [pars25, init, taus] = load_global_easy_2025();
    pars_keys25 = pars25.keys;
    pars_not_to_sens25 = load_pars_not_to_sens_2025();
    pars_free2move25 = setdiff(pars_keys25, pars_not_to_sens25);

    sens_matrix_to_2025 = zeros(size(pars_free2move25,1),size(sens_matrix2017,2)); 

    for par_index = 1:size(sens_matrix2017,1)
        par_string_on_index2017 = pars2sens2017(par_index);
        index_of_that_par_in2025 = find(pars_free2move25 == par_string_on_index2017);
        try
            sens_matrix_to_2025(index_of_that_par_in2025, :) = sens_matrix2017(par_index, :);  
        catch
            disp('no existence!')
        end
    end
end


