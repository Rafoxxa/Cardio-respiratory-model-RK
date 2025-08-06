vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");


%% SENSITIVITY ANALYSIS
%Set up
sensMatrices = {};
cnt = 0;
patients = [1,4,5,6];
for i = 1:length(patients)
    
    
    mode = "loading"; %"loading"

    setup = set_up("sens", patients(i), "hipoxia", "mix");
    basePath = sprintf('../Sens_analysis/%d/', patients(i));
    latestFullPath = getLatestFittingDateStr(basePath);
    setup.sensitivity_load_filename = latestFullPath;
    LSA_output = sens_functions(mode, "-", setup);
    sens_matrix = LSA_output{1};
    pars_to_sens = LSA_output{2}; 

    sensMatrices{i} = sens_matrix;
end

sens_matrix = simple_merge_S(sensMatrices);
subplot(2,2,1);
custom_plot("LSA-plot", {sens_matrix, pars_to_sens, setup.variables_of_interest, setup.idx_variable_of_interest});
title("A");
hold on;

% % Apply sensitivity threshold
sens_threshold = 0.33;
idx_variable_of_interest = [1:numel(setup.variables_of_interest)];
setup.idx_variable_of_interest = idx_variable_of_interest;
setup.sens_threshold = sens_threshold;
setup.sens_matrix = sens_matrix;
out = sens_functions("single", "sens_threshold", setup);
sens_reduced = out{1};
pars_to_sens_reduced = out{2};
subplot(2,2,2);
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
% 
% % %% IDENTIFIABILITY ANALYSIS
% % 
ident_args.sens_matrix = sens_matrix_filtered; %sens_reduced; %sens_matrix_filtered;
ident_args.pars_list = pars_to_sens_filtered; %pars_to_sens_reduced; %pars_to_sens_filtered;
ident_args.corr_threshold = 0.8; % Set the correlation threshold for the analysis
IDENT_output = ident_functions("compute-corr", "-", ident_args);
subplot(2,2,4);
custom_plot("ident-plot", IDENT_output);
title("D")



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
