vectorize_dicts("run_ode.m", "model_basic.m", "run_ode_vec_hipoxia.m", "model_vec_hipoxia.m");

%% SENSITIVITY ANALYSIS
%Set up
patient_idx = 5;
mode = "loading"; %"loading"
%mode = "saving";
%mode = "p-space";
%setup = set_up("sens", 1, "hipoxia", "mix", "requestedDate", "28-04-2025");
setup = set_up("sens", patient_idx, "hipoxia", "mix");
%setup = set_up("sens", 1, "hipoxia", "mix", "epsilon", 1, "params_sample_size",10);  
basePath = sprintf('../Sens_analysis/%d/', patient_idx);
latestFullPath = getLatestFittingDateStr(basePath);

setup.sensitivity_load_filename = latestFullPath;
% setup.p_space_files = ["../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_10.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_9.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_8.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_7.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_6.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_5.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_4.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_3.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_2.mat", ...
%                         "../Sens_analysis/5/SensMatrix_15-05-2025.mat_pspace_sampling_1.mat"];
%Run
LSA_output = sens_functions(mode, "-", setup);
if mode == "saving"
    sens_matrix = LSA_output{1};
    pars_to_sens = LSA_output{5};
    %sensitivities = LSA_output{2};
    %error = LSA_output{3};
    %perturbed = LSA_output{4};
elseif mode == "loading"    
    sens_matrix = LSA_output{1};
    pars_to_sens = LSA_output{2};   

end
%custom_plot("LSA-plot", {sens_matrix, pars_to_sens, setup.variables_of_interest, setup.idx_variable_of_interest});



% % Apply sensitivity threshold
sens_threshold = 0.5;
idx_variable_of_interest = [1:numel(setup.variables_of_interest)];
setup.idx_variable_of_interest = idx_variable_of_interest;
setup.sens_threshold = sens_threshold;
setup.sens_matrix = sens_matrix;
out = sens_functions("single", "sens_threshold", setup);
sens_reduced = out{1};
pars_to_sens_reduced = out{2};
%custom_plot("LSA-plot", {sens_reduced, pars_to_sens_reduced,setup.variables_of_interest, idx_variable_of_interest});

% %Apply fitler by class
setup.sens_final_time_matrix = sens_reduced; %sens_matrix;%
setup.pars_to_sens = pars_to_sens_reduced;%pars_to_sens;%
filtered_output = sens_functions("single", "filter_params_by_class", setup);
sens_matrix_filtered = filtered_output{1};
pars_to_sens_filtered = filtered_output{2};
%custom_plot("LSA-plot", {sens_matrix_filtered, pars_to_sens_filtered, setup.variables_of_interest, idx_variable_of_interest});
% 
% %% IDENTIFIABILITY ANALYSIS
% 
ident_args.sens_matrix = sens_matrix_filtered; %sens_reduced; %sens_matrix_filtered;
ident_args.pars_list = pars_to_sens_filtered; %pars_to_sens_reduced; %pars_to_sens_filtered;
ident_args.corr_threshold = 0.6; % Set the correlation threshold for the analysis
IDENT_output = ident_functions("compute-corr", "-", ident_args);
custom_plot("ident-plot", IDENT_output);

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
