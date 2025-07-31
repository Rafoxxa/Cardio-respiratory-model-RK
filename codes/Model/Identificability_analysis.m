
folder_name = "../Sens_analysis/";
file_load_name = "SensMatrix_07-04-2025.mat";
sens_name = strcat(folder_name, file_load_name); 
sens_out_struct = load(sens_name);
sens_final_time_matrix = sens_out_struct.sens_final_time_matrix;
pars_to_sens = sens_out_struct.pars_to_sens;


S = sens_final_time_matrix';
S = sens_final_time_matrix(mask,:)';
%S = S.*(S > 40);


%FIM computation with arbitrary noise level 
sigma = 0.05;   
W = diag(1 ./ sigma^2);  
FIM = S' * W * S;  

%FIM = FIM(logical(sum(FIM>0,2)),logical(sum(FIM>0,2)));



% Condition number (check for identifiability issues)
cond_FIM = cond(FIM);
disp(['Condition number: ', num2str(cond_FIM)]);



%Correlation matrix for parameters
correlation_matrix = FIM ./ sqrt(diag(FIM) * diag(FIM)');

FIM_ = FIM;
correlation_matrix_ = correlation_matrix;
correlation_matrix_cumulated = sum( abs(correlation_matrix_) > 1, 2);


%code for best parameters
while sum(correlation_matrix_cumulated) > 0
    correlation_matrix_cumulated = sum( abs(correlation_matrix_) > 1, 2);
    [~, idx_max] = max(correlation_matrix_cumulated);
    [rows, cols] = size(FIM_);
    indexes = 1:rows;
    new_idxs = indexes(indexes ~= idx_max);
    FIM_ = FIM_(new_idxs, new_idxs);
    correlation_matrix_ = FIM_ ./ sqrt(diag(FIM_) * diag(FIM_)');
end


%

%correlation_matrix = abs(correlation_matrix) > 0.9;
disp('Correlation matrix:');


%[pars, init, taus] = load_global_easy();
%pars_keys = pars.keys;
%%%
figure;
imagesc(abs(correlation_matrix));

colorbar; % Add a color scale
%caxis([-1, 1]); % Normalize between -1 and 1
colormap jet; % Use a color map for better visualization

% Set axis labels
xticks(1:length(pars_to_sens_red));
yticks(1:length(pars_to_sens_red));
xticklabels(pars_to_sens_red);
yticklabels(pars_to_sens_red);
xtickangle(45);

%Cramer Rao Intervals
param_variances = diag(inv(FIM));
param_stddevs = sqrt(param_variances);
disp('Parameter standard deviations:');
disp(param_stddevs);


%data = randn(1000,1); % Example data
%histogram(data, 'Normalization', 'probability'); % Normalize to probability
%yt = yticks; % Get current y-tick values
%yticklabels(string(yt * 100) + "%"); % Convert y-ticks to percentages
%ylabel('Percentage'); % Update y-axis label


zero_rows = all(S == 0, 2); 
find(zero_rows)
zero_cols = all(S == 0, 1); 
find(zero_cols)
S(zero_rows,:) = [];
S(:,zero_cols) = [];

data = correlation_matrix;
threshold = 0.95;
%function reduced_data = reduce_correlation(data, threshold)
    % Reduce dataset dimensions until no correlation exceeds a given threshold.
    % Inputs:
    %   - data: A matrix where rows are samples and columns are variables.
    %   - threshold: Maximum allowed absolute correlation.
    % Outputs:
    %   - reduced_data: The dataset with reduced dimensions.
    
    num_vars = size(data, 2);
    vars_to_keep = true(1, num_vars); % Track which variables to keep
    
    while true
        % Compute absolute correlation matrix
        R = abs(corrcoef(data(:, vars_to_keep))); 
        
        % Find the maximum correlation (excluding diagonal)
        R(logical(eye(sum(vars_to_keep)))) = 0; % Remove self-correlation (diagonal)
        [max_corr, idx] = max(R(:)); 
        
        if max_corr <= threshold
            break; % Stop if all correlations are below the threshold
        end
        
        % Get indices of the most correlated pair
        [i, j] = ind2sub(size(R), idx);
        
        % Compute mean absolute correlation of each variable
        avg_corr = mean(R, 2); 
        
        % Remove the variable with the highest average correlation
        all_indices = find(vars_to_keep);
        if avg_corr(i) > avg_corr(j)
            vars_to_keep(all_indices(i)) = false;
        else
            vars_to_keep(all_indices(j)) = false;
        end
    end

    % Return reduced dataset
    reduced_data = data(:, vars_to_keep);
    final_corr = corrcoef(reduced_data);
    pars_to_keep = pars_to_sens_red(vars_to_keep);
%end
figure;
imagesc(abs(final_corr));

colorbar; % Add a color scale
%caxis([-1, 1]); % Normalize between -1 and 1
colormap jet; % Use a color map for better visualization

% Set axis labels
xticks(1:length(pars_to_keep));
yticks(1:length(pars_to_keep));
xticklabels(pars_to_keep);
yticklabels(pars_to_keep);
xtickangle(45);
