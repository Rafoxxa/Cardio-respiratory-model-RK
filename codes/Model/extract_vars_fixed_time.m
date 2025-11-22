function rb = extract_vars_fixed_time(results_base_, vars_needed, extra_vars, setup)

[n_vars, n_timepoints] = size(results_base_);

% Interpolación a tamaño fijo (primero reducir)
n_target_points = 500;
idx = round(linspace(1, n_timepoints, n_target_points));
win = max(1, round(n_timepoints / n_target_points / 2));

%results_reduced = zeros(n_vars, n_target_points);
results_reduced = results_base_;

% for v = 1:n_vars
%     vec = results_base_(v, :);
% 
%     vec_sub = zeros(1, n_target_points);
% 
%     for k = 1:n_target_points
%         low  = max(1, idx(k) - win);
%         high = min(n_timepoints, idx(k) + win);
%         vec_sub(k) = mean(vec(low:high));
%     end
% 
%     results_reduced(v, :) = vec_sub;
% end

% Extraer variables
rb = zeros(8, size(results_reduced, 2));

for ki = 1:length(vars_needed)
    if ismember(vars_needed{ki}, extra_vars)
        extra_idx = find(strcmp(extra_vars, vars_needed{ki}));
        var_index = extra_idx + length(setup.init_keys);
    else
        var_index = find(strcmp(setup.init_keys, vars_needed{ki}));
    end
    rb(ki, :) = results_reduced(var_index, :);
end

% Poner dimensión -> [1 x 8 x T]
rb = reshape(rb, [1, size(rb,1), size(rb,2)]);

% Replicar 271 veces
rb = repmat(rb, [272, 1, 1]);
end