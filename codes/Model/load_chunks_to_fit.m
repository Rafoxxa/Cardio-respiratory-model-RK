function subset = load_chunks_to_fit(n)

    full_params = load_pars_to_fit();
    subset = get_param_subset(full_params, n);

    function subset = get_param_subset(full_set, index)
        % Returns the k-th 3-element combination from the given set.
        %
        % Usage:
        % subset = get_param_subset(full_set, index)
        %
        % Inputs:
        % - full_set: cell array of parameters (e.g. {"a","b","c",...})
        % - index: the index (1-based) of the combination to return
        %
        % Output:
        % - subset: 3-element cell array with the selected combination
        
            n = numel(full_set);
            
            if n < 3
                error('Need at least 3 elements to generate 3-element combinations.');
            end
        
            combos = nchoosek(1:n, 3);  % indices of all 3-combinations
        
            if index > size(combos,1)
                error('Index exceeds number of combinations (%d)', size(combos,1));
            end
        
            idx = combos(index, :);  % get indices for the desired combination
            subset = full_set(idx);  % extract the subset
        
        end
        
end