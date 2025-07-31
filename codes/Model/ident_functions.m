function  out = ident_functions(mode, fun, args)

  
    if mode == "single" 
        if fun == "compute_corr_matrix"
            corr_matrix = compute_corr_matrix(args.sens_matrix);         
            out = corr_matrix;
        elseif fun == "apply-threshold"   
            [corr_matrix, reduced_pars] = apply_threshold(args.corr_matrix, args.sens_matrix, args.pars_list, args.corr_threshold);
            out = {corr_matrix, reduced_pars};
        end

    elseif mode == "compute-corr"
        corr_matrix = compute_corr_matrix(args.sens_matrix);
        [corr_matrix, reduced_pars, best_idx] = apply_threshold(corr_matrix, args.pars_list, args.corr_threshold);
        [~, idxs] = ismember(reduced_pars, args.pars_list);
        sens_matrix_out = args.sens_matrix(idxs, :);
        out = {corr_matrix, reduced_pars, sens_matrix_out};
    end


function corr_matrix = compute_corr_matrix(sens_matrix)
    S = sens_matrix';

    sigma = 0.05;   
    W = diag(1 ./ sigma^2);  
    FIM = S' * W * S;  

    %intentar pearson sobre S si FIM singular.

    covariance_matrix = pinv(FIM);
    std_devs = sqrt(diag(covariance_matrix));
    corr_matrix = covariance_matrix ./ (std_devs * std_devs');



    %corr_matrix = FIM ./ sqrt(diag(FIM) * diag(FIM)');


end

function [C_final, par_list_final, best_param_idx] = apply_threshold(corr_matrix, par_list, threshold)
   
    C_final = corr_matrix;
    par_list_final = par_list;
    
    % Check condition: any off-diagonal element > threshold
    thresh = threshold;
    over_thresh = abs(C_final) > thresh;
    over_thresh(logical(eye(size(C_final)))) = false;
    
    while any(over_thresh(:))
        % Find all indices of over-threshold values in upper triangle
        [row_idx, col_idx] = find(triu(over_thresh, 1));
        involved_params = unique([row_idx; col_idx]);
    
        best_score = Inf;
        best_score2 = Inf;
        best_param_idx = -1;
    
        % Special case: only one over-threshold pair
        if length(row_idx) == 1
            i = row_idx(1);
            j = col_idx(1);
            total_corr_i = sum(abs(C_final(i, :)));
            total_corr_j = sum(abs(C_final(j, :)));
    
            if total_corr_i >= total_corr_j
                best_param_idx = j;
            else
                best_param_idx = i;
            end
        else
            % General case: evaluate each involved parameter
            for k = 1:length(involved_params)
                idx_remove = involved_params(k);
    
                % Simulate removing the parameter
                C_temp = C_final;
                C_temp(idx_remove, :) = [];
                C_temp(:, idx_remove) = [];
    
                % Count how many values exceed threshold after removal
                over_temp = abs(C_temp) > thresh;
                over_temp(logical(eye(size(C_temp)))) = false;
    
                score = sum(over_temp(:)); % number of over-threshold pairs
                score2 = sum(abs(C_temp(over_temp))); % total magnitude of those correlations
    
                % Choose the parameter that minimizes the score (and then magnitude)
                if score < best_score || (score == best_score && score2 < best_score2)
                    best_score = score;
                    best_score2 = score2;
                    best_param_idx = idx_remove;
                end
            end
        end
    
        % Remove the selected parameter
        C_final(best_param_idx, :) = [];
        C_final(:, best_param_idx) = [];
        par_list_final(best_param_idx) = [];

    
        % Update condition
        over_thresh = abs(C_final) > thresh;
        over_thresh(logical(eye(size(C_final)))) = false;
    end
    end
    
end