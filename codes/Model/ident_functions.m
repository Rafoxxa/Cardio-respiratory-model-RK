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
    p1m = [1.0870    2.9189    4.2307    9.3023  110.4937   48.3194    1/3.0849  176.6589  136.7549  155.4391];
    p4m = [1.0330    2.8000    1.9661    6.0127  111.0803   44.0947    2.7840 * 60   203.6085  132.7027  161.3660];
    p5m = [0.9996    3.1379    2.9470    7.0525  102.3333   50.9686    2.0473 * 60  229.6756  203.3467  212.4784];
    p6m = [0.9346    2.5861    3.6130    9.9826  103.8647   51.6685    2.6860 * 60  183.5269  145.2870  154.6818];
    pmax_averaged = (p1m + p4m + p5m + p6m)/4;
    p = pmax_averaged; %p1m;

    %"dVE", "VT", "TI", "Tresp", "PAO2", "PACO2", "HR", "PS", "PD", "PM"
    %datos
    %{'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal', 'dVE'};
    p1_maximums = [p(5) p(6) p(9) p(8) p(10) p(7) p(3) p(4) p(2) p(1)];
    I = eye(size(sens_matrix, 2));  
    W_different_sigmas = diag((1./(sigma.^2 * p1_maximums)).^1);
    W_single_sigma = I./sigma.^2;
    FIM_d = S' * W_different_sigmas * S;  
    FIM_s = S' * W_single_sigma * S;

    

    %intentar pearson sobre S si FIM singular.

    covariance_matrix_d = inv(FIM_d);
    std_devs_d = sqrt(diag(covariance_matrix_d));
    corr_matrix_d = covariance_matrix_d ./ (std_devs_d * std_devs_d');

    covariance_matrix_s = inv(FIM_s);
    std_devs_s = sqrt(diag(covariance_matrix_s));
    corr_matrix_s = covariance_matrix_s ./ (std_devs_s * std_devs_s');



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