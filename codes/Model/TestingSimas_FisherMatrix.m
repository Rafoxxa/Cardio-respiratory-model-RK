
%% Fisher Information Matrix con varias sensibilidades y sujetos

% % === Cargar datos ===
% load('../Sens_analysis/1/SensMatrix_hipoxia_26-09-2025.mat','S_tensor', 'pars_to_sens');  % (nVar x nPar x nTime)
% %load('S1_tensor.mat','S_tensor', 'pars_to_sens');  % (nPar x nVar x nTime)
% S_sens_1 = S_tensor;
% load('../Sens_analysis/4/SensMatrix_hipoxia_26-09-2025.mat','S_tensor');
% %load('S2_tensor.mat','S_tensor');
% S_sens_2 = S_tensor;
% load('../Sens_analysis/5/SensMatrix_hipoxia_26-09-2025.mat','S_tensor');
% %load('S3_tensor.mat','S_tensor');
% S_sens_3 = S_tensor;
% load('../Sens_analysis/6/SensMatrix_hipoxia_26-09-2025.mat','S_tensor');
% %load('S4_tensor.mat','S_tensor');
% S_sens_4 = S_tensor;

%loadedTensor = load('../Sens_analysis/STensor_new_16-11-2025.mat');
%STensor_new = loadedTensor.STensor_new.tensor;
%pars_to_sens = loadedTensor.STensor_new.parameters;
% 
% pars_keep = pars_to_sens_filtered; %pars_to_sens
% 
% 
% 
% % ---- Parámetros que quiero conservar ----
% %pars_keep = {'Kbg', 'P_n', 'PaO2_ac_n', 'Wb_h_s', 'Wb_p_s', 'Wb_v_s', 'Wp_p_s', 'Wp_v', 'Wt_v_s', 'aO2_n', 'f_ab_max', 'kev', 'x_h_s'};
% 
% 
% %pars_keep = pars_to_sens_reduced;%{'Tsys_0', 'gO2_e', 'gO2_p', 'gO2_s', 'kcc_p_s', 'tau_M', 'ub_a1', 'vO2_b_n', 'vO2_e_n', 'vO2_s_n', 'x_p_s'};
% 
% % ---- Encontrar índices de esos parámetros ----
% [found, idx] = ismember(pars_keep, pars_to_sens);
% 
% % Verificar si alguno NO se encontró
% missing = pars_keep(~found);
% if ~isempty(missing)
%     warning('Los siguientes parámetros no se encontraron en STensor_new:');
%     disp(missing);
% end
% 
% % Filtrar solo los encontrados
% idx = idx(found);
% filtered_pars   = pars_to_sens(idx);
% filtered_tensor = STensor_new(idx, :, :);  
% %filtered_tensor(:,[7,8],:) = [];
% 
% %STensor_new = filtered_tensor;
% pars_to_sens = filtered_pars;
pars_values_filtered_ = [];


for index = 1:size(pars_to_sens_filtered,2)
    key = pars_to_sens_filtered{index};    
    pars_values_filtered_(index) = setup.pars(key);
end


pars_to_sens = pars_to_sens_filtered;

filtered_tensor = STensor_new./pars_values_filtered_';






% p_max de cada sujeto (ejemplo)
p1m = [1.0870 2.9189 4.2307 9.3023 110.4937 48.3194 1/3.0849 176.6589 136.7549 155.4391]; %el dato de HR viene en Hz, lo pasé a tiempo (para los 4 sujetos)
p4m = [1.0330 2.8000 1.9661 6.0127 111.0803 44.0947 1/2.7840 203.6085 132.7027 161.3660];
%p5m = [0.9996 3.1379 2.9470 7.0525 102.3333 50.9686 1/2.0473 229.6756 203.3467 212.4784];
p5m = [0.9996 3.1379 2.9470 7.0525 102.3333 50.9686 0.8 229.6756 203.3467 212.4784];
p6m = [0.9346 2.5861 3.6130 9.9826 103.8647 51.6685 1/2.6860 183.5269 145.2870 154.6818];







% Reordenar parámetros si corresponde
%cc = @(p) [p(5) p(6) p(9) p(8) p(10) p(7) p(3) p(4) p(2) p(1)];
cc = @(p) [p(5) p(6) p(9) p(8) p(10) p(7) p(2) p(1)];
p1m = cc(p1m); 
p4m = cc(p4m); 
p5m = cc(p5m); 
p6m = cc(p6m);

p_all = [p1m; p4m; p5m; p6m];
common_max = p5m;%max(p_all);

%% === Stackear sensibilidades ===
%S_stack = cat(3, S_sens_1, S_sens_2, S_sens_3, S_sens_4);
S_stack = STensor_new;
S_stack = STensor_new./pars_values_filtered_';


% Ahora la prueba que no falla nunca:

%disp('=================================================================');
%%disp(['Valor MÁXIMO de sensibilidad RELATIVA (debe ser >0.1): ' num2str(max(abs(S_stack(:))))]);
%disp(['Valor MÍNIMO no cero: ' num2str(min(abs(S_stack(S_stack~=0))))]);
%disp('=================================================================');

%S_stack = S_stack./sqrt(rbt_fixed_filtered); 
S_stack([12:end], :, :) = [];
pars_to_sens([12:end])= [];

%S_stack(:, :, [500:end]) = [];

[nPar, nVar_tot, nTime] = size(S_stack);

%% === Construir pesos W stackeado ===
sigma = 0.05;%sqrt(0.05);

W1 = diag(1./(sigma .* p1m).^2);
W2 = diag(1./(sigma .* p4m).^2);
W3 = diag(1./(sigma .* p5m).^2);
W4 = diag(1./(sigma .* p6m).^2);

% Matriz en bloque diagonal
%W_different_sigmas = blkdiag(W1, W2, W3, W4);
W_different_sigmas =  diag(1./(sigma .* common_max).^2);

% Caso con un único sigma (misma varianza para todos)
W_single_sigma = eye(nVar_tot) ./ sigma.^2;

%
%W_time_vars = diag(1./(sigma .* squeeze(rbt_fixed_filtered(1, :, 1))'));
%% === Calcular Fisher Information Matrix ===
FIM_d = zeros(nPar, nPar);
FIM_s = zeros(nPar, nPar);
FIM_t = zeros(nPar, nPar);
FIM_ = zeros(nPar, nPar);

for t = 1:nTime
    S_t = S_stack(:,:,t);   % (nVar_tot x nPar) en tiempo t
    
    FIM_ = FIM_ + S_t * S_t';
    FIM_d = FIM_d + S_t * W_different_sigmas * S_t';
    FIM_s = FIM_s + S_t * W_single_sigma * S_t';
    FIM_t = FIM_t + S_t * diag(1./(sigma .* squeeze(rbt_fixed_filtered(1, :, t))')) * S_t';
    %disp('check');
end

%% === Calcular matrices de covarianza y correlación ===
cov_d = inv(FIM_d);
std_d = sqrt(diag(cov_d));
corr_d = cov_d ./ (std_d * std_d');

cov_s = inv(FIM_s);
std_s = sqrt(diag(cov_s));
corr_s = cov_s ./ (std_s * std_s');


cov_t = inv(FIM_t);
std_t = sqrt(diag(cov_t));
corr_t = cov_t ./ (std_t * std_t');

cov_ = inv(FIM_);
std_ = sqrt(diag(cov_));
corr_ = cov_ ./ (std_ * std_');

%%

%[corr_d, pars_to_sens, best_param_idx] = apply_threshold(corr_d, pars_to_sens, 0.2);

pars_to_sens_red = pars_to_sens;
figure;
%imagesc(FIM_d);
%subplot(2,2,4)
imagesc(corr_d);
colorbar; % Add a color scale
%caxis([-1, 1]); % Normalize between -1 and 1
colormap jet; % Use a color map for better visualization

% Set axis labels
set(gca, 'XTick', 1:length(pars_to_sens_red));
set(gca, 'YTick', 1:length(pars_to_sens_red));
if iscell(pars_to_sens_red)
    set(gca, 'XTickLabel', pars_to_sens_red);
    set(gca, 'YTickLabel', pars_to_sens_red);
else
    set(gca, 'XTickLabel', cellstr(pars_to_sens_red));
    set(gca, 'YTickLabel', cellstr(pars_to_sens_red));
end
% Rotate x-axis labels (manual approach for 2017)
hx = get(gca, 'XLabel');
set(hx, 'Units', 'data');
h_labels = get(gca, 'XTickLabel');
delete(get(gca, 'XLabel'));
for i = 1:length(h_labels)
    text(i, 0.5, h_labels{i}, 'HorizontalAlignment', 'right', 'Rotation', 45, ...
         'Units', 'data', 'VerticalAlignment', 'top');
end

%title('FIM');
title('Correlations');
%title("D")

%% === Output ===
disp('FIM con sigmas diferentes (stackeado por sujeto):');
disp(FIM_d);

disp('FIM con un sigma único:');
disp(FIM_s);


% loading_file = load('sens_matrix_testingSigmas.mat');
% sens_matrix = loading_file.sens_reduced;
% 
% 
% sens_matrix = S_stack;
% 
% S = sens_matrix';
% 
% sigma = 0.05;   
% %Valores maximos de normoxia de los 4 sujetos:
% p1m = [1.0870    2.9189    4.2307    9.3023  110.4937   48.3194    1/3.0849  176.6589  136.7549  155.4391];
% p4m = [1.0330    2.8000    1.9661    6.0127  111.0803   44.0947    1/2.7840   203.6085  132.7027  161.3660];
% p5m = [0.9996    3.1379    2.9470    7.0525  102.3333   50.9686    1/2.0473  229.6756  203.3467  212.4784];
% p6m = [0.9346    2.5861    3.6130    9.9826  103.8647   51.6685    1/2.6860  183.5269  145.2870  154.6818];
% cc = @(p) [p(5) p(6) p(9) p(8) p(10) p(7) p(3) p(4) p(2) p(1)];
% [cc(p1m) cc(p4m) cc(p5m) cc(p6m)];
% %pmax_averaged = (p1m + p4m + p5m + p6m)/4; %promedio de los valores máximos (para testear nomás)
% p = pmax_averaged; %p1m; acá seleccionar alguno de los 4, o el promedio
% 
% %"dVE", "VT", "TI", "Tresp", "PAO2", "PACO2", "HR", "PS", "PD", "PM" orden
% %de las variables en los datos experimentales
% 
% %{'PAO2', 'PACO2', 'pd', 'ps', 'pm', 'Theart', 'TI', 'BF', 'VTidal',
% %'dVE'}; order de las variables en la matriz de sensibilidad
% 
% %p_max = [p(5) p(6) p(9) p(8) p(10) p(7) p(3) p(4) p(2) p(1)]; %corrección para ordenar los datos
% 
% 
% I = eye(size(sens_matrix, 2));  
% W_different_sigmas = diag((1./(sigma * p_max).^2).^1);  %acá no sé si p_max también va al cuadrado o no
% W_single_sigma = I./sigma.^2;
% FIM_d = S' * W_different_sigmas * S;  
% FIM_s = S' * W_single_sigma * S;
% 
% %Corr matrix para 'sigmas' distintos
% covariance_matrix_d = inv(FIM_d);
% std_devs_d = sqrt(diag(covariance_matrix_d));
% corr_matrix_d = covariance_matrix_d ./ (std_devs_d * std_devs_d');
% 
% %Corr matrix para mismo sigma de 0.05
% covariance_matrix_s = inv(FIM_s);
% std_devs_s = sqrt(diag(covariance_matrix_s));
% corr_matrix_s = covariance_matrix_s ./ (std_devs_s * std_devs_s');






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