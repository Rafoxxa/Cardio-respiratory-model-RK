%% =========================================================================
%%  CÁLCULO DE LA MATRIZ DE INFORMACIÓN DE FISHER (FIM) Y CORRELACIONES
%%  Caso: Un único sujeto → STensor único (ya filtrado previamente)
%% =========================================================================

clear; clc; close all;

%% -------------------------------------------------------------------------
%% 1. CARGA DEL TENSOR DE SENSIBILIDAD Y PARÁMETROS (ya filtrados)
%% -------------------------------------------------------------------------
loaded = load('../Sens_analysis/STensor_new_16-11-2025.mat');  % Ajusta la ruta/nombre
STensor_fm          = loaded.STensor_new.tensor;      % [nPar × nVar × nTime]
pars_to_sens_fm     = loaded.STensor_new.parameters;  % cell con nombres
y_nominal_loaded = load('y_nominal.mat');
y_nominal = y_nominal_loaded.rbt_fixed;

% Si ya tienes una versión reducida (recomendado), úsala directamente:
%pars_to_sens = pars_to_sens_filtered;      % ← Aquí usas tu lista ya filtrada
%STensor      = STensor_new;                % ← Tu tensor ya filtrado (mismo orden!)
[STensor_new, ~, ~] = reduce_STensor_by_names(STensor, pars_to_sens, pars_to_sens_filtered);
STensor_fm(:, [7,8],:) = []; %remover las variables de tiempo;
pars_to_sens_fm = pars_to_sens;%pars_to_sens_filtered;

[nPar, nVar, nTime] = size(STensor_fm);
fprintf('Tensor cargado: %d parámetros × %d variables × %d puntos temporales\n', nPar, nVar, nTime);

%% -------------------------------------------------------------------------
%% 2. OBTENER VALORES MAXIMOS DEL SUJETO (para sensibilidad relativa)
%% -------------------------------------------------------------------------
% Valores maximos del sujeto de referencia (ej: sujeto 5)
v_max_p5 = [0.9996, 3.1379, 2.9470, 7.0525, 102.3333, 50.9686, 0.8, ...
             229.6756, 203.3467, 212.4784];

% Reordenar según el orden real de las variables en tu modelo/simulación
reorder = @(p) [p(5) p(6) p(9) p(8) p(10) p(7) p(2) p(1)];  % Ajustado a tu mapeo
v_max = reorder(v_max_p5);

p_nominal = cellfun(@(name) setup.pars(name), pars_to_sens_fm, 'UniformOutput', true);

S_abs = STensor_fm ./ p_nominal';  % Broadcasting: divide cada parámetro por su valor nominal



%% -------------------------------------------------------------------------
%% 3. DEFINIR MATRICES DE PONDERACIÓN (W) - Ruido de medición
%% -------------------------------------------------------------------------
sigma = 0.05;  % Desviación estándar relativa del ruido (5%)

% Caso 1: Mismo sigma para todas las variables
W_single = eye(nVar) / sigma^2;

% Caso 2: Sigma proporcional al valor nominal de cada variable (más realista)
W_prop = diag(1 ./ (sigma * v_max).^2);

% Caso 3: Sigma constante en valor absoluto (menos realista, pero común)
% W_const = eye(nVar) / sigma^2;

%% -------------------------------------------------------------------------
%% 4. CÁLCULO DE LA MATRIZ DE INFORMACIÓN DE FISHER (FIM)
%% -------------------------------------------------------------------------


%S_abs_ = S_abs;
%S_abs_([3, 16], :, :) = [];
%S_abs_([8:9, 15:20, 29:31, 34:38], :, :) = [];
[nPar, nVar, nTime] = size(S_abs);

FIM_single = zeros(nPar, nPar);
FIM_prop   = zeros(nPar, nPar);
max_over_time = zeros(1,nTime);
median_over_time = zeros(1,nTime);

for t = 1:nTime
    S_t = squeeze(S_abs_(:, :, t))';        % → [nVar × nPar]
    
    FIM_single = FIM_single + S_t' * W_single * S_t;
    FIM_prop   = FIM_prop   + S_t' * W_prop   * S_t;
    FIM_rel = FIM_rel + St' * diag(1./(sigma .* y_nominal(:, t)')) * St;
    max_over_time(t) = max(log(abs(FIM_prop(:)))/(log(10)));
    median_over_time(t) = median(log(abs(FIM_prop(:)))/(log(10)));

end

fprintf('FIM calculada correctamente (acumulación temporal completa)\n');

%% -------------------------------------------------------------------------
%% 5. MATRIZ DE COVARIANZA Y CORRELACIÓN
%% -------------------------------------------------------------------------
cov_single = inv(FIM_single);
cov_prop   = inv(FIM_prop);
cov_rel   = inv(FIM_rel);

std_single = sqrt(diag(cov_single));
std_prop   = sqrt(diag(cov_prop));
std_rel   = sqrt(diag(cov_rel));

corr_single = cov_single ./ (std_single * std_single');
corr_prop   = cov_prop   ./ (std_prop   * std_prop');
corr_rel   = cov_rel  ./ (std_rel   * std_rel');


%% 
% =================================================================
%  CÁLCULO ROBUSTO DE CORRELACIONES SIN INVERTIR LA FIM (SVD)
% =================================================================
[nPar, nVar, nTime] = size(STensor_fm);
S_stack = zeros(nVar*nTime, nPar);   % Matriz sensibilidad apilada
y_nominal = squeeze(rbt_fixed(1,:, :));

k = 1;
for t = 1:nTime
    St = squeeze(STensor_fm(:, :, t))';           % [nVar × nPar]
    %St = St ./ (sigma * y_nominal(:,t));           % W = 1/(σ·y(t))² → sqrt(W)
    St = St ./ (sigma);           % W = 1/(σ·y(t))² → sqrt(W)
    S_stack(k:k+nVar-1, :) = St;
    k = k + nVar;
end

% SVD (esto es numéricamente perfecto aunque cond(FIM) = 10²⁰)
[U, SingularValues, V] = svd(S_stack, 'econ');

% Correlaciones prácticos = correlaciones de los vectores V (right singular vectors)
corr_matrix_practical = corr(V);           % ¡Esta sí es la verdadera!
% o más estricto:
corr_matrix_practical = abs(corr(V));      % porque signo es arbitrario

% Ranking de identificabilidad (1/σ_i² proporcional a info)
info_per_param = sum(V.^2, 2);             % o diag(V'*V)
[~, ranking] = sort(info_per_param, 'descend');

%% -------------------------------------------------------------------------
%% 6. VISUALIZACIÓN DE LA MATRIZ DE CORRELACIÓN (la más usada)
%% -------------------------------------------------------------------------
figure('Position', [100, 100, 900, 800]);
imagesc(corr_prop_);
colorbar;
%caxis([-1 1]);
colormap(jet);
axis square;

% Etiquetas de parámetros
set(gca, 'XTick', 1:nPar, 'XTickLabel', pars_to_sens, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nPar, 'YTickLabel', pars_to_sens);

title('Matriz de Correlación de Parámetros (FIM con \sigma proporcional)', 'FontSize', 12);
xlabel('Parámetros');
ylabel('Parámetros');

% Mejorar legibilidad de etiquetas largas
set(gca, 'TickLabelInterpreter', 'none');

%% -------------------------------------------------------------------------
%% 7. (Opcional) ELIMINAR PARÁMETROS NO IDENTIFICABLES ITERATIVAMENTE
%% -------------------------------------------------------------------------
% Descomenta si quieres limpiar automáticamente parámetros muy correlacionados
% threshold = 0.95;
% [corr_clean, pars_clean, removed_idx] = remove_highly_correlated(corr_prop, pars_to_sens, threshold);
% 
% figure;
% imagesc(corr_clean); colorbar; caxis([-1 1]); colormap(jet);
% set(gca, 'XTick', 1:length(pars_clean), 'XTickLabel', pars_clean, 'XTickLabelRotation', 45);
% set(gca, 'YTick', 1:length(pars_clean), 'YTickLabel', pars_clean);
% title(sprintf('Correlaciones tras eliminar > %.2f', threshold));

%% -------------------------------------------------------------------------
%% 8. RESUMEN FINAL
%% -------------------------------------------------------------------------
fprintf('\n=== RESUMEN IDENTIFICABILIDAD ===\n');
fprintf('Parámetros totales: %d\n', nPar);
fprintf('Determinante de FIM (proporcional): %.3e\n', det(FIM_prop));
fprintf('Condición de FIM: %.3e\n', cond(FIM_prop));
fprintf('Parámetro con mayor varianza (peor estimado): %s (var = %.4f)\n', ...
    pars_to_sens{find(std_prop == max(std_prop),1)}, max(std_prop));

disp('¡Análisis de FIM completado con éxito!');

%% =========================================================================
%% FUNCIÓN AUXILIAR: Eliminar parámetros altamente correlacionados
%% =========================================================================
function [C_clean, par_clean, removed] = remove_highly_correlated(C, par_list, thresh)
    C_clean = C;
    par_clean = par_list;
    removed = [];
    
    while true
        over = abs(C_clean) > thresh;
        over = over & ~eye(size(C_clean));  % Excluir diagonal
        if ~any(over(:)), break; end
        
        [i, j] = find(triu(over,1), 1);
        % Eliminar el de mayor correlación total
        corr_sum = sum(abs(C_clean), 2);
        if corr_sum(i) >= corr_sum(j)
            idx_remove = i;
        else
            idx_remove = j;
        end
        
        removed = [removed; par_clean(idx_remove)];
        fprintf('Eliminado por alta correlación: %s\n', par_clean{idx_remove});
        
        C_clean(idx_remove, :) = [];
        C_clean(:, idx_remove) = [];
        par_clean(idx_remove) = [];
    end
end

%% =========================================================================
%%  FUNCIÓN: reduce_STensor_by_names
%%  Reduce un STensor usando listas explícitas de nombres
