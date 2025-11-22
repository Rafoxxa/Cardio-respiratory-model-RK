%% =========================================================================
%%  CÁLCULO DE LA MATRIZ DE INFORMACIÓN DE FISHER (FIM) Y CORRELACIONES
%%  Caso: Un único sujeto → STensor único (ya filtrado previamente)
%% =========================================================================

clear; clc; close all;

%% -------------------------------------------------------------------------
%% 1. CARGA DEL TENSOR DE SENSIBILIDAD Y PARÁMETROS (ya filtrados)
%% -------------------------------------------------------------------------
loaded = load('STensor_new_red.mat');  % Ajusta la ruta/nombre
STensor_fm          = loaded.STensor_new;      % [nPar × nVar × nTime]
%pars_to_sens_fm     = loaded.pars_to_sens_reduced;  % cell con nombres si es que se cara STensor_new
pars_to_sens_fm     = loaded.pars_to_sens_filtered;  % cell con nombres  si es que se carga STensor_new_red
y_nominal_loaded = load('y_nomnial.mat');
y_nominal = y_nominal_loaded.rbt_fixed;
pars_cont_loaded = load('pars_containers.mat');
pars_cont = pars_cont_loaded.pars_cont;

%% En caso de que las variables ya estén cargadas
pars_cont = setup.pars;
STensor_fm = STensor_red;
pars_to_sens_fm = pars_red;
y_nominal = rbt_fixed;

%% Submuestrear a la frecuencia de muestreo experimental (0.33 Hz o periodo de muestreo de 3s)
STensor_fm = STensor_fm(:, :, 1:3/0.1:end);
y_nominal = y_nominal(:, :, 1:3/0.1:end);

STensor_fm(:, [5,6], :) = [];
y_nominal(:, [5,6], :) = [];

%{'PAO2','PACO2','pd','ps','pm','Theart','VTidal','dVE'};

[nPar, nVar, nTime] = size(STensor_fm);
fprintf('Tensor cargado: %d parámetros × %d variables × %d puntos temporales\n', nPar, nVar, nTime);

%% -------------------------------------------------------------------------
%% 2. OBTENER VALORES MAXIMOS DEL SUJETO (para sensibilidad relativa)
%% -------------------------------------------------------------------------
% Valores maximos del sujeto de referencia (ej: sujeto 5)
v_max_p5 = [0.9996, 3.1379, 2.9470, 7.0525, 102.3333, 50.9686, 0.8, ...
             229.6756, 203.3467, 212.4784];

% Reordenar según el orden real de las variables en tu modelo/simulación
%reorder = @(p) [p(5) p(6) p(9) p(8) p(10) p(7) p(2) p(1)];  % Ajustado a tu mapeo
%reorder = @(p) [p(5) p(6) p(10) p(7) p(2) p(1)];  % Ajustado a tu mapeo
reorder = @(p) [p(5) p(6) p(9) p(8) p(2) p(1)];  % Ajustado a tu mapeo
v_max = reorder(v_max_p5);

p_nominal = cellfun(@(name) pars_cont(name), pars_to_sens_fm, 'UniformOutput', true);

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
FIM_rel   = zeros(nPar, nPar);
corr_time = zeros(nPar, nPar);
max_over_time = zeros(1,nTime);
median_over_time = zeros(1,nTime);

for t = 1:nTime
    S_t = squeeze(S_abs(:, :, t))';        % → [nVar × nPar]
    
    FIM_single = FIM_single + S_t' * W_single * S_t;
    FIM_prop   = FIM_prop   + S_t' * W_prop   * S_t;
    %FIM_rel = FIM_rel + S_t' * diag(1./(sigma .* squeeze(y_nominal(1,:, t)'))) * S_t;
    %max_over_time(t) = max(log(abs(FIM_prop(:)))/(log(10)));
    %median_over_time(t) = median(log(abs(FIM_prop(:)))/(log(10)));
   


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

%% 6. VISUALIZACIÓN DE LA MATRIZ DE CORRELACIÓN
%% -------------------------------------------------------------------------
figure('Position', [100, 100, 900, 800]);
imagesc(corr_single);
colorbar;
%caxis([-1 1]);
colormap(jet);
axis square;

% Etiquetas de parámetros
set(gca, 'XTick', 1:nPar, 'XTickLabel', pars_to_sens_fm, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nPar, 'YTickLabel', pars_to_sens_fm);

title('Matriz de Correlación de Parámetros (FIM single sigma)', 'FontSize', 12);
xlabel('Parámetros');
ylabel('Parámetros');

% Mejorar legibilidad de etiquetas largas
set(gca, 'TickLabelInterpreter', 'none');


figure('Position', [100, 100, 900, 800]);
imagesc(corr_prop);
colorbar;
%caxis([-1 1]);
colormap(jet);
axis square;

% Etiquetas de parámetros
set(gca, 'XTick', 1:nPar, 'XTickLabel', pars_to_sens_fm, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nPar, 'YTickLabel', pars_to_sens_fm);

title('Matriz de Correlación de Parámetros (FIM prop)', 'FontSize', 12);
xlabel('Parámetros');
ylabel('Parámetros');

% Mejorar legibilidad de etiquetas largas
set(gca, 'TickLabelInterpreter', 'none');


figure('Position', [100, 100, 900, 800]);
imagesc(corr_rel);
colorbar;
%caxis([-1 1]);
colormap(jet);
axis square;

% Etiquetas de parámetros
set(gca, 'XTick', 1:nPar, 'XTickLabel', pars_to_sens_fm, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:nPar, 'YTickLabel', pars_to_sens_fm);

title('Matriz de Correlación de Parámetros (FIM rel)', 'FontSize', 12);
xlabel('Parámetros');
ylabel('Parámetros');

% Mejorar legibilidad de etiquetas largas
set(gca, 'TickLabelInterpreter', 'none');