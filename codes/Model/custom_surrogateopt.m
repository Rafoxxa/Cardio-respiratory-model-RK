function [x_best, fval_best, exitflag, output, log_struct] = custom_surrogateopt(objective_fun, lb, ub, options)
% CUSTOM_SURROGATEOPT - Implementación propia de Surrogate Optimization
% Basada en algoritmos de Radial Basis Function (RBF) y Expected Improvement
%
% Algoritmo basado en:
% - Jones et al. (1998) "Efficient Global Optimization of Expensive Black-Box Functions"
% - Gutmann (2001) "A Radial Basis Function Method for Global Optimization"
% - Regis & Shoemaker (2007) "A stochastic radial basis function method"
%
% Sintaxis:
%   [x, fval, exitflag, output] = custom_surrogateopt(fun, lb, ub)
%   [x, fval, exitflag, output] = custom_surrogateopt(fun, lb, ub, options)
%
% Entradas:
%   objective_fun - Función objetivo a minimizar
%   lb            - Límites inferiores (vector)
%   ub            - Límites superiores (vector)
%   options       - Estructura de opciones (opcional)
%
% Salidas:
%   x_best        - Mejor punto encontrado
%   fval_best     - Valor de la función en x_best
%   exitflag      - Código de salida
%   output        - Información adicional
%   log_struct    - Estructura de logging para compatibilidad

    %% Validación de entradas
    if nargin < 3
        error('Se requieren al menos 3 argumentos: fun, lb, ub');
    end
    
    if nargin < 4 || isempty(options)
        options = struct();
    end
    
    % Dimensión del problema
    nvars = length(lb);
    
    if length(ub) ~= nvars
        error('lb y ub deben tener la misma longitud');
    end
    
    %% Configuración por defecto
    default_options = struct(...
        'MaxFunctionEvaluations', 200, ...
        'MinSampleDistance', 1e-6, ...
        'MinSurrogatePoints', max(5, nvars+1), ...
        'Display', 'iter', ...
        'FunctionTolerance', 1e-6, ...
        'MaxStallIterations', 20, ...
        'OutputFcn', [], ...
        'InitialPoints', [] ...
    );
    
    % Combinar opciones del usuario con las por defecto
    field_names = fieldnames(default_options);
    for i = 1:length(field_names)
        if ~isfield(options, field_names{i})
            options.(field_names{i}) = default_options.(field_names{i});
        end
    end
    
    %% Inicialización
    X_sampled = [];           % Puntos evaluados [n_samples x nvars]
    Y_sampled = [];           % Valores de función [n_samples x 1]
    n_evals = 0;              % Contador de evaluaciones
    stall_counter = 0;        % Contador de iteraciones sin mejora
    
    % Logging
    log_struct = struct();
    log_struct.iteration = [];
    log_struct.bestfval = [];
    log_struct.bestx = [];
    log_struct.fval_history = [];
    
    % Normalización de variables (importante para RBF)
    X_range = ub - lb;
    
    %% Fase 1: Muestreo inicial (Latin Hypercube Sampling)
    if isempty(options.InitialPoints)
        n_initial = options.MinSurrogatePoints;
        X_initial = lhsdesign(n_initial, nvars); % Latin Hypercube en [0,1]
        X_initial = repmat(lb, n_initial, 1) + X_initial .* repmat(X_range, n_initial, 1);
    else
        X_initial = options.InitialPoints;
        n_initial = size(X_initial, 1);
    end
    
    % Evaluar puntos iniciales
    Y_initial = zeros(n_initial, 1);
    for i = 1:n_initial
        Y_initial(i) = objective_fun(X_initial(i, :));
        n_evals = n_evals + 1;
        
        if strcmp(options.Display, 'iter')
            fprintf('Initial point %d/%d: f = %.6f\n', i, n_initial, Y_initial(i));
        end
    end
    
    X_sampled = X_initial;
    Y_sampled = Y_initial;
    
    % Mejor punto inicial
    [fval_best, best_idx] = min(Y_sampled);
    x_best = X_sampled(best_idx, :);
    
    %% Bucle principal de optimización
    iteration = 0;
    exitflag = 0;
    
    while n_evals < options.MaxFunctionEvaluations
        iteration = iteration + 1;
        
        % Logging
        log_struct.iteration(end+1) = iteration;
        log_struct.bestfval(end+1) = fval_best;
        log_struct.bestx(:, end+1) = x_best';
        log_struct.fval_history(end+1) = fval_best;
        
        % Llamar OutputFcn si existe
        if ~isempty(options.OutputFcn)
            optimValues = struct('iteration', iteration, 'fval', fval_best, ...
                'x', x_best, 'funccount', n_evals);
            try
                [stop, ~, ~] = options.OutputFcn(optimValues, options, 'iter');
                if stop
                    exitflag = -1;
                    break;
                end
            catch
                % Si hay error en OutputFcn, continuar
            end
        end
        
        %% Construir modelo surrogate (RBF)
        [surrogate_model, model_params] = buildRBFModel(X_sampled, Y_sampled, lb, ub);
        
        %% Optimizar función de adquisición (Expected Improvement)
        x_next = optimizeAcquisition(surrogate_model, model_params, X_sampled, Y_sampled, lb, ub);
        
        %% Verificar distancia mínima
        if size(X_sampled, 1) > 1
            distances = sqrt(sum((X_sampled - repmat(x_next, size(X_sampled, 1), 1)).^2, 2));
            if min(distances) < options.MinSampleDistance
                % Punto muy cercano, generar uno aleatorio
                x_next = lb + rand(1, nvars) .* X_range;
            end
        end
        
        %% Evaluar nuevo punto
        f_next = objective_fun(x_next);
        n_evals = n_evals + 1;
        
        % Actualizar conjunto de datos
        X_sampled = [X_sampled; x_next];
        Y_sampled = [Y_sampled; f_next];
        
        % Actualizar mejor solución
        if f_next < fval_best
            fval_best = f_next;
            x_best = x_next;
            stall_counter = 0;
        else
            stall_counter = stall_counter + 1;
        end
        
        % Display progress
        if strcmp(options.Display, 'iter')
            fprintf('Iter %3d | FuncEvals %3d | Best fval = %.6f | Current = %.6f\n', ...
                iteration, n_evals, fval_best, f_next);
        end
        
        %% Criterios de parada
        if stall_counter >= options.MaxStallIterations
            exitflag = 1;
            if strcmp(options.Display, 'iter')
                fprintf('Stopped: No improvement for %d iterations\n', options.MaxStallIterations);
            end
            break;
        end
        
        if abs(f_next - fval_best) < options.FunctionTolerance && fval_best < options.FunctionTolerance
            exitflag = 2;
            if strcmp(options.Display, 'iter')
                fprintf('Stopped: Function tolerance reached\n');
            end
            break;
        end
    end
    
    % Final OutputFcn call
    if ~isempty(options.OutputFcn)
        optimValues = struct('iteration', iteration, 'fval', fval_best, ...
            'x', x_best, 'funccount', n_evals);
        try
            options.OutputFcn(optimValues, options, 'done');
        catch
            % Continue if OutputFcn fails
        end
    end
    
    %% Preparar outputs
    output = struct();
    output.iterations = iteration;
    output.funcCount = n_evals;
    output.algorithm = 'Custom Surrogate Optimization (RBF + Expected Improvement)';
    output.message = sprintf('Optimization completed after %d function evaluations', n_evals);
    
    if exitflag == 0
        exitflag = 3; % Maximum function evaluations reached
    end
    
end

%% Función para construir modelo RBF
function [surrogate_model, model_params] = buildRBFModel(X, Y, lb, ub)
    % Normalizar datos a [0,1]
    X_norm = (X - repmat(lb, size(X, 1), 1)) ./ repmat(ub - lb, size(X, 1), 1);
    
    % Parámetros del modelo
    model_params = struct();
    model_params.lb = lb;
    model_params.ub = ub;
    model_params.X_train = X_norm;
    model_params.Y_train = Y;
    model_params.Y_mean = mean(Y);
    model_params.Y_std = std(Y);
    
    % Normalizar Y
    Y_norm = (Y - model_params.Y_mean) / max(model_params.Y_std, 1e-10);
    
    % Seleccionar parámetro de forma para RBF (basado en literatura)
    n_points = size(X_norm, 1);
    nvars = size(X_norm, 2);
    
    % Parámetro adaptativo basado en densidad de puntos
    avg_distance = sqrt(nvars) / (n_points^(1/nvars)); % Distancia promedio en espacio unitario
    model_params.sigma = max(0.1, min(2.0, avg_distance));
    
    % Función surrogate usando RBF Gaussiana
    surrogate_model = @(x_query) predictRBF(x_query, model_params, Y_norm);
end

%% Función de predicción RBF
function [y_pred, std_pred] = predictRBF(X_query, params, Y_norm)
    % Normalizar puntos de consulta
    X_query_norm = (X_query - repmat(params.lb, size(X_query, 1), 1)) ./ ...
                   repmat(params.ub - params.lb, size(X_query, 1), 1);
    
    n_query = size(X_query_norm, 1);
    n_train = size(params.X_train, 1);
    
    % Matriz de distancias
    y_pred = zeros(n_query, 1);
    std_pred = zeros(n_query, 1);
    
    for i = 1:n_query
        % Calcular distancias a todos los puntos de entrenamiento
        distances = sqrt(sum((params.X_train - repmat(X_query_norm(i, :), n_train, 1)).^2, 2));
        
        % RBF Gaussiana
        weights = exp(-(distances / params.sigma).^2);
        
        % Evitar división por cero
        if sum(weights) < 1e-10
            weights = ones(size(weights)) / n_train;
        else
            weights = weights / sum(weights);
        end
        
        % Predicción
        y_pred(i) = sum(weights .* Y_norm);
        
        % Estimación de incertidumbre (basada en distancia al punto más cercano)
        min_dist = min(distances);
        std_pred(i) = min_dist * std(Y_norm);
    end
    
    % Desnormalizar
    y_pred = y_pred * params.Y_std + params.Y_mean;
    std_pred = std_pred * params.Y_std;
end

%% Optimización de función de adquisición (Expected Improvement)
function x_next = optimizeAcquisition(surrogate_model, model_params, X_sampled, Y_sampled, lb, ub)
    nvars = length(lb);
    f_best = min(Y_sampled);
    
    % Función Expected Improvement
    acquisition_fun = @(x) -computeExpectedImprovement(x, surrogate_model, f_best);
    
    % Múltiples puntos de inicio para optimizar adquisición
    n_starts = min(10, 3*nvars);
    best_acq_val = Inf;
    x_next = [];
    
    for i = 1:n_starts
        % Punto inicial aleatorio
        x0 = lb + rand(1, nvars) .* (ub - lb);
        
        try
            % Usar fmincon para optimizar adquisición (disponible en 2017)
            options_fmincon = optimoptions('fmincon', ...
                'Algorithm', 'sqp', ...
                'Display', 'off', ...
                'MaxIterations', 50, ...
                'FiniteDifferenceStepSize', 0.01);
            
            [x_candidate, acq_val] = fmincon(acquisition_fun, x0, [], [], [], [], ...
                lb, ub, [], options_fmincon);
            
            % Mantener mejor candidato
            if acq_val < best_acq_val
                best_acq_val = acq_val;
                x_next = x_candidate;
            end
        catch
            % Si fmincon falla, usar el punto aleatorio
            if isempty(x_next)
                x_next = x0;
            end
        end
    end
    
    % Fallback si no se encontró nada
    if isempty(x_next)
        x_next = lb + rand(1, nvars) .* (ub - lb);
    end
end

%% Función Expected Improvement
function EI = computeExpectedImprovement(x, surrogate_model, f_best)
    % Predecir con modelo surrogate
    [mu, sigma] = surrogate_model(x);
    
    % Evitar división por cero
    sigma = max(sigma, 1e-10);
    
    % Expected Improvement formula
    if size(x, 1) == 1
        improvement = f_best - mu;
        z = improvement / sigma;
        
        % Usar approximaciones numéricas para evitar problemas con erf
        if z > 0
            EI = improvement * normcdf(z) + sigma * normpdf(z);
        else
            EI = sigma * normpdf(z);
        end
        
        % Penalizar predicciones muy inciertas
        if sigma > abs(mu)
            EI = EI * 0.5;
        end
    else
        % Manejar múltiples puntos
        EI = zeros(size(x, 1), 1);
        for i = 1:size(x, 1)
            EI(i) = computeExpectedImprovement(x(i, :), surrogate_model, f_best);
        end
    end
end

%% Función Latin Hypercube Sampling (para compatibilidad con 2017)
function X = lhsdesign(n, p)
    % Implementación simple de Latin Hypercube Sampling
    X = zeros(n, p);
    
    for i = 1:p
        % Crear permutación aleatoria
        perm = randperm(n);
        % Generar valores uniformes en cada intervalo
        X(:, i) = (perm' - rand(n, 1)) / n;
    end
end

%% Función de distribución normal estándar
function p = normcdf(x)
    % Aproximación de la CDF normal estándar
    p = 0.5 * (1 + erf(x / sqrt(2)));
end

function p = normpdf(x)
    % PDF normal estándar
    p = (1 / sqrt(2 * pi)) * exp(-0.5 * x.^2);
end