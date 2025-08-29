rng(2);  
%clc
clear -global delays_global
clear -global all_global
clear -global externals_global

vectorize_dicts('run_ode.m', 'model_basic.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m');

patient_idx = 1;
[setup] = set_up("simulation", patient_idx, "normoxia", "-", "dt", 0.1);%, "pars_from_fitting", 1, "fitting_mat_file", "Fitting_test.mat");
s = setup;

% Ejecutar simulación ODE
[t, x_dot, x_vars, x_keys, index] = s.run_ode_fun(s.model, s.pars, s.init, s.simulation_time, s.dt);

% Organizar resultados
struct_vars = arrange_results(x_dot, x_vars, x_keys, t);

% Guardar resultados
save(s.simulation_filename, 'struct_vars', 't');

% -----------------------------
% LINEARIZACIÓN LOCAL DESDE DATOS (robusta)
% -----------------------------
var_names = fieldnames(struct_vars);
n = numel(var_names);
N = length(t);

% Construir matriz de estados
X = zeros(n, N);
for i = 1:n
    X(i,:) = struct_vars.(var_names{i});
end

% Entradas: si no tienes, deja m=0 y U=[]
m = 0; 
U = [];

dt = s.dt;
win = 200;             % tamaño de ventana
ridge = 1e-3;          % regularización (sube a 1e-2 o 1e-1 si sigue mal condicionada)
var_thresh = 1e-8;     % umbral de varianza para descartar ventanas
rho_max = 1.05;        % límite opcional al radio espectral (discreto)
use_shrink = true;     % activa proyección de estabilidad
use_affine = true;     % incluye término afín

Ad = cell(1, N-1); 
Bd = cell(1, N-1); 
d_term = cell(1, N-1);

for k = 1:(N-1)
    % Ventana [k ... k+win-1]
    k1 = k; 
    k2 = min(k+win-1, N-1);
    idx = k1:k2; 
    T = numel(idx);
    
    Xk  = X(:, idx);
    Xkp = X(:, idx+1);
    if m > 0, Uk = U(:, idx); else, Uk = []; end
    
    % --- 1) Centrado + escalado (z-score por ventana) ---
    muX = mean(Xk, 2);
    sigX = std(Xk, 0, 2);  % std cols
    sigX(sigX < 1e-12) = 1e-12;  % evita división por cero
    Xk_s  = bsxfun(@rdivide, bsxfun(@minus, Xk,  muX), sigX);
    Xkp_s = bsxfun(@rdivide, bsxfun(@minus, Xkp, muX), sigX);
    
    if m > 0
        muU = mean(Uk, 2);
        sigU = std(Uk, 0, 2);
        sigU(sigU < 1e-12) = 1e-12;
        Uk_s = bsxfun(@rdivide, bsxfun(@minus, Uk, muU), sigU);
    else
        Uk_s = [];
    end
    
    % Descarta ventanas sin variación (mal condicionadas por construcción)
    if mean(var(Xk_s, 0, 2)) < var_thresh
        % Fallback: identidad y offset nulo
        Ad{k} = eye(n);
        Bd{k} = [];
        d_term{k} = zeros(n,1);
        continue;
    end
    
    % --- 2) Regressor estandarizado ---
    if m > 0
        if use_affine
            Phi = [Xk_s; Uk_s; ones(1,T)];
        else
            Phi = [Xk_s; Uk_s];
        end
    else
        if use_affine
            Phi = [Xk_s; ones(1,T)];
        else
            Phi = Xk_s;
        end
    end
    
    % --- 2b) Solución ridge vía SVD truncada (robusta) ---
    % Theta_s = Xkp_s * Phi' * inv(Phi*Phi' + ridge*I)
    [U_svd, S_svd, V_svd] = svd(Phi, 'econ');  % Phi = U S V'
    svals = diag(S_svd);
    tol = max(size(Phi)) * eps(max(svals));
    r = sum(svals > tol);                      % rango efectivo
    
    % Construye filtro ridge en singulares: s / (s^2 + ridge)
    filt = svals(1:r) ./ (svals(1:r).^2 + ridge);
    Theta_s = Xkp_s * V_svd(:,1:r) * diag(filt) * U_svd(:,1:r)';  % (n x rows(Phi))
    
    % --- 3) Des-estandarizar para recuperar A_d, B_d, d ---
    % Xkp_s ? A_s * Xk_s + B_s * Uk_s + d_s
    % En variables originales: Xkp ? A*Xk + B*Uk + d
    % Con z-score: A = D_x * A_s * D_x^{-1}, etc., y d ajusta con medias.
    Dx  = diag(sigX);
    iDx = diag(1./sigX);
    if m > 0
        Du  = diag(sigU);
        iDu = diag(1./sigU);
    end
    
    % Particiona Theta_s
    if m > 0
        if use_affine
            A_s = Theta_s(:, 1:n);
            B_s = Theta_s(:, n+(1:m));
            d_s = Theta_s(:, end);
        else
            A_s = Theta_s(:, 1:n);
            B_s = Theta_s(:, n+(1:m));
            d_s = zeros(n,1);
        end
    else
        if use_affine
            A_s = Theta_s(:, 1:n);
            B_s = [];
            d_s = Theta_s(:, end);
        else
            A_s = Theta_s(:, 1:n);
            B_s = [];
            d_s = zeros(n,1);
        end
    end
    
    % Mapea a escala original
    Ahat = Dx * A_s * iDx;
    if m > 0
        Bhat = Dx * B_s * iDu;
    else
        Bhat = [];
    end
    dhat = Dx * (d_s - A_s * ((-muX)./sigX) - ( (m>0)*B_s * ((- (m>0)*muU)./max(sigU,1e-12)) ));
    % La expresión de dhat equivale a corregir por las medias; más claro:
    % Xkp = Dx*A_s*iDx*Xk + Dx*B_s*iDu*Uk + Dx*(d_s - A_s*(muX./sigX) - B_s*(muU./sigU)) + muX
    
    % Alternativa explícita (más estable numéricamente):
    dhat = Dx*(d_s - A_s.*(muX./sigX)*0) + muX - Ahat*muX;
    if m > 0
        dhat = dhat - Bhat*muU;
    end
    
    % --- 4) Proyección opcional para evitar explosiones numéricas ---
    if use_shrink
        [V_eig, D_eig] = eig(Ahat);
        lam = diag(D_eig);
        rho = max(abs(lam));
        if isfinite(rho) && rho > rho_max
            Ahat = Ahat * (rho_max / rho);
        end
    end
    
    Ad{k}    = Ahat;
    Bd{k}    = Bhat;
    d_term{k}= dhat;
end

% -----------------------------
% SIMULACIÓN LINEAL CON A,B,d (LTV)
% -----------------------------
x_lin = zeros(n, N);
x_lin(:,1) = X(:,1);

for k = 1:(N-1)
    if isempty(Bd{k})
        x_lin(:,k+1) = Ad{k}*x_lin(:,k) + d_term{k};
    else
        x_lin(:,k+1) = Ad{k}*x_lin(:,k) + Bd{k}*U(:,k) + d_term{k};
    end
end

% -----------------------------
% GRAFICAR COMPARACIÓN
% -----------------------------
figure;
var_idx = 1; % elige variable a comparar
plot(t, X(var_idx,:), 'b', 'LineWidth',1.5); hold on;
plot(t, x_lin(var_idx,:), 'r--', 'LineWidth',1.5);
xlabel('Tiempo [s]');
ylabel(var_names{var_idx});
legend('No lineal','Lineal (robusto)');
title(['Comparación ', var_names{var_idx}]);