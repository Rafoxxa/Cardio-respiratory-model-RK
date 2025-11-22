function [t_filled, v_filled] = fill_dis(t, v, fs)
    % Función para rellenar discontinuidades en vectores de tiempo y datos
    % 
    % Entradas:
    %   t  - Vector de tiempos (puede tener discontinuidades)
    %   v  - Vector de valores asociados a cada tiempo
    %   fs - Frecuencia de muestreo (Hz)
    %
    % Salidas:
    %   t_filled - Vector de tiempos con discontinuidades rellenadas
    %   v_filled - Vector de valores interpolados
    
    % Validar entradas
    if length(t) ~= length(v)
        error('Los vectores t y v deben tener la misma longitud');
    end
    
    % Asegurar que sean vectores columna
    t = t(:);
    v = v(:);
    
    % Calcular el período de muestreo esperado
    dt = 1/fs;
    
    % Detectar discontinuidades
    % Una discontinuidad existe cuando la diferencia entre tiempos consecutivos
    % es mayor que 1.5 veces el período de muestreo (tolerancia para ruido)
    diferencias = diff(t);
    umbral = 1.5 * dt;
    indices_gaps = find(diferencias > umbral);
    
    if isempty(indices_gaps)
        % No hay discontinuidades
        t_filled = t;
        v_filled = v;
        fprintf('No se detectaron discontinuidades.\n');
        return;
    end
    
    fprintf('Se detectaron %d discontinuidades.\n', length(indices_gaps));
    
    % Crear vectores de salida
    t_filled = t(1);
    v_filled = v(1);
    
    % Procesar cada segmento
    for i = 1:length(t)
        if i == 1
            continue; % Ya agregamos el primer elemento
        end
        
        % Verificar si hay un gap antes de este índice
        if any(indices_gaps == i-1)
            % Hay una discontinuidad entre i-1 e i
            t_inicio = t(i-1);
            t_fin = t(i);
            v_inicio = v(i-1);
            v_fin = v(i);
            
            % Crear vector de tiempos interpolados
            t_interp = (t_inicio + dt : dt : t_fin - dt/2)';
            
            % Interpolar valores de v linealmente
            if ~isempty(t_interp)
                v_interp = interp1([t_inicio; t_fin], [v_inicio; v_fin], t_interp, 'linear');
                
                % Agregar puntos interpolados
                t_filled = [t_filled; t_interp];
                v_filled = [v_filled; v_interp];
                
                fprintf('Gap %d: Rellenado desde t=%.3f hasta t=%.3f con %d puntos.\n', ...
                    find(indices_gaps == i-1), t_inicio, t_fin, length(t_interp));
            end
        end
        
        % Agregar el punto actual
        t_filled = [t_filled; t(i)];
        v_filled = [v_filled; v(i)];
    end
    
    fprintf('Vector original: %d puntos -> Vector rellenado: %d puntos.\n', ...
        length(t), length(t_filled));
end

% ========== EJEMPLO DE USO ==========
% Crear datos de ejemplo con discontinuidades
fs = 100; % Hz
dt = 1/fs;

% Crear segmentos de tiempo con gaps
t1 = (0:dt:1)';           % Primer segmento: 0 a 1 segundo
t2 = (2:dt:3)';           % Segundo segmento: 2 a 3 segundos (gap de 1 seg)
t3 = (4.5:dt:5.5)';       % Tercer segmento: 4.5 a 5.5 segundos (gap de 1.5 seg)

% Concatenar
t = [t1; t2; t3];

% Crear valores de v (por ejemplo, una señal sinusoidal)
v = sin(2*pi*0.5*t) + 0.1*randn(size(t));

% Llamar a la función
[t_filled, v_filled] = rellenar_discontinuidades(t, v, fs);

% Visualizar resultados
figure('Position', [100, 100, 1200, 500]);

subplot(2,1,1);
plot(t, v, 'bo-', 'MarkerSize', 3, 'LineWidth', 1.5);
hold on;
plot(t_filled, v_filled, 'r.-', 'MarkerSize', 2);
xlabel('Tiempo (s)');
ylabel('Valor');
title('Comparación: Original vs Rellenado');
legend('Original con gaps', 'Rellenado', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t_filled, v_filled, 'k-', 'LineWidth', 1);
xlabel('Tiempo (s)');
ylabel('Valor');
title('Señal Rellenada (continua)');
grid on;