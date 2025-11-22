function [t_filled, v_filled] = rellenar_discontinuidades(t, v, fs, intervalos_eliminar)
    % Función para rellenar discontinuidades en vectores de tiempo y datos
    % 
    % Entradas:
    %   t  - Vector de tiempos (puede tener discontinuidades)
    %   v  - Vector de valores asociados a cada tiempo
    %   fs - Frecuencia de muestreo (Hz)
    %   intervalos_eliminar - (Opcional) Matriz Nx2 con [t_inicio, t_fin] de intervalos
    %                         a eliminar de v y rellenar con interpolación
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
    
    %% PASO 1: Eliminar intervalos especificados de v
    if nargin >= 4 && ~isempty(intervalos_eliminar)
        fprintf('=== ELIMINACIÓN DE INTERVALOS ESPECIFICADOS ===\n');
        
        if size(intervalos_eliminar, 2) ~= 2
            error('intervalos_eliminar debe ser una matriz Nx2 con [t_inicio, t_fin]');
        end
        
        mask_eliminar = false(size(v));
        
        for i = 1:size(intervalos_eliminar, 1)
            t_ini = intervalos_eliminar(i, 1);
            t_fin = intervalos_eliminar(i, 2);
            
            % Encontrar índices en ese rango de tiempo
            idx_intervalo = find(t >= t_ini & t <= t_fin);
            
            if ~isempty(idx_intervalo)
                mask_eliminar(idx_intervalo) = true;
                fprintf('  Intervalo %d: t=[%.3f, %.3f]s -> %d puntos marcados para eliminar\n', ...
                    i, t_ini, t_fin, length(idx_intervalo));
            else
                fprintf('  Intervalo %d: t=[%.3f, %.3f]s -> No se encontraron puntos\n', ...
                    i, t_ini, t_fin);
            end
        end
        
        % Eliminar los puntos marcados
        fprintf('Total de puntos eliminados: %d\n', sum(mask_eliminar));
        t(mask_eliminar) = [];
        v(mask_eliminar) = [];
        fprintf('\n');
    end
    
    %% PASO 2: Detectar discontinuidades en tiempo
    fprintf('=== DETECCIÓN DE GAPS TEMPORALES ===\n');
    
    % Detectar discontinuidades
    diferencias = diff(t);
    umbral_tiempo = 1.5 * dt;
    indices_gaps = find(diferencias > umbral_tiempo);
    
    if isempty(indices_gaps)
        t_filled = t;
        v_filled = v;
        fprintf('No se detectaron discontinuidades temporales.\n');
        return;
    end
    
    fprintf('Se detectaron %d discontinuidades temporales.\n', length(indices_gaps));
    
    %% PASO 3: Rellenar gaps
    t_filled = t(1);
    v_filled = v(1);
    
    % Procesar cada segmento
    for i = 1:length(t)
        if i == 1
            continue;
        end
        
        % Verificar si hay un gap antes de este índice
        if any(indices_gaps == i-1)
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
                
                fprintf('  Gap %d: t=[%.3f, %.3f]s -> %d puntos interpolados.\n', ...
                    find(indices_gaps == i-1), t_inicio, t_fin, length(t_interp));
            end
        end
        
        % Agregar el punto actual
        t_filled = [t_filled; t(i)];
        v_filled = [v_filled; v(i)];
    end
    
    fprintf('\n=== RESUMEN ===\n');
    fprintf('Puntos finales: %d\n', length(t_filled));
end

% % ========== EJEMPLO DE USO ==========
% % Crear datos de ejemplo
% fs = 100; % Hz
% dt = 1/fs;
% 
% % Crear segmentos de tiempo con gaps temporales
% t1 = (0:dt:2)';
% t2 = (3:dt:5)';      % Gap temporal de 1 seg
% t3 = (6.5:dt:8.5)';  % Gap temporal de 1.5 seg
% 
% t = [t1; t2; t3];
% 
% % Crear señal base
% v = 2*sin(2*pi*0.3*t) + 0.2*randn(size(t));
% 
% % Agregar intervalos con valores anómalos (outliers)
% idx_outlier1 = find(t >= 1.0 & t <= 1.5);
% v(idx_outlier1) = v(idx_outlier1) + 10;
% 
% idx_outlier2 = find(t >= 4.0 & t <= 4.8);
% v(idx_outlier2) = v(idx_outlier2) - 12;
% 
% % Guardar datos originales
% t_original = t;
% v_original = v;
% 
% %% CASO 1: Solo rellenar gaps temporales (sin eliminar outliers)
% fprintf('========== CASO 1: Solo gaps temporales ==========\n');
% [t_filled1, v_filled1] = rellenar_discontinuidades(t_original, v_original, fs);
% 
% fprintf('\n\n');
% 
% %% CASO 2: Especificar manualmente intervalos a eliminar
% fprintf('========== CASO 2: Eliminar intervalos específicos + gaps temporales ==========\n');
% 
% % Definir intervalos problemáticos que quieres eliminar
% % Cada fila es [t_inicio, t_fin]
% intervalos_eliminar = [
%     1.0, 1.5;    % Primer intervalo anómalo
%     4.0, 4.8     % Segundo intervalo anómalo
% ];
% 
% [t_filled2, v_filled2] = rellenar_discontinuidades(t_original, v_original, fs, intervalos_eliminar);
% 
% %% Visualizar resultados
% figure('Position', [100, 100, 1400, 800]);
% 
% % Gráfico 1: Señal original
% subplot(3,1,1);
% plot(t_original, v_original, 'b.-', 'MarkerSize', 4);
% hold on;
% % Marcar los intervalos a eliminar
% for i = 1:size(intervalos_eliminar, 1)
%     idx = find(t_original >= intervalos_eliminar(i,1) & t_original <= intervalos_eliminar(i,2));
%     if ~isempty(idx)
%         plot(t_original(idx), v_original(idx), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
%     end
% end
% xlabel('Tiempo (s)');
% ylabel('Valor');
% title('Señal Original (círculos rojos = intervalos a eliminar)');
% legend('Señal original', 'Intervalos a eliminar', 'Location', 'best');
% grid on;
% 
% % Gráfico 2: Solo gaps temporales
% subplot(3,1,2);
% plot(t_filled1, v_filled1, 'g.-', 'MarkerSize', 3);
% xlabel('Tiempo (s)');
% ylabel('Valor');
% title('Caso 1: Solo gaps temporales rellenados (outliers se mantienen)');
% grid on;
% 
% % Gráfico 3: Intervalos eliminados + gaps rellenados
% subplot(3,1,3);
% plot(t_filled2, v_filled2, 'r.-', 'MarkerSize', 3);
% xlabel('Tiempo (s)');
% ylabel('Valor');
% title('Caso 2: Intervalos eliminados + gaps temporales rellenados');
% grid on;