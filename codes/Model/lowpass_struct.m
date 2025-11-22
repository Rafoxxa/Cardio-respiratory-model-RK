function struct_filt = lowpass_struct(struct_vars, t, cutoff, order, fieldsToPlot)
% LOWPASS_STRUCT aplica un filtro pasabajo Butterworth a los campos de struct_vars
% y grafica solo los campos especificados.
%
%   struct_filt = LOWPASS_STRUCT(struct_vars, t, cutoff, order, fieldsToPlot)
%
%   INPUTS:
%       struct_vars   : estructura con campos que contienen vectores de señales
%       t             : vector de tiempo (s)
%       cutoff        : frecuencia de corte (Hz)
%       order         : orden del filtro (opcional, por defecto 4)
%       fieldsToPlot  : celda con nombres de campos a graficar, ej. {'fsh','hr'}
%
%   OUTPUT:
%       struct_filt   : estructura filtrada con los mismos campos que struct_vars
%
%   Ejemplo:
%       struct_filt = lowpass_struct(struct_vars, t, 0.2, 4, {'fsh','hr'});

    % -------------------------
    % Validación de argumentos
    % -------------------------
    if nargin < 4 || isempty(order)
        order = 4;
    end
    if nargin < 5
        fieldsToPlot = {};
    end

    if ~isvector(t)
        error('El argumento "t" debe ser un vector de tiempo.');
    end

    % -------------------------
    % Calcular frecuencia de muestreo
    % -------------------------
    t = t(:);
    dt = mean(diff(t));
    fs = 1 / dt;
    fprintf('Frecuencia de muestreo estimada: %.4f Hz\n', fs);

    % -------------------------
    % Diseño del filtro pasabajo
    % -------------------------
    [b, a] = butter(order, cutoff / (fs/2), 'low');

    % Inicializa estructura de salida
    struct_filt = struct_vars;

    % -------------------------
    % Aplicar filtro a cada campo
    % -------------------------
    fnames = fieldnames(struct_vars);
    for i = 1:numel(fnames)
        fname = fnames{i};
        data = struct_vars.(fname);

        if isvector(data)
            filt_data = filtfilt(b, a, data(:));
        else
            filt_data = filtfilt(b, a, data);
        end

        struct_filt.(fname) = filt_data;

        % -------------------------
        % Graficar solo campos seleccionados
        % -------------------------
        % if ismember(fname, fieldsToPlot)
        %     figure('Name', ['Filtro pasabajo: ', fname], 'NumberTitle', 'off');
        %     plot(t, data, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Original');
        %     hold on;
        %     plot(t, filt_data, 'LineWidth', 1.5, 'DisplayName', 'Filtrada');
        %     xlabel('Tiempo (s)');
        %     ylabel(fname);
        %     title(sprintf('Filtro pasabajo Butterworth - %s (fc = %.3f Hz)', fname, cutoff));
        %     legend('show');
        %     grid on;
        % end
    end
end
