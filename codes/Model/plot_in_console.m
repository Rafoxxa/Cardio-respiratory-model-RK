function plot_in_console(data, max_width, mask)
    % plot_in_console(data, max_width, mask)
    %   Visualiza un gráfico de barras vertical con eje Y y permite resaltar
    %   zonas especiales usando una máscara lógica (con '*').
    %
    % Inputs:
    %   data      : vector de datos numéricos
    %   max_width : (opcional) número máximo de columnas (default: 80)
    %   mask      : (opcional) vector lógico del mismo tamaño que data

    if nargin < 2 || isempty(max_width)
        max_width = 80;
    end
    if nargin < 3 || isempty(mask)
        mask = false(size(data));
    end

    max_height = 10;

    % Asegurar que la máscara tenga la misma longitud que data
    mask = logical(mask(:));
    data = data(:);  % asegurar columna

    % Interpolar si es más largo que el ancho deseado
    if length(data) > max_width
        x = linspace(1, length(data), max_width);
        data_interp = interp1(1:length(data), data, x);
        mask_interp = interp1(1:length(mask), double(mask), x) > 0.5;
    else
        data_interp = data;
        mask_interp = mask;
    end

    % En R2017, rescale no existe - usar implementación manual
    data_min = min(data_interp);
    data_max = max(data_interp);
    if data_max == data_min
        scaled = ones(size(data_interp));  % Evitar división por cero
    else
        scaled = round(1 + (data_interp - data_min) * (max_height - 1) / (data_max - data_min));
    end
    
    real_max = max(data_interp);
    real_min = min(data_interp);
    real_range = linspace(real_min, real_max, max_height);
    label_width = 10;

    for row = max_height:-1:1
        label = sprintf('%*.2f |', label_width, real_range(row));
        line = '';
        for i = 1:length(scaled)
            if scaled(i) >= row
                if mask_interp(i)
                    line = [line, '*'];  % zona especial
                else
                    line = [line, '|'];
                end
            else
                line = [line, ' '];
            end
        end
        disp([label, line])
    end
end