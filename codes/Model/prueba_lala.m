function archivoPatch = obtenerArchivoReciente(carpeta, tipoCondicion)
    % OBTENERARCHIVORECIENTE Obtiene el path del archivo más reciente según la condición
    %
    % Sintaxis:
    %   archivoPatch = obtenerArchivoReciente(carpeta, tipoCondicion)
    %
    % Parámetros:
    %   carpeta - string/char con la ruta de la carpeta (ej: 'C:/datos/' o 'C:\datos\')
    %   tipoCondicion - string/char: 'hipoxia' o 'normoxia'
    %
    % Salida:
    %   archivoPatch - char con el path completo del archivo más reciente
    %                 vacío si no se encuentran archivos
    %
    % Ejemplo:
    %   archivo = obtenerArchivoReciente('C:/datos/', 'hipoxia');
    %   archivo = obtenerArchivoReciente('C:\mis_datos\', 'normoxia');
    %   archivo = obtenerArchivoReciente('/home/user/experimentos/', 'hipoxia');
    
    % Validar parámetros de entrada
    if nargin < 2
        error('Se requieren ambos parámetros: carpeta y tipoCondicion');
    end
    
    % Convertir a char si es necesario (compatibilidad MATLAB 2017)
    if iscell(carpeta)
        carpeta = carpeta{1};
    end
    if iscell(tipoCondicion)
        tipoCondicion = tipoCondicion{1};
    end
    
    % Convertir a char (compatible con MATLAB 2017)
    carpeta = char(carpeta);
    tipoCondicion = char(tipoCondicion);
    
    % Validar tipo de condición
    tipoCondicion_lower = lower(tipoCondicion);
    if ~(strcmp(tipoCondicion_lower, 'hipoxia') || strcmp(tipoCondicion_lower, 'normoxia'))
        error('tipoCondicion debe ser "hipoxia" o "normoxia"');
    end
    
    % Verificar que la carpeta existe (compatibilidad MATLAB 2017)
    if exist(carpeta, 'dir') ~= 7
        error('La carpeta especificada no existe: %s', carpeta);
    end
    
    % Buscar archivos .mat que contengan la condición especificada
    patronBusqueda = fullfile(carpeta, ['*' tipoCondicion_lower '*.mat']);
    archivos = dir(patronBusqueda);
    
    % Si no se encuentran archivos, retornar vacío
    if isempty(archivos)
        archivoPatch = '';
        warning('No se encontraron archivos para la condición: %s', tipoCondicion);
        return;
    end
    
    % Extraer fechas de los nombres de archivos
    fechasNum = [];
    nombresValidos = {};
    contador = 0;
    
    for i = 1:length(archivos)
        nombre = archivos(i).name;
        
        % Buscar patrón de fecha DD-MM-YYYY al final del nombre
        patron = '(\d{2})-(\d{2})-(\d{4})\.mat$';
        coincidencia = regexp(nombre, patron, 'tokens');
        
        if ~isempty(coincidencia)
            % Extraer día, mes y año
            dia = str2double(coincidencia{1}{1});
            mes = str2double(coincidencia{1}{2});
            anio = str2double(coincidencia{1}{3});
            
            % Validar fecha manualmente
            if validarFecha(dia, mes, anio)
                % Convertir fecha a número (formato YYYYMMDD para comparación)
                fechaNum = anio * 10000 + mes * 100 + dia;
                contador = contador + 1;
                fechasNum(contador) = fechaNum;
                nombresValidos{contador} = nombre;
            else
                warning('Fecha inválida encontrada en archivo: %s', nombre);
            end
        else
            warning('Formato de fecha no reconocido en archivo: %s', nombre);
        end
    end
    
    % Si no hay fechas válidas, retornar vacío
    if isempty(fechasNum)
        archivoPatch = '';
        warning('No se encontraron archivos con formato de fecha válido');
        return;
    end
    
    % Encontrar el índice del archivo con fecha más reciente
    [~, indiceReciente] = max(fechasNum);
    
    % Construir el path completo del archivo más reciente
    archivoReciente = nombresValidos{indiceReciente};
    archivoPatch = fullfile(carpeta, archivoReciente);
    
    fprintf('Archivo más reciente encontrado: %s\n', archivoReciente);
end

function esValida = validarFecha(dia, mes, anio)
    % Función auxiliar para validar fechas manualmente
    esValida = false;
    
    % Verificar rangos básicos
    if anio < 1900 || anio > 2100 || mes < 1 || mes > 12 || dia < 1 || dia > 31
        return;
    end
    
    % Días por mes
    diasPorMes = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    % Verificar año bisiesto
    if mod(anio, 4) == 0 && (mod(anio, 100) ~= 0 || mod(anio, 400) == 0)
        diasPorMes(2) = 29;
    end
    
    % Verificar si el día es válido para el mes
    if dia <= diasPorMes(mes)
        esValida = true;
    end
end