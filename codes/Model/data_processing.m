function out = data_processing(mode, var, time)
    if strcmp(mode, 'pressure')
        [pm, ps, pd] = compute_pressure(var, time);
        out = {pm, ps, pd};
    elseif strcmp(mode, 'volume')  
        out = compute_VT(var, time);
    elseif strcmp(mode, 'filtering')
        out = filtering(var); 
    elseif strcmp(mode, 'add-desired')
        out = add_desired_variables(var{1}, var{2}, time);
    elseif strcmp(mode, 'match-size')
        [arr1,arr2] = matchSize(var{1}, var{2});
        out = {arr1, arr2};
    elseif strcmp(mode, 'downsample')  
        out = downsampleArray(var{1}, var{2});
    end

    function [pm, ps, pd] = compute_pressure(presion, t)
        % Crear un vector de tiempo basado en el tamaño de 'presion'
        n = length(presion);  % Número de puntos en la curva
        tiempo = t;  % Suponiendo un paso de tiempo uniforme
        
        % Derivadas de la curva de presión
        dp = diff(presion)./diff(tiempo);  % Primera derivada
        d2p = diff(dp)./diff(tiempo(1:end-1));  % Segunda derivada
        
        % Encontrar picos sistólicos (máximos locales) y diastólicos (mínimos locales)
        sistole_indices = find(dp(1:end-1) > 0 & dp(2:end) <= 0 & d2p < 0) + 1; % Máximos locales
        diastole_indices = find(dp(1:end-1) < 0 & dp(2:end) >= 0 & d2p > 0) + 1; % Mínimos locales
        
        % Extraer valores de presión sistólica y diastólica
        presion_sistolica = presion(sistole_indices);
        presion_diastolica = presion(diastole_indices);
        
        % Interpolación de los valores encontrados para tener un vector continuo
        ps = interp1(tiempo(sistole_indices), presion_sistolica, tiempo, 'linear', 'extrap');
        pd = interp1(tiempo(diastole_indices), presion_diastolica, tiempo, 'linear', 'extrap');
        
        % Calcular la presión media
        pm = trapz(tiempo, presion) / n;  % Integral dividida por la longitud del vector
        pm = ones(size(pd)) * pm;
        %[pm, ps, pd] = compute_presion(presion);
    
    end
    
    function [VT] = compute_VT(V, t)
        % Crear un vector de tiempo basado en el tamaño de 'presion'
        n = length(V);  % Número de puntos en la curva
        tiempo = t;  % Suponiendo un paso de tiempo uniforme
        
        % Derivadas de la curva de presión
        dp = diff(V)./diff(tiempo);  % Primera derivada
        d2p = diff(dp)./diff(tiempo(1:end-1));  % Segunda derivada
        
        % Encontrar picos sistólicos (máximos locales) y diastólicos (mínimos locales)
        VT_indices = find(dp(1:end-1) > 0 & dp(2:end) <= 0 & d2p < 0) + 1; % Máximos locales
        volumen_tidal = V(VT_indices);

        VT = interp1(tiempo(VT_indices), volumen_tidal, tiempo, 'linear', 'extrap');
    
    end

    function [BF] = compute_BF(TI, TE)
        BF = 1./(TI + TE);
    end



    function X = filtering(X)
        %Filtering
        fs = 1/0.995;
        cutoff = fs/80;
        %for i = 1:size(X,2)
        X(:, 5) = lowpass(X(:, 5),cutoff, fs, 'Steepness', 0.95);
        X(:, 6) = lowpass(X(:, 6),cutoff, fs, 'Steepness', 0.95);                                          
        X(:, 7) = lowpass(X(:, 7),cutoff, fs, 'Steepness', 0.95);      
        X(:, 8) = lowpass(X(:, 8),cutoff, fs, 'Steepness', 0.95);                                          
        X(:, 9) = lowpass(X(:, 9),cutoff, fs, 'Steepness', 0.95);                 
        X(:, 10) = lowpass(X(:, 10),cutoff, fs, 'Steepness', 0.95);                 
        %end 
    end

    function [Xout] = add_desired_variables(Xin, init, t)
        
        key_press = 'P_sa';
        key_V = 'V';
        key_TI = 'TI';
        key_TE = 'TE';
        
        press_idx = find(strcmp(init.keys, key_press));
        V_idx = find(strcmp(init.keys, key_V));
        TI_idx = find(strcmp(init.keys, key_TI));
        TE_idx = find(strcmp(init.keys, key_TE));

        pressure = Xin(press_idx, :);
        V = Xin(V_idx, :);
        TI = Xin(TI_idx, :);
        TE = Xin(TE_idx, :);        

        [pm, ps, pd] = compute_pressure(pressure, t);
        [VT] = compute_VT(V, t);
        [BF] = compute_BF(TI, TE);

        Xout = [Xin; pm; ps; pd; VT; BF];
    end

    function [array1, array2] = matchSize(array1, array2)
        sz1 = size(array1, 2);
        sz2 = size(array2, 2);

        if sz1 == sz2
            return;
        end

        targetLength = min(sz1, sz2);
        array1 = array1(:, 1:targetLength);
        array2 = array2(:, 1:targetLength);

        %if sz1 > sz2            
            %array1 = downsampleArray(array1, targetLength);
        %else            
            %array2 = downsampleArray(array2, targetLength);
        %end
    end

    function outArray = downsampleArray(inArray, newSize)
        indices = round(linspace(1, length(inArray), newSize));
        outArray = inArray(indices);
    end



end