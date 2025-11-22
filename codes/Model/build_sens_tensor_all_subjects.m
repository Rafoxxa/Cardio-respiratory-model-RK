% function STensor = build_sens_tensor_all_subjects()
% 
% subjects = [1,4,5,6];
% base_folder = fullfile('..','Sens_analysis');
% 
% all_subject_params = {};     % Lista de parámetros que tiene cada sujeto
% subject_data      = struct;  % Guardaremos matrices temporalmente por sujeto
% 
% %% ---- 1) Leer parámetros disponibles por sujeto ----
% for s = subjects
%     folder = fullfile(base_folder, num2str(s), 'IndependentFilesServer');
%     files  = dir(fullfile(folder, 'SensTensor_*.mat'));
% 
%     params_s = {};
%     temp_mats = struct();
% 
%     for i = 1:length(files)
%         S = load(fullfile(folder, files(i).name));
% 
%         % Nombre del parámetro
%         if isfield(S, 'pars_to_sens')
%             pname = S.pars_to_sens{1};
%         else
%             error('Archivo %s no tiene pars_to_sens', files(i).name);
%         end
% 
%         % Extraer matriz VxT
%         varNames = fieldnames(S);
%         varNames = varNames(~cellfun(@(x) iscell(S.(x)), varNames));
%         mat      = S.(varNames{1});
% 
%         params_s{end+1} = pname;
%         temp_mats.(pname) = mat;
%     end
% 
%     all_subject_params{end+1} = params_s;
%     subject_data.(sprintf('S%d', s)) = temp_mats;
% end
% 
% %% ---- 2) Encontrar parámetros comunes ----
% common_params = all_subject_params{1};
% for i = 2:length(all_subject_params)
%     common_params = intersect(common_params, all_subject_params{i}, 'stable');
% end
% 
% fprintf('\nSe usarán %d parámetros comunes:\n', length(common_params));
% disp(common_params')
% 
% if isempty(common_params)
%     error('No hay parámetros en común entre los sujetos.');
% end
% 
% %% ---- 3) Construir el tensor concatenando solo los comunes ----
% ST_all = [];
% 
% for s = subjects
%     subj_name = sprintf('S%d', s);
%     mats = subject_data.(subj_name);
% 
%     P = length(common_params);
%     subj_tensor = [];
% 
%     for i = 1:P
%         pname = common_params{i};
% 
%         if ~isfield(mats, pname)
%             error('Parámetro %s faltante en sujeto %d (esto no debería pasar)', pname, s);
%         end
% 
%         subj_tensor(i,:,:) = mats.(pname);
%     end
% 
%     % Concatenar en la dimensión del tiempo
%     ST_all = cat(3, ST_all, subj_tensor);
% end
% 
% %% ---- 4) Guardar el tensor final ----
% STensor.tensor     = ST_all;
% STensor.parameters = common_params;
% 
% save(fullfile(base_folder, 'STensor.mat'), 'STensor');
% 
% fprintf('\nTensor final creado.\n');
% fprintf('Dimensiones: P = %d, V = %d, T = %d\n', size(ST_all,1), size(ST_all,2), size(ST_all,3));
% fprintf('Archivo guardado en: %s\n', fullfile(base_folder, 'STensor.mat'));
% 
% end

function STensor = build_sens_tensor_all_subjects()

subjects = [5];
base_folder = fullfile('..','Sens_analysis');
target_date = '20-11-2025';   % <-- fecha que quieres filtrar

Ttarget = 45053;  % <---- NUEVO: tamaño fijo al que interpolar

all_subject_params = {};     % Lista de parámetros que tiene cada sujeto
subject_data      = struct;  % Guardaremos matrices temporalmente por sujeto

%% ---- 1) Leer parámetros disponibles por sujeto ----
for s = subjects
    folder = fullfile(base_folder, num2str(s), 'IndependentFIlesServerEps-5Corrections');

    % Buscar solo archivos con la fecha específica
    pattern = sprintf('SensTensor_*_%s.mat', target_date);
    files  = dir(fullfile(folder, pattern));

    if isempty(files)
        warning('No se encontraron archivos con fecha %s en sujeto %d.', target_date, s);
        continue;
    end

    params_s = {};
    temp_mats = struct();

    for i = 1:length(files)
        S = load(fullfile(folder, files(i).name));

        % Nombre del parámetro
        if isfield(S, 'pars_to_sens')
            pname = S.pars_to_sens{1};
        else
            error('Archivo %s no tiene pars_to_sens', files(i).name);
        end

        % Extraer matriz VxT (busca la primera variable numérica)
        varNames = fieldnames(S);
        varNames = varNames(~cellfun(@(x) iscell(S.(x)), varNames));
        mat      = S.(varNames{1});

        params_s{end+1} = pname;
        temp_mats.(pname) = mat;
    end

    all_subject_params{end+1} = params_s;
    subject_data.(sprintf('S%d', s)) = temp_mats;
end

%% ---- 2) Encontrar parámetros comunes ----
common_params = all_subject_params{1};
for i = 2:length(all_subject_params)
    common_params = intersect(common_params, all_subject_params{i}, 'stable');
end

fprintf('\nSe usarán %d parámetros comunes:\n', length(common_params));
disp(common_params')

if isempty(common_params)
    error('No hay parámetros en común entre los sujetos.');
end

%% ---- 3) Construir el tensor concatenando solo los comunes ----
ST_all = [];
torigss = [];

for s = subjects
    subj_name = sprintf('S%d', s);
    if ~isfield(subject_data, subj_name)
        warning('Sujeto %d no tiene datos cargados. Se omite.', s);
        continue;
    end

    mats = subject_data.(subj_name);
    P = length(common_params);
    subj_tensor = [];

    for i = 1:P
        pname = common_params{i};

        if ~isfield(mats, pname)
            warning('Parámetro %s faltante en sujeto %d (omitido).', pname, s);
            continue;
        end

        mat = mats.(pname);   % V x T_orig
        [~,V, T_orig] = size(mat);
        torigss(i) = T_orig;

        % ----------------------------
        %     INTERPOLACIÓN AQUÍ
        % ----------------------------
        t_orig = linspace(0, 1, T_orig);

        t_new  = linspace(0, 1, Ttarget);
        mat_interp = zeros(V, Ttarget);

        for v = 1:V
            mat_interp(v, :) = interp1(t_orig, squeeze(mat(1,v,:))', t_new, 'linear');
        end
        % ----------------------------

        subj_tensor(i,:,:) = mat_interp;   % (P) x V x Ttarget
    end

    % Concatenar en la dimensión del tiempo: 3
    ST_all = cat(3, ST_all, subj_tensor);
end

%% ---- 4) Guardar el tensor final ----
STensor.tensor     = ST_all;
STensor.parameters = common_params;

save(fullfile(base_folder, sprintf('STensorEps-5C_%s.mat', target_date)), 'STensor');

fprintf('\nTensor final creado.\n');
fprintf('Dimensiones: P = %d, V = %d, T = %d\n', size(ST_all,1), size(ST_all,2), size(ST_all,3));
fprintf('Archivo guardado en: %s\n', fullfile(base_folder, sprintf('STensor_%s.mat', target_date)));

end



