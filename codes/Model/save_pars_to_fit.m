function save_pars_to_fit(patient, tag, pars_to_fit, folder)

   if nargin < 4
    folder = 'last';  % Usar char array en lugar de string
    end
    
    
    if strcmp(folder, 'historic')  % Usar strcmp en lugar de ==                 
        path = sprintf('../Fitting/pars2fit/%d/historic/', patient);    
        if ~exist(path, 'dir')
            mkdir(path);
        end
        currentDate = now;  % Usar now() en lugar de datetime('today')
        formattedDate = datestr(currentDate, 'dd-mm-yyyy');
        name = sprintf('%s_%s.mat', tag, formattedDate);
        filename = sprintf('%s%s', path, name);    
        save(filename, 'pars_to_fit', 'filename');
    elseif strcmp(folder, 'last')  % Usar strcmp en lugar de ==
        last_path = sprintf('../Fitting/pars2fit/%d/last/', patient);
        hist_path = sprintf('../Fitting/pars2fit/%d/historic/', patient);
        
        if ~exist(last_path, 'dir')
            mkdir(last_path);
        end
        if ~exist(hist_path, 'dir')
            mkdir(hist_path);
        end

        files = dir(fullfile(last_path, '*.mat'));

        if ~isempty(files)
            old_file = files(1).name;
            src_file = fullfile(last_path, old_file);
            dst_file = fullfile(hist_path, old_file);
        
            if exist(dst_file, 'file') == 2  % En R2017, isfile podría no existir
                delete(dst_file);  % Remove existing file in historic
            end
            movefile(src_file, dst_file);  % Move from last to historic
        end

        % Save new file in "last"
        currentDate = now;  % Usar now() en lugar de datetime('today')
        formattedDate = datestr(currentDate, 'dd-mm-yyyy');
        name = sprintf('%s_%s.mat', tag, formattedDate);
        filename = fullfile(last_path, sprintf('%s', name));
        save(filename, 'pars_to_fit');
    end
end