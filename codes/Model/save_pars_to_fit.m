function save_pars_to_fit(patient, tag, pars_to_fit, folder)

   if nargin < 4
    folder = "last";
    end
    
    
    if folder  == "historic"                     
        path = sprintf("../Fitting/pars2fit/%d/historic/", patient);    
        if ~exist(path, 'dir')
            mkdir(path);
        end
        currentDate = datetime('today');  
        formattedDate = datestr(currentDate, 'dd-mm-yyyy');
        name = sprintf("%s_%s.mat", tag, formattedDate);
        filename = sprintf("%s%s", path, name);    
        save(filename, "pars_to_fit", "filename");
    elseif folder == "last"
        last_path = sprintf("../Fitting/pars2fit/%d/last/", patient);
        hist_path = sprintf("../Fitting/pars2fit/%d/historic/", patient);
        
        if ~exist(last_path, 'dir')
            mkdir(last_path);
        end
        if ~exist(hist_path, 'dir')
            mkdir(hist_path);
        end

        files = dir(fullfile(last_path, "*.mat"));

        if ~isempty(files)
            old_file = files(1).name;
            src_file = fullfile(last_path, old_file);
            dst_file = fullfile(hist_path, old_file);
        
            if isfile(dst_file)
                delete(dst_file);  % Remove existing file in historic
            end
            movefile(src_file, dst_file);  % Move from last to historic
        end

        % Save new file in "last"
        currentDate = datetime('today');  
        formattedDate = datestr(currentDate, 'dd-mm-yyyy');
        name = sprintf("%s_%s.mat", tag, formattedDate);
        filename = fullfile(last_path, sprintf("%s", name));
        save(filename, "pars_to_fit");
    end
                 

end