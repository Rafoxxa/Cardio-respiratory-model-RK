function cell_of_pars = load_pars_to_fit(patient, test_folder)
    if nargin < 2
        path = sprintf('../Fitting/pars2fit/%d/last/', patient);  % Usar char array en lugar de string
    else
        if strcmp(test_folder, '')
            path = sprintf('../Fitting/pars2fit/%d/last/', patient);  % Usar char array en lugar de string
        else
            path = sprintf('../Fitting/pars2fit/%d/%s/last/', patient, test_folder);  % Usar char array en lugar de string
        end
    end
    files = dir(fullfile(path, '*.mat'));
    filename_path = fullfile(path, files(1).name);
    pars_to_fit_array = load(filename_path);
    cell_of_pars = cellstr(pars_to_fit_array.pars_to_fit);
end