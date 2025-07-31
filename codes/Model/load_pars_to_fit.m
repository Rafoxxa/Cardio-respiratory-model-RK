function cell_of_pars = load_pars_to_fit(patient)
    path = sprintf("../Fitting/pars2fit/%d/last/", patient);
    files = dir(fullfile(path, "*.mat"));
    filename_path = fullfile(path, files(1).name);
    pars_to_fit_array = load(filename_path);
    cell_of_pars = cellstr(pars_to_fit_array.pars_to_fit);

     
end