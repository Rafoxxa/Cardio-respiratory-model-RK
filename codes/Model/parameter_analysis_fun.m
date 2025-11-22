function parameter_analysis_fun(p, pindex, load_rb)

    % Convertir los inputs (que llegan como string desde el .sh)

    load_rb_ = logical(str2double(load_rb));
    if ischar(p) || isstring(p)
        p = str2double(p);
    end

    if ischar(pindex) || isstring(pindex)
        pindex = str2double(pindex);
    end

    patients = [1, 4, 5, 6];

    % Validación opcional para evitar errores silenciosos
    if isnan(p) || isnan(pindex)
        error('Los argumentos p o pindex no son números válidos.');
    end

    %% SENSITIVITY ANALYSIS
    patient_idx = patients(p);    
    setup_n = set_up('sens', patient_idx, 'normoxia', '.', 'pars_from_fitting', 1, 'time_from_data', 1,  'param_index_sens', pindex, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'load_rb', load_rb_);
    setup_h = set_up('sens', patient_idx, 'hipoxia', 'mix', 'pars_from_fitting', 1, 'time_from_data', 1, 'param_index_sens', pindex, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'load_rb', load_rb_);
    setup = {setup_n, setup_h};
    sens_functions('saving', '-', setup);
end
