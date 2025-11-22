%%%%%%%%%%%%%%%%%%%%%%%
function setup = ForwardFittingModel(mode, patient_idx, old_date, new_date)
    if nargin < 3
        old_date = 'error';
        new_date = 'error';
    end
 
    
    if ischar(patient_idx)
        idx = str2double(patient_idx);
    else
        idx = patient_idx;
    end
    iterations = 'None';
    initial_noise = 0.8;   
    patient_list = [1,4,5,6];
    patient_number = patient_list(idx);
    %parsfolder = 'custom';
    %parsfolder = 'test-cardiovasc';
    parsfolder = 'custom';
    previous = {'07-09-2025'};
    %patient keys: number: [1,4,5,6]  idx: [1,2,3,4]
    % Setting up and run
    vectorize_dicts('run_ode.m', 'model_basic.m', 'run_ode_vec_hipoxia.m', 'model_vec_hipoxia.m'); 
    %try 
    %    parallel_fitting(patient_idx, 'find-best-solution', 'last');
    %catch
    %    disp('in the catch of generating best.mat: it was an error or was previosuly done by one node')
    %end
    %pause(5);

    if strcmp(mode, 'only-compute-J')
        bestPars_simulation_folder = sprintf('../Simulations/simulation_after_fitting/%d/', patient_number);
        [actual_n, old_n, actual_h, old_h] = simulation_files_to_load(new_date, old_date);
        actual_best_normoxia_simulation_filename = sprintf('%s%s',bestPars_simulation_folder, actual_n);
        actual_best_hipoxia_simulation_filename = sprintf('%s%s',bestPars_simulation_folder, actual_h);

        lfN = load(actual_best_normoxia_simulation_filename);
        lfH = load(actual_best_hipoxia_simulation_filename);

        struct_vars_normoxia = lfN.struct_vars_normoxia;
        t_normoxia = lfN.t_normoxia;
        struct_vars_hipoxia = lfH.struct_vars_hipoxia;
        t_hipoxia = lfH.t_hipoxia;

        fast_data_filename_hipoxia = sprintf('../fast_data/%d/%s_data_preprocessed.mat', patient_number, 'hipoxia');
        fast_data_filename_normoxia = sprintf('../fast_data/%d/%s_data_preprocessed.mat', patient_number, 'normoxia');
        lffd_normoxia = load(fast_data_filename_normoxia);
        lffd_hipoxia = load(fast_data_filename_hipoxia);
        texp_normoxia = lffd_normoxia.texp;
        yexp_normoxia = lffd_normoxia.yexp;

        texp_hipoxia = lffd_hipoxia.texp;
        yexp_hipoxia = lffd_hipoxia.yexp;

        [J_normoxia, ~] = compute_J_direct(texp_normoxia, yexp_normoxia, t_normoxia, struct_vars_normoxia);
        [J_hipoxia, ~] = compute_J_direct(texp_hipoxia, yexp_hipoxia, t_hipoxia, struct_vars_hipoxia);

        setup.J_normoxia = J_normoxia;
        setup.J_hipoxia = J_hipoxia;
        
        return



    end
    
    if strcmp(mode, 'optimize')
        %save_pars_to_fit(patient_number, "LSA_cardiovascular_parameters", next_to_optim);
        setup = parallel_fitting(patient_idx, 'pattern_search', {'07-09-2025', '30-09-2025', '20-10-2025', '02-11-2025'}, initial_noise, iterations, parsfolder);
        %setup = parallel_fitting(patient_idx, 'pattern_search', {}, initial_noise, iterations, parsfolder);
        return;
    end

   



    %%%%%%%%%%%%%%%%%%%%%%%

    % Take the best solution   
    %solver_output = exec_solver('find-best-solution', setup);
    global JH JN
    %%setup = parallel_fitting(patient_idx, 'find-best-solution-local', {'07-09-2025', '30-09-2025', '06-10-2025'}, initial_noise, iterations, parsfolder);
    %setup = parallel_fitting(patient_idx, 'find-best-solution-local', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025'}, initial_noise, iterations, parsfolder);
    %setup = parallel_fitting(patient_idx, 'find-best-solution-local', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025'}, initial_noise, iterations, parsfolder);
    setup = parallel_fitting(patient_idx, 'find-best-solution-local', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025' , '02-11-2025'}, initial_noise, iterations, parsfolder);
    out_solver = setup.out_solver;
    setup.JH = JH;
    setup.JN = JN;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    

    % Run simulation
    
        %parsfitted_filename = setup.best_fitting_filename;
        %parts = strsplit(parsfitted_filename, '/');
        %lastTwo = parts(end-1:end);
        %date = lastTwo{1}(9:end);
        %date = [previous, date];
        setup_out_normoxia_sim = set_up('simulation', patient_number, 'normoxia', '-', 'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'});
        sn = setup_out_normoxia_sim;
        setup_out_hipoxia_sim = set_up('simulation', patient_number, 'hipoxia', 'mix', 'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'});
        sh = setup_out_hipoxia_sim;
        if strcmp(mode, 'hard:only-loading')
            return
        end
   
        

        
    if ~strcmp(mode, 'only-plot')
    
        global all_global;
        all_global = zeros(15, round(10 * sn.simulation_time/sn.dt) + 1) + 0;  %This array saves all the data used for delays and for integration 
        [t_normoxia, x_dot, x_vars, x_keys, ~] = sn.run_ode_fun(sn.model, sn.pars, sn.init, sn.simulation_time, sn.dt);
        struct_vars_normoxia = arrange_results(x_dot, x_vars, x_keys, t_normoxia);
        all_global = zeros(15, round(10 * sh.simulation_time/sh.dt) + 1) + 0;  %This array saves all the data used for delays and for integration 
        [t_hipoxia, x_dot, x_vars, x_keys, ~] = sh.run_ode_fun(sh.model, sh.pars, sh.init, sh.simulation_time, sh.dt);
        struct_vars_hipoxia = arrange_results(x_dot, x_vars, x_keys, t_hipoxia);
    
        %%%%%%%%%%%%%%%%%%%%%%%%
    
        % Save simulation
        
        sn_saving_name_file = sprintf('%s.mat', sn.simulation_filename);
        sh_saving_name_file = sprintf('%s.mat', sh.simulation_filename);
        save(sn_saving_name_file, 'struct_vars_normoxia', 't_normoxia');
        save(sh_saving_name_file, 'struct_vars_hipoxia', 't_hipoxia');

        if strcmp(mode, 'only-loading')
            return
        end
    
        %%%%%%%%%%%%%%%
    end
    
        % Plot
    if ~strcmp(mode, 'optimize') && ~strcmp(mode, 'not-plot')
        
        figure;    
        old_mode = 'off';  %that means, the struct_vars have the old version       
        
        bestPars_simulation_folder = sprintf('../Simulations/simulation_after_fitting/%d/', patient_number);
        Jtot = JH(end) + JN(end);
        [actual_n, old_n, actual_h, old_h] = simulation_files_to_load(new_date, old_date);
        actual_best_normoxia_simulation_filename = sprintf('%s%s',bestPars_simulation_folder, actual_n);%getMostRecentFile(bestPars_simulation_folder, 'normoxia', 1);
        last_best_normoxia_simulation_filename = sprintf('%s%s',bestPars_simulation_folder, old_n); %getMostRecentFile(bestPars_simulation_folder, 'normoxia', 2);
        actual_best_hipoxia_simulation_filename = sprintf('%s%s',bestPars_simulation_folder, actual_h);%getMostRecentFile(bestPars_simulation_folder, 'hipoxia', 1);
        last_best_hipoxia_simulation_filename = sprintf('%s%s',bestPars_simulation_folder, old_h); %getMostRecentFile(bestPars_simulation_folder, 'hipoxia', 2);
        %time_sim = t_normoxia;
        %X_sim = struct_vars_normoxia;
        %custom_plot('sim_vs_exp', {time_sim, sn.texp, X_sim, sn.yexp, setup.xnames_fitting, setup.units_table, 5, 2, sn.simulation_filename, old_mode});
        custom_plot('sim_vs_exp', {'-', sn.texp, '-', sn.yexp, setup.xnames_fitting, setup.units_table, 5, 2, actual_best_normoxia_simulation_filename, old_mode, patient_number, 'actual', JN, 'Normoxia'});
        %original_normoxia_simulation_filename = sprintf('../Simulations/only_simulation/%d/1200_sec_normoxia-28-07-2025.mat', patient_number);
        
        
        %custom_plot('sim_vs_exp', {'-', sn.texp, '-', sn.yexp, setup.xnames_fitting, setup.units_table, 5, 2, last_best_normoxia_simulation_filename, old_mode, patient_number, 'old', JN, 'Normoxia'});
        
        figure;
        %time_sim = t_hipoxia; 
        %X_sim = struct_vars_hipoxia;
        custom_plot('sim_vs_exp', {'-', sh.texp, '-', sh.yexp, setup.xnames_fitting, setup.units_table, 5, 2, actual_best_hipoxia_simulation_filename, old_mode, patient_number, 'actual', JH, 'Hypoxia'});
        %original_hipoxia_simulation_filename = sprintf('../Simulations/only_simulation/%d/3300_sec_hipoxia-28-07-2025.mat', patient_number);
        
        %custom_plot('sim_vs_exp', {'-', sh.texp, '-', sh.yexp, setup.xnames_fitting, setup.units_table, 5, 2, last_best_hipoxia_simulation_filename, old_mode, patient_number, 'old', JH, 'Hipoxia'});

    
        
%         sh_filename = sh.simulation_filename;
%         sn_filename = sn.simulation_filename;
%         name = sprintf('../plot_data/%d/logJ_is:%.4f.mat', patient_number, out_solver.fval);
%         if ~isfolder(sprintf('../plot_data/%d', patient_idx))
%             mkdir(sprintf('../plot_data/%d', patient_idx));
%         end
%         save(name, 'sh_filename', 'sn_filename', 'original_hipoxia_simulation_filename', 'original_normoxia_simulation_filename');
    end
    
    function filePath = getMostRecentFile(folder, condition, index)
    % GETMOSTRECENTFILE Gets the path of the nth most recent file based on condition
    %
    % Syntax:
    %   filePath = getMostRecentFile(folder, condition, index)
    %
    % Parameters:
    %   folder - char with folder path
    %   condition - char: 'hipoxia' or 'normoxia'
    %   index - integer: 1 for most recent, 2 for second most recent, etc.
    %
    % Output:
    %   filePath - char with full path of nth most recent file
    
    % Convert to char
    folder = char(folder);
    condition = char(condition);
    
    % Default to most recent if index not provided
    if nargin < 3
        index = 1;
    end
    
    % Search for .mat files containing the specified condition
    searchPattern = fullfile(folder, ['*' lower(condition) '*.mat']);
    files = dir(searchPattern);
    
    % Extract dates from file names
    dateNumbers = [];
    validNames = {};
    counter = 0;
    
    for i = 1:length(files)
        name = files(i).name;
        
        % Find date pattern DD-MM-YYYY at end of name
        pattern = '(\d{2})-(\d{2})-(\d{4})\.mat$';
        match = regexp(name, pattern, 'tokens');
        
        if ~isempty(match)
            day = str2double(match{1}{1});
            month = str2double(match{1}{2});
            year = str2double(match{1}{3});
            
            % Convert date to number (YYYYMMDD format for comparison)
            dateNum = year * 10000 + month * 100 + day;
            counter = counter + 1;
            dateNumbers(counter) = dateNum;
            validNames{counter} = name;
        end
    end
    
    % Sort dates in descending order (most recent first)
    [~, sortedIndices] = sort(dateNumbers, 'descend');
    targetIndex = sortedIndices(index);
    
    % Build full path of nth most recent file
    targetFile = validNames{targetIndex};
    filePath = fullfile(folder, targetFile);
end
 
function [actual_n, old_n, actual_h, old_h] = simulation_files_to_load(new_date, old_date)    
    actual_n = sprintf('1200_sec_normoxia-%s.mat', new_date);
    old_n =  sprintf('1200_sec_normoxia-%s.mat', old_date);
    actual_h = sprintf('3300_sec_hipoxia-%s.mat', new_date);
    old_h = sprintf('3300_sec_hipoxia-%s.mat', old_date);
end

function [J, J_total, R] = compute_J_direct(t_exp, y_exp, t_sim, struct_vars)
% compute_J_direct  Compute per-variable cost J between experimental and simulated data
%   [J, J_total, R] = compute_J_direct(t_exp, y_exp, t_sim, struct_vars)
%
% Inputs:
%   t_exp       - vector tiempo experimental (Nx1 or 1xN)
%   y_exp       - matriz experimental (N x M) con columnas como en tu formato original
%   t_sim       - vector tiempo simulación (Kx1 or 1xK)
%   struct_vars - struct con campos: PAO2, PACO2, dVE, VT, TI, Tresp, HR, pm, ps, pd
%
% Outputs:
%   J       - vector (10x1) con coste por variable (orden: PAO2,PACO2,dVE,VT,TI,Tresp,HR,PM,PS,PD)
%   J_total - scalar sum(J)
%   R       - residuos normalizados (10 x T) después de sincronizar/interpolar

    % --- Experimental variables (orden usado en tu función original) ---
    PAO2_exp  = y_exp(:,5);
    PACO2_exp = y_exp(:,6);
    dVE_exp   = y_exp(:,1);
    VT_exp    = y_exp(:,2);
    TI_exp    = y_exp(:,3);
    Tresp_exp = y_exp(:,4);
    HR_exp    = y_exp(:,7) * 60;    % pass to bpm
    PS_exp    = y_exp(:,8);
    PD_exp    = y_exp(:,9);
    PM_exp    = y_exp(:,10);

    exp_vars = [PAO2_exp, PACO2_exp, dVE_exp, VT_exp, TI_exp, Tresp_exp, HR_exp, PM_exp, PS_exp, PD_exp]';

    % --- Finapres mask (penultimate column in y_exp in your pipeline) ---
    % Keep as row vector for later broadcasting
    finapres_notnan_mask = y_exp(:, end-1)';  % 1 x N_exp

    % --- Simulation vars from struct (assume column vectors of same length as t_sim) ---
    PAO2_sim  = struct_vars.PAO2;
    PACO2_sim = struct_vars.PACO2;
    dVE_sim   = struct_vars.dVE;
    VT_sim    = struct_vars.VT;
    TI_sim    = struct_vars.TI;
    Tresp_sim = struct_vars.Tresp;
    HR_sim    = struct_vars.HR;
    PM_sim    = struct_vars.pm;
    PS_sim    = struct_vars.ps;
    PD_sim    = struct_vars.pd;

    % Ensure row/column orientation: sim_vars should be variables x time (10 x K)
    sim_vars = [PAO2_sim(:)'; PACO2_sim(:)'; dVE_sim(:)'; VT_sim(:)'; TI_sim(:)'; ...
                Tresp_sim(:)'; HR_sim(:)'; PM_sim(:)'; PS_sim(:)'; PD_sim(:)'];

    % --- Ensure time vectors are row vectors ---
    t_sim = t_sim(:)';
    t_exp = t_exp(:)';

    % --- Synchronize times: crop to common interval [t_start, t_end] ---
    t_start = max(t_sim(1), t_exp(1));
    t_end   = min(t_sim(end), t_exp(end));

    if t_start >= t_end
        error('No overlapping time interval between simulation and experiment.');
    end

    % Indices inside common interval
    sim_idx = find(t_sim >= t_start & t_sim <= t_end);
    exp_idx = find(t_exp >= t_start & t_exp <= t_end);

    % Crop
    t_sim_crop = t_sim(sim_idx);
    sim_vars_crop = sim_vars(:, sim_idx);

    t_exp_crop = t_exp(exp_idx);
    exp_vars_crop = exp_vars(:, exp_idx);
    finapres_mask_crop = finapres_notnan_mask(exp_idx);

    % Remove duplicated experimental times (keep stable order)
    if ~isempty(t_exp_crop)
        [t_exp_crop_u, ia, ~] = unique(t_exp_crop, 'stable');
        exp_vars_crop = exp_vars_crop(:, ia);
        finapres_mask_crop = finapres_mask_crop(ia);
        t_exp_crop = t_exp_crop_u;
    end

    % Safety check: if after cropping arrays empty
    if isempty(t_sim_crop) || isempty(t_exp_crop)
        error('After cropping to the overlapping interval, one of the time arrays is empty.');
    end

    % --- Interpolate so both series have same time base ---
    if size(sim_vars_crop,2) > size(exp_vars_crop,2)
        % Interp experimental onto sim times
        exp_i = interp1(t_exp_crop, exp_vars_crop', t_sim_crop, 'linear', 'extrap')';  % 10 x Tsim
        sim_i = sim_vars_crop;                                                        % 10 x Tsim
        mask_i = interp1(t_exp_crop, finapres_mask_crop', t_sim_crop, 'linear', 'extrap'); % 1 x Tsim
        mask_i = mask_i; % keep as row
    else
        % Interp simulation onto exp times
        sim_i = interp1(t_sim_crop, sim_vars_crop', t_exp_crop, 'linear', 'extrap')';   % 10 x Texp
        exp_i = exp_vars_crop;                                                        % 10 x Texp
        mask_i = finapres_mask_crop;                                                  % 1 x Texp
    end

    % --- Compute normalized squared residuals R (10 x T) ---
    % protect against zeros in exp_i by using eps + exp_i.^2 as denominator (same as original)
    R = ((exp_i - sim_i).^2) ./ (eps + exp_i.^2);

    % --- Apply finapres mask to cardiovascular vars (indices 7..10 in your original) ---
    % mask_i is 1 x T, so broadcast multiply
    R(7,:)  = R(7,:)  .* mask_i;
    R(8,:)  = R(8,:)  .* mask_i;
    R(9,:)  = R(9,:)  .* mask_i;
    R(10,:) = R(10,:) .* mask_i;

    % --- (Optional) exercise weighting & NaN handling left out (you can re-add if needed) ---

    % --- Compute per-variable J as in original ---
    [sz_vars, sz_time] = size(R);
    J = sqrt( sum(R,2) ./ sz_time );   % sqrt of time-mean per variable
    J = J / sz_vars;                   % normalize by number of variables

    % --- Ignore TI and Tresp as in original ---
    J(5) = 0;
    J(6) = 0;

    % --- Total cost ---
    J_total = sum(J);

end


 

end

