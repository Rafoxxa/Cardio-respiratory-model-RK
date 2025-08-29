function logJ = obj_fun(optpars_iter,objFunParams) %agregar pa

    % Extract variables from params struct
    t_exp_list = objFunParams.texp_list;
    y_exp_list = objFunParams.yexp_list;
    dt = objFunParams.dt;
    settling_time = objFunParams.settling_time;
    init_dict = objFunParams.init;
    pars_list = objFunParams.pars_list;
    idx_optpars = objFunParams.idx_optpars;
    percentages = objFunParams.percentages;  
    simulation_time_list = objFunParams.simulation_time_list;
    initial_point = objFunParams.initial_point;
    lb = objFunParams.lb;
    ub = objFunParams.ub;
    xnames_fitting = objFunParams.xnames_fitting;

    x_keys = init_dict.keys;
    init = cell2mat(init_dict.values);
    
    simulation_time_normoxia = simulation_time_list{1} + settling_time;  
    simulation_time_hipoxia = simulation_time_list{2} + settling_time;  
    global all_global;
    model = @(varargin) model_vec_hipoxia(varargin{:});
    run_ode_fun = @(varargin) run_ode_vec_hipoxia(varargin{:});

    init_time = tic;
    if length(pars_list) < 2        
        [texp, yexp, pars] = extract_input_args(pars_list{1}, t_exp_list{1}, y_exp_list{1}, percentages, idx_optpars, optpars_iter);
    else        
        [texp_normoxia, yexp_normoxia, pars_normoxia] = extract_input_args(pars_list{1}, t_exp_list{1}, y_exp_list{1}, percentages, idx_optpars, optpars_iter);
        [texp_hipoxia, yexp_hipoxia, pars_hipoxia] = extract_input_args(pars_list{2}, t_exp_list{2}, y_exp_list{2}, percentages, idx_optpars, optpars_iter);
    end
   %[t, x_dot, x_vars, x_keys, index] = run_ode_fun(model, pars, init, taus, simulation_time, dt, control_on);
   %isp("fin");
     try
        if length(pars_list) < 2
            J = compute_J(texp, yexp, model, pars, init, simulation_time, dt);
        elseif length(pars_list) == 2
            J_normoxia = compute_J(texp_normoxia, yexp_normoxia, model, pars_normoxia, init, simulation_time_normoxia, dt);
            J_hipoxia = compute_J(texp_hipoxia, yexp_hipoxia, model, pars_hipoxia, init, simulation_time_hipoxia, dt);
            
            Jp = J_normoxia + J_hipoxia;

            show_each_J(J_normoxia, 'INDIVIDUAL JN: ');
            show_each_J(J_hipoxia, 'INDIVIDUAL JH: ');
            J = sum(Jp);
            logJ = log(J);
            %J = 10*Jp(7)/numel(Jp) + J;

            real_total_time = toc(init_time);
            sim_total_time = simulation_time_normoxia + simulation_time_hipoxia;
            disp(sprintf('EVAL COMPLETED | REAL TIME: %.4f | RATIO: %.4f | JN: %.4f | JH: %.4f | J: %.4f', real_total_time, real_total_time / sim_total_time, sum(J_normoxia), sum(J_hipoxia) , J));
            show_parameter_variation(initial_point, optpars_iter, pars_hipoxia, idx_optpars, lb, ub);

            
            

          
        end      
    catch ME
        J = 10^10;
        logJ = log(J);
        disp(sprintf('EVAL INCOMPLETED | Error: %s | J: %.4f', ME.message, J));

        %disp(ME.message);
        %disp(ME.identifier);
        %disp(ME.stack(1));
    end

    function Jp = compute_J(texp, yexp, model, pars, init, simulation_time, dt )
        
        all_global = zeros(15, round(10 * simulation_time/dt) + 1) + 0;
        finapres_notnan_mask = yexp(:, end-1);
        on_exercise_mask = yexp(:, end);

        [t, x_dot, x_vars, not_x_keys, not_index] = run_ode_fun(model, pars, init, simulation_time, dt);  
        %% fit experimental and simulated vectors  
        %Simulation        
        [sim_vars, time_sim] = find_sim_vars(x_vars', x_keys', x_dot',t');       
        %Experimental
        [exp_vars, time_exp] = find_exp_vars(yexp, texp);
        %Adjust based on settling time
        [~, sim_vars_crop, ~, exp_vars_crop, finapres_notnan_mask_crop, on_exercise_mask_crop] = adjust_sizes(time_sim, time_exp, sim_vars, exp_vars, finapres_notnan_mask, on_exercise_mask);

        %% Residuals

        R_all = ((exp_vars_crop - sim_vars_crop).^2)./(eps + exp_vars_crop).^2; %the normalization is related to the experimental variable
        
        %exercise weight over rest weight
        ex_rest_weight_ratio = 5;  %5 is the maximum value, 2 is too low
        R_all = R_all .* on_exercise_mask_crop * ex_rest_weight_ratio + R_all .* (1 - on_exercise_mask_crop); 
        %remove "nan" values
        
        
        R_all(7,:)  = R_all(7,:)  .* finapres_notnan_mask_crop; 
        R_all(8,:)  = R_all(8,:)  .* finapres_notnan_mask_crop; 
        R_all(9,:)  = R_all(9,:)  .* finapres_notnan_mask_crop;
        R_all(10,:) = R_all(10,:) .* finapres_notnan_mask_crop;

        nan_values = sum(~finapres_notnan_mask_crop(:));
        [sz_vars, sz_time] = size(R_all);

        sum_only_for_plotting = sum(R_all, 1);
        plot_in_console(sum_only_for_plotting, 100, on_exercise_mask_crop); %just testing

        sum_without_nan_vars = sum(R_all(1:6,:),2)/sz_time;
        sum_with_nan_vars = sum(R_all(7:10,:),2)/(sz_time - nan_values + eps);
        sum_along_time = [sum_without_nan_vars; sum_with_nan_vars];


        Jp  = sqrt(sum_along_time);   
        Jp = Jp/sz_vars;
        
        %J = 10*Jp(7)/numel(Jp) + J;

        %Ignored data points
        %Jp(1) = 0;
        Jp(5) = 0;
        Jp(6) = 0;
        
    end

    function [texp, yexp, pars] = extract_input_args(pars, t_exp, y_exp, percentages, idx_optpars, optpars_iter)
        p_values = cell2mat(values(pars));
        p_keys = keys(pars);

        p_values(idx_optpars) = optpars_iter;   
        texp = t_exp;
        yexp = y_exp;      
        pars = estimate_newton_ohm(percentages, containers.Map(p_keys, p_values));
    end

    function [sim_vars, time_sim] = find_sim_vars(x_vars, x_keys, x_dot, t)

       PAO2_sim = x_vars(:, find(strcmp(x_keys, 'PAO2')));%
       PACO2_sim = x_vars(:, find(strcmp(x_keys, 'PACO2')));%
       dVE_sim = x_vars(:, find(strcmp(x_keys, 'dVE')));%
       V_sim = x_vars(:, find(strcmp(x_keys, 'V')));
       TI_sim = x_dot(:, find(strcmp(x_keys, 'fake_TI')));%
       Tresp_sim = x_dot(:, find(strcmp(x_keys, 'fake_Tresp')));%
       Theart_sim = x_vars(:, find(strcmp(x_keys, 'Theart')));
       P_sa_sim = x_vars(:, find(strcmp(x_keys, 'P_sa')));
       PM_sim = x_vars(:, find(strcmp(x_keys, 'mean_P_sa')));
       
       %Theart_sim = (Theart_sim <= 0).*(Theart_sim);
       HR_sim = Theart_sim.^1 * 60;
       VT_sim = data_processing('volume', V_sim, t);
       out_pressure = data_processing('pressure', P_sa_sim, t);
       PS_sim = out_pressure{2};
       PD_sim = out_pressure{3};

       sim_vars = vertcat(PAO2_sim',  PACO2_sim',  dVE_sim',  VT_sim',  TI_sim',  Tresp_sim',  HR_sim',  PM_sim',  PS_sim', PD_sim');
       

       time_sim = t;
    end
    function [exp_vars, time_exp] = find_exp_vars(yexp, texp)
       PAO2_exp = yexp(:,5);
       PACO2_exp = yexp(:,6);
       dVE_exp = yexp(:,1);
       VT_exp = yexp(:,2);
       TI_exp = yexp(:,3);
       Tresp_exp = yexp(:,4);
       HR_exp = yexp(:,7);    %this has to be multiplied by 60 in datapreprocessing
       PM_exp = yexp(:,10);
       PS_exp = yexp(:,8);
       PD_exp = yexp(:,9);

       HR_exp = HR_exp * 60;

       exp_vars = [PAO2_exp'; PACO2_exp'; dVE_exp'; VT_exp'; TI_exp'; Tresp_exp'; HR_exp'; PM_exp'; PS_exp'; PD_exp'];
       
       time_exp = texp - texp(1);
    end
    function [time_sim_crop, sim_vars_crop, time_exp_crop, exp_vars_crop, finapres_notnan_mask_crop, on_exercise_mask_crop] = adjust_sizes(time_sim, time_exp, sim_vars, exp_vars, finapres_notnan_mask, on_exercise_mask)
       try
           time_sim_crop = time_sim(time_sim > settling_time);
       sim_vars_crop = sim_vars(:, time_sim > settling_time);
       

       time_exp_crop = time_exp(time_exp < time_sim(end));
       exp_vars_crop = exp_vars(:, time_exp < time_sim(end));
       finapres_notnan_mask_crop = finapres_notnan_mask(time_exp < time_sim(end));
       on_exercise_mask_crop = on_exercise_mask(time_exp < time_sim(end));

       [time_exp_crop, indexes, ~] = unique(time_exp_crop, 'stable');
       exp_vars_crop = exp_vars_crop(:, indexes); 
       finapres_notnan_mask_crop = finapres_notnan_mask_crop(indexes);
       on_exercise_mask_crop = on_exercise_mask_crop(indexes);

       if size(sim_vars_crop, 2) >  size(exp_vars_crop, 2)
        exp_vars_crop = interp1(time_exp_crop, exp_vars_crop', time_sim_crop, 'linear', 'extrap'); 
        exp_vars_crop = exp_vars_crop';
        finapres_notnan_mask_crop = interp1(time_exp_crop, finapres_notnan_mask_crop', time_sim_crop, 'linear', 'extrap');
        finapres_notnan_mask_crop = finapres_notnan_mask_crop'; 
        on_exercise_mask_crop = interp1(time_exp_crop, on_exercise_mask_crop', time_sim_crop, 'linear', 'extrap');
        on_exercise_mask_crop = on_exercise_mask_crop'; 
        
       else
        sim_vars_crop = interp1(time_sim_crop, sim_vars_crop', time_exp_crop, 'linear', 'extrap'); 
        sim_vars_crop = sim_vars_crop';
       end
       catch
           disp('probably: interpolation error');
       end
    end

function show_parameter_variation(initial_point, optpars_iter, pars, idx_optpars, lb, ub)
    str = ['PARAMETERS:', sprintf('\n')];
    
    % Define column widths
    param_width = 15; % Width for the parameter name
    num_width = 8;    % Width for each numeric value
    sep_width = 2;    % Spaces between numbers (==)

    for i = 1:length(initial_point)
        value = optpars_iter(i);
        delta_percent = 100 * (value - initial_point(i)) / initial_point(i);
        param_name = keys(pars);
        param_name = param_name{idx_optpars(i)};
        upper_boundary = ub(idx_optpars(i)); 
        lower_boundary = lb(idx_optpars(i));
        
        % Build formatted number section
        numbers_part = sprintf('%*.2f%s==%*.2f%s==%*.2f', ...
            num_width, lower_boundary, repmat(' ', 1, sep_width), ...
            num_width, value, repmat(' ', 1, sep_width), ...
            num_width, upper_boundary);
        
        % Build whole line with fixed width until "DELTA"
        line = sprintf('%-*s %s   DELTA: %6.1f%%', ...
            param_width, param_name, ...
            numbers_part, ...
            delta_percent);
        
        str = [str, line, sprintf('\n')];
    end

    % Display
    disp(str);
end


    function show_each_J(J, msg)
        str = msg;
        xnames = {'PAO2',  'PACO2',  'dVE',  'VT',  'TI',  'Tresp',  'HR',  'PM',  'PS', 'PD'};
        for i = 1:length(xnames)           
            
            str = [str, sprintf('%s: %.4f | ', xnames{i}, J(i))];
        end

        % Remove last separator
        str = strtrim(str);
        if length(str) > 0 && str(end) == '|'
            str = str(1:end-2);
        end

        disp(str);
    end
   
   
  end