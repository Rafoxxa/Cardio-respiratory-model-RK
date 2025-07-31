function [t_exp, y_exp, VO2_poly, VCO2_poly, fiO2_poly, basal] = data_preprocessing(patient_idx, hipoxia_state, ascend_state, plot_state, VO2_type)
        
        

        if nargin < 4  % Check if 'n' is missing
            plot_state = 0;  % Default value
            VO2_type = "poly";
        elseif nargin < 5
            VO2_type = "poly";
        end
    
        testing = 0;
        %hipoxia = 1;
        %normoxia = 1-hipoxia;
        %ascend = 1;
        %exercise = 1-ascend;
        %patient_idx = 6;
         
        hipoxia = hipoxia_state == "hipoxia";
        normoxia = hipoxia_state == "normoxia";

        ascend = ascend_state == "ascend";
        exercise = ascend_state == "exercise";
        mix = ascend_state == "mix";
        
        

        
        if testing
            testing_stuff();
        end
        
        if hipoxia
            if ~mix
                [time_cpet, HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2, AT, ~, finapres_nan_intervals, on_exercise, ladder_points, Pow] = hipoxia_preprocessing(patient_idx, ascend, exercise);           
            else
                [time_cpet_ascend, HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2, AT, ~, finapres_nan_intervals_asc_, on_exercise_asc, ~, Pow_asc] = hipoxia_preprocessing(patient_idx, 1, 0);           
                ascend_ = [time_cpet_ascend, HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2];
                [time_cpet_exercise, HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2, AT,  ~, finapres_nan_intervals_ex_, on_exercise_ex, ladder_points, Pow_ex] = hipoxia_preprocessing(patient_idx, 0, 1);          
                exercise = [time_cpet_exercise + time_cpet_ascend(end), HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2];
                mixed = [ascend_; exercise];
                finapres_nan_intervals = [finapres_nan_intervals_asc_; finapres_nan_intervals_ex_];
                on_exercise = [on_exercise_asc; on_exercise_ex];
                Pow = [Pow_asc; Pow_ex];
                [time_cpet, HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2] = deal(mixed(:,1), mixed(:,2), mixed(:,3), mixed(:,4), mixed(:,5), mixed(:,6), mixed(:,7), mixed(:,8), mixed(:,9), mixed(:,10), mixed(:,11), mixed(:,12), mixed(:,13), mixed(:,14));
                

            end 

            
        
        elseif normoxia
            [time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR, PS, PD, PM, AT, finapres_nan_intervals, on_exercise, ladder_points, Pow] = normoxia_preprocessing(patient_idx);
        end

        
        

        
        %unit corrections
        dVE = dVE/60;
        HR = HR/60;

        
        t_exp = time_cpet;
        x_exp = [dVE, VT, TI, Tresp, PAO2, PACO2, HR, PS, PD, PM];

        x_names = {"dVE", "VT", "TI", "Tresp", "PAO2", "PACO2", "HR", "PS", "PD", "PM"};


        
        % window_size = 6;
        % x_exp(:,3) = movmin(x_exp(:,3), window_size, 'Endpoints', 'fill'); %TI

        %filtering
        fs = 1/(mean(diff(t_exp)));
        %cutoff = fs/1000;
        x_exp = filloutliers(x_exp,"linear", "movmedian",5);
        cutoff = fs/8; %10*fs/2;
        %cutoff = 10*fs/2;
        
        for i = 1:size(x_exp,2)
            x_exp(:, i) = lowpass(x_exp(:, i),cutoff, fs, "Steepness", 0.95);                     
        end

        y_exp = x_exp(4:end, :);
        t_exp = t_exp(4:end, :);

        %VO2 / VCO2 management: polynomial fitting of high order

        %VO2_poly = 1;
        %VCO2_poly = 1;

        if VO2_type == "poly"
            [~, VO2_poly, ~, ~, ~] = bestPolynomialFit(time_cpet, VO2/1000, 8, "VO2 (l/min)", plot_state);
            [~, VCO2_poly, ~, ~, ~] = bestPolynomialFit(time_cpet, VCO2/1000, 8, "VCO2 (l/min)" , plot_state);             
            [~, fiO2_poly, ~, ~, ~] = bestPolynomialFit(time_cpet, fiO2, 4, "fiO2 (%)", plot_state);
        elseif VO2_type == "ori"
            VO2_poly = VO2(1:end-3,:)/1000;
            VCO2_poly = VCO2(1:end-3,:)/1000;
            fiO2_poly = fiO2(1:end-3,:);

        end

        
        finapres_nan_mask = zeros(size(t_exp));
        for interval = finapres_nan_intervals'  
            if interval(2) == -1
                finapres_nan_mask = finapres_nan_mask + (t_exp > interval(1)) .* (t_exp <= t_exp(end));
            else
                finapres_nan_mask = finapres_nan_mask + (t_exp > interval(1)) .* (t_exp <= interval(2));
            end
            
        end
        finapres_notnan_mask = ~finapres_nan_mask;
        % Pow_ = 75*(Pow > 75) + Pow.*(Pow <= 75);
        % indexes = find(Pow~=0);
        % begin_on = indexes(1);
        % ending_on = indexes(end)
        % Pow__ = Pow;
        % Pow__(1:begin_on) = 1;
        % Pow__(begin_on:end) = 0;

        if VO2_type ~= "ori"
            y_exp = [y_exp, Pow(1:end-3, :), finapres_notnan_mask, on_exercise(4:end)];
        else
            y_exp = [y_exp, VO2_poly, VCO2_poly, Pow(1:end-3, :), finapres_notnan_mask, on_exercise(4:end)];
        end
        basal = [mean(abs(VO2(1:30))), mean(abs(VCO2(1:30)))]/1000;
        basal = [basal, AT/1000, ladder_points(:,1)', ladder_points(:,2)'];

        

         
        
        %units_table = readtable("variables_units.xlsx");
        %plot_vars(t_exp, y_exp, x_names, units_table, 5, 2);
        %plot_vars(t_exp, [VO2, VCO2], {"VO2", "VCO2"}, 2, 1);

        
    function [time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR, PS, PD, PM, AT, finapres_nan_intervals, on_exercise, ladder_points, Pow] = normoxia_preprocessing(patient_idx)
        switch patient_idx
            case 1
                finapres_nan_intervals = [0, 846; 2690, -1];
            case 2
                finapres_nan_intervals = [];
            case 3
                finapres_nan_intervals = [];
            case 4
                finapres_nan_intervals = [];
            case 5
                finapres_nan_intervals = [0, 350; 1846, -1];
            case 6    
                finapres_nan_intervals = [2164,2240, 2310, -1];

        end
        %finapres data is 800 seconds over cpet's
        finapres_delay = [800, 0, 0, 0, 0, 0]; %[0, 0, 0, 0, 0, 0];
        cpet_delay = [0,500,100, 0, 0, 0];
    
        filename_finapres = ["../../data/NORMOXIA/AG/FINAPRES/08.30.2024_10.30.07.nsc", "../../data/NORMOXIA/RK/FINAPRES/09.02.2024_09.53.50_1.nsc", "../../data/NORMOXIA/AM/FINAPRES/09.13.2024_08.53.11_1.nsc", "../../data/NORMOXIA/CV/FINAPRES/01.31.2025_13.03.17.nsc", "../../data/NORMOXIA/VC/FINAPRES/01.17.2025_15.32.39.nsc", "../../data/NORMOXIA/EC/FINAPRES/04.08.2025_08.42.03.nsc"];
        
        [time_finapres, HR, PS, PD, PM] = read_finapress(char(filename_finapres(patient_idx)), cpet_delay(patient_idx));

        

        filename_cpet = ["../../data/NORMOXIA/AG/CPET/CPET-1_AG_cleaned.xlsx", "../../data/NORMOXIA/RK/CPET/CPET-2_RK_cleaned.xlsx", "../../data/NORMOXIA/AM/CPET/CPET-3_AM_cleaned.xlsx", "../../data/NORMOXIA/CV/CPET/CPET-4_CV_cleaned.xlsx", "../../data/NORMOXIA/VC/CPET/CPET-5_VC_cleaned.xlsx", "../../data/NORMOXIA/EC/CPET/CPET-6_EC_cleaned.xlsx"];
        [time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_CPET, AT, on_exercise, ladder_points, Pow] = read_cpet(filename_cpet(patient_idx),hipoxia, ascend); 
        
        
       

        
        [~, HR] = finapres_merge_cpet(time_cpet, dVE, time_finapres', HR', finapres_delay(patient_idx));
        [~, PS] = finapres_merge_cpet(time_cpet, dVE, time_finapres', PS', finapres_delay(patient_idx));
        [~, PD] = finapres_merge_cpet(time_cpet, dVE, time_finapres', PD', finapres_delay(patient_idx));
        [~, PM] = finapres_merge_cpet(time_cpet, dVE, time_finapres', PM', finapres_delay(patient_idx));

        if HR_CPET ~= 0
            HR = HR_CPET;
        end
    end
        
    function [time_cpet, HR, PS, PD, PM, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2, AT,  HR_CPET, finapres_nan_intervals,on_exercise, ladder_points, Pow] = hipoxia_preprocessing(patient_idx, ascend, exercise)
        
            switch patient_idx
            case 4
                subject_filenames = {"../../data/HIPOXIA/CV/12.12.2024_09.36.27.nsc", "../../data/HIPOXIA/CV/12.12.2024_10.22.26.nsc", "../../data/HIPOXIA/CV/12.12.2024_10.35.08.nsc"};
                split_instructions = {...
                 {'split1', [70, 2132-70]}, ...
                 {'split2', [176]},...
                 {'split3', [1]}
                 };
                merge_instructions = {
                 {'split1_2'},
                 {'split2_2', 82, 'split3_2'}...
                 };
                finapres_delay_ascend = -70;
                finapres_delay_exercise = -150;
                finapres_nan_intervals_asc = [1900, 2180];
                finapres_nan_intervals_ex = [];
                ascend_cpet = "../../data/HIPOXIA/CV/ASCENSO/CPET_HIPOXIA_ASCENSO_CV.xlsx";
                exercise_cpet = ["../../data/HIPOXIA/CV/NEW_CPET_HIPOXIA_EJERCICIO_CV.xlsx"];

            case 1    
                subject_filenames = {"../../data/HIPOXIA/AG/12.21.2024_09.27.47_1.nsc", "../../data/HIPOXIA/AG/12.21.2024_10.38.12.nsc"};
                split_instructions = {...
                 {'split1', [2500]}, ...
                 {'split2', [1]},...                 
                 };
                merge_instructions = {
                 {'split1_1'}, ...
                 {'split1_2', 10, 'split2_2'}...
                 };
                finapres_delay_ascend = 0;
                finapres_delay_exercise = -2500;
                finapres_nan_intervals_asc = [];
                finapres_nan_intervals_ex = [];
                ascend_cpet = "../../data/HIPOXIA/AG/ASCENSO/CPET_HIPOXIA_ASCENSO_AG.xlsx";
                exercise_cpet = "../../data/HIPOXIA/AG/NEW_CPET_HIPOXIA_EJERCICIO_AG.xlsx";


            case 5
                subject_filenames = {"../../data/HIPOXIA/VC/12.20.2024_10.12.27.nsc", "../../data/HIPOXIA/VC/12.20.2024_11.15.58.nsc"};
                split_instructions = {...
                 {'split1', [1]}, ...
                 {'split2', [1]}                     
                 };
                merge_instructions = {
                 {'split1_2'},
                 {'split2_2'}...
                 };
                finapres_delay_ascend = -20;
                finapres_delay_exercise = -15;
                finapres_nan_intervals_asc = [2170, 2303];
                finapres_nan_intervals_ex = [];
                ascend_cpet = "../../data/HIPOXIA/VC/ASCENSO/CPET_HIPOXIA_ASCENSO_VC.xlsx";
                exercise_cpet = ["../../data/HIPOXIA/VC/NEW_CPET_HIPOXIA_EJERCICIO_VC.xlsx"];
            
            case 6
                subject_filenames = {"../../data/HIPOXIA/EC/12.11.2024_13.14.50.nsc", "../../data/HIPOXIA/EC/12.11.2024_14.13.45.nsc", "../../data/HIPOXIA/EC/12.11.2024_14.32.55.nsc",  "../../data/HIPOXIA/EC/12.11.2024_14.53.33.nsc"};
                split_instructions = {...
                 {'split1', [2040]}, ...
                 {'split2', [1]},...
                 {'split3', [1]},...
                 {'split4', [1]}
                 };
                merge_instructions = {
                 {'split1_2', 240, 'split2_2'},
                 {'split3_2', 400, 'split4_2'}...
                 };
                finapres_delay_ascend = -2040;
                finapres_delay_exercise = 0;
                finapres_nan_intervals_asc = [2360, 2530];
                finapres_nan_intervals_ex = [];
                ascend_cpet = "../../data/HIPOXIA/EC/ASCENSO/CPET_HIPOXIA_ASCENSO_EC.xlsx";
                exercise_cpet = ["../../data/HIPOXIA/EC/NEW_CPET_HIPOXIA_EJERCICIO_EC.xlsx"];
            end
                                
            subject_vars = cell(1, length(subject_filenames)); % Preallocate cell array
            
            for filename_index = 1:length(subject_filenames)
                % Read the data from each file
                [time_finapres, HR, PS, PD, PM] = read_finapress(char(subject_filenames{filename_index}), 0);
            
                % Concatenate all variables into a single matrix (rows = observations, columns = variables)
                finapres_vars = horzcat(time_finapres, HR, PS, PD, PM);
            
                % Store the concatenated matrix in the cell array
                subject_vars{filename_index} = finapres_vars'; % Directly assign to cell
            end

         
            vars_list = subject_vars;  
            struct_out = finapres_split_merge(vars_list, split_instructions, merge_instructions);
            ascend_finapres = struct_out.merge1;
            exercise_finapres = struct_out.merge2;

            ascend_finapres = interpolate_nan(ascend_finapres);
            exercise_finapres = interpolate_nan(exercise_finapres);


            if ascend 
                
                
                time_finapres = ascend_finapres(1,:);
                HR = ascend_finapres(2,:);
                PS = ascend_finapres(3,:);
                PD = ascend_finapres(4,:);
                PM = ascend_finapres(5,:);

                
                
                [time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_CPET, AT, on_exercise, ladder_points, Pow] = read_cpet(ascend_cpet(1), hipoxia, ascend);            

                finapres_delay = finapres_delay_ascend;
                finapres_nan_intervals = finapres_nan_intervals_asc; 
                

            elseif exercise 

                time_finapres = exercise_finapres(1,:);
                HR = exercise_finapres(2,:);
                PS = exercise_finapres(3,:);
                PD = exercise_finapres(4,:);
                PM = exercise_finapres(5,:);

                

                [time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2,fiO2, HR_CPET, AT, on_exercise, ladder_points, Pow] = read_cpet(exercise_cpet(1), hipoxia, ascend); 

                finapres_delay = finapres_delay_exercise;
                finapres_nan_intervals = finapres_nan_intervals_ex;
            end
        
        
            [time_finapres, dependant] = remove_duplicates(time_finapres, {HR, PS, PD, PM});
            HR = dependant{1};
            PS = dependant{2};
            PD = dependant{3};
            PM = dependant{4};
            
            
            [~, HR] = finapres_merge_cpet(time_cpet, dVE, time_finapres, HR, finapres_delay);
            [~, PS] = finapres_merge_cpet(time_cpet, dVE, time_finapres, PS, finapres_delay);
            [~, PD] = finapres_merge_cpet(time_cpet, dVE, time_finapres, PD, finapres_delay);
            [~, PM] = finapres_merge_cpet(time_cpet, dVE, time_finapres, PM, finapres_delay);

            if HR_CPET ~= 0
                HR = HR_CPET;
            end

            

    end

        
    function [output, dependant_] = remove_duplicates(input, dependant)
    [~, unique_indices] = unique(input, 'first');
    sorted_indices = sort(unique_indices); 

    output = input(sorted_indices);          

    dependant_ = cell(1, numel(dependant));  
    
    for idx = 1:numel(dependant)
        dependant_var = dependant{idx};
        dependant_{idx} = dependant_var(sorted_indices);
    end
end


    function data_interpolated = interpolate_nan(data)
        % Initialize the output matrix
        data_interpolated = data;
        
        % Loop over each row of the matrix
        for row = 1:size(data, 1)
            % Find indices of NaN values in the row
            nan_indices = isnan(data(row, :));
            
            % Interpolate only the NaN values
            data_interpolated(row, nan_indices) = interp1(find(~nan_indices), data(row, ~nan_indices), ...
                find(nan_indices), 'linear', 'extrap');
        end
    end
    
        
    function [time_finapres, HR, PS, PD, PM] = read_finapress(filename, cpet_delay)
    
    [output] = readBeatscopeData20(filename, {'all'});
    HR_struct = output.HR_AP;
    PS_struct = output.reSYS;
    PD_struct = output.reDIA;
    PM_struct = output.reMAP;

    time_finapres_ = PS_struct.time - cpet_delay;    %finapres takes a lot of measurements with different sample rates
    time_finapres = time_finapres_(time_finapres_ >= 0);
    PS = PS_struct.values;
    PS = PS(time_finapres_ >= 0);
    PD = PD_struct.values;
    PD = PD(time_finapres_ >= 0);
    PM = PM_struct.values;
    PM = PM(time_finapres_ >= 0);

    PS = remove_outliers_with_interpolation(PS, 2);
    PD = remove_outliers_with_interpolation(PD, 2);
    PM = remove_outliers_with_interpolation(PM, 2);
    
    HR_values_resized = interp1(HR_struct.time, HR_struct.values, PS_struct.time, 'linear', 'extrap');
    HR = HR_values_resized(time_finapres_ >= 0);

    
        

    end

    function [t_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_CPET, AT, on_exercise, ladder_points, Pow] = read_cpet(filename_cpet, hipoxia,  ascend)
        
        

        X_cpet = readtable(filename_cpet);

        %time
        days2seconds = @(x) x*24*60*60; 
        t_cpet = days2seconds(X_cpet.t); 
        
        %direct vars
        VO2 = X_cpet.VO2;
        VCO2 = X_cpet.VCO2;
        


        
        VT = X_cpet.VT;
        dVE = X_cpet.VE;
        TI = X_cpet.Ti;
        Tresp = X_cpet.Ttot;

        fiO2 = X_cpet.FiO2/100;
        FeO2 = X_cpet.FeO2/100;
        FeCO2 = X_cpet.FeCO2/100;
  

        if hipoxia            
 
            final_value_fio2 = 12.3;  %this will depend on the experiment
            step_fio2 = 0.5;  %each minute
            step_duration = 2; %in minutes
            ascend_duration = 20; %in minutes
            inital_value_fio2 = final_value_fio2 + step_fio2 * ascend_duration/step_duration;
            if ascend
                fiO2 = vo2_correction(t_cpet, inital_value_fio2,  final_value_fio2, step_fio2, step_duration * 60);
                fiO2 = (fiO2)/100;
                %VO2 = 1000 * ((dVE .* fiO2 .* (1 - (FeO2 + FeCO2)))./(1.2 - 2*fiO2) - (dVE .* FeO2));  %Manual correction, derived from OMNIA outputs
                %VO2 = 1000 * ((dVE .* fiO2 .* (1 - (FeO2 + FeCO2)))./(1 - FeO2 + 0.3*FeCO2) - (dVE .* FeO2));  %Manual correction, derived from OMNIA outputs
                %VO2 = 1000 * 0.826 * (dVE .* fiO2 - (dVE .* FeO2));
                %VO2 = 1000 * 0.7694 * (dVE .* fiO2 - (dVE .* FeO2));
                STPD = X_cpet.STPD(1);
                VO2 = 1000 * STPD * (dVE .* fiO2 .*(1 - FeO2 - FeCO2)./(1 - fiO2) - dVE .* FeO2);
                
                %Alternatives to factor: 0.79, (1 - fiO2), (1.2 - 2*fiO2)
            else
                fiO2 = 12.3 + 0*(fiO2);
                fiO2 = fiO2/100;
                VO2 = X_cpet.VO2;   %using NEW_CPET files, which have the correction from OMNIA
            end            
            
            
        else
            fiO2 = X_cpet.FiO2/100;
            VO2 = X_cpet.VO2;

            
            
        end
        VCO2 = X_cpet.VCO2;
        
        %non direct vars
        PETCO2 = X_cpet.PetCO2;        
        PACO2 = 5.5 + 0.9 * PETCO2 - 0.0021*VT; 
        Patm = 719; 
         
        PAO2 = fiO2 .* (Patm - 47) - PACO2 .* VO2./VCO2;  

        %HR
        
        try
          HR_CPET = X_cpet.HR;
        catch 
          HR_CPET = 0;  
        end
        
        
        %filters
        %PAO2 = lowpass(PAO2,0.05 , 1/(t_cpet(2) - t_cpet(1)));         
        %PACO2 = lowpass(PACO2, 0.05, 1/(t_cpet(2) - t_cpet(1)));

        %VO2 and VCO2 values at the end of each step
        Pow = X_cpet.Potencia;
        on_exercise = Pow > 0;
        pre_selector = diff([Pow(1);Pow]) ~= 0;
        selector = find(pre_selector == 1);
        VCO2_final_points = VCO2(selector);
        VO2_final_points = VO2(selector);
        ladder_points = [VO2_final_points, VCO2_final_points];

        %Anaerobic thresholds based on VCO2 value at VT1
        try
            AT = X_cpet.VCO2_at_VT1(1);
        catch
            AT = 1000;
        end

  


        
        

    
    
    end

    function o2_corrected = vo2_correction(time, o2_start, o2_final, step_size, step_interval)
        % Generates an O2 vector with ascending values based on time and parameters.
        
        
        % Calculate the number of steps for each time point
        steps = floor(time / step_interval);
        
        % Generate the corrected O2 vector
        o2_corrected = o2_start - steps * step_size;
        o2_corrected = (o2_corrected >= o2_final) .* (o2_corrected) + (o2_corrected < o2_final) * o2_final;

    end


    function [t_merged, xfin_merged] = finapres_merge_cpet(t_cpet, x_cpet, t_fin, x_fin, time_delay)
        
        
        shifted_common_time = t_fin + time_delay;
        unified_time = t_cpet;

        x_fin_unified = zeros(size(unified_time));

        overlap_start = find(unified_time >= shifted_common_time(1), 1, 'first');
        overlap_end = find(unified_time <= shifted_common_time(end), 1, 'last');

        x_fin_unified(overlap_start: overlap_end) = interp1(shifted_common_time, x_fin, unified_time(overlap_start:overlap_end), 'linear', 0);

        %xcpet_merged = interp1(t_cept, x_cpet, unifed_time, 'linear', 'extrap');
        
        xfin_merged = x_fin_unified;
        t_merged = unified_time;
        



    end 

    function plot_vars(t, X, Xnames, Xunits, width, height)
        t_exp = t;
        x_exp = X;
        variable_names = Xnames;
        units_table = Xunits;
        % Create the figure and 5x2 grid of subplots
        %figure;
        
        % Loop through each variable
        for i = 1:size(x_exp, 2)
            subplot(width, height, i); % Define the subplot position
            hold on;
            plot(t_exp, x_exp(:, i)); % Plot the variable against time
            title(variable_names{i}); % Set the title to the variable name
            xlabel('Time (s)'); % Label the x-axis
            unit = find_unit(units_table, variable_names{i});
            yl = strcat(variable_names{i},"(", unit, ")");
            ylabel(yl); % Label the y-axis with the variable name
            grid on; % Add grid for better visualization
        end
        
        % Adjust the layout to prevent overlap
        sgtitle('Time Series Data'); % Add a super title for the entire figure 
    end

   function [output, merge_times] = finapres_split_merge(vars_list, split_instructions, merge_instructions)
    % Output structure to hold split and merged variables
    output = struct();
    merge_times = struct(); % Structure to hold the times for each merge

    % Step 1: Split Variables (handles one or multiple splits)
    splits = struct();
    for i = 1:length(vars_list)
        vars = vars_list{i}; % Current set of variables (matrix form)
        time = vars(1, :);  % First row is the time vector

        % Check for split instruction for this set of variables
        file_splits = split_instructions{i}; % Example: {'split1', [234, 543]}
        split_name = file_splits{1};          % Base name for the splits
        split_times = file_splits{2};         % Times at which to split

        % Define split points (first point is the start time, last point is the end time)
        split_points = [time(1), split_times, time(end)];

        % Split the data into chunks
        for j = 1:(length(split_points) - 1)
            % Define split indices
            indices.start = find(time >= split_points(j), 1, 'first');
            indices.end = find(time <= split_points(j + 1), 1, 'last');

            % Extract splits and store them
            chunk_name = sprintf('%s_%d', split_name, j);
            splits.(chunk_name) = vars(:, indices.start:indices.end);
        end
    end

    % Step 2: Merge Variables (preserving the sizes and creating merged vectors)
    for i = 1:length(merge_instructions)
        % Example: {'split1_2'} or {'split2_2', 82, 'split3_2'}
        merge_instr = merge_instructions{i};
        merged_data = [];  % Initialize matrix to store merged variables
        merged_times = []; % Initialize vector to store merged times

        last_time = []; % Variable to store the last time value of the merged data

        for j = 1:2:length(merge_instr)
            name = merge_instr{j};
            padding = 0; % Default no padding

            % Check if there is padding specified
            if j + 1 <= length(merge_instr) && isnumeric(merge_instr{j + 1})
                padding = merge_instr{j + 1};
            end

            % Check if the split data is available
            if isfield(splits, name)
                % Append the current split data
                split_data = splits.(name);
                split_time = split_data(1, :); % Time vector of current split data

                % If it's not the first chunk, adjust the time of the next chunk
                if ~isempty(last_time)
                    % Adjust the time by adding padding and the last time value
                    
                    disp(last_time(end))
                    disp(padding)
                    disp(split_time(1))
                    split_time = split_time + (last_time(end) + padding - split_time(1)) + 1;

                end

                % Append the current time and data
                merged_times = [merged_times, split_time];
                merged_data = [merged_data, split_data(2:end, :)];  % Only append data rows (exclude time row)

                % Update the last time for the next chunk
                last_time = split_time; 

                % If padding is required, insert padding with NaN values
                if padding > 0
                    % Calculate the padding time range
                    pad_time = (last_time(end) + 1):(last_time(end) + padding); % Padding time

                    % Ensure same number of rows in padded data (matching original var size)
                    pad_matrix = nan(size(split_data, 1) - 1, length(pad_time)); % Create NaN padding (skip time row)
                    merged_data = [merged_data, pad_matrix];  % Append padding to the data
                    merged_times = [merged_times, pad_time];   % Append padding time
                end
            else
                warning('Split "%s" not found in splits structure.', name);
            end
        end

        % Store merged variables in the output structure with names like finapres_merge_1, finapres_merge_2, etc.
        merge_name = sprintf('merge%d', i); % Create dynamic name like finapres_merge_1, finapres_merge_2, etc.
        
        % Store the merged time and data
        output.(merge_name) = [merged_times; merged_data];
        merge_times.(merge_name) = merged_times;
    end
end

function clean_data = remove_outliers_with_interpolation(data, threshold, method)
    if nargin < 2 || isempty(threshold)
        threshold = 3;
    end
    if nargin < 3 || isempty(method)
        method = 'linear';
    end
    mu = mean(data, 'omitnan');
    sigma = std(data, 'omitnan');
    Z = (data - mu) / sigma;
    data_with_nans = data;
    data_with_nans(abs(Z) > threshold) = NaN;
    clean_data = fillmissing(data_with_nans, method);
end











    


    function unit = find_unit(table, var)
        row_index = strcmp(table.Variable, var);
        unit = table.MeasureUnit{row_index};    
    end

   function testing_stuff()
    disp("pass");
    % You can cross the data using HR as a common variable, with
            % that you can find the splitting times and merging portions

            % filename = "../../data/HIPOXIA/VC/12.20.2024_11.13.26.nsc";
            % [time_finapres, HR, PS, PD, PM] = read_finapress(char(filename), 0);
            % plot(time_finapres, HR);
            % title("1")
            % figure;
            % 
            % 
            %filename = "../../data/HIPOXIA/EC/12.11.2024_13.14.50.nsc";
            %[time_finapres, HR, PS, PD, PM] = read_finapress(char(filename), 0);
            %plot(time_finapres, HR);
            %title("2")
            %figure;

            %filename = "../../data/HIPOXIA/EC/12.11.2024_14.13.45.nsc";
            %[time_finapres, HR, PS, PD, PM] = read_finapress(char(filename), 0);
            %plot(time_finapres, HR);
            %title("3")
            %figure;
            %% 
            %% 
            %filename = "../../data/HIPOXIA/EC/ASCENSO/CPET_HIPOXIA_ASCENSO_EC.xlsx";
            %[time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_CPET] = read_cpet(filename, hipoxia, ascend); 
            %plot(time_cpet, HR_CPET)
            %title("Ascenso");
            %figure;
            % 
            % filename = "../../data/HIPOXIA/EC/CPET_HIPOXIA_EJERCICIO_EC.xlsx";
            % [time_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, HR_CPET] = read_cpet(filename); 
            % plot(time_cpet, HR_CPET)
            % title("Exercise");
   end   

end

