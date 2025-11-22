function [t_exp, y_exp, VO2_poly, VCO2_poly, fiO2_poly, TI_poly, Tresp_poly, basal] = data_preprocessing2(patient_idx, hipoxia_state, ascend_state, plot_state, VO2_type)
% %%Generating excel files for finapres
% % Suponiendo que tu estructura principal se llama "data"
%   % obtiene los nombres de los substructs% 
if nargin < 4  % Check if 'n' is missing
    plot_state = 0;  % Default value
    VO2_type = "poly";
elseif nargin < 5
    VO2_type = "poly";
end

if strcmp(hipoxia_state, 'hipoxia')
    hipoxia_value = 1;
elseif strcmp(hipoxia_state, 'normoxia')
    hipoxia_value = 0;
end

general_structures.hipoxia = hipoxia_value;
general_structures.normoxia = 1 - general_structures.hipoxia;
general_structures = setting_up_general_structures(general_structures);


%for i = 1:4
    %if i > 0 
        %dsd = [1 1 2 2 3 3 4 4];
        %CONDICIONES INICIALES
        ies = [1, 0, 0, 2, 3, 4];
        i = ies(patient_idx);
        patient_number = [1 2 3 4];       
        idx = patient_number(i);
        subject_number = [1,4,5,6];
        
        general_structures = sel_general_structures(idx, general_structures);
        ga = general_structures;

        
        %ITERACION SOBRE LOS FINAPRES
        for i = 1:size(ga.finapres_filenames,2)
            [output] = readBeatscopeData20(char(ga.finapres_filenames(i)), {'all'});
            data = output;
            fields = fieldnames(data);            
        
            HR_struct = output.HR_AP;
            PS_struct = output.reSYS;
            PD_struct = output.reDIA;
            PM_struct = output.reMAP;
            time_finapres = PS_struct.time;
            PM = PM_struct.values;
            PD = PD_struct.values;
            PS = PS_struct.values;
        
            HR = interp1(HR_struct.time, HR_struct.values, PS_struct.time, 'linear', 'extrap');
            f_all = [HR(4:end), PM(4:end), PS(4:end), PD(4:end)];

       
            
            %ga.HR_finapres_all = [ga.HR_finapres_all; HR(4:end)];
            %ga.PM = [ga.PM; PM(4:end)];   
            ga.finapres_all = [ga.finapres_all;f_all];
            
            ga.time_finapres_all  = [ga.time_finapres_all; time_finapres(4:end) + ga.finapres_inits(i)];
            
            
        end

        general_structures = finapres_specifcs(idx, ga);
        %figure;
        %plot(general_structures.time_finapres_all, general_structures.finapres_all(:, 1));
        
        %figure;
        %tcpet_antiguo = 0;
        %ITERACION SOBRE LOS CPETS
      
        for i = 1:size(general_structures.cpet_filenames,2)   

            [t_cpet, VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_cpet, AT, on_exercise, ladder_points, Pow] = read_cpet(general_structures.cpet_filenames(i), general_structures.hipoxia, 2-i);
            if size(HR_cpet, 1) == 1   %for subject 1, there is no cpet HR.
                HR_cpet = VT * 0;
            end
            
            %plot(t_cpet + tcpet_antiguo, HR_cpet);
            %hold on;
            %tcpet_antiguo = t_cpet(end);
            
    

    
            X_cpet_.all = [VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_cpet, on_exercise, Pow];

            X_cpet_.AT = AT;
            X_cpet_.ladder_points = ladder_points;

            
            %ladder_points is pendant;
            

            data_structures.X = X_cpet_;
            data_structures.t = t_cpet;
            general_structures = cpet_specifcs(idx, i, general_structures, data_structures);
            
            
        end    
        general_structures = finapres2cpet_corrections(idx, general_structures);

        %plotting_data(idx, general_structures, subject_number);

        general_structures = align_and_resample(idx, general_structures);

    

        [t_exp, y_exp, VO2_poly, VCO2_poly, fiO2_poly, TI_poly, Tresp_poly, basal] = polish_and_filter(idx, general_structures, plot_state, VO2_type);
        

        


        %filters

        %saving
    

        

function [general_structures] =  cpet_specifcs(idx, iteration, general_structures, data_structures)

            i = iteration;
            %HR_cpet_all = general_structures.HR_cpet_all;
            %Power_all = general_structures.Power_all;
            cpet_all = general_structures.cpet_all;
            AT_all = general_structures.AT_all;            
            cpet_inits = general_structures.cpet_inits;
            finapres_inits = general_structures.finapres_inits;
            hipoxia = general_structures.hipoxia;

            time_cpet_all = general_structures.time_cpet_all;
            last_time_cpet = general_structures.last_time_cpet;

            X_cpet = data_structures.X;
            t_cpet = data_structures.t;


           
            %Variable Managment
            % Defaults
            %HR_cpet_all_new = [HR_cpet_all; X_cpet.HR];    
            %Power_all_new = [Power_all; X_cpet.Potencia];
            cpet_all_new = [cpet_all; X_cpet.all];
            AT_all_new = [AT_all; X_cpet.AT];

            switch idx
                case 3  
                    if hipoxia 
                        if i == 2
                            patient5_delay = 151;
                            mini_mask = (t_cpet + cpet_inits(i)) > (finapres_inits(4) + patient5_delay);
                            %mini_mask_big = repmat(mini_mask', size(X_cpet.all, 2));
                            %HR_cpet_all_new = [HR_cpet_all; X_cpet.HR(mini_mask)];
                            %Power_all_new = [Power_all; X_cpet.Potencia(mini_mask)];   
                            cpet_all_new = [cpet_all; X_cpet.all(mini_mask,:)];
                        end
                    end

                case  1
                    if hipoxia == 0
                        disp('in normoxia, s1 does not have HR');
                        cpet_all_new(:, 10) = 0;
                        
                    end   

            end

            %Time managment

            switch hipoxia 

                case 1
                    if i == 1
                        time_cpet_all = [time_cpet_all; t_cpet + cpet_inits(i)];
                        last_time_cpet = t_cpet(end) + last_time_cpet + cpet_inits(i);
                    elseif i == 2
                        if idx == 3
                            t_cpet = t_cpet(mini_mask) - patient5_delay;   
                        end
                        time_cpet_all = [time_cpet_all; t_cpet + last_time_cpet];
                        time_interval_to_erase = [last_time_cpet, cpet_inits(i)];
                    end


                case 0
                    time_cpet_all = [time_cpet_all; t_cpet + cpet_inits(i)];
                    
            end
        general_structures.cpet_all = cpet_all_new;
        %general_structures.Power_all = Power_all_new;
        general_structures.time_cpet_all = time_cpet_all;      
        general_structures.AT_all = AT_all_new;
        
        if hipoxia
            general_structures.last_time_cpet = last_time_cpet;
            if i == 2 
                general_structures.time_interval_to_erase = time_interval_to_erase;
                general_structures.ladder_points = X_cpet.ladder_points;
            end
        else
            general_structures.ladder_points = X_cpet.ladder_points;
        end
end

function general_structures = finapres_specifcs(idx, general_structures)
        
        %PM_all = general_structures.PM;
        finapres_all = general_structures.finapres_all;

        if general_structures.hipoxia == 1

            if (idx == 1) 
                sample_rate = 2;
                minute2seconds = 60;            
                final = 1.5;
                %PM_all_ = PM_all;
                %PM_all = mean(PM_all(2:final*minute2seconds*sample_rate)) * ones(size(PM_all));
                %finapres_all = mean(finapres_all(2:final*minute2seconds*sample_rate,:)) .* ones(size(finapres_all));
                finapres_all_new = repmat(mean(finapres_all(2:final*minute2seconds*sample_rate,:)),size(finapres_all, 1) , 1);
                finapres_all_new(:, 1) = finapres_all(:, 1);  %the transformation above is not for HR (first variable)
                finapres_all = finapres_all_new; 
            end
    
            if (idx == 0) %this was for subject 4 
                sample_rate = 2;
                minute2seconds = 60;  
                init = 4;
                final = 10;
                %PM_all_ = PM_all;
                %finapres_all = mean(finapres_all(init * minute2seconds * sample_rate:final*minute2seconds*sample_rate)) * ones(size(finapres_all));            
                %PM_all = PM_all_;
                finapres_all_new = repmat(mean(finapres_all(init * minute2seconds * sample_rate:final*minute2seconds*sample_rate,:)), size(finapres_all, 1), 1);
                finapres_all_new(:, 1) = finapres_all(:, 1);  %the transformation above is not for HR (first variable)
                finapres_all = finapres_all_new; 
            end
    
            general_structures.finapres_all = finapres_all;
        end
end

function general_structures = sel_general_structures(idx, general_structures)
        
        general_structures.finapres_filenames = general_structures.finapres_filenames_cell{idx};
        general_structures.cpet_filenames = general_structures.cpet_filenames_cell{idx};
        general_structures.finapres_inits = general_structures.finapres_inits_cell{idx};        
        general_structures.cpet_inits = general_structures.cpet_inits_cell{idx};       
        
        general_structures.last_time = 0;
        general_structures.last_time_cpet = 0;
        general_structures.HR_finapres_all = [];
        general_structures.HR_cpet_all = [];
        general_structures.PM = [];
        general_structures.PM_all_ = [];
        general_structures.Power_all = [];
        general_structures.VO2_all = [];
        
        general_structures.time_finapres_all = [];
        general_structures.time_cpet_all = [];
        general_structures.cpet_all = [];
        general_structures.AT_all = [];
        general_structures.finapres_all = [];
        general_structures.ladder_points = [];

        general_structures.finapres_notnan_interval = general_structures.finapres_notnan_intervals{idx};  

end

function general_structures = finapres2cpet_corrections(idx, general_structures)
        ga = general_structures;
        time_cpet_all = ga.time_cpet_all - ga.cpet_inits(1);
        time_finapres_all = ga.time_finapres_all - ga.cpet_inits(1);
        hipoxia = ga.hipoxia;
        
        cpet_all = ga.cpet_all;
        finapres_all = ga.finapres_all;
        %HR_finapres_all = ga.HR_finapres_all;
        %PM_all = ga.PM;
        
        if hipoxia == 0
            %if idx == 1
            %    HR_cpet_all = zeros(size(time_cpet_all));    
            %end
        else
            time_interval_to_erase = ga.time_interval_to_erase - ga.cpet_inits(1); 
            if (idx == 1 || idx == 2)
                finapres_mask1 = (time_finapres_all < time_interval_to_erase(1));
                finapres_mask1_big = repmat(finapres_mask1,1, size(finapres_all, 2));
                finapres_all = finapres_all .* finapres_mask1_big;
                %PM_all = PM_all .* (time_finapres_all < time_interval_to_erase(1));
            end
            finapres_mask2 = (time_finapres_all < time_interval_to_erase(1)) + (time_finapres_all >  time_interval_to_erase(end));
            finapres_mask2 = logical(finapres_mask2);
            
            time_finapres_all = time_finapres_all(finapres_mask2);
            time_finapres_all( time_finapres_all > time_interval_to_erase(end)) =  time_finapres_all( time_finapres_all > time_interval_to_erase(end)) - (time_interval_to_erase(end) - time_interval_to_erase(1)); 
            
            finapres_all = finapres_all(finapres_mask2, :);
            %HR_finapres_all = HR_finapres_all(mask);
            %PM_all = PM_all(mask);

        end            
        
        
        if ga.finapres_crazy_intervals{idx} ~= 0
            finapres_all_new = [];
            for finapres_var_index = 1:size(finapres_all, 2)
                [tf, var] = fill_dis(time_finapres_all, finapres_all(:, finapres_var_index), 2, ga.finapres_crazy_intervals{idx} + time_finapres_all(1));
                finapres_all_new = [finapres_all_new, var];
            end
            


            %[tf, PM_all] = fill_dis(time_finapres_all, PM_all, 2, ga.finapres_crazy_intervals{idx} + time_finapres_all(1));

        else
            finapres_all_new = finapres_all;
            tf = time_finapres_all;
        end

        time_finapres_all = tf;
        time_cpet_all = time_cpet_all + ga.corrections_delay(idx);

        ga.cpet_all = cpet_all;
        %ga.HR_finapres_all = HR_finapres_all;
        ga.finapres_all = finapres_all_new;
        %ga.PM = PM_all;
        ga.time_cpet_all = time_cpet_all;
        ga.time_finapres_all = time_finapres_all;
        ga.AT_all = ga.AT_all(end);
        
        general_structures = ga;      



end

function plotting_data(idx, general_structures, subject_number)
        
        ga = general_structures;
        if ga.hipoxia
            title_text = sprintf('Sujeto %d Hipoxia', subject_number(idx));
        else
            title_text = sprintf('Sujeto %d Normoxia', subject_number(idx));
        end
    
        
        figure;
        subplot(2,1,1);
        plot((ga.time_finapres_all)/60, ga.finapres_all(:, 1));
        hold on;
    
        %figure;
        plot((ga.time_cpet_all)/60, ga.cpet_all(:,10));
    
    
        %xline((1200 + corrections_delay(idx))/60);
        xlabel('tiempo (min)'); ylabel('HR (bpm)');
        yyaxis right;
        plot((ga.time_cpet_all)/60, ga.cpet_all(:, end));
        legend('finapres', 'cpet', 'Power');
        ylabel('Power');
        
        title(title_text);
        subplot(2,1,2);        
        
        %plot((time_finapres_all)/60, PM_all);cl
        plot((ga.time_finapres_all)/60, ga.finapres_all(:, 2));
        %xline((1200 + corrections_delay(idx))/60);
        legend('PM finapres');
        xlabel('tiempo (min)'); ylabel('PM (mmHg)');       
        xlim([ga.time_cpet_all(1)/60,ga.time_cpet_all(end)/60])
        title(title_text);
end


function general_structures = setting_up_general_structures(general_structures)
    ga = general_structures;
    if ga.hipoxia
        finapres_filenames4 = ["../../data/HIPOXIA/CV/12.12.2024_09.36.27.nsc", "../../data/HIPOXIA/CV/12.12.2024_10.22.26.nsc", "../../data/HIPOXIA/CV/12.12.2024_10.35.08.nsc"];
        finapres_inits4 = [9*3600 + 36*60 + 27, 10*3600 + 22*60 + 26, 10*3600 + 35*60 + 8];
        finapres_not_nan_intevals4 = [35 * 60, -1];
        cpet_filenames4 = ["../../data/HIPOXIA/CV/ASCENSO/CPET_HIPOXIA_ASCENSO_CV.xlsx", "../../data/HIPOXIA/CV/NEW_CPET_HIPOXIA_EJERCICIO_CV.xlsx"];
        cpet_inits4 = [34787, 37448];  
        hipoxia_finapres_crazy4 = [90  170];
        
        
        finapres_filenames1 = ["../../data/HIPOXIA/AG/12.21.2024_09.27.47_1.nsc", "../../data/HIPOXIA/AG/12.21.2024_10.38.12.nsc"];
        %finapres_inits1 = [9*3600 + 27*60 + 47, 10*3600 + 38*60 + 12];
        finapres_inits1 = [9*3600 + 27*60 + 47 + 8.2*60, 10*3600 + 38*60 + 12 - 0.4*60];  %adding 8.5 minutes into the first and substracting 1 minute
        finapres_not_nan_intevals1 = [39 * 60, -1]; 
        cpet_filenames1 = ["../../data/HIPOXIA/AG/ASCENSO/CPET_HIPOXIA_ASCENSO_AG.xlsx", "../../data/HIPOXIA/AG/NEW_CPET_HIPOXIA_EJERCICIO_AG.xlsx"];
        cpet_inits1 = [34678, 37344];    
        
        finapres_filenames5 = ["../../data/HIPOXIA/VC/12.20.2024_09.59.39.nsc", "../../data/HIPOXIA/VC/12.20.2024_10.12.27.nsc","../../data/HIPOXIA/VC/12.20.2024_10.49.38.nsc" , "../../data/HIPOXIA/VC/12.20.2024_11.15.58.nsc"];
        finapres_inits5 = [9*3600 + 59*60 + 39, 10*3600 + 12*60 + 27, 10*3600 + 49*60 + 38, 11*3600 + 15*60 + 58 - 151];
        finapres_not_nan_intevals5 = [];
        cpet_filenames5 = ["../../data/HIPOXIA/VC/ASCENSO/CPET_HIPOXIA_ASCENSO_VC.xlsx", "../../data/HIPOXIA/VC/NEW_CPET_HIPOXIA_EJERCICIO_VC.xlsx"];
        cpet_inits5 = [36564, 40407];
        hipoxia_finapres_crazy5 = [101, 176; 336, 770];
        
        finapres_filenames6 = ["../../data/HIPOXIA/EC/12.11.2024_13.14.50.nsc", "../../data/HIPOXIA/EC/12.11.2024_14.13.45.nsc", "../../data/HIPOXIA/EC/12.11.2024_14.32.55.nsc",  "../../data/HIPOXIA/EC/12.11.2024_14.53.33.nsc"];
        finapres_inits6 = [13*3600 + 14*60 + 50, 14*3600 + 13*60 + 45, 14*3600 + 32*60 + 55, 14*3600 + 53*60 + 33];
        finapres_not_nan_intevals6 = [];
        cpet_filenames6 = ["../../data/HIPOXIA/EC/ASCENSO/CPET_HIPOXIA_ASCENSO_EC.xlsx", "../../data/HIPOXIA/EC/NEW_CPET_HIPOXIA_EJERCICIO_EC.xlsx"];
        cpet_inits6 = [49726, 52595];
        hipoxia_finapres_crazy6 = [104, 180];
    
        corrections_delay = [0,-2.2*60,30,1*60];
    
    elseif ga.normoxia
        finapres_filenames4 = ["../../data/NORMOXIA/CV/FINAPRES/01.31.2025_13.03.17.nsc"];
        finapres_inits4 = [13*3600 + 3*60 + 17];
        finapres_not_nan_intevals4 = [];
        cpet_filenames4 = ["../../data/NORMOXIA/CV/CPET/CPET-4_CV_cleaned.xlsx"];
        cpet_inits4 = [47570];
        
        finapres_filenames1 = ["../../data/NORMOXIA/AG/FINAPRES/08.30.2024_10.30.07.nsc"];
        finapres_inits1 = [10*3600 + 30*60 + 7];
        finapres_not_nan_intevals1 = [1, 14.64 * 60];
        cpet_filenames1 = ["../../data/NORMOXIA/AG/CPET/CPET-1_AG_cleaned.xlsx"];
        cpet_inits1 = [36958];
        
        
        finapres_filenames5 = ["../../data/NORMOXIA/VC/FINAPRES/01.17.2025_15.08.09.nsc","../../data/NORMOXIA/VC/FINAPRES/01.17.2025_15.30.00.nsc", "../../data/NORMOXIA/VC/FINAPRES/01.17.2025_15.32.39.nsc"];
        finapres_inits5 = [15*3600 + 8*60 + 9, 15*3600 + 30*60 + 0, 15*3600 + 32*60 + 39];
        finapres_not_nan_intevals5 = [];
        cpet_filenames5 = ["../../data/NORMOXIA/VC/CPET/CPET-5_VC_cleaned.xlsx"];
        cpet_inits5 = [55275];
        normoxia_finapres_crazy5 = [90, 209; 1617, 1817; 2758, 2810];
        
        finapres_filenames6 = ["../../data/NORMOXIA/EC/FINAPRES/04.08.2025_08.42.03.nsc"];
        cpet_filenames6 = ["../../data/NORMOXIA/EC/CPET/CPET-6_EC_cleaned.xlsx"];
        finapres_inits6 = [8*3600 + 42*60 + 3];
        finapres_not_nan_intevals6 = [1775,-1];
        cpet_inits6 = [31877];
        normoxia_finapres_crazy6 = [99, 175; 325, 340; 375, 410];
    
        corrections_delay = [0, -1.5*60,0,0];
    end
    
    general_structures.finapres_filenames_cell = {finapres_filenames1, finapres_filenames4, finapres_filenames5, finapres_filenames6};
    general_structures.cpet_filenames_cell = {cpet_filenames1, cpet_filenames4, cpet_filenames5, cpet_filenames6};
    general_structures.finapres_inits_cell = {finapres_inits1, finapres_inits4, finapres_inits5, finapres_inits6};
    general_structures.cpet_inits_cell = {cpet_inits1, cpet_inits4, cpet_inits5, cpet_inits6};
    general_structures.corrections_delay = corrections_delay;
    
    if ga.hipoxia
        general_structures.finapres_crazy_intervals = {0, hipoxia_finapres_crazy4, hipoxia_finapres_crazy5, hipoxia_finapres_crazy6};
    else
        general_structures.finapres_crazy_intervals = {0,0, normoxia_finapres_crazy5, normoxia_finapres_crazy6};
    end

    general_structures.finapres_notnan_intervals = {finapres_not_nan_intevals1, finapres_not_nan_intevals4, finapres_not_nan_intevals5, finapres_not_nan_intevals6};

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

function general_structures = align_and_resample(idx, general_structures)
% ALIGN_AND_RESAMPLE_STRUCTS
% 
%   [time2_new, vars2_new] = align_and_resample_structs(time1, time2, vars1, vars2)
%
%   Esta función alinea las series de tiempo `time2`, `vars2` con respecto a 
%   `time1`, `vars1`, eliminando los tiempos negativos y generando nuevas
%   versiones (`time2_new`, `vars2_new`) con el mismo tamaño que `time1`.
%
%   Entradas:
%       time1  - vector temporal (Nx1 o 1xN)
%       time2  - vector temporal (Mx1 o 1xM)
%       vars1  - matriz asociada a time1 (N x P)
%       vars2  - matriz asociada a time2 (M x Q)
%
%   Salidas:
%       time2_new - vector temporal con la misma longitud que time1
%       vars2_new - matriz interpolada con la misma cantidad de filas que vars1
%
%   Comportamiento:
%     1. Elimina cualquier valor de time2 < 0 y sus correspondientes en vars2.
%     2. Ajusta el rango de time2 al rango cubierto por time1.
%     3. Interpola (linealmente) vars2 sobre los valores de time1.
%     4. Si vars2 tiene NaN fuera del rango, se extienden los valores más cercanos.

    

    ga = general_structures;
    time1 = ga.time_cpet_all;
    time2 = ga.time_finapres_all;
    vars1 = ga.cpet_all;
    vars2 = ga.finapres_all;

    if idx == 4 && general_structures.hipoxia
        time1 = time1 - 400;   %to remove the initial 400 seconds of no adquisition
    end

    % --- Validaciones ---
    if numel(time1) ~= size(vars1,1)
        error('time1 y vars1 deben tener la misma cantidad de filas.');
    end
    if numel(time2) ~= size(vars2,1)
        error('time2 y vars2 deben tener la misma cantidad de filas.');
    end

    % --- Eliminar tiempos negativos ---
    mask_valid = time2 >= 0;
    time2 = time2(mask_valid);
    vars2 = vars2(mask_valid, :);

    if isempty(time2)
        % Si todo fue eliminado, devolver matriz de ceros
        vars2_new = zeros(size(vars1,1), size(vars2,2));
        return;
    end

    % --- Interpolar dentro del rango ---
    vars2_interp = interp1(time2, vars2, time1, 'linear', NaN);

    % --- Rellenar con ceros fuera del rango ---
    out_of_range = (time1 < min(time2)) | (time1 > max(time2));
    vars2_interp(out_of_range, :) = 0;

    % --- Output ---
    vars2_new = vars2_interp;


    %Incorporing finapres not nan mask

    f_n_i = general_structures.finapres_notnan_interval;    
    finapres_not_nan_mask = ones(size(vars2_new,1));
    if ~isempty(f_n_i)
        init_ = f_n_i(1);
        end_  = f_n_i(end);

        
        if end_ == -1
            mask = (time1 > init_) .* (time1 < time1(end));
            mask = logical(mask);
            finapres_not_nan_mask(mask) = 0;
        else
            mask = (time1 > init_) .* (time1 < end_);
            mask = logical(mask);
            finapres_not_nan_mask(mask) = 0;
        end
        
    end 

    %remove negative values
    non_negative_time_mask = (time1 >= 0);
    time1 = time1(non_negative_time_mask);
    vars1 = vars1(non_negative_time_mask,:);
    vars2_new = vars2_new(non_negative_time_mask, :);
    finapres_not_nan_mask = finapres_not_nan_mask(non_negative_time_mask);

    %obtain the time at 3rd ladder value
    Power = vars1(:,end);
    max_power = max(Power);
    max_power_time = time1(Power == max_power);
    last_time_max_power = max_power_time(end);
     
    last_level_time = time1(logical((Power == 75) .* (time1 < last_time_max_power)));
    last_level_time_value = last_level_time(end);


    general_structures.finapres_all = vars2_new;
    general_structures.time_cpet_all = time1;
    general_structures.cpet_all = vars1;
    general_structures.finapres_not_nan_mask = finapres_not_nan_mask;
    general_structures.simulation_time = last_level_time_value;
    
end

function [t_exp, y_exp, VO2_poly, VCO2_poly, fiO2_poly, TI_poly, Tresp_poly, basal] = polish_and_filter(idx, general_structures, plot_state, VO2_type)

        %plot_state = 0;
        %VO2_type = "poly";

        time_cpet = general_structures.time_cpet_all;       
        finapres_all = general_structures.finapres_all;
        cpet_all = general_structures.cpet_all;
        ladder_points = general_structures.ladder_points;
        AT = general_structures.AT_all;
        hipoxia = general_structures.hipoxia;
        finapres_notnan_mask = general_structures.finapres_not_nan_mask;

        simulation_time = general_structures.simulation_time;
        disp(simulation_time);

        PM = finapres_all(:, 2);
        PS = finapres_all(:, 3);
        PD = finapres_all(:, 4);

        VT = cpet_all(:, 1);
        dVE = cpet_all(:, 2);
        TI = cpet_all(:, 3);
        Tresp = cpet_all(:, 4);
        PAO2 = cpet_all(:, 5);
        PACO2 = cpet_all(:, 6);
        VO2 = cpet_all(:, 7);
        VCO2 = cpet_all(:, 8);
        fiO2 = cpet_all(:, 9);
        HR = cpet_all(:, 10);
        on_exercise = cpet_all(:, 11);
        Pow = cpet_all(:, 12);

        if idx == 1 && ~hipoxia
            HR = finapres_all(:, 1);
        end 


        %f_all = [HR(4:end), PM(4:end), PS(4:end), PD(4:end)];
        %X_cpet.all = [VT, dVE, TI, Tresp, PAO2, PACO2, VO2, VCO2, fiO2, HR_cpet, on_exercise, Pow];


        
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
            [~, TI_poly, ~, ~, ~] = bestPolynomialFit(time_cpet, TI, 15, "TI (s)" , plot_state);             
            [~, Tresp_poly, ~, ~, ~] = bestPolynomialFit(time_cpet, Tresp, 15, "Tresp (s)" , plot_state); 
            [~, fiO2_poly, ~, ~, ~] = bestPolynomialFit(time_cpet, fiO2, 4, "fiO2 (%)", plot_state);
        elseif VO2_type == "ori"
            VO2_poly = VO2(1:end-3,:)/1000;
            VCO2_poly = VCO2(1:end-3,:)/1000;
            fiO2_poly = fiO2(1:end-3,:);
            TI_poly = TI(1:end-3,:);
            Tresp_poly = Tresp(1:end-3,:);

        end

       

        if VO2_type ~= "ori"
            y_exp = [y_exp, Pow(1:end-3, :), finapres_notnan_mask(1:end-3), on_exercise(4:end)];
            
        else
            disp('wait');
            y_exp = [y_exp, VO2_poly, VCO2_poly, Pow(1:end-3, :), finapres_notnan_mask(1:end-3), on_exercise(4:end)];
            
        end
        basal = [mean(abs(VO2(1:30))), mean(abs(VCO2(1:30)))]/1000;
        basal = [basal, AT/1000, ladder_points(:,1)', ladder_points(:,2)', simulation_time];


end

end
        

        
        
        