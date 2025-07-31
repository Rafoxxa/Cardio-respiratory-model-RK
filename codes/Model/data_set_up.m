function data_set_up(num_all_patients, normoxia_patients, hipoxia_patients)
    
    for patient_idx = 1:num_all_patients
        if ismember(patient_idx, normoxia_patients)
            [texp, yexp, VO2_poly, VCO2_poly, fO2_poly, basal] = data_preprocessing(patient_idx, "normoxia", "-");
            filename_n = sprintf("../fast_data/%d/normoxia_data_preprocessed.mat", patient_idx); 
            ladder_points = basal(4:end);
            AT = basal(3);
            basal = basal(1:2);
            text = sprintf("basal normoxia para paciente %d", patient_idx);
            disp(text);
            disp(basal);
            
            len_ladder_points = length(ladder_points);
            VO2_ladder_points = ladder_points(1:len_ladder_points/2);
            VCO2_ladder_points = ladder_points((len_ladder_points/2 + 1):end);
            save(filename_n, "texp", "yexp", "VO2_poly", "VCO2_poly", "fO2_poly", "basal", "VO2_ladder_points", "VCO2_ladder_points", "AT");
        end

        if ismember(patient_idx, hipoxia_patients)
            [texp, yexp, VO2_poly, VCO2_poly, fO2_poly, basal] = data_preprocessing(patient_idx, "hipoxia", "mix");         
            filename_h = sprintf("../fast_data/%d/hipoxia_data_preprocessed.mat", patient_idx);
            ladder_points = basal(4:end);
            %AT = basal(3); for the moment just use normoxia's AT.
            basal = basal(1:2);        
            text = sprintf("basal hipoxia para paciente %d", patient_idx);
            disp(text);
            disp(basal);
            disp("______________________")
            len_ladder_points = length(ladder_points);
            VO2_ladder_points = ladder_points(1:len_ladder_points/2);
            VCO2_ladder_points = ladder_points((len_ladder_points/2 + 1):end);
            save(filename_h, "texp", "yexp", "VO2_poly", "VCO2_poly", "fO2_poly", "basal", "VO2_ladder_points", "VCO2_ladder_points", "AT");
        end
    end



end