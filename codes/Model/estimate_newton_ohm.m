function new_pars = estimate_newton_ohm(x, pars, patient_idx, format_output)

    %structura : [BW, Hgt, Gender]
    


    if nargin < 4
        format_output = 'dict';
        if nargin < 3
            patient_idx = 0;
        end
    end
    % Compute estimated parameters based on x multipliers

    % Take data from patient measurements
    
    BW        = pars('BW');%70; %data.BW;                         % Body weight (kg)                
    Hgt       = pars('Hgt');%172;%data.Hgt;                        % Height (cm)                     
    Gender    = pars('Gender');%2;%data.Gender;                     % Gender (1=female, 2=male)       
    
    
    %% Total blood volume (mL) from general measurements
    if (Gender == 2)  
           ToTBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * BW) + 0.6041) * 1000; % Total blood volume (mL)
           BSA   = 0.000579479 * BW^0.38 * Hgt^1.24;                          % Body surface area (m^2)
       else
           ToTBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * BW) + 0.1833) * 1000; % Total blood volume (mL)
           BSA   = 0.000975482 * BW^0.46 * Hgt^1.08;                          % Body surface area (m^2)
    end
    
    %% Cardiac Output (mL/s)

    Ave_HR        = 60;%data.Ave_HR;                 % Average heart rate (beats/min)    
    %TotFlow       = TotBV/60;                    % Total flow
    %CO            = TotFlow;                     % Cardiac output 

    HI            = Ave_HR/60;                   % Heart rate (s)
    %HI              =100;                   % Heart rate

    PM = 100; %1/3 * 120 + 2/3*75;
    
    %fprintf("Vtot: %f\n", pars('V_tot'));
    %fprintf("ToTBV: %f\n", ToTBV);
    %ToTBV = 5027.6;  
    %ToTBV = pars('V_tot');
    TotFlow = ToTBV/60;  
    CO = TotFlow;
    %% Blood Pressure
    P_m_p = x(1) * PM;
    P_h_p = x(2) * PM;
    P_p_p = x(3) * PM;
    P_s_p = x(4) * PM;
    P_e_p = x(5) * PM;
    P_b_p = x(6) * PM - 10; %Intracraneal pressure has to be substracted
    

    P_m_v = 3;
    P_h_v = 4.5;
    P_p_v = 5;
    P_s_v = 6;
    P_e_v = 7.5;
    P_b_v = 5;
    P_vc_v = 2.5;   %This value has to be correctly acquired,.

    %% Volume (Delta Volume = VTotal - Vunestressed)
    V_m_Total =  x(7)  * ToTBV;  % Muscle Circulation (Active + Resting, Arterial + Venous): 15%
    V_h_Total =  x(8)  * ToTBV;   % Heart: 7%
    V_p_Total =  x(9)  * ToTBV;   % Pulmonary Circulation: 9%
    V_s_Total =  x(10) * ToTBV;   % Splanchnic Circulation (Arterial + Venous): 20%
    V_e_Total =  x(11) * ToTBV;   % Extrasplanchnic Circulation (Arterial + Venous): 15%
    V_b_Total =  x(12) * ToTBV;   % Brain Circulation (Arterial + Venous): 8%
    V_vc_Total = x(13) * ToTBV;  % Vena Cava (Superior & Inferior): 5%   
   
    V_m_v = x(14) * V_m_Total;  
    V_h_v =  x(15) * V_h_Total;  
    V_p_v =  x(16) * V_p_Total;  
    V_s_v =  x(17) * V_s_Total;  
    V_e_v =  x(18) * V_e_Total;  
    V_b_v =  x(19) * V_b_Total;  
    V_vc_v =  V_vc_Total;  
    
    V_m_p =  V_m_Total  - V_m_v   ;
    V_h_p =   V_h_Total  - V_h_v  ; 
    V_p_p =   V_p_Total  - V_p_v  ; 
    V_s_p =   V_s_Total  - V_s_v  ; 
    V_e_p =   V_e_Total  - V_e_v  ; 
    V_b_p =   V_b_Total  - V_b_v  ; 
    V_vc_p = 0; 

    V_m_v_un = x(20) * V_m_v  ;
    V_h_v_un = x(21) * V_h_v  ;
    V_p_v_un = x(22) * V_p_v  ;
    V_s_v_un = x(23) * V_s_v  ;
    V_e_v_un = x(24) * V_e_v  ;
    V_b_v_un = x(25) * V_b_v  ;
    V_vc_v_un =  V_vc_Total ;

    V_m_p_un = x(26) * V_m_p  ;
    V_h_p_un = x(27) * V_h_p  ;
    V_p_p_un = x(28) * V_p_p  ;
    V_s_p_un = x(29) * V_s_p  ;
    V_e_p_un = x(30) * V_e_p  ;
    V_b_p_un = x(31) * V_b_p  ;
    %V_vc_p_un =  V_vc_Total ;

    V_m_v_delta = V_m_v - V_m_v_un;
    V_h_v_delta =   V_h_v -  V_h_v_un;
    V_p_v_delta =   V_p_v -  V_p_v_un;
    V_s_v_delta =  V_s_v - V_s_v_un;
    V_e_v_delta =   V_e_v -  V_e_v_un;
    V_b_v_delta =   V_b_v -  V_b_v_un;


    V_m_p_delta = V_m_p - V_m_p_un;
    V_h_p_delta =   V_h_p -  V_h_p_un;
    V_p_p_delta =   V_p_p -  V_p_p_un;
    V_s_p_delta =  V_s_p - V_s_p_un;
    V_e_p_delta =   V_e_p -  V_e_p_un;
    V_b_p_delta =   V_b_p -  V_b_p_un;
    
    
    % V_rm_p = 0.45 * ToTBV * 0.01;
    % V_h_p = 0.103 * ToTBV * 0.01;
    % V_p_p = 0.459 * ToTBV * 0.01;
    % V_s_p = 1.19 * ToTBV * 0.01;
    % V_e_p = 0.57 * ToTBV * 0.01;
    % V_b_p = 0.3 * ToTBV * 0.01;
    % 
    % 
    % V_rm_v = 2.18 * ToTBV * 0.01;
    % V_h_v = 0.42 * ToTBV * 0.01;
    % V_p_v = 0.4556 * ToTBV * 0.01;
    % V_s_v = 6.2 * ToTBV * 0.01;
    % V_e_v = 2.77 * ToTBV * 0.01;
    % V_b_v = 1.27 * ToTBV * 0.01;
    % V_vc_v = 0.534 * ToTBV * 0.01;

    %% Flows
    Q_m =  x(32) * TotFlow;
    Q_h =  x(33) * TotFlow;
    Q_p =  x(34) * TotFlow;
    Q_s =  x(35) * TotFlow;
    Q_e =  x(36) * TotFlow;
    Q_b =  x(37) * TotFlow;
    Q_vc = x(38) * TotFlow;

    %% Resistances
    
    R_m_p = (P_m_p - P_m_v )/Q_m;
    R_h_p = (P_h_p - P_h_v )/Q_h;
    R_p_p = (P_p_p - P_p_v )/Q_p;
    R_s_p = (P_s_p - P_s_v )/Q_s;
    R_e_p = (P_e_p - P_e_v )/Q_e;
    R_b_p = (P_b_p - P_b_v )/Q_b;
    %R_vc = (P_vc_p - P_vc_v )/Q_vc;
    R_m_v = (P_m_v - P_vc_v )/Q_m;
    R_h_v = (P_h_v - P_vc_v )/Q_h;
    R_p_v = (P_p_v - P_vc_v )/Q_p;
    R_s_v = (P_s_v - P_vc_v )/Q_s;
    R_e_v = (P_e_v - P_vc_v )/Q_e;
    R_b_v = (P_b_v - P_vc_v )/Q_b;


    %% Capacitances

    %C_m_p = V_m_p_delta / P_m_p ;
    %C_h_p =  V_h_p_delta / P_h_p;
    %C_p_p =  V_p_p_delta / P_p_p;
    %C_s_p =  V_s_p_delta / P_s_p;
    %C_e_p =  V_e_p_delta / P_e_p;
    %C_b_p =  V_b_p_delta / P_b_p;
    %
    %C_m_v = V_m_v_delta / P_m_v;
    %C_h_v =  V_h_v_delta / P_h_v;
    %C_p_v =  V_p_v_delta / P_p_v;
    %C_s_v =  V_s_v_delta / P_s_v;
    %C_e_v =  V_e_v_delta / P_e_v;
    %C_b_v =  V_b_v_delta / P_b_v;
    %C_vc_v = V_vc_v / P_vc_v;
%
    % Store old values before updating
    old_pars = containers.Map(keys(pars), values(pars));
    
    %Renaming
    pars('V_unstressed_am_p') = V_m_p_un * 0.6;
    %V_unstressed_am_v = V_m_v_un;
    pars('V_unstressed_b_p') = V_b_p_un;
    pars('V_unstressed_b_v') = V_b_v_un;
    pars('V_unstressed_e_p') = V_e_p_un;
    %V_unstressed_e_v = V_e_v_un;
    pars('V_unstressed_h_p') = V_h_p_un;
    pars('V_unstressed_h_v') = V_h_v_un;
    %V_unstressed_la = 
    %V_unstressed_lv
    %V_unstressed_pa
    pars('V_unstressed_pp') = V_p_p_un;
    pars('V_unstressed_pv') = V_p_v_un;
    %V_unstressed_ra
    pars('V_unstressed_rm_p') = V_m_p_un * 0.4;
    %V_unstressed_rm_v = V_m_v_un;
    %V_unstressed_rv 
    pars('V_unstressed_s_p') = V_s_p_un;
    %V_unstressed_s_v = V_s_v_un;
    %V_unstressed_sa = 
    pars('V_unstressed_vc') = V_vc_v_un;

    pars('V_u_am_v_0') = V_m_v_un * 0.6;
    pars('V_u_e_v_0')= V_e_v_un;
    pars('V_u_rm_v_0') = V_m_v_un * 0.4;
    pars('V_u_s_v_0') = V_s_v_un;

    pars('R_am_n') = R_m_v;
    pars('R_am_p_0') = R_m_p;
    pars('R_b_n') = R_b_v;
    pars('R_e_n') = R_e_v;
    pars('R_e_p_0') = R_e_p;
    pars('R_h_n') = R_h_v;
    pars('R_h_p_n') = R_h_p;
    %R_la = ;
    %R_pa = ;
    pars('R_p_p_n') = R_p_p;
    pars('R_pv') = R_p_v;
    %R_ra = ;
    pars('R_rm_n') = R_m_v;
    pars('R_rm_p_0') = R_m_p;
    pars('R_s_n') = R_s_v;
    pars('R_s_p_0') = R_s_p;

    %C_am_p = C_m_p;
    %C_rm_p = C_m_p;
    %C_am_v = C_m_v;
    %C_rm_v = C_m_v;
    %C_pp = C_p_p;
    %C_pv = C_p_v;

    %new_pars = pars;

    if strcmp(format_output, 'dict')
        new_pars = pars;
    elseif strcmp(format_output, 'array')
        new_pars = values(pars);
    end

    % %% Generate comparison table
    % fprintf('\n=== TABLA COMPARATIVA DE PARÁMETROS ===\n');
    % fprintf('%-25s %15s %15s %15s\n', 'Parámetro', 'Valor Antiguo', 'Valor Nuevo', 'Diferencia');
    % fprintf('%s\n', repmat('-', 73, 1));
    % 
    % % List of parameters that were updated
    % updated_params = {'V_unstressed_am_p', 'V_unstressed_b_p', 'V_unstressed_b_v', ...
    %                   'V_unstressed_e_p', 'V_unstressed_h_p', 'V_unstressed_h_v', ...
    %                   'V_unstressed_pp', 'V_unstressed_pv', 'V_unstressed_rm_p', ...
    %                   'V_unstressed_s_p', 'V_unstressed_vc', 'V_u_am_v_0', ...
    %                   'V_u_e_v_0', 'V_u_rm_v_0', 'V_u_s_v_0', 'R_am_n', ...
    %                   'R_am_p_0', 'R_b_n', 'R_e_n', 'R_e_p_0', 'R_h_n', ...
    %                   'R_h_p_n', 'R_p_p_n', 'R_pv', 'R_rm_n', 'R_rm_p_0', ...
    %                   'R_s_n', 'R_s_p_0'};
    % 
    % % Prepare data for table
    % param_names = {};
    % old_vals = [];
    % new_vals = [];
    % diff_vals = [];
    % 
    % for i = 1:length(updated_params)
    %     param_name = updated_params{i};
    %     if isKey(old_pars, param_name) && isKey(pars, param_name)
    %         old_val = old_pars(param_name);
    %         new_val = pars(param_name);
    %         diff_val = new_val - old_val;
    % 
    %         param_names{end+1} = param_name;
    %         old_vals(end+1) = old_val;
    %         new_vals(end+1) = new_val;
    %         diff_vals(end+1) = diff_val;
    % 
    %         fprintf('%-25s %15.6f %15.6f %15.6f\n', param_name, old_val, new_val, diff_val);
    %     elseif isKey(pars, param_name)
    %         new_val = pars(param_name);
    % 
    %         param_names{end+1} = param_name;
    %         old_vals(end+1) = NaN;
    %         new_vals(end+1) = new_val;
    %         diff_vals(end+1) = NaN;
    % 
    %         fprintf('%-25s %15s %15.6f %15s\n', param_name, 'N/A', new_val, 'Nuevo');
    %     end
    % end
    % 
    % fprintf('%s\n', repmat('-', 73, 1));
    % fprintf('Paciente: %d\n', patient_idx);
    % fprintf('\n');
    % 
    % %% Save comparison table to Excel
    % % Create table data
    % tableData = cell(length(param_names) + 1, 4);
    % 
    % % Add headers
    % tableData{1, 1} = 'Parámetro';
    % tableData{1, 2} = 'Valor Antiguo';
    % tableData{1, 3} = 'Valor Nuevo';
    % tableData{1, 4} = 'Diferencia';
    % 
    % % Add data rows
    % for i = 1:length(param_names)
    %     tableData{i+1, 1} = param_names{i};
    %     if isnan(old_vals(i))
    %         tableData{i+1, 2} = 'N/A';
    %     else
    %         tableData{i+1, 2} = old_vals(i);
    %     end
    %     tableData{i+1, 3} = new_vals(i);
    %     if isnan(diff_vals(i))
    %         tableData{i+1, 4} = 'Nuevo';
    %     else
    %         tableData{i+1, 4} = diff_vals(i);
    %     end
    % end
    % 
    % % Create summary statistics
    % valid_diffs = diff_vals(~isnan(diff_vals));
    % summary_data = cell(5, 2);
    % summary_data{1, 1} = 'RESUMEN';
    % summary_data{1, 2} = '';
    % summary_data{2, 1} = 'Paciente';
    % summary_data{2, 2} = patient_idx;
    % summary_data{3, 1} = 'Total parámetros';
    % summary_data{3, 2} = length(param_names);
    % summary_data{4, 1} = 'Parámetros modificados';
    % summary_data{4, 2} = length(valid_diffs);
    % summary_data{5, 1} = 'Parámetros nuevos';
    % summary_data{5, 2} = sum(isnan(old_vals));
    % 
    % % Generate filename
    % filename = sprintf('estimate_ohmnewton_parametros_paciente_%d.xlsx', patient_idx);
    % 
    % % Write to Excel
    % try
    %     % Write main table
    %     writecell(tableData, filename, 'Sheet', 'Comparación', 'Range', 'A1');
    % 
    %     % Write summary
    %     writecell(summary_data, filename, 'Sheet', 'Comparación', 'Range', 'F1');
    % 
    %     fprintf('Tabla guardada exitosamente en: %s\n', filename);
    % catch ME
    %     warning('Error al guardar el archivo Excel: %s', ME.message);
    %     fprintf('Intentando guardar en formato CSV...\n');
    % 
    %     % Fallback to CSV if Excel fails
    %     csv_filename = sprintf('parametros_paciente_%d.csv', patient_idx);
    %     writecell(tableData, csv_filename);
    %     fprintf('Tabla guardada en formato CSV: %s\n', csv_filename);
    % end
   
end