    clear all
    BW        = 70; %data.BW;                         % Body weight (kg)                
    Hgt       = 172;%data.Hgt;                        % Height (cm)                     
    Gender    = 2;%data.Gender;                     % Gender (1=female, 2=male)       

    %% Total blood volume (mL)
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
    ToTBV = 5027.6;  
    TotFlow = 5027.6/60;  
    CO = TotFlow;
    %% Blood Pressure
    P_m_p = 0.95*PM;
    P_h_p = 1*PM;
    P_p_p = 0.25 * PM;
    P_s_p = 0.7*PM;
    P_e_p = 0.8*PM;
    P_b_p = 1*PM - 10; %Intracraneal pressure has to be substracted
    

    P_m_v = 3;
    P_h_v = 4.5;
    P_p_v = 5;
    P_s_v = 6;
    P_e_v = 7.5;
    P_b_v = 5;
    P_vc_v = 2.5;   %This value has to be correctly acquired,.

    %% Volume (Delta Volume = VTotal - Vunestressed)
    V_m_Total = 0.15 * ToTBV;  % Muscle Circulation (Active + Resting, Arterial + Venous): 15%
    V_h_Total = 0.07 * ToTBV;   % Heart: 7%
    V_p_Total = 0.09 * ToTBV;   % Pulmonary Circulation: 9%
    V_s_Total = 0.20 * ToTBV;   % Splanchnic Circulation (Arterial + Venous): 20%
    V_e_Total = 0.15 * ToTBV;   % Extrasplanchnic Circulation (Arterial + Venous): 15%
    V_b_Total = 0.08 * ToTBV;   % Brain Circulation (Arterial + Venous): 8%
    V_vc_Total = 0.05 * ToTBV;  % Vena Cava (Superior & Inferior): 5%   
   
    V_m_v = 0.85 * V_m_Total;  
    V_h_v =  0.5 * V_h_Total;  
    V_p_v =  0.7 * V_p_Total;  
    V_s_v =  0.9 * V_s_Total;  
    V_e_v =  0.8 * V_e_Total;  
    V_b_v =  0.7 * V_b_Total;  
    V_vc_v =  V_vc_Total;  
    
    V_m_p =  V_m_Total  - V_m_v   ;
    V_h_p =   V_h_Total  - V_h_v  ; 
    V_p_p =   V_p_Total  - V_p_v  ; 
    V_s_p =   V_s_Total  - V_s_v  ; 
    V_e_p =   V_e_Total  - V_e_v  ; 
    V_b_p =   V_b_Total  - V_b_v  ; 
    V_vc_p = 0; 

    V_m_v_un = 0.65 * V_m_v  ;
    V_h_v_un =  0.25 * V_h_v  ;
    V_p_v_un =  0.55 * V_p_v  ;
    V_s_v_un =  0.8 * V_s_v  ;
    V_e_v_un =  0.65 * V_e_v  ;
    V_b_v_un =  0.45 * V_b_v  ;
    V_vc_v_un =  V_vc_Total ;

    V_m_p_un = 0.05 * V_m_p  ;
    V_h_p_un =  0.1 * V_h_p  ;
    V_p_p_un =  0.05 * V_p_p  ;
    V_s_p_un =  0.05 * V_s_p  ;
    V_e_p_un =  0.05 * V_e_p  ;
    V_b_p_un =  0.1 * V_b_p  ;
    %V_vc_p_un =  V_vc_Total ;

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
    Q_rm = 0.2 * TotFlow;
    Q_h = 0.05 * TotFlow;
    Q_p = 1 * TotFlow;
    Q_s = 0.225 * TotFlow;
    Q_e = 0.235 * TotFlow;
    Q_b = 0.15 * TotFlow;
    Q_vc = 1 * TotFlow;

    %% Resistances
    
    R_m_p = (P_m_p - P_m_v )/Q_rm;
    R_h_p = (P_h_p - P_h_v )/Q_h;
    R_p_p = (P_p_p - P_p_v )/Q_p;
    R_s_p = (P_s_p - P_s_v )/Q_s;
    R_e_p = (P_e_p - P_e_v )/Q_e;
    R_b_p = (P_b_p - P_b_v )/Q_b;
    %R_vc = (P_vc_p - P_vc_v )/Q_vc;
    R_m_v = (P_m_v - P_vc_v )/Q_rm;
    R_h_v = (P_h_v - P_vc_v )/Q_h;
    R_p_v = (P_p_v - P_vc_v )/Q_p;
    R_s_v = (P_s_v - P_vc_v )/Q_s;
    R_e_v = (P_e_v - P_vc_v )/Q_e;
    R_b_v = (P_b_v - P_vc_v )/Q_b;


    %% Capacitances

    C_m_p = V_m_p / P_m_p ;
    C_h_p =  V_h_p / P_h_p;
    C_p_p =  V_p_p / P_p_p;
    C_s_p =  V_s_p / P_s_p;
    C_e_p =  V_e_p / P_e_p;
    C_b_p =  V_b_p / P_b_p;
    C_m_v = V_m_v / P_m_v;
    C_h_v =  V_h_v / P_h_v;
    C_p_v =  V_p_v / P_p_v;
    C_s_v =  V_s_v / P_s_v;
    C_e_v =  V_e_v / P_e_v;
    C_b_v =  V_b_v / P_b_v;
    C_vc_v = V_vc_v / P_vc_v;

    %Renaming
    V_unstressed_am_p = V_m_p_un;
    V_unstressed_am_v = V_m_v_un;
    V_unstressed_b_p = V_b_p_un;
    V_unstressed_b_v = V_b_v_un;
    V_unstressed_e_p = V_e_p_un;
    V_unstressed_e_v = V_e_v_un;
    V_unstressed_h_p = V_h_p_un;
    V_unstressed_h_v = V_h_v_un;
    %V_unstressed_la = 
    %V_unstressed_lv
    %V_unstressed_pa
    V_unstressed_pp = V_p_p_un;
    V_unstressed_pv = V_p_v_un;
    %V_unstressed_ra
    V_unstressed_rm_p = V_m_p_un;
    V_unstressed_rm_v = V_m_v_un;
    %V_unstressed_rv 
    V_unstressed_s_p = V_s_p_un;
    V_unstressed_s_v = V_s_v_un;
    %V_unstressed_sa = 
    V_unstressed_vc = V_vc_v_un;

    V_u_am_v_0 = V_m_v_un;
    V_u_e_v_0 = V_e_v_un;
    V_u_rm_v_0 = V_m_v_un;
    V_u_s_v_0 = V_s_v_un;

    R_am_n = R_m_v;
    R_am_p_0 = R_m_p;
    R_b_n = R_b_v;
    R_e_n = R_e_v;
    R_e_p_0 = R_e_p;
    R_h_n = R_h_v;
    R_h_p_n = R_h_p;
    %R_la = ;
    %R_pa = ;
    R_pp = R_p_p;
    R_pv = R_p_v;
    %R_ra = ;
    R_rm_n = R_m_v;
    R_rm_p_0 = R_m_p;
    R_s_n = R_s_v;
    R_s_p_0 = R_s_p;

    C_am_p = C_m_p;
    C_rm_p = C_m_p;
    C_am_v = C_m_v;
    C_rm_v = C_m_v;
    C_pp = C_p_p;
    C_pv = C_p_v;

    % Load parameter names from the Excel file
    file_name = 'excel_para_ori.xlsx'; % Replace with your file name
    param_data = readtable(file_name);
    
    % Assume the first column contains the parameter names
    param_names = param_data{:, 1};
    
    % Create a dictionary (struct) from your script variables
    vars = whos;  % Get all variables in the workspace
    data_dict = struct();
    
    for i = 1:length(vars)
        var_name = vars(i).name;
        data_dict.(var_name) = eval(var_name); % Save all variables into struct
    end
    
    % Match parameter names and retrieve values
    n_params = length(param_names);
    estimation_results = cell(n_params, 1); % Placeholder for estimated results
    
    for i = 1:n_params
        param_name = param_names{i};
        if isfield(data_dict, param_name)
            estimation_results{i} = data_dict.(param_name); % Fetch value if exists
        else
            estimation_results{i} = NaN; % Mark as missing if not found
        end
    end
    
    % Add a new column to the table
    param_data.EstimationResults = estimation_results;
    
    % Save the updated table to a new Excel file
    output_file_name = 'updated_parameters.xlsx';
    writetable(param_data, output_file_name);
    
    disp(['Results have been saved to ', output_file_name]);








    

    % vars = whos;  % Get all variables in workspace
    % data_dict = struct();
    % 
    % for i = 1:length(vars)
    %     var_name = vars(i).name;
    %     data_dict.(var_name) = eval(var_name);  % Store variable in struct
    % end
    % 
    % %% Save to Excel
    % var_names = fieldnames(data_dict);
    % var_values = struct2cell(data_dict);
    % T = table(var_names, var_values);
    % writetable(T, 'estimated_OhmNewton_parameters.xlsx', 'WriteVariableNames', false);
    % 
    % disp('All variables have been saved to variables.xlsx');


    %% Blood Pressure
    % P_SAsys = 120;
    % P_SAdia = 75;
    % % Systemic Arteries (SA)
    % P_TAsys  = P_SAsys;                   % Thoracic systolic SA pressure  (mmHg)    
    % P_TAdia  = P_SAdia;                   % Thoracic diastolic SA pressure (mmHg)   
    % P_TAm    = P_TAsys/3 + 2*P_TAdia/3;        % Thoracic SA mean pressure      (mmHg)
    % pm= P_TAm;

    % P_ABsys = P_TAsys*0.98;                    % Abdomen systolic SA pressure (mmHg)
    % P_ABdia = P_TAdia*0.98;                    % Abdomen diastolic SA pressure (mmHg)
    % P_ABm   = P_ABsys/3 + 2*P_ABdia/3;         % Abdomen SA mean pressure (mmHg)

    % P_LAsys  = P_SAsys*0.95;              % Lower body systolic SA pressure  (mmHg)    
    % P_LAdia  = P_SAdia*0.95;              % Lower body diastolic SA pressure (mmHg)   
    % P_LAm    = P_LAsys/3 + 2*P_LAdia/3;        % Lower body SA mean pressure      (mmHg)

    % % Systemic Veins (SV)
    % P_TVU      = 2.5;                          % Thoracic SV pressure (mmHg)
    % P_AB       = 2.75;                         % Abdomen SV pressure (mmHg) 
    % P_VL       = 3;                            % Lower body SV pressure 
    % pt         = P_TVU;                         % Superior cava vei

    % % Pulmonary Arteries (PA)
    % P_PAsys   = 24;                            % Systolic  PA pressure (mmHg)  
    % P_PAdia   = 8;                             % Diastolic PA pressure (mmHg)  
    % P_PAmean   =P_PAsys/3 + 2*P_PAdia/3;       % PA mean pressure (mmHg)  

    % % Pulmonary Veins (PU)
    % P_PV       = 5;                            % Pulmonary vein pressure (mmHg) 

    % % Left ventricle (LV)
    % P_LVdia       = 2.5;                       % Diastolic  LV pressure (mmHg)      
    % P_LVsyst  =  1.01 * P_TAsys;               % Systolic  LV pressure (mmHg)



    % % Right Ventricle (RV)
    % P_RVsys   = P_PAsys*1.01;                 % Systolic  RV pressure (mmHg) BC Lampert,  1.05
    % P_RVdia    = 2.5;                          %Recordar que esto es por ahora          1.05


    % V_d_lvf    = 10;                % LV end systolic zero pressure volume (mL) from Smith & Andreassen
    % V_d_rvf    = 0.9*V_d_lvf;       % RV end systolic zero pressure volume (mL) from Smith & Andreassen

    % V_LVm      = 50-V_d_lvf;        % Minimum left ventricular volume minus dead volume 50 average diastolic left ventricular volume
    % V_LVM      = 110-V_d_lvf;       % Minimum left ventricular volume minus dead volume 50 average diastolic left ventricular volume
    % V_RVM      = 1.01*V_LVM;        % Max RV volume (10% higher than the  V_LV)  
    % V_RVm      = 1.01*V_LVm;        % Min RV volume (10% higher than the  V_LV)

    % % Distribution of volume outside the heart (fractions add to 1)

    % % Circulating (total volume)
    % Circ_pa   = 0.02*TotBV;    % PA volume from Benekin  2% of total volume 
    % pa_UN     = Circ_pa*0.42;  % Pulmonary unestressed volume in arteries

    % Circ_pu   = 0.08*TotBV;    % PU volume from Benekin  8% of total volume
    % pu_UN     = Circ_pu*0.89;  % Pulmonary unestressed volume in veins

    % Circ_saT =  TotBV*0.0595;  % Arteries thorax volume TVol*SA(0.2)*Upper body(0.85)*thorax(0.35)  from Benekin0.27 
    % saT_UN   =  Circ_saT*0.7;  % Thoracic unestressed volume in arteries

    % Circ_saAB =  TotBV*0.1105; % Arteries Abdomen volume TVol*SA(0.2)*Upper body(0.85)*thorax(0.65)  from Benekin0.27
    % saAB_UN   =  Circ_saAB*0.7;% Abdomen unestressed volume in arteries

    % Circ_svAB =  TotBV*0.3867; % Veins Abdomen volume TVol*SA(0.7)*Upper body(0.85)*thorax(0.65)  from Benekin
    % svAB_UN   =  Circ_svAB*0.92; % Abdomen unestressed volume in veins

    % Circ_svT =  TotBV*0.2082;  % Veins Thorax volume TVol*SV(0.7)*Upper body(0.85)*thorax(0.35)  from Benekin
    % svT_UN   =  Circ_svT*0.92; % Thorax unestressed volume in veins

    % Circ_saL  =  TotBV*0.03;   % Arteries lower body volume TVol*SA(0.2)*Lower body(0.15)  from Benekin 
    % saL_UN    =  Circ_saL*0.7; % Lower body unestressed volume in arteries

    % Circ_svL  =  TotBV*0.105;  % Veins lower body volume TVol*SV(0.7)*Lower body(0.15)  from Benekin   % SV volume from Benekin 0.08
    % svL_UN    =  Circ_svL*0.92;% Lower body unestressed volume in veins

    % % Total circulating (stressed volume)
    % CircBV =  Circ_pa + Circ_pu + Circ_saT+Circ_saAB+Circ_saL+Circ_svT+Circ_svAB+Circ_svL;
    % % TotBV
    % % pause

    % global V_un
    % V_un = [saT_UN,pa_UN,pu_UN,saAB_UN,svAB_UN,svT_UN,saL_UN,svL_UN];


    % %%
    % % Flows
    % qtu  = (TotFlow*0.85)*0.35;      % Thoracic flow, arteries --> veins (Upper)
    % ql   = TotFlow*0.15;             % Lower body flow,, arteries--> veins
    % qab  = ((TotFlow*0.85)*0.65)-ql; % Abdominal flow, arteries--> veins
    % qtul  = qab;                     % Thoracic to abdominal flow
    % qul  = ql;                       % Upper body arteries --> lower body arteries
    % qlu  = ql;                       % Lower body veins --> upper body veins
    % qtlu = qtul;                     % Abdominal to thoracic flow

    % %% Resistences

    % % Peripheral resistances (Ohm's law)
    % R_t  = (P_TAdia-P_TVU)/qtu;  % Thoracic resistance
    % R_tul = (P_TAm-P_ABm)/qtul;  % Thoracic to abdomen resistance
    % R_ab  = (P_ABm-P_AB)/qab;    % Abdomen resistance
    % R_l   = (P_LAdia-P_VL)/ql;   % lower body resistance 
    % R_ul  = (P_TAm-P_LAm)/qul;   %Upper body to lower body 
    % R_lu  = (P_VL-P_AB)/qlu;     %Lower body to upper boddy
    % R_tlu = (P_AB-P_TVU)/qtlu;   %Abdomen to thorax resistance

    % % Valve resistances
    % R_mt  = 0.0025;        % Mitral valve resistance (mmHg*s/mL)
    % R_av  = 0.0025;        % Aortic valve resistance (mmHg*s/mL)
    % R_tc  = 0.0025;        % Tricuspid valve resistance (mmHg*s/mL)
    % R_pv  = 0.0025;        % Pulmonary valve resistance (mmHg*s/mL)

    % % Pumonary resistances
    % R_pul   = (P_PAdia - P_PV)/TotFlow ;        % Pulmonary vascular resistance (mmHg*s/mL) 
    % %% Elastances

    % % Pulmonary artery and vein parameters
    % E_pa = P_PAsys/(Circ_pa-pa_UN);                % PA artery elastance (mmHg/mL)  
    % E_pu = P_PAmean/(Circ_pu-pu_UN);               % PU elastance (mmHg/mL)  


    % % Systemic artery and vein parameters
    % E_ta    = P_TAsys/(Circ_saT-saT_UN);     % Thoracic arteries elastance (mmHg/mL)

    % E_tv    = P_TVU/(Circ_svT-svT_UN);       % Thoracic veins elastance

    % E_la    = P_LAsys/(Circ_saL-saL_UN);     % Lower body arteries elastance
    % E_vl    = P_VL/(Circ_svL-svL_UN);        % Lower body veins elastance

    % E_aba   = P_ABsys/(Circ_saAB-saAB_UN);   % Abdominal arteries elastance
    % E_abv   = P_AB/(Circ_svAB-svAB_UN);      % Abdominal veins elastance


    % %%
    % % Left ventricle free wall parameters

    % Ed     = P_LVdia/V_LVM;                     % Maximum elastance of the heart during diastole
    % Es     = P_LVsyst/(V_LVm-V_d_lvf);          % Maximum elastance of the heart during systole


    % % Right ventricle free wall parameters
    % Edr     = P_RVdia/V_RVM;                    % Maximum elastance of the heart during diastole
    % Esr     = P_RVsys/(V_RVm-V_d_rvf);          % Maximum elastance of the heart during systole




    % %% Respiration parameters
    % A_u= 3;
    % B_u=  0;

    % B_l= 0;
    % A_l=  -2;

    % Tinsp= 1.5;
    % Tc=8/3*(Tinsp);

    % %% Baroreceptor parameters
    % tauPm = 2.5; 
    % tauZ  = tauPm/10;

    % %% Autonomic parameters
    % kP  = 10;
    % TP0 = 0.8;
    % p2P = (pm^kP*(1-TP0)/TP0)^(1/kP);  
    % tauP = 2.5;

    % kS  = 7;
    % kSv= 10; 
    % TS0 = 0.2;
    % p2S = (pm^kS*TS0/(1-TS0))^(1/kS);
    % p2Sv = (pt^kSv*TS0/(1-TS0))^(1/kSv);

    % tauS = 12.5;
    % %left
    % Edm= 0.01*Ed;
    % EdX= Ed/(0.8*(TP0+Edm)+0.2*(TS0+Edm));
    % Eds= 0.2*EdX;
    % Edp=0.8*EdX;
    % %rigth
    % Edmr= 0.01*Edr;
    % EdXr= Edr/(0.8*(TP0+Edmr)+0.2*(TS0+Edmr));
    % Edsr= 0.2*EdXr;
    % Edpr=0.8*EdX;

    % %Resistances modulated by ANS

    % %Resistance thorax
    % Rtps = 0.8/TS0*R_tc;
    % Rtpm = 0.2*R_tc;

    % %Resistance abdomen
    % Rabps = 0.8/TS0*R_ab;
    % Rabpm = 0.2*R_ab;

    % %Resistance Lower Body Pars
    % Rlps = 0.8/TS0*R_l;
    % Rlbpm = 0.2*R_l;

    % %% Hear rate [Beregne HI fra p2H givet]
    % % delta=2.5;
    % tauH= 0.5;%delta*2.5; % probar 0.5
    % H0= 1.6; % (s)
    % Hp=0.5;
    % Hs= 0.3;


    % %% Vein lower
    % Vvl= Circ_svL-svL_UN; 
    % VMvl = 2*Vvl; %Choosen based on BP drop without controls 
    % mvl  = log(VMvl/(VMvl - Vvl))/P_VL;  %Cvl/(VM_vl-Vvl)

    % Vvab= Circ_svAB-svAB_UN; 
    % VMvab = 2*Vvab; %Choosen based on BP drop without controls 
    % mvab  = log(VMvab/(VMvab - Vvab))/P_AB;  %Cvl/(VM_vl-Vvl)