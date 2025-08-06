function [pars, init, taus] = load_global_easy()
  
init = containers.Map();  
pars = containers.Map();  
taus = containers.Map();
local_vars = containers.Map();  %this is just for the metaprogramming feature

%Notes:
% - pars('tauMR') and pars('tauMRv') defines how fast the input reaches its defined values, is important to change this value if the simulation scale changes.
% - pars('L_sa') and pars('L_pa') are related to inertance of flow, but it is too caotic, the authors had this values marked with red in supplementary material.
% - It's really important to check the measure units of flows, volumes and consumption. Brain consumption (MRbCO2: l/s) has different units than tissue consumption (MRtCO2: l/min) and different from Qpp, Qbp and Qla (ml/s)
%   Volumes of cardiovascular are in ml, Volumes in respiratory are in L.
%   respiratory pressures and resistances are in cmH20, for cardiovascular it is in mmHg (that includes Ptor and Pabd).

  

%Regarding volume and flows, everything should be in terms of ml/s, even though MRt are defined in l/min


%% ------------PARAMETERS:    There will be a description of every parameter in the form of Name |Value | Units |References | Description 
%Fitting hyperparameters
pars('settling_time') = 0;
pars('dt') = 0.01;
pars('tiny_dt') = 0.1;
%Time-related  

pars('T1') = 1;             % s              | [2]    | Time constant for cardiovascular mixing
pars('T2') = 2;             % s              | [2]    | Time constant for cardiovascular mixing
pars('tauMR') = 50;%50;%5;%50;         % s              | [0]    | Metabolic rate time constant     %maybe if we change tauMR we should change tauMRv too
pars('tauMRv') = 50;%5;%50;        % s              | [8]    | Metabolic rate time constant
pars('Ta') = 0.8;             % s              | [0]    | Time delay for gas exchange
% Pressure-related          
pars('gabd') = 3.39;        % mmHg/l         | [5]    | Constant gain factor linking tidal volume changes to abdominal pressure variations
pars('gthor') = 6.8;        % mmHg/l         | [5]    | Constant gain factor linking tidal volume changes to intrathoracic pressure variations
pars('Pabdmax') = 0;        % mmHg           | [5]    | Basal value of abdominal pressure at the end of expiration
pars('Pabdmin') = -2.5;     % mmHg           | [5]    | Basal value of abdominal pressure at the end of inspiration
pars('Pthormax') = -4;      % mmHg           | [5]    | Basal value of intrathoracic pressure at the end of expiration
pars('Pthormin') = -9;      % mmHg           | [5]    | Basal value of intrathoracic pressure at the end of inspiration
pars('Pcrit_min') = -40;    % cmH2O          | [2]    | Critical UPPER AIRWAY pressure
pars('Patm') = 719;%640;         % mmHg           | [0]    | Atmospheric pressure
pars('Pws') = 47;           % mmHg           | [1]    | Water vapor pressure
pars('Pao') = 0;            % cmH2O          | [1]    | Airway pressure
% Volume-related          
pars('VTn') = 0.73;         % l              | [5]    | Basal value of tidal volume
pars('VLO2') = 2.5;         % l              | [2]    | Lungs storage volume for O2
pars('VLCO2') = 3;          % l              | [2]    | Lungs storage volume for CO2
pars('Vtissue_CO2') = 15;   % l              | [8]    | Tissue storage volume for CO2
pars('Vtissue_O2') = 6;     % l              | [2]    | Tissue storage volume for O2

pars('Vdead') = 0.1587;     % l              | [3]    | Dead space volume
%Resistance-related          
pars('Rcw') = 0.8326;       % cmH2Os/l       | [0]    | Chest wall resistance
pars('Rl') = 1.366;         % cmH2Os/l       | [0]    | Lung transmural resistance
pars('Raw') = 0.82128;      % cmH2Os/l       | [0]    | Airway wall resistance
pars('Rrs') = 3.02;         % cmH2Os/l       | [6][7] | Overall resistance
pars('Rtrachea') = 10^6;    % cmH2Os/l       | [5]    | Upper airway wall resistance  
      
      
%Elastance-related         
pars('Cua') = 0.001;        % l/cmH2O        | [2]    | upper airways compliance
pars('Ecw') = 10.545;       % cmH2O/l        | [0]    | chest wall elastance
pars('El') = 10.545;        % cmH2O/l        | [0]    | lung transmural elastance
      
%Gas-related         
pars('fO2') = 21.0379;      %    %           | [1]    | Inspired fraction of O2
pars('fCO2') = 0.0421;      %    %           | [1]    | Inspired fraction of CO2
      
%Constants         
%UPPER AIRWAYS         
pars('bua') = 1;            % l              | [2]    | upper airways mechanic constant
pars('A0ua') = 1;           % --------       | [2]    | Maximum area of oppening in airway
pars('Kua') = 1;            % l/cmH2O/s      | [2]    | Proportionality coefficient
pars('kaw1')  = 1.85;       % cmH2O s/l      | [2]    | constant for upper airway pressure
pars('kaw2') = 0.43;        % cmH2O s^2/l^2      | [2]    | constant for upper airway pressure

%GAS         
pars('Z') = 0.0227;         % l/mmol         | [-]    | molar conversion factor 
pars('K1') = 13;            % mmHg           | [2]    | Parameter in O2 dissociation equation
pars('K2') = 194.4;         % mmHg           | [2][1] | Parameter in CO2 dissociation equation
pars('A1') = 0.3836;        % ----           | [2]    | Parameter in O2 dissociation equation (Bohr and Haldane effect)
pars('A2') = 1.819;         % ----           | [2]    | Parameter in CO2 dissociation equation (Bohr and Haldane effect)
pars('alpha1') = 0.03198;   % 1/mmHg         | [1][9] | Parameter in O2 dissociation equation  (Bohr and Haldane effect)
pars('alpha2') = 0.05591;   % 1/mmHg         | [1][2] | Parameter in CO2 dissociation equation (Bohr and Haldane effect)
pars('beta1') = 0.008275;   % 1/mmHg         | [9]    | Parameter in O2 dissociation equation (Bohr and Haldane effect)
pars('beta2') = 0.03255;    % 1/mmHg         | [2][1] | Parameter in CO2 dissociation equation (Bohr and Haldane effect)
pars('C1') = 9;             % mmol/l         | [2][1] | Max concentration of Hb bound O2
pars('C2') = 87;            % mmol/l         | [2][1] | Max concentration of CO2 
pars('LCTV') = 0.588;       % l              | [2]    | Lung to chemorreceptor transportation vascular volume constant

%Consumption-related         
pars('AT') = 1; %1.2 %1;             % l/min          | [0]    | Anaerobic treshold
pars('MRtCO2_basal') = 0.3;%0.05;%0.3; % l/min STPD     | [0]    | Basal production rate for CO2
pars('MRtO2_basal') = 0.12;%0.05;%0.3;%0.28;%0.33; % l/min STPD     | [0]    | Basal production rate for O2
pars('MRCO2') = pars('MRtCO2_basal');        % l/min STPD     | [0]    | Input metabolic production rate for CO2, starting with basal values
pars('MRO2') = pars('MRtO2_basal');        % l/min STPD     | [0]    | Input metabolic production rate for O2,  starting with basal values


%Brain-related
pars('KCSFCO2') = 320;      % s              | [3]    | CO2 difussion time constant of cerebrospinal fluid
pars('KCCO2') = 346000;     % s/cm^2/l       | [3]    | CO2 central receptor constant
pars('dc') = 0.015;         % cm             | [3]    | Depth of central receptor
pars('h') = 0.0183;         % ml/100g/s      | [2]    | Cerebral blood flow constant
pars('SCO2') = 0.0043;      % 1/mmHg         | [2]    | Dissociation slope for CO2 in brain
pars('SbCO2') = 0.36;       % ml 100/g/mmHg  | [2]    | Dissociation slope for CO2 in blood
pars('MRbCO2') = 0.0009;    % l/s STPD       | [3]    | Metabolic production rate in brain of CO2
pars('MRbO2') = 0.000925;   % l/s STPD       | [3]    | Metabolic production rate in brain of O2

%Ventilatory control-related
pars('KpCO2') = 0.225;%0.225;%0.255;%0.205;%0.35;%0.31;%0.2025;%0.35; %0.2025;%0.35;%0.2025; %TESTING 0.2832;%0.2025;       % mmHg^-1
pars('KcCO2') = 0.2432;%0.2432;%0.2832;%0.2332;%0.35;%0.31;%0.2332;%0.35; %0.2332;  %0.35;%0.2332; %TESTING 0.2255;%0.2332;       % mmHg^-1
pars('KpO2') = 2 * 10^(-9); %2 * 10^(-10);%4.72 %1.72 * 10^(-10);%4.72 * 10-9 %* 10^-8.75; %* 10^-8.2 %10^-9;   % mmHg^(-4.9) 
pars('KcMRv') = 0.8;%1;%0.8;%0.35;%.35; %0.35; % 1           %--
pars('Kbg') = 17.4;%17.4;           %--
pars('GVdead') = 0.1698;      %s  [3]
pars('dV_rest') = 0.0673;%0.0673;%0.0773; %0.|0673;     %lt/s

%Optimizer control - related
pars('n') = 1.101;            %--
pars('lambda1') = 0.86;       %-- 
pars('lambda2') = 0.489;      %-- 
pars('Pmax') = 50;%100;           %cmH2O
pars('dPmax') = 1000;         %cmH2O/s

%cycle pars
pars('t0') = 0;


%Cardiovascular pars
pars('C_sa') = 1.5;%1.2;%0.28;%0.9;%0.28;               % ml/mmHg      | [2][5] | Systemic arterial compliance
pars('L_sa') = 0.22;%0.22; %* 10^-3;       % mmHgs^2/ml   | [2][5] | Systemic arterial inertance
pars('R_sa') = 0.15; %0.3;%0.06;%1;%0.1;%0.06;               % mmHgs^2/ml   | [2][5] | Systemic arterial hydraulic resistance
pars('V_unstressed_sa') = 0;      % ml           | [2][5] | Systemic arterial unstressed volume

pars('k_r_am') = 24.17;            % s/ml         | [5]    | Constant parameter
pars('P0') = 3.93;                 % mmHg         | [5]    | Constant parameter
pars('C_e_p') = 0.668;             % ml/mmHg      | [2][5] | Extra-splanchnic compliance
pars('C_s_p;') = 2.05;             % ml/mmHg      | [2][5] | Splanchnic peripheral compliance
pars('C_b_p') = 0.668;             % ml/mmHg      | [2][5] | Brain peripheral compliance
pars('C_h_p') = 0.119;             % ml/mmHg      | [2][5] | Coronary peripheral compliance
pars('C_rm_p') = 0.21;             % ml/mmHg      | [5]    | Resting skeletal muscle peripheral compliance
pars('C_am_p') = 0.315;            % ml/mmHg      | [5]    | Active skeletal muscle peripheral compliance
pars('C_e_v') = 20;                % ml/mmHg      | [2]    | Extra-splachnic venous compliance
pars('C_s_v') = 61.1;              % ml/mmHg      | [2]    | Splachnic venous compliance
pars('C_b_v') = 10.71;             % ml/mmHg      | [2][5] | Brain venous compliance
pars('C_h_v') = 3.57;              % ml/mmHg      | [2]    | Coronary venous compliance
pars('C_rm_v') = 6.28;             % ml/mmHg      | --     | Resting skeletal muscle venous compliance
pars('C_am_v') = 9.4;              % ml/mmHg      | --     | Active skeletal muscle venous compliance
pars('R_e_n') = 0.24623; %0.04;              % mmHgs/ml     | [2]    | Extra-splanchnic venous resistance
pars('R_s_n') = 0.14072; %0.038;             % mmHgs/ml     | [2]    | Splanchnic venous resistance
pars('R_b_n') = 0.20629;%0.075;             % mmHgs/ml     | [2][5] | Brain venous resistance
pars('R_h_n') = 0.5030; %0.224;             % mmHgs/ml     | [2][5] | Coronary venous resistance
pars('R_rm_n') = 0.0305; %0.125;            % mmHgs/ml     | [5]    | Resting skeletal muscle venous resistance
pars('R_am_n') = 0.03054; %0.0833;           % mmHgs/ml     | [5]    | Active skeletal muscle venous resistance
pars('V_unstressed_e_p') = 128.931; %127.72; % ml           | --     | Extra-splachnic peripheral unstressed volume
pars('V_unstressed_s_p') = 109.394; %260.3; % ml           | --     | Splachnic peripheral unstressed volume
pars('V_unstressed_b_p') = 68.480; %68.42;  % ml           | --     | Brain peripheral unstressed volume
pars('V_unstressed_h_p') = 22.17444;%23;     % ml           | --     | Coronary peripheral unstressed volume
pars('V_unstressed_rm_p') = 38.6584; %40.1;  % ml           | --     | Resting skeletal muscle peripheral unstressed volume
pars('V_unstressed_am_p') = 57.9877; %60.22; % ml           | --     | Active skeletal muscle peripheral unstressed volume

pars('V_unstressed_e_v') = 50;    %this is actually a var from carfiac control
pars('V_unstressed_s_v') = 50;    %this is actually a var from carfiac control
pars('V_unstressed_b_v') = 275.005; %279.49; % ml           | --     | Brain venous unstressed volume
pars('V_unstressed_h_v') = 87.3010; %93.16;  % ml           | --     | Coronary venous unstressed volume
pars('V_unstressed_rm_v') = 50;  %this is actually a var from carfiac contro = l
pars('V_unstressed_am_v') = 7;  %this is actually a var from carfiac contro = l
pars('V_tot') = 5027.6;            % ml           | --    | Total blood volume

pars('D1') = 0.3855;                % mmHg                 | [2]        | Parameter for P-V curve of vena cava
pars('D2') = -5;                    % mmHg                 | [2]        | Parameter for P-V curve of vena cava
pars('K1_vc') = 0.15;               % mmHg                 | [2]        | Parameter for P-V curve of vena cava
pars('K2_vc') = 0.4;                % mmHg                 | [2]        | Parameter for P-V curve of vena cava
pars('K_r_vc') = 0.001;             % mmHgs/ml             | [2]        | Gain for vena cava flow resistance
pars('R_vc_n') = 0.025;             % mmHgs/ml             | [2]        | Nominal vena cava flow resistance
pars('V_unstressed_vc') = 154.279; %123;      % ml                   | [2]        | Vena cava unstressed volume
pars('V_vc_max') = 350;             % ml                   | [2]        | Maximum volume of vena cava
pars('V_vc_min') = 50;              % ml                   | [2]        | Minimum volume of vena cava

pars('C_pa') = 0.76;                % ml/mmHg              | [2][5]     | Pulmonary arterial compliances
pars('C_pp') = 5.8;                 % ml/mmHg              | [2][5]     | Pulmonary peripheral compliances
pars('C_pv') = 25.37;               % ml/mmHg              | [2][5]     | Pulmonary venous compliances
pars('L_pa') = 0.18;%0.18;%0.00018;             % mmHg^2/ml            | [2][5]     | Pulmonary arterial intertance
pars('R_pa') = 0.023;               % mmHgs/ml             | [2][5]     | Pulmonary arterial flow resistance
pars('R_pp') = 0.24266; %0.0894;              % mmHgs/ml             | [2][5]     | Pulmonary peripheral flow resistance
pars('R_pv') = 0.0298; %0.0056;              % mmHgs/ml             | [2][5]     | Pulmonary venous flow resistance
pars('V_unstressed_pa') = 0;        % ml                   | [2][5]     | Pulmonary arterial unstressed volume
pars('V_unstressed_pp') = 94.9390; %116.6775; % ml                   | [2]        | Pulmonary peripheral unstressed volume
pars('V_unstressed_pv') = 114.393; %114;      % ml                   | [2]        | Pulmonary venous unstressed volume
    



pars('C_la') = 19.3;%5.23; %19.23;               % ml/mmHg              | [2]        | Left atrial compliances
pars('C_ra') = 31.25;               % ml/mmHg              | [2]        | Right atrial compliances
pars('K_E_lv') = 0.014;             % ml^-1                | [10]       | Parameters that describe the end-diastolic pressure-volume relationship in the left ventricle
pars('K_E_rv') = 0.011;             % ml^-1                | [10]       | Parameters that describe the end-diastolic pressure-volume relationship in the left ventricle
pars('KR_lv') = 0.000375;           % s/mmHg               | [10]       | Parameters that describe the end-diastolic pressure-volume relationship in the left ventricle
pars('KR_rv') = 0.0014;             % s/mmHg               | [10]       | Parameters that describe the end-diastolic pressure-volume relationship in the left ventricle
pars('ksys') = 0.075;               % s^2                  | [10]       | Constant parameter which describes the duration of systole as a function of heart rate
pars('P0_lv') = 1.5;                % mmHg                 | [10]       | Parameters that describe the end-diastolic pressure-volume relationship in the left ventricle
pars('P0_rv') = 1.5;                % mmHg                 | [10]       | Parameters that describe the end-diastolic pressure-volume relationship in the left ventricle
pars('R_la') = 0.0025;              % mmHgs/ml             | [2]        | Left atrial flow resistance
pars('R_ra') = 0.0025;              % mmHgs/ml             | [2]        | Right atrial flow resistance
pars('Tsys_0') = 0.5;               % s                    | [10]       | Constant parameter which describes the duration of systole as a function of heart rate
pars('V_unstressed_la') = 24;       % ml                   | --         | Left atrial unstressed volume
pars('V_unstressed_lv') = 15.908;   % ml                   | --         | Left ventricular unstressed volume
pars('V_unstressed_rv') = 24;       % ml                   | --         | Right atrial unstressed volume
pars('V_unstressed_ra') = 38.703;   % ml                   | --         | Right ventricular unstressed volume


pars('Aim') = 50;                   % mmHg                 | [5]        | Peak value of intramuscular pressure
pars('Tc') = 0.75;                  % s                    | [5]        | The overall duration of muscular contraction
pars('Tim') = 1;                    % s                    | [5]        | Duration of the muscular contraction-relaxation cycle
    
%Autonomic control pars

%aferent barorreflex
pars('f_ab_min') = 2.52;            % spikes/s             | [2][11]    | Lower saturation level of the frequency discharge in the baroreceptor afferent fibers
pars('f_ab_max') = 47.78;           % spikes/s             | [2][11]    | Upper saturation level of the frequency discharge in the baroreceptor afferent fibers
pars('kab') = 11.76;                % mmHg                 | [2][11]    | Parameter related to the slope of the static function at the central point
pars('P_n') = 92;                   % mmHg                 | [2][11]    | Value of baroreceptor pressure at the central point of the sigmoidal function
pars('tau_p') = 2.076;              % s                    | [2][11]    | Time constant for the real pole
pars('tau_z') = 6.37;               % s                    | [2][11]    | Time constant for the real zero

%aferent chemorreflex
pars('f_ac_CO2_n') = 1.4;           % --                   | [2][4]     | Constant parameter tuned to reproduce the CO2 static response
pars('f_ac_max') = 12.3;            % spikes/s             | [2][4]     | Upper saturation level of the frequency discharge in the chemoreceptor afferent fibers
pars('f_ac_min') = 0.835;           % spikes/s             | [2][4]     | Lower saturation level of the frequency discharge in the chemoreceptor afferent fibers 
pars('kac') = 29.27;                % mmHg                 | [2][4]     | Parameter related to the slope of the sigmoid at the central point
pars('KH') = 3;                     % --                   | [2][4]     | Constant parameter tuned to reproduce the CO2 static response
pars('PaO2_ac_n') = 45;             % mmHg                 | [4]        | Arterial PO2 at the central point of the sigmoid
pars('PaCO2_n') = 40;               % mmHg                 | [2][4]     | PaCO2 basal value
pars('tau_ac') = 2;                 % s                    | [2][4]     | Time constant of the chemoreceptor mechanism

%stretch receptors
pars('G_ap') = 11.76;               % spikes/s/l           | [1]        | Constant gain factor
pars('tau_ap') = 2;                 % s                    | [2]        | Time constant of the lung inflation afferent response

%blood flow local
pars('A') = 20.9;                   % --                   | [2][4]     | Constant parameter
pars('B') = 92.8;                   % --                   | [2][4]     | Constant parameter
pars('C') = 10570;                  % --                   | [2][4]     | Constant parameter
pars('D') = -5.251;                 % --                   | [4]        | Constant parameter
pars('vO2_b_n') =  0.14;            % --                   | [2][4]     | O2 concentration in venous blood leaving the brain under normal conditions
pars('gO2_b') = 10;                 % --                   | [2][4]     | Constant gain factor
pars('MO2_bp') = 0.7917;            % ml/s                 | [11]       | Oxygen consumption rate in the brain compartment
pars('R_bmp') = 6.57;               % mmHg s/ml            | [2]        | Constant parameter denoting the basal value of peripheral cerebrovascular conductance
pars('tau_CO2') = 20;               % s                    | [2][4]     | Time constants of the effect of CO2 on cerebral circulation
pars('tau_O2') = 10;                % s                    | [2][4]     | Time constants of the effect of O2 on cerebral circulation

%coronoary and resting muscle blood flow
pars('vO2_h_n') = 0.11;            % --                   | [2][4]     | O2 concentration in venous blood leaving the heart under normal conditions
pars('vO2_rm_n') = 0.155;          % --                   | [2][4]     | O2 concentration in venous blood leaving the skeletal resting muscle under normal conditions        
pars('gO2_h') = 35;                % --                   | [2][4]     | Constant gain factor
pars('gO2_rm') = 30;               % --                   | [2][4]     | Constant gain factor
pars('KCO2_h') = 11.11;            % mmHg                 | [4]        | Parameter related to the slope of the sigmoidal function at the central point
pars('KCO2_rm') = 142.8;           % mmHg                 | [4]        | Parameter related to the slope of the sigmoidal function at the central point
pars('MO2_h_p_n') = 0.4;           % ml/s                 | [11]       | Nominal value of O2 consumption rate in the heart
pars('MO2_rm_p') = 0.86;           % ml/s                 | [11]       | Consumption rate in the resting muscle
pars('R_h_p_n') = 24.0225; %19.71;           % mmHgs/ml             | [2][4]     | Normal peripheral resistance in coronary compartment
pars('tau_w') = 5;                 % s                    | [11]       | Time constant of the filter
pars('Whn') = 12660;               % mmHg/sml             | [11]       | Nominal value of the average power of the cardiac pump

%active muscle flow
pars('vO2_am_n') = 0.1555;         % --                   | [5]        | O2 concentration in venous blood leaving the heart under normal conditions
pars('delay_met') = 4;             % s                    | [5]        | Pure delay
pars('gO2_am') = 30;               % --                   | [5]        | Constant gain factor
pars('g_M') = 40;                  % --                   | [5]        | Static gain
pars('I0_met') = 0.4266;           % --                   | [5]        | Is I at the central point of the sigmoid
pars('kmet') = 0.18;               % --                   | [5]        | Parameter related to the slope of the sigmoid at the central point
pars('MO2_am_p_n') = 0.516;        % ml/s                 | [5]        | Nominal oxygen consumption rate
pars('phi_max') = 20;              % --                   | [5]        | Upper saturation of the static sigmoidal characteristic
pars('phi_min') = -1.87;           % --                   | [5]        | Lower saturation of the static sigmoidal characteristic 
pars('tau_M') = 40;                % s                    | [5]        | Time constant
pars('tau_met') = 10;                % s                    | [5]        | Time constant

%CNS Ischemic response (integration)
pars('gcc_h_s') = 1;               % mmHg^-1s^-1          | [4]        | Constant gain factor tuned to reproduce experimental results
pars('gcc_p_s') = 1.5;             % mmHg^-1s^-1          | [4]        | Constant gain factor tuned to reproduce experimental results
pars('gcc_v_s') = 0;               % mmHg^-1s^-1          | [4]        | Constant gain factor tuned to reproduce experimental results  
pars('k_isc_h_s') = 6;             % mmHg                 | [4]        | Parameter related to the slope of the static function at the central point for heart
pars('k_isc_p_s') = 2;             % mmHg                 | [2][4]     | Parameter related to the slope of the static function at the central point for peripheral resistance
pars('k_isc_v_s') = 2;             % mmHg                 | [2][4]     | Parameter related to the slope of the static function at the central point for unstressed volume of veins
pars('PO2_ref_h_s') = 45;          % mmHg                 | [2][4]     | Value of PO2 at the central point of the sigmoidal function for heart
pars('PO2_ref_p_s') = 30;          % mmHg                 | [2][4]     | Value of PO2 at the central point of the sigmoidal function for peripheral resistance
pars('PO2_ref_v_s') = 30;          % mmHg                 | [2][4]     | Value of PO2 at the central point of the sigmoidal function for unstressed volume of veins
pars('tau_cc') = 20;               % s                    | [2][4]     | Time constant
pars('tau_isc') = 30;              % s                    | [2][4]     | Time constant of the mechanism
pars('Theta_h_s_n') = 3.6;         % spike/s              | [2][4]     | Offset term in basal condition for heart
pars('Theta_p_s_n') = 13.32;       % spike/s              | [2][4]     | Offset term in basal condition for peripheral resistance
pars('Theta_v_s_n') = 13.32;       % spike/s              | [2][4]     | Offset term in basal condition for unstressed volume of veins
pars('x_h_s') = 20; %20 %53;                % spike/s              | [4]        | Saturation of the hypoxic response for heart 
pars('x_p_s') = 6;                 % spike/s              | [2][4]     | Saturation of the hypoxic response for peripheral resistance
pars('x_v_s') = 6;                 % spike/s              | [2][4]     | Saturation of the hypoxic response for unstressed volume of veins

%Efferent pathways (symp and vagal)     
pars('fab_0') = 25;                % spike/s              | [4]        | Central value in the curve of fab
pars('fes_0') = 16.11;             % spike/s              | [5]        | Constant parameter
pars('fes_inf') = 2.1;             % spike/s              | [5]        | Constant parameter
pars('fes_max') = 60;              % spike/s              | [5]        | Saturation level above which the sympathetic activity cannot increase
pars('fev_0') = 3.2;               % spike/s              | [5]        | Constant parameter
pars('fev_inf') = 6.3;             % spike/s              | [5]        | Constant parameter
pars('kes') = 0.0675;              % spike/s              | [5]        | Constant parameter
pars('kev') = 7.06;                % spike/s              | [5]        | Constant parameter
pars('I_0_h_s') = 0.658;           % spike/s              | [5]        | Value of exercise intensity at the central point of the sigmoid
pars('I_0_p_s') = 0.65;            % spike/s              | [5]        | Value of exercise intensity at the central point of the sigmoid
pars('I_0_v_s') = 0.45;            % spike/s              | [5]        | Value of exercise intensity at the central point of the sigmoid
pars('I_0_v') = 0.126;             % spike/s              | [5]        | Value of exercise intensity at the central point of the sigmoid
pars('kcc_h_s') = 0.114;           % spike/s              | [5]        | Parameter related to the slope of the characteristic at the central point
pars('kcc_p_s') = 0.13;            % spike/s              | [5]        | Parameter related to the slope of the characteristic at the central point
pars('kcc_v_s') = 0.09;            % spike/s              | [5]        | Parameter related to the slope of the characteristic at the central point
pars('kcc_v') = 0.0162;            % spike/s              | [5]        | Parameter related to the slope of the characteristic at the central point
pars('gamma_h_s_max') = 9;         % spike/s              | [5]        | Upper saturation of the central command response
pars('gamma_p_s_max') = 5.5;       % spike/s              | [5]        | Upper saturation of the central command response
pars('gamma_v_s_max') = 64.9;      % spike/s              | [5]        | Upper saturation of the central command response
pars('gamma_v_max') = 1.9;         % spike/s              | [5]        | Upper saturation of the central command response
pars('gamma_h_s_min') = -0.0283;   % spike/s              | [5]        | Lower saturation of the central command response
pars('gamma_p_s_min') = -0.037;    % spike/s              | [5]        | Lower saturation of the central command response
pars('gamma_v_s_min') = -0.437;    % spike/s              | [5]        | Lower saturation of the central command response
pars('gamma_v_min') = -0.0008;     % spike/s              | [5]        | Lower saturation of the central command response
pars('Theta_v') = -0.68;           % spike/s              | [11]       | Offset term
pars('Wb_h_s') = -1.2671; %-1.75;            % spike/s              | [1]        | Synaptic weight tuned to reproduce physiological results |change by my analysis of HR and Psa (01-02-25)|
pars('Wb_p_s') = -1.1375;          % spike/s              | [1]        | Synaptic weight tuned to reproduce physiological results
pars('Wb_v_s') = -1.1375;          % spike/s              | [0]        | Synaptic weight tuned to reproduce physiological results
pars('Wc_h_s') = 1;                % spike/s              | [4]        | Synaptic weight tuned to reproduce physiological results
pars('Wc_p_s') = 1.716;            % spike/s              | [1]        | Synaptic weight tuned to reproduce physiological results
pars('Wc_v_s') = 1.716;            % spike/s              | [1]        | Synaptic weight tuned to reproduce physiological results
pars('Wc_v') = 0.2;                % spike/s              | [11]       | Synaptic weight tuned to reproduce physiological results
pars('Wp_h_s') = 0;                % spike/s              | [5]        | Synaptic weight tuned to reproduce physiological results
pars('Wp_p_s') = -0.3997;          % spike/s              | [1]        | Synaptic weight tuned to reproduce physiological results
pars('Wp_v_s') = -0.3997;          % spike/s              | [0]        | Synaptic weight tuned to reproduce physiological results
pars('Wp_v') = -0.103;             % spike/s              | [5]        | Synaptic weight tuned to reproduce physiological results
pars('Wt_h_s') = 0.4;              % spike/s              | [2]        | Synaptic weight tuned to reproduce physiological results
pars('Wt_p_s') = 0.4;              % spike/s              | [2]        | Synaptic weight tuned to reproduce physiological results
pars('Wt_v_s') = 0.4;              % spike/s              | [2]        | Synaptic weight tuned to reproduce physiological results
pars('Wt_v') = 0.4;                % spike/s              | [2]        | Synaptic weight tuned to reproduce physiological results

%Effectors: Rp, Vu, E 
pars('delay_Emax_lv') = 2;         % s                    | [2][11]    | Pure latency of the mechanism
pars('delay_Emax_rv') = 2;         % s                    | [2]        | Pure latency of the mechanism
pars('delay_R_am_p') = 2;          % s                    | [2]        | Pure latency of the mechanism
pars('delay_R_e_p') = 2;           % s                    | [2]        | Pure latency of the mechanism
pars('delay_R_rm_p') = 2;          % s                    | [2]        | Pure latency of the mechanism
pars('delay_R_s_p') = 2;           % s                    | [2]        | Pure latency of the mechanism
pars('delay_V_u_am_v') = 5;        % s                    | [2]        | Pure latency of the mechanism
pars('delay_V_u_e_v') = 5;         % s                    | [2]        | Pure latency of the mechanism
pars('delay_V_u_rm_v') = 5;        % s                    | [2]        | Pure latency of the mechanism
pars('delay_V_u_s_v') = 5;         % s                    | [2]        | Pure latency of the mechanism
pars('Emax_lv_0') = 2.392;         % mmHg/ml                 | [11]       | Basal level of maximum end-systolic elastance of the left ventricle
pars('Emax_rv_0') = 1.412;         % mmHg/ml                 | [11]       | Basal level of maximum end-systolic elastance of the right ventricle
pars('R_am_p_0') = 4.7088; %3.510;          % mmHgs/ml             | [0]        | Basal level of active skeletal peripheral resistance
pars('R_e_p_0') = 3.09159; %1.655;           % mmHgs/ml             | [11]       | Basal level of extra-splanchnic peripheral resistance
pars('R_rm_p_0') = 4.708827; %5.270;          % mmHgs/ml             | [0]        | Basal level of resting skeletal peripheral resistance
pars('R_s_p_0') = 2.4986; %2.49;            % mmHgs/ml             | [11]       | Basal level of splanchnic skeletal peripheral resistance 
pars('V_u_am_v_0') = 258.668;%286.4;        % ml                   | [0]        | Basal level of active skeletal muscle venous unstressed volume
pars('V_u_rm_v_0') = 172.445; %190.95;       % ml                   | [0]        | Basal level of resting skeletal muscle venous unstressed volume
pars('V_u_e_v_0') = 618.31; %607.8;         % ml                   | [0]        | Basal level of extra-splanchnic venous unstressed volume
pars('V_u_s_v_0') = 606.888; %1361.6;        % ml                   | [0]        | Basal level of splanchnic venous unstressed volume
pars('G_Emax_lv') = 0.475;         %mmHg/ml/v             | [11]       | Constant gain factor
pars('G_Emax_rv') = 0.282;         %mmHg/ml/v             | [11]       | Constant gain factor
pars('G_R_am_p') = 2.47;           %mmHg/ml/v             | [11]       | Constant gain factor
pars('G_R_e_p') = 1.94;            %mmHg/ml/v             | [11]       | Constant gain factor
pars('G_R_rm_p') = 2.47;           %mmHg/ml/v             | [11]       | Constant gain factor
pars('G_R_s_p') = 0.695;           %mmHg/ml/v             | [11]       | Constant gain factor
pars('G_V_u_am_v') = -58.29;       %ml/v                  | [11]       | Constant gain factor
pars('G_V_u_e_v') = -74.21;        %ml/v                  | [11]       | Constant gain factor
pars('G_V_u_rm_v') = -58.29;       %ml/v                  | [11]       | Constant gain factor
pars('G_V_u_s_v') = -265.4;        %ml/v                  | [11]       | Constant gain factor
pars('tau_Emax_lv') = 8;           %s                     | [11]       | Time constant
pars('tau_Emax_rv') = 8;           %s                     | [11]       | Time constant
pars('tau_R_am_p') = 2;            %s                     | [2]        | Time constant
pars('tau_R_e_p') = 2;             %s                     | [2]        | Time constant
pars('tau_R_rm_p') = 2;            %s                     | [2]        | Time constant
pars('tau_R_s_p') = 2;             %s                     | [2]        | Time constant
pars('tau_V_u_am_v') = 20;         %s                     | [11]       | Time constant
pars('tau_V_u_e_v') = 20;          %s                     | [11]       | Time constant
pars('tau_V_u_rm_v') = 20;         %s                     | [11]       | Time constant
pars('tau_V_u_s_v') = 20;          %s                     | [11]       | Time constant
pars('fes_min') = 2.66;            %spike/s               | [11]       | Threshold for sympathetic stimulation

pars('delay_Tsym') = 2;            %s                     | [2][11]    | Pure latency of the mechanism
pars('delay_Tvagal') = 0.2;        %s                     | [2][11]    | Pure latency of the mechanism
pars('GTsym') = -0.18; %-0.13; %-0.18;%-0.13;             %--                    | [2][11]    | Constant gain factor  %this parameter is too sensitive!!
pars('GTvagal') = 0.0652; %0.09;            %--                    | [2][11]    | Constant gain factor  |change by my analysis of HR and Psa (01-02-25)|
pars('T0') = 0.4;%0.425; %0.58;%0.4;%0.58;                 %s                     | [2][11]    | Heart period in the absence of cardiac inervation, |change by my analysis of HR and Psa (01-02-25)|
pars('tau_Tsym') = 2;              %s                     | [2][11]    | Time constant
pars('tau_Tvagal') = 1.5;          %s                     | [2][11]    | Time constant

%For fitting
pars('MRO2_poly_0') = 0;
pars('MRO2_poly_1') = 0;
pars('MRO2_poly_2') = 0;
pars('MRO2_poly_3') = 0;
pars('MRO2_poly_4') = 0;
pars('MRO2_poly_5') = 0;
pars('MRO2_poly_6') = 0;
pars('MRO2_poly_7') = 0;
pars('MRO2_poly_8') = 0;

pars('MRCO2_poly_0') = 0;
pars('MRCO2_poly_1') = 0;
pars('MRCO2_poly_2') = 0;
pars('MRCO2_poly_3') = 0;
pars('MRCO2_poly_4') = 0;
pars('MRCO2_poly_5') = 0;
pars('MRCO2_poly_6') = 0;
pars('MRCO2_poly_7') = 0;
pars('MRCO2_poly_8') = 0;

pars('fiO2_poly_0') = 0;
pars('fiO2_poly_1') = 0;
pars('fiO2_poly_2') = 0;
pars('fiO2_poly_3') = 0;
pars('fiO2_poly_4') = 0;


pars('vO2_e_n') = 0.13;
pars('vO2_s_n') = 0.13;
pars('aO2_n') = 0.9;

pars('gO2_e') = 1*pars('gO2_rm'); %0.5*
pars('gO2_s') = 1*pars('gO2_rm'); %0.5 * 
pars('gO2_p') = 0.1 * pars('gO2_rm');%

pars('MO2_e') = 0.6;
pars('MO2_s') = 0.6;
pars('MO2_p') = 0;

pars('R_p_p_n') = 0.24266;
pars('Hgt') = 170;
pars('BW') = 70;
pars('Gender') = 2;
 






%% --------INITIAL CONDITIONS: There will be a description of every initial condition in the form of Name |Value | Units |References | Description 

%Blood-related
init('Qpp') = 100; %1000;           % ml/s| []  | blood flow in lungs, but in the end it is really ml/s because it depends on Rpp
init('Qbp') = 10;%1000*0.1;        % ml/s  (it seems that this quantity has to be in lt/s in another part)| []  | blood flow in brain (the same here, compared with what happens above: ml/s)
init('Qla') = 255; %0.225;%0.01;             % ml/s | []  | blood flow in left atrium


% Time-related
init('TI') = 1.6;%1;%1;%3;             % s    | []  | Inspiration time ,its important to recall that TI of 1 at least in their research was for maximum excercise.
init('Tresp') = 3.5;%3;%3;%5.6;        % s    | []  | Respiration time
init('TE') = init('Tresp') - init('TI');

init('fake_TI') = init('TI');
init('fake_Tresp') = init('Tresp');

% Volume-related
init('VT') = 0;           % l    | []  | Tidal volume
init('dV') = 0;            % l/s  | []  | Flow
init('dVua') = 0;          % l    | []  | Upper airways flow volume
init('V') = 0;             % l    | []  | Total volume
init('dVE') = 0;           % l/s  | []  | Minute ventilation, inspired volume in a minute (it must be a positive value)



%Pressure-related
init('Ppl') = 0;           % mmHg | []  | Pleural pressure
init('Pmusc') = 0;         % mmHg | []  | diafragmatic pressure
init('Pua') = 0;

%Gains
init('Gaw') = 0;           %   ?  | []  | Airflow's gain factor
init('Nt') = 0;            %   ?  | []  | Central respiratory neuromuscular drive response
%init('RR') = 0;            %   ?  | []  | Respiratory rythm

%Control optimizer
init('a0')  = 0;           %   ?  | []  | Control parameter to optimize (asociated with pmusc curve shape)
init('a1') = 30; %8;%40;%12;          %   ?  | []  | Control parameter to optimize (asociated with pmusc curve shape)
init('a2') = -3;%-1;            %   ?  | []  | Control parameter to optimize (asociated with pmusc curve shape)
init('t1') = init('TI');   %   ?  | []  | Control parameter to optimize (inspiration time)
init('t2') = init('TE');   %   ?  | []  | Control parameter to optimize (expiration time)
init('tau') = 0.3;           %   ?  | []  | Control parameter to optimize (decaying time constant)
init('insp_integrand') = 0;
init('exp_integrand') = 0;
init('J') = 0;
init('insp_work_power') = 0;
init('exp_work_power') = 0;

%Partial pressure - related
init('vO2') = 0.08; %0.1339; %0.1639;      %   mmol/l  | []  | venous O2 concentration
init('vCO2') = 0.5247;     %   mmol/l  | []  | venous CO2 concentration
init('aO2') = 0.2;%0.9;         %   ?  | []  | arterial O2 concentration
init('aCO2') = 0.1;        %   ?  | []  | arterial CO2 concentration
init('PaO2') = 103.1435;   %   ?  | []  | arterial O2 partial pressure
init('PaCO2') = 40.3928;   %   ?  | []  | aeterial CO2 partial pressure
init('dPaO2') = 0;         %   ?  | []  | arterial O2 partial pressure derivative
init('dPaCO2') = 0;        %   ?  | []  | arterial CO2 partial pressure derivative
init('PvbCO2') = 40;        %   ?  | []  | venous brain CO2 partial pressure
init('PCSFCO2') = 40;       %   ?  | []  | CO2 cerebroespinal fluid partial pressure
init('PbCO2') = 40; %36; %40;   %   ?  | []  | total brain CO2 partial pressure
init('P_1O2') = 104.3637;  %   ?  | []  | number 1 compartiment O2 partial pressure
init('P_1CO2') = 39.5616;  %   ?  | []  | number 1 compartiment CO2 partial pressure
init('P_2O2') = 104.2258;  %   ?  | []  | number 2 compartiment O2 partial pressure
init('P_2CO2') = 39.6736;  %   ?  | []  | number 2 compartiment CO2 partial pressure
init('P_3O2') = 104.0505;  %   ?  | []  | number 3 compartiment O2 partial pressure
init('P_3CO2') = 39.8127;  %   ?  | []  | number 3 compartiment CO2 partial pressure
init('P_4O2') = 103.8005;  %   ?  | []  | number 4 compartiment O2 partial pressure
init('P_4CO2') = 40.0061;  %   ?  | []  | number 4 compartiment CO2 partial pressure
init('P_5O2') = 103.3579;  %   ?  | []  | number 5 compartiment O2 partial pressure
init('P_5CO2') = 40.3359;  %   ?  | []  | number 5 compartiment CO2 partial pressure
init('PAO2') = 92;%110;%88;%102.5153;%88; %;   %   ?  | []  | alveolar O2 partial pressure
init('PACO2') = 40.9432;   %   ?  | []  | alveolar CO2 partial pressure
init('mean_PaO2') = init('PaO2');     %   ?  | []  | mean value of PaO2
init('mean_PaCO2') = init('PaCO2');    %   ?  | []  | mean value of PaCO2
init('mean_PbCO2') = init('PbCO2');    %   ?  | []  | mean value of PbCO2



%Consumption-related
init('MRtO2') = pars('MRtO2_basal');%0.9;       %   ?  | []  | metabolic O2 production/consumption rate in tissues
init('MRtCO2') = pars('MRtCO2_basal');%0.9;      %   ?  | []  | metabolic CO2 production/consumption rate in tissues
init('M_Rv') = 0;          %   ?  | []  | fundamental part of MRv variable
init('MRv') = 1.2;           %   ?  | []  | Metabolic rate neural drive to ventilation
init('MRR') = 0;           %   ?  | []  | Metabolic rate ratio
init('I') = 0;             %   ?  | []  | Anaerobic excersice relative intensity

%Systemic arteries

%Pressure      %the same as below, we will have to set it depending on if its systole or diastole
init('P_sa') = 80;  % systemic arteries pressure 
init('P_sp') = 75;%0.85 * init('P_sa');   % systemic peripheric pressure (capilar pressure?)

init('mean_P_sa') = init('P_sa');

%flows  %we will have to set it depending if simulation starts at systole or diastole.
local_vars('TotBV') = 5000; %blood volume (ml)
local_vars('TotFlow') = local_vars('TotBV')/60; %ml/s
init('Q_sa') = local_vars('TotFlow')*0.85*0.35;   % flow of systemic arteries  
init('Q_pa') = local_vars('TotFlow')/2;   % flow of pulmonary arteries

%volumes   define if I will begin with systole or dyastole!. check Vt - Vu = V!.

% Use the proportion between Total venous volume and unstressed blood volume
init('V_total_e_v') = pars('V_u_e_v_0') * 1.3;%pars('V_unstressed_e_p') * 1.3;%/0.97;
init('V_total_s_v') = pars('V_u_s_v_0') * 1.5;%pars('V_unstressed_s_p') * 1.3;%/0.97;
init('V_total_b_v') = pars('V_unstressed_b_v') * 1.1;%/0.97;
init('V_total_h_v') = pars('V_unstressed_h_v') * 1.1;%/0.97;
init('V_total_rm_v') = pars('V_u_rm_v_0') * 1.5;%pars('V_unstressed_rm_p') * 1.3;%/0.97;
init('V_total_am_v') = pars('V_u_am_v_0') * 1.5;%pars('V_unstressed_am_v') * 2;%/0.97;

init('V_total_vc')  = 280;%pars('V_unstressed_vc') * 1.5;%1/0.97;

init('V_total_lv') = pars('V_unstressed_lv') * 3; %1.8;%1.1189;%1/0.97;
init('V_total_la') = pars('V_unstressed_la') * 1.1;%1/0.97;
init('V_total_rv') = pars('V_unstressed_rv') *2;%1/0.97;
init('V_total_ra') = pars('V_unstressed_ra') * 1.1;%1/0.97;   

init('V_total_pa') = 20; %pars('V_unstressed_pa') * 1.1;%1/0.97;  precaution!!! V_unstressed_pa is 0!!
init('V_total_pp') = pars('V_unstressed_pp') * 1;%1/0.97;
init('V_total_pv') = pars('V_unstressed_pv') * 1.1;%1/0.97;






%Peripheric resistances  (This are all arbitrary values) use ohms law to compute it (this is the best approach)
 %factor para hacer dPs_p no sea tan brusco%

init('R_e_p') = 0.08;
init('R_s_p') = 0.8;
init('R_b_p') = 0.8;
init('R_h_p') = 0.8;
init('R_rm_p') = 0.8;
init('R_am_p') = 0.8;

init('R_rm_p_n') = 0.8;
init('R_am_p_n') = 0.8;

%Elastances
init('E_max_lv') = 26;%50;%26;%8;
init('E_max_rv') = 4;%4;

%Theart
init('Theart') = 1;                   % s 
init('HR') = 1;%1/init('Theart');
init('t0_heart') = 0;     %we will begin with dyastole's end
init('zheta_heart') = 0;
init('u_t0') = 0;

%Control associated variables
init('P_mean') = 80;
init('f_ac') = 8.0807;
init('f_ap') = 4.4492;   %[spikes/s]
init('xO2_b') = 0;
init('xCO2_b') = 0;
init('xO2_h') = 0;
init('xO2_rm') = 0;
init('xCO2_h') = 0;
init('xCO2_rm') = 0;
init('Wh') = 0;
init('xO2_am') = 0;
init('x_met') = 0;
init('x_M') = 0;
init('DThetaO2_h_s') = 0;
init('DThetaO2_p_s') = 0;
init('DThetaO2_v_s') = 0;
init('DThetaCO2_h_s') = 0;
init('DThetaCO2_p_s') = 0;
init('DThetaCO2_v_s') = 0;
init('DTsym') = 0;
init('DTvagal') = 0;
init('DTheta_R_e_p') = 0;
init('DTheta_R_s_p') = 0;
init('DTheta_R_rm_p_n') = 0;
init('DTheta_R_am_p_n') = 0;
init('DTheta_V_unstressed_e_v') = 0;
init('DTheta_V_unstressed_s_v') = 0;
init('DTheta_V_unstressed_rm_v') = 0;
init('DTheta_V_unstressed_am_v') = 0;
init('DTheta_Emax_lv') = 0;
init('DTheta_Emax_rv') = 0;

init('phi_met') = 0;

init('fh_s') = 0;
init('fp_s') = 0;
init('fv_s') = 0;

init('fv') = 0;

%Testing_variables
init('var1') = 0;
init('var2') = 0;
init('var3') = 0;
init('var4') = 0;
init('var5') = 0;
init('var6') = 0;
init('var7') = 0;
init('var8') = 0;
init('var9') = 0;
init('R_bp') = 0;
init('Q_e') = 0;
init('Q_s') = 0;
init('Q_h') = 0;
init('Q_rm') = 0;
init('Q_am') = 0;
init('P_v_e') = 0;
init('P_v_s') = 0;
init('P_v_h') = 0;
init('P_v_rm') = 0;
init('P_v_am') = 0;
init('R_e_p') = 0;
init('R_s_p') = 0;
init('R_b_p') = 0;
init('R_h_p') = 0;
init('R_rm_p') = 0;
init('R_am_p') = 0;

%HIPOXIA
init('xO2_e') = 0;
init('xCO2_e') = 0;
init('xO2_s') = 0;
init('xCO2_s') = 0;
init('xO2_p') = 0;
init('xCO2_p') = 0;









%%%%%% DELAYS: delays for control equations and some of open-loop equations


taus('tau_gases') =  pars('LCTV') / init('Qla');




pars('lb_TI')  = 1;%1.2;  %we could change them to be 1 the minimum value
pars('lb_TE') = 1.2;%1.2;   %cambiar a mas bajo
pars('lb_a1')  = 20; %TESTING 11.8; %15 %11; %36.2;%10;   %36.2  %bajar este también, en base a los papers este debería estar en 20
pars('lb_a2')  = -20;%TESTING -2.5;%-20;  %-5.5
pars('lb_tau') = 0.25;

pars('ub_TI') = 3;%1.8;%2;    %this could be 2, its minimum value, cambiar a más alto
pars('ub_TE') = 5;%1.95;%2.2;
pars('ub_a1') = 40; %60;  
pars('ub_a2') = -2.5;%TESTING 0;
pars('ub_tau') = 0.45;

pars('TI') = 0;
pars('TE') = 0;
pars('a0') = 0;
pars('a1') = 0;
pars('a2') = 0;
pars('tau') = 0;
pars('Tresp') = 0;

%% Metasimulation parameters

pars('type_of_input') = 0;













    
    




end 


%% References

 % 1.Albanese A, Cheng L, Ursino M, Chbat NW. An integrated mathematical model of the human cardiopulmonary system: model
 % development. Am J Physiol Circ Physiol 310: H899–H921, 2016.
 % 2. Cheng L, Ivanova O, Fan H-H, Khoo MCK. An integrative model of respiratory and cardiovascular control in sleep-disordered
 % breathing. Respir Physiol Neurobiol 174: 4–28, 2010.
 % 3. Fincham WF, Tehrani FT. A mathematical model of the human respiratory system. J Biomed Eng 5: 125–33, 1983.
 % 4. Magosso E, Ursino M. A mathematical model of CO2 effect on cardiovascular regulation. Am J Physiol Heart Circ Physiol 281:
 % H2036-52, 2001.
 % 5. Magosso E, Ursino M. Cardiovascular response to dynamic aerobic exercise: a mathematical model. Med Biol Eng Comput 40: 660–
 % 74, 2002.
 % 6. Poon CS, Lin SL, Knudson OB. Optimization character of inspiratory neural drive. J Appl Physiol 72: 2005–17, 1992.
 % 7. Serna Higuita LY, Mananas MA, Mauricio Hernandez A, Marina Sanchez J, Benito S. Novel Modeling of Work of Breathing
 % for Its Optimization During Increased Respiratory Efforts. IEEE Syst J 10: 1003–1013, 2016.
 % 8. Serna LY, Mañanas MA, Hernández AM, Rabinovich RA. An Improved Dynamic Model for the Respiratory Response to
 % Exercise. Front Physiol 9: 69, 2018.
 % 9. Spencer JL, Firouztale E, Mellins RB. Computational expressions for blood oxygen and carbon dioxide concentrations. Ann Biomed
 % Eng 7: 59–66, 1979.
 % 10. Ursino M. Interaction between carotid baroregulation and the pulsating heart: a mathematical model. Am J Physiol 275: H1733–
 % H1747, 1998.
 % 11. Ursino M, Magosso E. Acute cardiovascular response to isocapnic hypoxia. I. A mathematical model. Am J Physiol Heart Circ
 % Physiol 279: H149-65, 2000.