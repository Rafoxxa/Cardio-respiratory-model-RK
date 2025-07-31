function [pars, init] = load_global_easy()

%BW        = data.BW;                         % Body weight (kg)                
%Hgt       = data.Hgt;                        % Height (cm)                     
%Gender    = data.Gender;                     % Gender (1=female, 2=male)     

%General parameters
a_0   =  0.2;        % Positive constants, represent changes in the O2 affinity to Hb in adults, Revow et al. (1989)
b_0   =  0.05;        % Positive constants, represent changes in the O2 affinity to Hb in adults, Revow et al. (1989)
c_0   = 2;       % Positive constants, standard value
k_CO2 = 0.224;   % (liters)
S=99;
%% % Body parameters (are constant, based in BTPS)
T       = 37;                      % Temperature (ºC)
P_amb   = 760;                     % Ambient pressure at sea level (mmHg)
W_vapor = 863;                     % Water vapor in terms of BTPS
K_CO2 = 0.0057;   % (L/mmHg) 

%% Frecuency
F_V   = 15;          % Frequency of breaths per minute (L/min)
FR    = 1/(60/15);          % Frequency of breaths per minute (breaths/min)



%% Volume

% Lung compartment

Va_CO2 = 3.2;        % Effective CO2 volume in arterial (liters)
Va_O2  = 2.5;        % Effective O2 volume in arterial (liters)
V_D    = 0.150;      % Dead space ventilation (liters)
V_AO2  = Va_O2;      
V_ACO2 = Va_CO2;

% Tissue compartment
V_TCO2 = 15.0;       % Effective CO2 volume in the tissue compartment (liters)
V_TO2  = 6.0;        % Effective O2 volume in the tissue compartment (liters)

% Brain compartment
V_bCO2 = 0.9;        % Effective CO2 volume in the brain compartment (liters)


%% Pressure

% Arterial pressure (simplifying assumptions)
P_aO2  =  log([(S / 100) - 1] / -2.4) / -0.05;    % Arterial O2 partial pressure (mmHg)
P_aCO2 = 40;     % Arterial CO2 partial pressure (mmHg)

% Alveolar pressure 
P_ACO2 = P_aCO2;  % Alveolar pressure of CO2 (mmHg) 
P_AO2  = P_aO2;                    % Alveolar pressure of O2 (mmHg) or F_AO2 * (P_amb - 47)

% Inspired air pressure 
P_iCO2 = 0.0;                    % Inspired air pessure of CO2 (mmHg) or F_iCO2 * (P_amb - 47)
P_iO2  = 765; %150.0 ;                 % Inspired air pessure of O2 (mmHg)
P_eO2 = 755;%100.0;                   % solo una suposición
P_air_O2 = (P_iO2 + P_eO2)/2;    % mean O2 pressure in the air.
% Venous gas pressure (are constant)
P_vO2   = 40.9;                  % Venous pressure of CO2 (mmHg) 
P_vCO2  = 46.0;                  % Venous pressure of O2 (mmHg)


% Interstitial partial pressure
P_tCO2 = P_vCO2;%45;                 % Interstitial CO2 partial pressure

% Concentration-partial pressure 
Ca_O2  = a_0 * (1 - exp(-b_0 * P_aO2))^c_0;       % O2 concentration
Ca_CO2 = K_CO2 * P_aCO2 + k_CO2;                 % CO2 concentration 

Cv_O2  = a_0 * (1 - exp(-b_0 * P_vO2))^c_0;
Cv_CO2 = K_CO2 * P_vCO2 + k_CO2;





%% Flows

%Ventilatory
qVE = 6.0;                                    % Minute ventilation (L/min) 
qVD = 2.2;                                    % Total dead space ventilation (L/min)
qVA = qVE - qVD;                              % Alveolar ventilation (L/min)

%Perfusion (cardiovascular submodel is not included, are equated)
qCO =  6.0;                                   % Cardiac output (l/min)
qS =  6.0;                                   % Tissue blood flow (l/min)
qFP =  6.0;                                   % Pulmonary blood flows (l/min)
qB =  0.8;                                   % Brain blood flows (liter/min)

%% Respiration parameters
G_PO2 = 60;      % (mmHg), Initial alveolar-capillary O2 partial pressure gradients
G_PCO2 = 6;      % (mmHg), Initial alveolar-capillary CO2 partial pressure gradients





%Cv_O2  = B_v + (m_v * P_vO2);                  % Venous blood gas levels/concentration of O2 
%Ca_CO2 = B_a + (m_a * P_aCO2);                 % Arterial blood gas levels/concentration of CO2 
                                       % Arterial blood gas levels/concentration of CO2 



% Tissue compartment parameters
Ct_CO2 = Cv_CO2;        % Tissue blood gas levels/concentration of CO2
Ct_O2  = Cv_O2;         % Tissue blood gas levels/concentration of O2
M_CO2 = 1 * 0.220/1000;         % Metabolic production of CO2 (L/min)
M_O2  = -0.270/1000;         % metabolic consumption of CO2 (L/min)

% Brain compartment parameters
Cvb_CO2 = 0.9/15 * Cv_CO2; %Solo haciendo proporciones de volumen de CO2 
Cb_CO2  = Cvb_CO2;      % Brain tissue blood gas levels/concentration of CO2
%Cvb_CO2 =               % Brain venous blood gas levels/concentration of CO2
M_brain_CO2 = 0.04;         % Brain CO2 metabolic production (L/min)


% Pressures in lungs
C0 = 10;
%C1 = 2;
R0 = 1;
%R1 = 1;
R1 = 0.0023/1000; %mmHg/Ls
C1 = 0.33*1000; % L/mmHg
Patm = 760;   %mmHg

P0 = 1013.25/0.0065;
P1 = 760;%1013.25/0.0065;
Vref0 = 0;%500;
Vref1 = 0;%500;
V0 = C0 * P0 + Vref0;
V1 = C1 * P1 + Vref1;


PA = 762;
VA = 5;
fO2 = 0.12;
fCO2 = 0.04;
PO2_capilar_lung = 90;
PCO2_capilar_lung = 35;
PO2_capilar_tissue = 80;
PCO2_capilar_tissue = 40; 


%init_values=[P_aO2,P_aCO2,Cv_O2,Cv_CO2,Cb_CO2, P0, P1, V0, V1];
init_values = [PA, VA, fO2, fCO2, PO2_capilar_lung, PCO2_capilar_lung, PO2_capilar_tissue, PCO2_capilar_tissue];
init_keys = ["PA", "VA", "fO2", "fCO2", "PO2_capilar_lung", "PCO2_capilar_lung", "PO2_capilar_tissue", "PCO2_capilar_tissue"];

init = dictionary(init_keys, init_values);
disp(init);

pars_values=[P_iO2,FR,qVA,qFP,qS,qB,...
    V_AO2,V_ACO2,V_TCO2,V_TO2,...
    a_0,b_0,c_0                 ,...
    K_CO2,k_CO2,                      ...
    M_CO2,M_O2,M_brain_CO2, P_air_O2, R0, R1, C0, C1, Patm];

pars_keys = ["P_iO2","FR","qVA","qFP","qS","qB",...
    "V_AO2","V_ACO2","V_TCO2","V_TO2",...
    "a_0","b_0","c_0"                 ,...
    "K_CO2","k_CO2",                      ...
    "M_CO2","M_O2","M_brain_CO2", "P_air_O2", "R0", "R1", "C0", "C1", "Patm"];

pars = dictionary(pars_keys, pars_values);


end 