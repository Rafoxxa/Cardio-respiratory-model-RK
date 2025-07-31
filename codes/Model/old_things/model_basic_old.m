function xdot= model_basic(time,init_values,pars, init_keys)

    
%disp(time);
    y = dictionary(init_keys, init_values);
   % Unpack the state variable vector   
    PA = y('PA');  
    VA = y('VA');   
   
    f_O2 = y('fO2');
    f_CO2 = y('fCO2'); 

    PO2_capilar_lung  = y('PO2_capilar_lung');       
    PCO2_capilar_lung  = y('PCO2_capilar_lung');   

    PO2_capilar_tissue = y('PO2_capilar_tissue');
    PCO2_capilar_tissue = y('PCO2_capilar_tissue');
    
    %cCO2_brain  = y('Cb_CO2');      
    
    % disp('PA')
    % disp(PA)
    % disp('VA')
    % disp(VA)
    % 
    % disp('fO2')
    % disp(f_alvO2)
    % disp('fCO2')
    % disp(f_alvCO2)
    % 
    % disp('PO2_capilar_lung');       
    % disp(PO2_capilar_lung);
    % disp('PCO2_capilar_lung');   
    % disp(PCO2_capilar_lung);
    % 
    % disp('PO2_capilar_tissue');
    % disp(PO2_capilar_tissue);
    % disp('PCO2_capilar_tissue');
    % disp(PCO2_capilar_tissue); 

    P_capilar_lung = [PO2_capilar_lung, PCO2_capilar_lung];
    P_capilar_tissue = [PO2_capilar_tissue, PCO2_capilar_tissue];
    f_lungs = [f_O2, f_CO2];
    
    
         
    


    %% Parameters

    P_iO2 = pars('P_iO2');
    FR = pars('FR');
    qVA = pars('qVA');
    qFP = pars('qFP');
    qS = pars('qS');
    qB = pars('qB');
    V_AO2 = pars('V_AO2');
    V_ACO2 = pars('V_ACO2');
    V_TCO2 = pars('V_TCO2');
    V_TO2 = pars('V_TO2');
    a0 = pars('a_0');
    b0 = pars('b_0');
    c0 = pars('c_0');
    K_CO2 = pars('K_CO2');
    k_CO2 = pars('k_CO2');
    M_CO2 = pars('M_CO2');
    M_O2 = pars('M_O2');
    M_brain_CO2 = pars('M_brain_CO2');
    P_air_O2 = pars('P_air_O2');

    pars("kPa") = 7.5;
    kPa = pars("kPa");
    pars("alpha_CO2_ery") = 0.195/(1000 * kPa);   %these are in mmol, so because pH function takes co2 concentrations as mol, we had to change it to divide it by 1000.
    pars("alpha_CO2_plasma") = 0.23/(1000 * kPa);
    pars("cHb") = 9.3/1000;   %9.3 mmol/L
    pars("cHb_ery") = 21/1000; % 21 mmol/L
    pars("co20first") = 0.008;   % 1mmol/L
    pars("co20last") = 0.04;     % 50 mmol/L
    
    pars("KaCO2") = 10^(-6.1);  %mol/L
    pars("NaOH_0") = 0.0462;  %mmol/L  % we have to check this thing, because is extremely rare that sodium travels as NaOH in blood.
    pars("HPr_0") = 0.0398;   %mmol/L
    pars("KaPr") = 10^(-7.3);
    
    pars("kmet") = 0;%1%0.06; %it is too tiny
    pars("Temp") = 36;
    pars("alpha_O2") = 9.83/(1000*1000*7.5);
    
    pars("Rgases") = 62.36367; %mmHg L /(mol K)
    pars("f_extO2") = 0.21;  % percentage of O2 in air
    pars("f_extCO2") = 0.04; % percentage of CO2 in air
    pars("kCO2") = 30.46/10000; %this values have to be changed
    pars("kO2") = 30.46/10000;  %this values have to be changed
    pars("Patm") = 760;        % 
    pars("V0_alv") = 0.15;     % Residual alveolar volume (L)
    %Some blood fluxes
    pars("Q_lungs") = 6;%4.5;
    pars("Q_tissue") = 6;%4.6;

    %pars('Tinsp') = 0.23;
    %pars('')


%%  Exterior compartment (Oscillator)
    h = 0;%time*(10000/60);
    altitude_factor = (1013.25 * (1 - 0.0065 * h / 288.15)^(9.81 * 0.029 / (8.314 * 0.0065)))/(1013.25); % altitude factor in mmHg using barometric formula in the troposphere
    
    
    

    
  

   
  %% Physiology
     
    Q_lungs = pars("Q_lungs");
    Q_tissue = pars("Q_tissue");
    Patm = pars("Patm");
   
    Ptor = Patm - 2 + pulmonary_osc(time, pars);
    [dPA, dVA] = alveolar(Patm, PA, Ptor, pars);

    
    P_venous = P_capilar_tissue;
    
    [dfO2, dfCO2] = lung_air_exchange(f_lungs, PA, Ptor, PO2_capilar_lung, PCO2_capilar_lung,  pars);
    [dPO2_capilar_lung, dPCO2_capilar_lung] = gas_transport_lungs(P_venous, f_lungs, Q_lungs, PA, P_capilar_lung, pars);
    
    PO2_arterial = PO2_capilar_lung;
    PCO2_arterial = PCO2_capilar_lung;
    P_arterial = [PO2_arterial, PCO2_arterial];
    
    M = [M_O2, M_CO2];
    
    [dPO2_capilar_tissue, dPCO2_capilar_tissue] = gas_transport_tissue(P_arterial, Q_tissue, P_capilar_tissue, M, pars);
    
    


 


 %% dxdt   
    xdot=[dPA, dVA, dfO2, dfCO2, dPO2_capilar_lung, dPCO2_capilar_lung, dPO2_capilar_tissue, dPCO2_capilar_tissue]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55



%% Functions

%Testing

model_parameters = define_parameters();
%dVt = @(t) sin(2*pi*t*0.15)*0.2;
D = @(t) sin(2*pi*t*0.15)*0.2;
dQ = 5;
dt = 0.05;


dP1v = zeros(2,300);
dP2v = zeros(2,300);
dP3v = zeros(2,300);
dP4v = zeros(2,300);
dP5v = zeros(2,300);

P1v = zeros(2,300);
P2v = zeros(2,300);
P3v = zeros(2,300);
P4v = zeros(2,300);
P5v = zeros(2,300);

PAlvgasv = zeros(2, 300);
Pav = zeros(2, 300);
vv = zeros(2, 300);

Vtv = zeros(1, 300);
dVtv = zeros(1, 300);
PAlvv = zeros(1,300);


%init_cond

P_1 = [150,5];
P_2 = [150,5];
P_3 = [150,5];
P_4 = [150,5];
P_5 = [150,5];
PAlvgas = [90, 80];
Pa = [100, 20];
dPa = [0, 0];
Vt = 0;
dVt = 0;

G = 1;
v = [0.1, 0.023];
a = [0.1, 0.023];

for t = 1:300
a = dissociation(Pa, model_parameters);   

Pao = expiratory_pressure(dVt, Vt, model_parameters);
[dVt, PAlv] = ventilatory_mechanics(Vt, Pao,D(t*dt),G,model_parameters);

[dP1, dP2, dP3, dP4, dP5] = deadspace(dVt, P_1,P_2, P_3, P_4, P_5,PAlvgas, model_parameters);


[dPAlvgas] =  alveolar_exchange(PAlvgas, Pa, P_5, dQ, dVt,v, model_parameters);
[ddPa]  = mixing(dQ, PAlvgas, dPa, Pa, model_parameters);
[dv]    = tissue(a,v,dQ, model_parameters);

PAlvgas = PAlvgas + dPAlvgas * dt;
dPa = dPa + ddPa * dt;
Pa = Pa + dPa * dt;
v = v + dv * dt;

dVt = D(t*dt+1);  %We have to change this
Vt = Vt + dVt*dt;


P_1 = P_1 +  dP1 * dt;
P_2 = P_2 + dP2 * dt;
P_3 = P_3 + dP3 * dt;
P_4 = P_4 + dP4 * dt;
P_5 = P_5 + dP5 * dt;

dP1v(:,t) = dP1;
dP2v(:,t) = dP2;
dP3v(:,t) = dP3;
dP4v(:,t) = dP4;
dP5v(:,t) = dP5;

P1v(:,t) = P_1;
P2v(:,t) = P_2;
P3v(:,t) = P_3;
P4v(:,t) = P_4;
P5v(:,t) = P_5;

Pav(:,t) = Pa;
PAlvgasv(:,t) = PAlvgas;
vv(:,t) = v;
Vtv(:, t) = Vt;
dVtv(:, t) = dVt;
PAlvv(:,t) = PAlv;
end

%model parameters

function model_parameters = define_parameters()
    model_parameters = containers.Map;
    model_parameters('Vdead1') = 0.15; % dead space volume
    model_parameters('Vdead2') = 0.15; % dead space volume
    model_parameters('Vdead3') = 0.15; % dead space volume
    model_parameters('Vdead4') = 0.15; % dead space volume
    model_parameters('Vdead5') = 0.15; % dead space volume
    
    model_parameters('PIO2') = 150;
    model_parameters('PICO2') = 5;

    model_parameters('a1') = -1.6; % parameter in O2 dissociation equation
    model_parameters('a2') = -0.06; % parameter in CO3 dissociation equation
    model_parameters('alpha1') = 2.7; % parameter in O3 dissociation equation
    model_parameters('alpha2') = 0.03; % parameter in CO3 dissociation equation
    model_parameters('K1') = 7.7e-4; % parameter in O3 dissociation equation
    model_parameters('K2') = 4.8e-4; % parameter in CO3 dissociation equation
    model_parameters('beta1') = -1.4; % parameter in O3 dissociation equation
    model_parameters('beta2') = -0.02; % parameter in CO3 dissociation equation
    model_parameters('CO2a_') = 0.1; % parameter in O3 dissociation equation
    model_parameters('O2a_') = 0.2; % parameter in CO3 dissociation equation

    model_parameters('VCO2') = 0.25; % lung storage volume for CO2
    model_parameters('VO2') = 0.25; % lung storage volume for O2

    model_parameters('T1') = 0.5; % time constant for cardiovascular mixing
    model_parameters('T2') = 0.5; % time constant for cardiovascular mixing
    model_parameters('Ta') = 0.5; % lung to chemoreceptor circulation delay
    model_parameters('LCTV') = 1.5; % lung to chemoreceptor transportation vascular volume constant

    model_parameters('MRCO2') = 0.003; % metabolic production rate for CO2
    model_parameters('MRO2') = 0.003; % metabolic consumption rate for O2

    model_parameters('VtissueO2') = 0.3; % metabolic production rate for CO2
    model_parameters('VtissueCO2') = 0.3; % metabolic consumption rate for O2
    model_parameters('SI') = 0; % state variable showing awake(0) or sleep(1)

    model_parameters('Yua') = 1;  %1/Rua Rua = 1;
    model_parameters('Raw') = 1.016;
    model_parameters('Rlt') = 1.69;
    model_parameters('Rcw') = 1.03;
    model_parameters('VC') =  8; %it is vital capacity
    model_parameters('Ers') = 1;
    model_parameters('b') = 1;
    



    model_parameters('C1') = 1.34; % maximum concentration of hemoglobin-bound oxygen
    model_parameters('C3') = 48; % maximum carbon dioxide concentration
    
    model_parameters('MRbCO2') = 0.0005; % metabolic production rate for CO2 in the brain tissue
    model_parameters('SCO2') = 0.03; % dissociation slope for CO2 in the blood
    model_parameters('SbCO2') = 0.03; % dissociation slope for CO2 in the brain tissue
    model_parameters('RBPN') = 10; % cerebral peripheral flow resistance
    model_parameters('CvO2n_b') = 15; % nominal venous O2 concentration in cerebral peripheral circulation
    model_parameters('PaCO2_n') = 40; % nominal arterial CO2 partial pressure
    model_parameters('Tau_CO2') = 10; % time constant for peripheral CO2 response
    model_parameters('Tau_O2') = 10; % time constant for peripheral O2 response
    model_parameters('A') = 10000; % parameter for flow regulation equation
    model_parameters('B') = 10000; % parameter for flow regulation equation
    model_parameters('C') = 10000; % parameter for flow regulation equation
    model_parameters('GO2_b') = 10000; % gain of local O2 response on cerebral vascular bed
    model_parameters('Ic') = 40; % central apneic threshold
    model_parameters('Ip_CO3') = 50; % peripheral apneic threshold for CO3
    model_parameters('Ip_O3') = 50; % peripheral apneic threshold for O3
    model_parameters('Gc') = -10/40;% gain for central chemical drive 
    % (negative sign is used to make it consistent with the sign convention used in the paper)
    % The value of Gc is calculated by dividing the change in ventilation by the change in PaCO2 
    % when PaO2 is held constant at a normal level of around 100 mmHg.
    % The value of -10 L/min/mmHg is obtained from Table II of the paper.
    % The value of Ic is taken from Table I of the paper.
    % The value of Swake is taken from Table III of the paper.
    % The values of Ip_CO3 and Ip_O3 are taken from Table IV of the paper.


end

%EQUATIONS (CHENG-SARMIENTO-SERNA)

% Gas transport through the dead space
function [dP_1, dP_2,dP_3,dP_4,dP_5] = deadspace(dVt, P_1,P_2, P_3, P_4, P_5,PAlv, model_parameters)
    Vdead1 = model_parameters('Vdead1'); % dead space volume
    Vdead2 = model_parameters('Vdead2'); 
    Vdead3 = model_parameters('Vdead3'); 
    Vdead4 = model_parameters('Vdead4'); 
    Vdead5 = model_parameters('Vdead5'); 
    PIO2 = model_parameters('PIO2');
    PICO2 = model_parameters('PICO2');
    
    %P = [PO2, PCO2];
    PI = [PIO2, PICO2];
    
    if dVt >= 0    %Inspiration
        dP_1 = (PI - P_1) * dVt/Vdead1;
        dP_2 = (P_1 - P_2) * dVt/Vdead2;
        dP_3 = (P_2 - P_3) * dVt/Vdead3;
        dP_4 = (P_3 - P_4) * dVt/Vdead4;
        dP_5 = (P_4 - P_5) * dVt/Vdead5;
    else        %Expiration
        dP_1 = (P_2 - P_1) * dVt/Vdead1;
        dP_2 = (P_3 - P_2) * dVt/Vdead2;
        dP_3 = (P_4 - P_3) * dVt/Vdead3;
        dP_4 = (P_5 - P_4) * dVt/Vdead4;
        dP_5 = (PAlv - P_5) * dVt/Vdead5;
    end
  
    
end

% Blood–gas dissociation curves
function gas = dissociation(P, model_parameters)
    
    
    a1 = model_parameters('a1'); % parameter in O2 dissociation equation
    a2 = model_parameters('a2'); % parameter in CO3 dissociation equation
    alpha1 = model_parameters('alpha1'); % parameter in O3 dissociation equation
    alpha2 = model_parameters('alpha2'); % parameter in CO3 dissociation equation
    K1 = model_parameters('K1'); % parameter in O3 dissociation equation
    K2 = model_parameters('K2'); % parameter in CO3 dissociation equation
    beta1 = model_parameters('beta1'); % parameter in O3 dissociation equation
    beta2 = model_parameters('beta2'); % parameter in CO3 dissociation equation
    CO2a_ = model_parameters('CO2a_');
    O2a_ = model_parameters('O2a_');
    
    PO2 = P(1);
    PCO2 = P(2);

    FCO2 = PCO2 * (1 + beta2 * PO2)/(K2 * (1 + alpha2 * PO2));
    CO2 = CO2a_ * FCO2^(1/a2)/(1 + FCO2^(1/a2));
    
    FO2 = PO2 * (1 + beta1 * PCO2)/(K1 * (1 + alpha1 * PCO2));
    O2 = O2a_ * FO2^(1/a1)/(1 + FO2^(1/a1));
    
    gas = [O2, CO2];
    
end

% Alveolar gas exchange
function [dPAlv] = alveolar_exchange(PAlv, Pa, P_5, dQ, dVt,v, model_parameters)
    VLCO2 = model_parameters('VCO2'); % lung storage volume for CO2
    VLO2 = model_parameters('VO2'); % lung storage volume for O2[^2^][2]
    VL = [VLO2, VLCO2];

   

    
    a = dissociation(Pa, model_parameters);
    
    if dVt >= 0
        dPAlv = (863 * dQ * (v - a) + dVt * (P_5 - PAlv)) * (eye(2) *  1./VL); % ODE for alveolar CO2 partial pressure
        
    else
        
        dPAlv = (863 * dQ * (v - a)) * (eye(2) *  1./VL);
    end
   
    
    
end

% Cardiovascular mixing (heart and vasculature)
function [ddPa] = mixing(dQ, PAlv, dPa, Pa, model_parameters)
    T1 = model_parameters('T1'); % time constant for cardiovascular mixing
    T2 = model_parameters('T2'); % time constant for cardiovascular mixing
    LCTV = model_parameters('LCTV'); % lung to chemoreceptor transportation vascular volume constant
    %Ta = LCTV/dQ;   %for now lets supose Ta is 0.

    ddPa = 1/(T1 * T2) * (PAlv - (T1 + T2)*dPa - Pa);

end
    


% Brain Compartment CO2 exchange
function [dPbCO2] = brain(PbCO2, PaCO2, model_parameters)
    MRbCO2 = model_parameters('MRbCO2'); % metabolic production rate for CO2 in the brain tissue[^4^][4]
    SCO2 = model_parameters('SCO2'); % dissociation slope for CO2 in the blood
    SbCO2 = model_parameters('SbCO2'); % dissociation slope for CO2 in the brain tissue
    dPbCO2 = (MRbCO2 - SCO2 * (PbCO2 - PaCO2))/SbCO2; % ODE for partial CO2 pressure from the brain
end

% Cerebral blood flow
function [dRBPN] = cerebral(PaO21, PaO22, PaCO21, PaCO22, model_parameters)
    RBPN = model_parameters('RBPN'); % cerebral peripheral flow resistance
    CvO2n_b = model_parameters('CvO2n_b'); % nominal venous O2 concentration in cerebral peripheral circulation
    PaCO2_n = model_parameters('PaCO2_n'); % nominal arterial CO2 partial pressure
    Tau_CO2 = model_parameters('Tau_CO2'); % time constant for peripheral CO2 response[^5^][5]
    Tau_O2 = model_parameters('Tau_O2'); % time constant for peripheral O2 response
    A = model_parameters('A'); % parameter for flow regulation equation
    B = model_parameters('B'); % parameter for flow regulation equation
    C = model_parameters('C'); % parameter for flow regulation equation
    GO2_b = model_parameters('GO2_b'); % gain of local O2 response on cerebral vascular bed[^6^][6]
    
    dRBPN = -RBPN * (A * exp(-B * (PaO21 - PaO22)/Tau_O2) * (GO2_b * (CvO21 - CvO22) - CvO21_n) ...
        + A * exp(-B * (PaCO21 - PaCO22)/Tau_CO3) * (C - 1) * (PaCO21 - PaCO22))/C; % ODE for cerebral peripheral flow resistance
end

% Gas exchange in the body tissues compartments
function [dv] = tissue(a,v,dQ, model_parameters)
    VtissueCO2 = model_parameters('VtissueCO2'); % body tissue storage volume for CO3
    VtissueO2 = model_parameters('VtissueO2'); % body tissue storage volume for O3
    MRCO2 = model_parameters('MRCO2'); % metabolic production rate for CO3
    MRO2 = model_parameters('MRO2'); % metabolic consumption rate for O3
    SI = model_parameters('SI');
    M_ = [MRO2, MRCO2];
    Vtissue = [VtissueO2, VtissueCO2];
    M = M_ * (0.375 * (1 - 0.4*SI) + 0.625);
    dv = (M * eye(2).*[-1,1] + dQ*(a-v)) * eye(2)./Vtissue;
end

% Ventilatory controller (For normal breathing)
function [D] = ventilatory_oscillator(t, model_parameters)
    Ic = model_parameters('Ic'); % central apneic threshold
    Ip_CO3 = model_parameters('Ip_CO3'); % peripheral apneic threshold for CO3
    Ip_O3 = model_parameters('Ip_O3'); % peripheral apneic threshold for O3
    Gc = model_parameters('Gc'); % gain for central chemical drive
    Gp = model_parameters('Gp'); % gain for peripheral chemical drive
    Swake = model_parameters('Swake'); % factor of wakefulness to sleep
    
    FR = 0.15;
    D = sin(2 * pi*FR*t);
end

function Pao = expiratory_pressure(dVt, Vt, model_parameters)
    Raw = model_parameters('Raw');
    Rlt = model_parameters('Rlt');
    Rcw = model_parameters('Rcw');
    Ers = model_parameters('Ers');
    Yua = model_parameters('Yua');

    Yrs = Yua/(1 + (Raw + Rlt + Rcw) * Yua);  %Remember that Y is 1/R

    Pao_ = Vt * Ers + dVt/Yrs;
    if Pao_ >= 0
        Pao = Pao_;
    else 
        Pao = 0;
    end

end
function [dVt, PAlv] = ventilatory_mechanics(Vt,Pao,D,G,model_parameters)
    Yua = model_parameters('Yua');
    Raw = model_parameters('Raw');
    Rlt = model_parameters('Rlt');
    Rcw = model_parameters('Rcw');
    VC = model_parameters('VC');
    Ers = model_parameters('Ers');
    b = model_parameters('b');
    %Ecw = model_parameters('Ecw');
    %Elt = model_parameters('Elt');
    Ecw = 1;
    Elt = 1;
    
    
    Pisom = G * D;
    Yrs = Yua/(1 + (Raw + Rlt + Rcw) * Yua);
    term1 = Pisom*exp(-Vt/0.28*VC);
    dVt = sqrt((0.25 ^ (term1) * Yrs + Vt * Ers*Yrs - Pao*Yrs)^2 + 4*b^(Vt)*(term1*Yrs - Vt*Ers*Yrs + Pao*Yrs) ) - (0.25 ^ (term1) * Yrs + Vt * Ers*Yrs - Pao*Yrs)/2;
    
    Pmus = (term1 * (b^Vt - 0.25*Vt))/(Vt + b^(Vt));
    
    Ppl = Rcw*dVt + Ecw*Vt - Pmus;
    PAlv = Rlt*dVt + Elt*Vt + Ppl;
end



%OLD EQUATIONS (BRAZIL)

function cO2 = testing_o2(PCO2, PO2, pars)
        cO2 = zeros(1, 100);
        for iter = 1:100
            cO2(iter) = 100 * O2_dissociation(PO2(iter),PCO2, pars);
        end
end



function [dP0, dV0] = airways(P0, P1, Ppl, pars)
    Patm = pars('Patm');
    R0 = pars('R0');
    R1 = pars('R1');
    C0 = pars('C0');
    

    dP0 = 1/(R0 * C0) * (Patm - P0 - R0/R1 * (P0 - P1 - Ppl));
    dV0 = C0 * dP0;
end

function Evar = pulmonary_osc(time, pars)
        Tinsp = 3; %pars("Tinsp");
        Texp  = 1; %pars("Texp");
        Tciclo = Tinsp + Texp; %pars("Tciclo");
        Emax = 1;%pars("Emax");
        Emin = -1;%pars("Emin");
        A = Emax;
        B = Emin;
        tresp= mod(time,Tciclo);
        for i=1:length(tresp)
            if tresp(i)<= Tinsp
                 el(i)= (A).*(1-cos((pi*tresp(i))/(Tinsp)))+B;
            elseif (tresp(i)>= Tinsp) && (tresp(i)<= time)
                el(i)= (A).*(1+cos(pi*((tresp(i)-Tinsp)./((Texp)))))+B;
            else
                el(i)= B;
            end
        end
        Evar = el;
    end

function [dP1, dV1] = alveolar(Patm, PA, Ptor, pars)
    R1 = pars('R1');
    C1 = pars('C1');
    
    P0_ = Patm - 760;
    P1_ = PA - 760;
    Ppl_ = Ptor  - 760;

    dP1 = 1/(R1 * C1) * (P0_ - P1_ - Ppl_);    
    dV1 = C1 * dP1;
end

function [df_O2, df_CO2] = lung_air_exchange(f_alv, P1, Ppl, PvO2, PvCO2,  pars)

    
    f_extCO2 = pars('f_extCO2');
    f_extO2 = pars('f_extO2');

    T = pars('Temp') + 273.15; %transform to Kelvin
    Rgases = pars('Rgases');
    kO2 = pars('kO2');
    kCO2 = pars('kCO2');
    Patm = pars('Patm');
    V0_alv = pars('V0_alv'); %is the alveolar volume at atmosferic pressure, I think it referes to residual volume or minimal volume at atmosferic pressure.
    R1 = pars('R1');
    C1 = pars('C1');
    %R1 = 0.0023/1000; %mmHg/Ls
    %C1 = 0.33*1000; % L/mmHg

    K = [kO2, 0; 0 kCO2];
    P1_ = P1 - 760;
    P0_ = Patm - 760;
    Ppl_ = Patm - 760;
    Pcp = [PvO2, PvCO2];   %In theory this should be Partial pressures of gases in alveolar capillars, which we consider to be the same as pulmonary arterial partial pressures which is venous pressure in systemic circulation
    
    %P0 = Patm;
    
    
    f_ext = [f_extO2, f_extCO2]; 
    factor = Rgases*T/(P1*(C1 * P1_ + V0_alv));
    alveolar_gas_fluxes = (((P0_ - P1_ - Ppl_) >  0) * (P0_ - P1_ - Ppl_)) * (f_ext - f_alv)/R1;
    
    alveolar_pressures = (K * (Pcp - (P1 + Ppl_)*f_alv)')';
    %factor = 0.001;
    df_alv = factor * (alveolar_gas_fluxes + alveolar_pressures);
    df_O2 = df_alv(1);
    df_CO2 = df_alv(2);
end

function [dPO2_cap, dPCO2_cap] = gas_transport_lungs(P_in, f_lungs, Q_lungs, PA, P_cap, pars)
    
    %´´´This allows computing partial presures in alveoli which correspond to those in systemic arterial
    %blood´´´
    
    %PARAMETERS and INPUT VARIABLES
    kCO2 = 30.46/10000;
    kO2 = 30.46/10000; 
    k = [kO2 0; 0 kCO2];                   %hematoalveolar constant.
    Vb = 4.6;                              %Blood volume in lungs 
    P_lungs = PA .* f_lungs;               %This converts gas fractions to partial pressures.
    
    PO2_in = P_in(1);                      %Input blood to alveoli (equivalent to systemic venous partial pressures)
    PCO2_in = P_in(2);

    PO2_cap = P_cap(1);                    %This correspond to systemic arterial partial pressures, or pulmonary venous partial pressures or capilar partial pressures in lungs
    PCO2_cap = P_cap(2);
    
    %Non linear functions for COMPUTING JACOBIAN
    O2fun  = @(P) O2_dissociation(P(1), P(2), pars);
    CO2fun = @(P) CO2_dissociation(P(1), P(2), pars);
    fun = @(P) [O2fun(P), CO2fun(P)];
    jacobian = jacobianest(fun, [PO2_cap, PCO2_cap]);  %maybe the relation is from pressures to mmols, and we computed it using mols so we will have to multiply by 1000
    
    %BALANCES
    %1. Fluid balance
    fluid_balanceO2 = Q_lungs * (O2_dissociation(PO2_in, PCO2_in, pars) - O2_dissociation(PO2_cap, PCO2_cap, pars));
    fluid_balanceCO2 = Q_lungs * (CO2_dissociation(PO2_in, PCO2_in, pars) - CO2_dissociation(PO2_cap, PCO2_cap, pars));
    fluid_balance = [fluid_balanceO2, fluid_balanceCO2];
    %2. Gas balance
    gas_balance = (k * (P_lungs - P_cap)')';    
    
    %ODE
    dP_cap =    (fluid_balance + gas_balance)/(jacobian * Vb);
    dPO2_cap =  dP_cap(1)/60;  %to seconds
    dPCO2_cap = dP_cap(2)/60;  %to seconds

end

function [dPO2_cap, dPCO2_cap] =  gas_transport_tissue(P_in, Q_tissue, P_cap, M, pars)
%´´´This functions allows the computation for variations of partial pressures in each tissue,
% if there's only one tissue in the model, the output pressure (P_cap) would be equivalent to systemic venous partial pressures´´´


%PARAMETERS and INPUT VARIABLES    
Vt = 20;                               %Tissue volume (L)
Vb = 0.12 * Vt;%20;                              %Blood volume in tissue (L)

%disp('linear change due tu metabolism');
%disp(Q_tissue/(Vb + Vt) * M(1));

PO2_in = P_in(1);                      %Input blood to tissue (equivalent to systemic arterial partial pressures)
PCO2_in = P_in(2);       

PO2_cap = P_cap(1);                    %This correspond to venous partial pressures, or pulmonary artery partial pressures or capilar partial pressures in tissue
PCO2_cap = P_cap(2);

%Non linear functions for COMPUTING JACOBIAN
O2blood_fun  = @(P) O2_dissociation(P(1), P(2), pars);
CO2blood_fun = @(P) CO2_dissociation(P(1), P(2), pars);

O2tissue_fun  = @(P) O2_tissue(P(1), P(2), pars);
CO2tissue_fun = @(P) CO2_tissue(P(1), P(2), pars);

blood_fun = @(P) [O2blood_fun(P), CO2blood_fun(P)];
tissue_fun = @(P) [O2tissue_fun(P), CO2tissue_fun(P)];
jacobian_b = jacobianest(blood_fun, [PO2_cap, PCO2_cap]);  %maybe the relation is from pressures to mmols, and we computed it using mols so we will have to multiply by 1000
jacobian_t = jacobianest(tissue_fun, [PO2_cap, PCO2_cap]);
jacobian_sum = Vt * jacobian_t + jacobian_b * Vb; 

%BALANCES
%1. Fluid balance
fluid_balanceO2 =  Q_tissue * (O2tissue_fun(P_in) - O2tissue_fun(P_cap));


fluid_balanceCO2 = Q_tissue * (CO2tissue_fun(P_in) - CO2tissue_fun(P_cap));
fluid_balance = [fluid_balanceO2, fluid_balanceCO2];
  

%ODE
dP_cap =    (fluid_balance + M)/(jacobian_sum);

dPO2_cap =  dP_cap(1)/60;  %to seconds
dPCO2_cap = dP_cap(2)/60;  %to seconds


end

function sO2 = satO2(PO2, PCO2, cCO2, pars)
    
    kmet = pars('kmet'); %This should be dynamic but its effect is too little, we will use mean values
    T = pars('Temp');

    kPa = 7.5; %mmHg

    acid_base_balance = -0.72 * (ph_model(cCO2, pars) - 7.4) + 0.09 * log(PCO2/(5.33 *kPa) ) + kmet -0.07*5;
    DeltaO2 = log(PO2/kPa) - acid_base_balance - 1.946 - 0.0055 * (T-37);  %This is the available oxygen for hemoglobin to bind, we have to substract other oxygen-binding-sustrate from total amount of oxygen 


   

    f = 1.875 + DeltaO2 + (acid_base_balance + 3.5) * tanh(0.5342 * DeltaO2);  %The acid-base balance is the bohr effect 
    sO2 = 1/(1 + exp(-f));

end

function cO2 = O2_dissociation(PO2,PCO2, pars)  
    alpha_O2 = pars("alpha_O2");
    cHb = pars("cHb");
    
    co2 =  CO2_dissociation(0, PCO2, pars);    

    henrys_law = alpha_O2 * PO2;  %plasma concentration
    sO2 = satO2(PO2,PCO2, co2, pars);  %hemoglobin concentration
    %cO2 = sO2;
    cO2 = henrys_law + cHb * sO2;
 end

function cCO2 = CO2_dissociation(PO2, PCO2, pars)  
    
    
    %Non linear solver
    %disp("Entrando")
    nonlinear_fun = @(co2) co2 - f(co2, PO2, PCO2, pars);
    co20first = pars('co20first');
    co20last = pars('co20last');
    
    
    co20 = [co20first,co20last];
    cCO2 = (1/5000 * PCO2) + 0.0154;
    % try
    %     cCO2 = fzero(nonlinear_fun, co20);
    % 
    % catch
    %     %because we know that if the previous line doesn´t work it means
    %     %that there's a problem with the simulation, so the best
    %     %approximation is to throw the linear function approx.
    %     cCO2 = (1/5000 * PCO2) + 0.0154;
    % end
    %disp(PCO2);
    function co2 = f(co2, PO2, PCO2, pars)
        %disp("Saliendo1")    
        cHb = pars('cHb');
        cHb_ery = pars('cHb_ery');
        alpha_CO2_ery = pars('alpha_CO2_ery');
        alpha_CO2_plasma = pars('alpha_CO2_plasma');
        
        sO2 = 0;%0;%0 * satO2(PO2, PCO2, cCO2, pars); %we will make the assumption that Haldane effect is almost null in comparision with the other effects. 
        
        
        pH = ph_model(co2, pars);
        
        pH_ery = 7.19 + 0.77 * (pH - 7.4) + 0.035 * (1 - sO2);
        pK_ery = 6.125 - log10(1 + 10^(pH_ery - 7.84 - 0.06 * sO2));
        pK_plasma = 6.125 - log10(1 + 10^(pH - 8.7));
    
        HCO3_ery = alpha_CO2_ery * PCO2 * 10^(pH_ery - pK_ery);
        HCO3_plasma = alpha_CO2_plasma * PCO2 * 10^(pH - pK_plasma);
    
        cCO2_ery = HCO3_ery + alpha_CO2_ery * PCO2;
        cCO2_plasma = HCO3_plasma + alpha_CO2_plasma * PCO2;

        co2 = cCO2_ery * cHb/cHb_ery + cCO2_plasma * (1 - cHb/cHb_ery);   
       
        % We will have to use a non linear method to obtain the cCO2 values
    end

    

    

end

function cCO2t = CO2_tissue(PO2,PCO2, pars)
    %Non linear solver
    nonlinear_fun = @(co2) co2 - f(co2, PCO2, pars);
    co20first = pars('co20first');
    co20last = pars('co20last');
    co20 = [co20first, co20last];
    cCO2t = (1/5000 * PCO2) + 0.0154;
    %cCO2t = fzero(nonlinear_fun, co20);
    function co2 = f(co2, PCO2, pars)
        alpha_CO2_plasma = pars('alpha_CO2_plasma');     
        pH = ph_model(co2, pars);
        pK_plasma = 6.125 - log10(1 + 10^(pH - 8.7));
        HCO3_plasma = alpha_CO2_plasma * PCO2 * 10^(pH - pK_plasma);
        cCO2_plasma = HCO3_plasma + alpha_CO2_plasma * PCO2;
        co2 = cCO2_plasma;
      % We will have to use a non linear method to obtain the cCO2 values
    end

end
function cO2t = O2_tissue(PO2,PCO2, pars)
    alpha_O2 = pars('alpha_O2');
    cO2t = alpha_O2 * PO2;
end
function pH = ph_model(cCO2, pars) 
    %disp("Saliendo3")     
    KaCO2 = pars('KaCO2');
    NaOH_0 = pars('NaOH_0');
    HPr_0 = pars('HPr_0');
    KaPr = pars('KaPr');
    %disp("Saliendo4") 
    g0 = KaCO2*KaPr *(NaOH_0 - HPr_0 - cCO2);
    g1 = KaCO2 * (NaOH_0 - cCO2) + KaPr*(KaCO2 + NaOH_0 - HPr_0);
    g2 = KaPr + NaOH_0 + KaCO2;
    
    protons = roots([1, g2, g1, g0]);
    
    
    filter = (protons > 10^(-9)) .* (protons < 10^(-6)) ;   %This ranges are veryyy important!!! if pH tolerances change the ph(cCO2) relations varies in important manner
    protons_ = protons' * (filter);
    
    
    %we have to select appropiate values for H+ values.
   
    
    pH = -log10(protons_);
        
        
end
%%



end

