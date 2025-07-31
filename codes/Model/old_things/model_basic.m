function xdot= model_basic(time, init_values, delays, pars, init_keys, taus_keys)

    %global hyperpars;
    global delays_global;
    global all_global;
    global externals_global;
    global dV_out;

    
    %delay driver
    index = @(keys, key) find(strcmp(keys, key));

    % if time > 1.9 && time < 2.1
    %     disp('en tiempo ');
    %     disp(time);
    %     disp('al principio')
    %     disp(externals_global(index(init_keys, 'Nt'), round(time/pars('dt'))));
    % end
    
    dt = pars('dt');
    tau = pars('tau');
    t0 = pars('t0');

    t = time;
    y = dictionary(init_keys, init_values);

    % if t == t0  %this is to reset the all_global in the new cycle when we want to re run the model basic
    %     all_global(:,round(t/dt) + 1:end) = 0;
    %     delays_global(:,round(t/dt) + 1:end) = 0; 
    % end
    


    if ((time - t0) < tau)  && ((time - t0) > 0) && ((time - tau ) > 0)
        delays = all_global(:,round((time-tau)/dt) + 1);  %this could be optimized but we are just using one tau, so it doesn't matter too much
        
    end
    
    
    delays_global(:,:,round(time/dt) + 1) = delays;
    externals_global(:,round(time/dt) + 1) = all_global(:,round(time/dt) + 1);
    all_global(:,round(time/dt) + 1) = init_values;

    prev_delays_global = delays_global(:,1:round(time/dt) + 1);
    delays_global(:,1:round(time/dt) + 1) = prev_delays_global + (prev_delays_global == 0).* delays;
    
    prev_all_global = all_global(:,1:round(time/dt) + 1);
    all_global(:,1:round(time/dt) + 1) = prev_all_global + (prev_all_global == 0).* init_values;
   
     
    %at the end we are going to interpolate, because otherwise it will take too much time
    
    
    
    
      
    



   
    
    

    

  %% Physiology
     
    
    [dPabd, dPtor] = respiratory_pump(t, y, pars);

    [y, ddVua, Gaw, dPua, dVua] = upper_airways(t, y, pars);
    [y, dPmusc, dPpl, dV] = pulmonary_mechanics(t, y, pars);
    

    [I] = metabolic_regulation(t, y, pars);
    
    [y, dPAgas, ddPa, dP_1, dP_2, dP_3, dP_4, dP_5, a] = exchange_mixing(t, y, pars, delays,index);

    [dv, dMRtgas] = tissue(t, y, pars);
    [dPvbCO2, dPCSFCO2, dPbCO2] = brain(t, y, pars);
    [dMRR, dM_Rv, dMRv] = metabolism_dynamics(t, y, pars);

    %Mean values calculator
    [mean_dPaO2, mean_dPaCO2, mean_dPbCO2] = mean_values_computer(t, y, pars);

    %Control
    [y] = ventilation_control(t, y, pars);    %VE computation, through 'y' I pass dVE to the next ones
    [y] = neuromuscular_drive(t, y, pars, index, init_keys);   %it integrates all_global after ventilation_control updated it correctly
    
    [y] = respiratory_mechanical_work(t, y, pars, index, init_keys);
    
    %Not ode-computable variables
    %disp(dV);
    %hyperpars.dV = [hyperpars.dV; dV];
    
    
    % if t > 1.8 && t < 2.4
    %     disp('index:');
    %     disp(round(t/dt) + 1);
    %     disp('Nt:')
    %     disp(externals_global(index(init_keys, 'Nt'), round(t/dt) + 1));
    %     externals_global(externals_global ~= 0)'
    % end
    
 %% dxdt   
    xdot_dict = dictionary();
    xdot_dict("Pabd") = dPabd;
    xdot_dict("Ptor") = dPtor;
    xdot_dict("Pmusc") = dPmusc;
    xdot_dict("Ppl") = dPpl;
    xdot_dict("V") = dV;
    xdot_dict("dVua") = ddVua;
    xdot_dict("Vua") = dVua;
    xdot_dict("Gaw") = Gaw;   %fake derivative
    xdot_dict("Pua") = dPua;   
    xdot_dict("Nt") = y('Nt'); %fake derivative
    xdot_dict("dVE") = y('dVE'); %fake derivative

    xdot_dict("I") = I;      %fake derivative
    xdot_dict("vO2") = dv(1);
    xdot_dict("vCO2") = dv(2);
    xdot_dict("MRtO2") = dMRtgas(1);
    xdot_dict("MRtCO2") = dMRtgas(2);
    xdot_dict("PvbCO2") = dPvbCO2;
    xdot_dict("PCSFCO2") = dPCSFCO2;
    xdot_dict("PbCO2") = dPbCO2;
    xdot_dict("PAO2") = dPAgas(1);
    xdot_dict("PACO2") = dPAgas(2);
    xdot_dict("dPaO2") = ddPa(1);
    xdot_dict("dPaCO2") = ddPa(2);
    xdot_dict("PaO2") = y('dPaO2');
    xdot_dict("PaCO2") = y('dPaCO2');
    xdot_dict("aO2") = a(1);   %fake derivative
    xdot_dict("aCO2") = a(2);  %fake derivative
    xdot_dict("P_1O2") = dP_1(1);
    xdot_dict("P_1CO2") = dP_1(2);
    xdot_dict("P_2O2") = dP_2(1);
    xdot_dict("P_2CO2") = dP_2(2);
    xdot_dict("P_3O2") = dP_3(1);
    xdot_dict("P_3CO2") = dP_3(2);
    xdot_dict("P_4O2") = dP_4(1);
    xdot_dict("P_4CO2") = dP_4(2);
    xdot_dict("P_5O2") = dP_5(1);
    xdot_dict("P_5CO2") = dP_5(2);
    xdot_dict("MRR") = dMRR;
    xdot_dict("M_Rv") = dM_Rv;
    xdot_dict("MRv") = dMRv;
    xdot_dict("J") = y('J'); %fake derivative;
    xdot_dict('mean_PaO2') = mean_dPaO2; 
    xdot_dict('mean_PaCO2') = mean_dPaCO2; 
    xdot_dict('mean_PbCO2') = mean_dPbCO2; 
    xdot_dict('fake_TI') = y('TI'); %fake derivative
    xdot_dict('fake_Tresp') = y('Tresp'); %fake derivative

    
    

    xdot = zeros(1, length(init_keys));
    for i = 1:length(init_keys)
        try
            xdot(i) = xdot_dict(init_keys(i));
            
        catch
            xdot(i) = 0;
        end
    end
    
    
    
    xdot = xdot';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55



%% Functions

%EQUATIONS (CHENG-SARMIENTO-SERNA)

%Respiratory mechanics
        %For oscillators, as the oscillates depending on a variable time, we cannot
        %use modular operations 


function  [dPabd, dPtor] = respiratory_pump(t, y, pars)
        
        %pars definition
        Pthormax_n = pars('Pthormax');
        Pthormin_n = pars('Pthormin');
        Pabdmax_n = pars('Pabdmax');
        Pabdmin_n = pars('Pabdmin');
        VTn = pars('VTn');
        gthor = pars('gthor');
        gabd = pars('gabd');
        t0 = pars('t0');

        %var definition
        TI = y('TI');
        Tresp = y('Tresp');
        VT = y('VT');
        
        
            
        
        
        
        %Equations

        %s = mod(t, Tresp)/Tresp;
        s = (t - t0)/Tresp;
        TE = Tresp - TI;
        DVT = VT - VTn;
        Pabdmax = Pabdmax_n + gabd * DVT;
        Pabdmin = Pabdmin_n - gabd * DVT;
        Pthormax = Pthormax_n + gthor * DVT;
        Pthormin = Pthormin_n - gthor * DVT;

        %Pthor
        if s >= 0 && s < TI/Tresp
            Ptor = Pthormax - (Pthormax - Pthormin) * Tresp/TI * s;
        elseif s >= TI/Tresp && s < (TI + TE)/Tresp
            Ptor = Pthormax - (Pthormax - Pthormin) * (TI + TE - Tresp * s)/TE;
        elseif s >= (TI + TE)/Tresp %&& s <= 1
            Ptor = Pthormax;
        end
        

        %Pabd
        if s >= 0 && s < TI/2 * 1/Tresp
            Pabd = Pabdmax - (Pabdmax - Pabdmin) * Tresp/(TI/2) * s;
        elseif s >= TI/2 * 1/Tresp && s < TI/Tresp
            Pabd = Pabdmin;
        elseif s >= TI/Tresp && s < (TI + TE)/Tresp
            Pabd = Pabdmax - (Pabdmax - Pabdmin) * (TI + TE - Tresp*s)/TE;
        elseif s >= (TI + TE)/Tresp %&& s <= 1
            Pabd = Pabdmax;
        end
       

        
        dPabd = Pabd - y('Pabd');
        
        dPtor = Ptor - y('Ptor');
end

function [y_, ddVua, Gaw, dPua, dVua] = upper_airways(t, y, pars)
    
    %pars
    R_trachea = pars('Rtrachea');    
    Rl = pars('Rl');
    Rcw = pars('Rcw');  
    Raw = pars('Raw');
    Cua = pars('Cua');
    bua = pars('bua');
    Pcrit_min = pars('Pcrit_min');
    A0ua = pars('A0ua');
    Kua = pars('Kua');
    
    %vars
    Ppl = y('Ppl');
    dV = y('dV');
    dVua = y('dVua');
    
    Rrs = Raw + Rl + Rcw;
    dVla = dVua + dV;

    Pua = Ppl + dVla * Rrs;
    dPua = (Pua - y('Pua'))/dt;
    

    ddVua = -1/R_trachea * (dPua + dVua/Cua);

    if 10 <= -1/(Cua * bua)
        Pcrit = 10;
    elseif Pcrit_min < -1/(Cua * bua) && -1/(Cua * bua) < 10  
        Pcrit = -1/(Cua * bua);
    elseif Pcrit_min >= -1/(Cua * bua)
        Pcrit = Pcrit_min;
    end

    if Pua <= Pcrit
        Gaw = 0;
    elseif (Pcrit < Pua && Pua <= 0) && ((1 - Pua/Pcrit) >= 0)
        Gaw = A0ua * (1 - Pua/Pcrit) * Kua;

    elseif Pua > 0
        Gaw = A0ua * Kua;
    end

    y('Gaw') = Gaw;
    y_ = y;
end
function [y_, dPmusc, dPpl, dV] = pulmonary_mechanics(t, y, pars)
    
    

    %pars definition
    Ecw = pars('Ecw');
    El = pars('El');
    kaw1  = pars('kaw1');
    kaw2 = pars('kaw2');
    Rcw = pars('Rcw');
    Rrs = pars('Rrs');
    Pao = pars('Pao');

    KpCO2 =pars('KpCO2');
    KcCO2 = pars('KcCO2');
    KpO2 = pars('KpO2');
    KcMRv = pars('KcMRv');
    Kbg = pars('Kbg');  
    dt = pars('dt');

    t0 = pars('t0');


    %vars definition
    
    V = y('V');
    %dV = y('dV'); %hyperpars.dV(end);%y('dV');
    Gaw = y('Gaw');
    
    a0 = y('a0');
    a1 = y('a1');
    a2 = y('a2');
    TI = y('TI');
    Tresp = y('Tresp');
    tau = y('tau');
    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');
    PbCO2 = y('PbCO2');
    MRv = y('MRv');



    Ers = Ecw + El;
    
    %Pmusc
    %t = mod(t, Tresp);   %we are using this only when controllers are off
    t = t - t0;
    if 0 <= t && t <= TI
        Pmusc = a0 + a1*t + a2*t^2;
        

    elseif  TI < t 
        PmuscTI = a0 + a1*TI + a2*TI^2;
        Pmusc = PmuscTI * exp(-(t - TI)/tau);
        
    end

    
    % dV
    %control_gain = 1 + 0*abs(KpCO2 * PaCO2 + KcCO2 * PbCO2 + (KpO2 * (104 - PaO2)^4.9) * (PaO2 > 104) + KcMRv * MRv - Kbg);
    %disp(control_gain);
    dV =  Gaw/Rrs * ((Pmusc - Pao) - Ers * V); %Gaw/Rrs * ((Pmusc - Pao) - Ers * V);
    
    dV_out.values = [dV_out.values dV];
    dV_out.time_points = [dV_out.time_points t];
    %Pcw
    if dV < 0
        Pcw = Ecw * V - 1;
        Pa_ = Pao;
    else
        Pcw = Ecw * V - 1 + Rcw * dV;
        Pa_ = Pao - kaw1 * dV - kaw2 * abs(dV)^2;
        %Pa = Pa_ * (Pa_ > 0);
        
    end
    Pa = Pa_ * (Pa_ > 0);
    
    %Ppl
    Ppl = Pcw + Pa - Pmusc;                  

    dPpl = (Ppl - y('Ppl'))/dt;
    dPmusc = (Pmusc - y('Pmusc'))/dt;

    y('dV') = dV;
    y_ = y;
    
    
    



end

function [y_] = neuromuscular_drive(t, y, pars, index, init_keys)
    
    %pars
    dt = pars('dt');
    t0 = pars('t0');
    %---

    %vars
    Tresp = y('Tresp');
    TI = y('TI');
    %dVE = y('dVE');   

    
    t_cycle = t - t0; 
    
    % if t0 > 0
    %     disp(t_cycle);
    %     disp(t);
    % end
    Nt = eps;
    
    if t_cycle < TI
        dVE_historic = all_global(index(init_keys, 'dVE'), round(t0/dt) + 1: round(t/dt) + 1);
        Nt = trapz(dVE_historic, 2) * dt;
        %Nt = 1;
    %Nt = trapz(dVE_historic, 2) * dt;
    %Nt = 0;   
    % elseif t_cycle >= TI
    %     Nt = 2;
    end 
    
    %disp(Nt);   
    y('Nt') = Nt;  %dejar derivada en 0 y guardar todas las variables externas en una global
    y_ = y;
    %externals_global(index(init_keys, 'Nt'), round(t/dt) + 1) = Nt;

     
    
    
    
end
function [I] = metabolic_regulation(t, y, pars)
    %pars
    MRtCO2_basal = pars("MRtCO2_basal");
    AT = pars("AT");

    %vars
    MRtCO2 = y("MRtCO2");

    I = (MRtCO2 - MRtCO2_basal)/(AT - MRtCO2_basal);
    %dI = I - y("I");

end

%Gas exchange
function gas = dissociation(P, pars)
    
    
    a1 = pars('a1'); % parameter in O2 dissociation equation
    a2 = pars('a2'); % parameter in CO3 dissociation equation
    alpha1 = pars('alpha1'); % parameter in O3 dissociation equation
    alpha2 = pars('alpha2'); % parameter in CO3 dissociation equation
    K1 = pars('K1'); % parameter in O3 dissociation equation
    K2 = pars('K2'); % parameter in CO3 dissociation equation
    beta1 = pars('beta1'); % parameter in O3 dissociation equation
    beta2 = pars('beta2'); % parameter in CO3 dissociation equation
    C1 = pars('C1');
    C2 = pars('C2');
    Z = pars('Z');
    

    CO2a_ = C2 * Z;
    O2a_ = C1 * Z;
    
    PO2 = P(1);
    PCO2 = P(2);

    FCO2 = PCO2 * (1 + beta2 * PO2)/(K2 * (1 + alpha2 * PO2));
    CO2 = CO2a_ * FCO2^(1/a2)/(1 + FCO2^(1/a2));
    
    FO2 = PO2 * (1 + beta1 * PCO2)/(K1 * (1 + alpha1 * PCO2));
    O2 = O2a_ * FO2^(1/a1)/(1 + FO2^(1/a1));
    
    gas = [O2, CO2];
    
end
function [y_, dPAgas, ddPa, dP_1, dP_2, dP_3, dP_4, dP_5, a] = exchange_mixing(t, y, pars, delays, index)
    

    %taus
    y_taugas = delays(:, index(taus_keys, 'tau_gases')); 
    PAO2_delayed = y_taugas(index(init_keys, 'PAO2'));
    PACO2_delayed = y_taugas(index(init_keys, 'PACO2'));
    
    %pars
    fO2 = pars('fO2');
    fCO2 = pars('fCO2');
    Patm = pars('Patm');
    Pws = pars('Pws');
    Vdead = pars('Vdead');
    VLO2 = pars('VLO2');
    VLCO2 = pars('VLCO2');
    T1 = pars('T1');
    T2 = pars('T2');
    Ta = pars('Ta');

    
    %vars
    
    P_1O2 = y('P_1O2');
    P_1CO2 = y('P_1CO2');
    P_2O2 = y('P_2O2');
    P_2CO2 = y('P_2CO2');
    P_3O2 = y('P_3O2');
    P_3CO2 = y('P_3CO2');
    P_4O2 = y('P_4O2');
    P_4CO2 = y('P_4CO2');
    P_5O2 = y('P_5O2');
    P_5CO2 = y('P_5CO2');

    PAO2 = y('PAO2');
    PACO2 = y('PACO2');

    dV = y('dV');
    V = y('V');

    vO2 = y('vO2');
    vCO2 = y('vCO2');    

    Qpp = y('Qpp');

    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');
    dPaO2 = y('dPaO2');
    dPaCO2 = y('dPaCO2');

    fgas = [fO2, fCO2];
   

    P_1 = [P_1O2, P_1CO2];
    P_2 = [P_2O2, P_2CO2];
    P_3 = [P_3O2, P_3CO2];
    P_4 = [P_4O2, P_4CO2];
    P_5 = [P_5O2, P_5CO2];
    
    PAgas = [PAO2, PACO2];
    PAgas_delayed = [PAO2_delayed, PACO2_delayed];
    %PAgas_delayed = PAgas;

    v = [vO2, vCO2];
    

    Pa = [PaO2, PaCO2];   
    dPa = [dPaO2, dPaCO2];   
    VL = [VLO2, VLCO2]; 

    % Inspired Air
    PI = fgas * (Patm - Pws)/100;
    
    
    % Deadspace flow
    if dV >= 0    %Inspiration
    
        dP_1 = (PI - P_1) * abs(dV)/Vdead * 1/0.2;
        dP_2 = (P_1 - P_2) * abs(dV)/Vdead * 1/0.2;
        dP_3 = (P_2 - P_3) * abs(dV)/Vdead * 1/0.2;
        dP_4 = (P_3 - P_4) * abs(dV)/Vdead * 1/0.2;
        dP_5 = (P_4 - P_5) * abs(dV)/Vdead * 1/0.2;
    else        %Expiration
    
        dP_1 = (P_2 - P_1) * abs(dV)/Vdead * 1/0.2;
        dP_2 = (P_3 - P_2) * abs(dV)/Vdead * 1/0.2;
        dP_3 = (P_4 - P_3) * abs(dV)/Vdead * 1/0.2;
        dP_4 = (P_5 - P_4) * abs(dV)/Vdead * 1/0.2;
        dP_5 = (PAgas - P_5) * abs(dV)/Vdead * 1/0.2;
    end
    
    
    % dissociation
    a = dissociation(Pa, pars); %we are going to asume that arterial gas partial pressures are alveolar gas partial pressure
    %da = a - [y('aO2'), y('aCO2')];
    
    % alveolar exchange 
    
    if dV >= 0
        dPAgas = (863 * Qpp/1000 * (v - a) + dV * (P_5 - PAgas)) * (eye(2) *  1./(VL + V)); % ODE for alveolar CO2 partial pressure %Here we have a 1/1000 to make units corrections (from ml to L)
        
    else
        
        dPAgas = (863 * Qpp/1000 * (v - a)) * (eye(2) *  1./(VL + V));    %Here we have a 1/1000 to make units corrections (from ml to L)
    end
    
    
    % mixing
    
    ddPa = 1/(T1 * T2) * (PAgas_delayed - (T1 + T2)*dPa - Pa);  %here we should take PAgas(t-Ta)
    


    y('aO2') = a(1);
    y('aCO2') = a(2);
    y_ = y;
end
function [dPvbCO2, dPCSFCO2, dPbCO2] = brain(t, y, pars)
    
    %pars
    KCSFCO2 = pars('KCSFCO2');
    KCCO2 = pars('KCCO2');
    dc = pars('dc');
    h = pars('h');
    SCO2 = pars('SCO2');
    SbCO2 = pars('SbCO2');
    MRbCO2 = pars('MRbCO2');

    %vars
    PvbCO2 = y('PvbCO2');
    PCSFCO2 = y('PCSFCO2');
    PbCO2 = y('PbCO2');
    PaCO2 = y('PaCO2');    
    Qbp = y('Qbp');
    
    dPvbCO2 = (MRbCO2 + Qbp * SCO2 * (PaCO2 - PvbCO2) - h)/SbCO2;
    dPCSFCO2 = (PvbCO2 - PCSFCO2)/KCSFCO2;
    PbCO2 = PvbCO2 + (PCSFCO2 - PvbCO2) * exp(-dc * (Qbp * KCCO2)^0.5);
    dPbCO2 = PbCO2 - y("PbCO2");
end

function [dv, dMRtgas] = tissue(t, y, pars)
    
    %Input variables- consumption-production rates
    pars = input_consumption(t, y, pars);

    %pars
    tauMR = pars('tauMR');    
    Vtissue_CO2 = pars('Vtissue_CO2');
    Vtissue_O2 = pars('Vtissue_O2');
    MRO2 = pars('MRO2');   %This parameters correspond to the input values.  
    MRCO2 = pars('MRCO2'); %This parameters correspond to the input values.

    %vars
    MRtO2 = y('MRtO2');
    MRtCO2 = y('MRtCO2');    
    vO2 = y('vO2');
    vCO2 = y('vCO2');
    aO2 = y('aO2');
    aCO2 = y('aCO2');
    Qpp = y('Qpp');
    Qbp = y('Qbp');
    
    MRtgas = [MRtO2, MRtCO2];
    MRgas = [MRO2, MRCO2];
    a = [aO2, aCO2];
    v = [vO2, vCO2];
    Vtissue = [Vtissue_O2, Vtissue_CO2];

    Qt = Qpp - Qbp;   
    dMRtgas = (MRgas - MRtgas)/tauMR;  
    dv = (MRtgas * eye(2).*[-1,1] + Qt*(a-v)) * eye(2)./Vtissue;

    function pars_ = input_consumption(t, y, pars)
        %pars
        MRO2 = pars('MRO2');
        MRCO2 = pars('MRCO2');
        MRtO2_basal = pars('MRtO2_basal');
        MRtCO2_basal = pars('MRtCO2_basal');
        
        if t < 6
            MRO2 = 2*MRtO2_basal;
            MRCO2 = 2 * MRtCO2_basal;
        else
            MRO2 = MRtO2_basal;
            MRCO2 = MRtCO2_basal;
        end
        
        pars('MRO2') = MRO2;
        pars('MRCO2') = MRCO2;

        pars_ = pars;      
        
        
    end

end
function [dMRR, dM_Rv, dMRv] = metabolism_dynamics(t, y, pars)
    %pars
    MRbCO2 = pars('MRbCO2');
    MRbO2 = pars('MRbO2');
    MRtCO2_basal = pars('MRtCO2_basal');
    MRtO2_basal = pars('MRtO2_basal');
    tauMRv = pars('tauMRv');

    dt = pars('dt');
    %vars
    
    
    M_Rv = y('M_Rv');
    MRtO2 = y('MRtO2');
    MRtCO2 = y('MRtCO2');

    
    M_RR = (MRbCO2 + MRbO2 + MRtCO2 + MRtO2)/(MRbCO2 + MRbO2 + MRtCO2_basal + MRtO2_basal);
    if M_RR >= 1
        MRR = M_RR;
    else
        MRR = 1;
    end

    dMRR = MRR - y('MRR');

    dM_Rv = ((MRR - 1) - M_Rv)/tauMRv;

    if M_Rv >= 0 && MRR > 1
        MRv = M_Rv;
    else
        MRv = 0;
    end

    dMRv = (MRv - y('MRv'))/dt;
end







%Ventilatory controller

%alveolar ventilation
function  [y_] = ventilation_control(t, y, pars)
    
    %pars
    
    KpCO2 = pars('KpCO2');
    KpO2 = pars('KpO2');
    KcMRv = pars('KcMRv');
    KcCO2 = pars('KcCO2');
    Kbg = pars('Kbg');

    dV_rest = pars('dV_rest');
    V0dead = pars('Vdead');
    GVdead = pars('GVdead');

    %vars
    Tresp = y('Tresp');
    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');
    PbCO2 = y('PbCO2');
    MRv = y('MRv');
    mean_PaO2 = y('mean_PaO2');
    mean_PaCO2 = y('mean_PaCO2');
    mean_PbCO2 = y('mean_PbCO2');
    
    dVA_ = dV_rest * (KpCO2 * mean_PaCO2 + KcCO2 * mean_PbCO2 + (KpO2 * (104 - mean_PaO2)^4.9) * (mean_PaO2 < 104) + KcMRv * MRv - Kbg); %This should be alveolar minute volume, and it's part of the minute ventilation we want the lungs to adquire
    dVA = dVA_ * (dVA_ > 0); % as minute ventilation is a positive value, we must take the absolute value, we will never be able to remove air from the lungs, the minimum volume value will always be dead space volume 
    Vd = GVdead * dVA + V0dead; %the first part refers to dead space volume that changes because respiration (expansion of airway channels?)
    dVd = 1/Tresp * Vd;     %the amount of volume that is exchanged during a minute due to dead space, will be the breathing frequency times dead space volume
    
    dVE = dVA + dVd;    %dVE corresponds to minute ventilation (how much volume do I want to exchange over a minute, it's a indicator of respiration flow and its different than instant flow)
    y("dVE") = dVE;
    %all_global(index(init_keys, 'dVE'), round(t0/dt) + 1: round(t/dt) + 1) = dVE;
    all_global(index(init_keys, 'dVE'), round(t/dt) + 1) = dVE;
    y_ = y;

    %In the simulation is correct kind of the values we got (0.11 l/s, which are indeed 7 l/min, reported as common minute ventilation values)

end

function [y_] = respiratory_mechanical_work(t, y, pars, index, init_keys)
    
    
    % Parameters
    n = pars('n');
    lambda1 = pars('lambda1');
    lambda2 = pars('lambda2');   
    Pmax = pars('Pmax');
    dPmax = pars('dPmax');

    Rrs = pars('Rrs');
    Ecw = pars('Ecw');
    El = pars('El');
    Ers = Ecw + El;

    dt = pars('dt');
    t0 = pars('t0');
    
    %vars 
    Tresp = y('Tresp');
    TI = y('TI');
    a0 = y('a0');
    a1 = y('a1');
    a2 = y('a2');
    tau = y('tau');

    Pmusc = y('Pmusc');
    dVE = y('dVE');

    t_cycle = t - t0;
    
    try  %avoiding out of range problem
        Pmusc_prev = all_global(index(init_keys, 'Pmusc'), round((t-dt)/dt) + 1);
        dVE_prev = all_global(index(init_keys, 'dVE'), round((t-dt)/dt) + 1);

        dPmusc_dt = (Pmusc - Pmusc_prev)/dt;    
        ddVE_dt = (dVE - dVE_prev)/dt;
    catch
        dPmusc_dt = 0;
        ddVE_dt = 0;
    end
    
    
    
    
    

    

    zheta1 = 1 - Pmusc/Pmax;
    zheta2 = 1 - dPmusc_dt/dPmax;

    insp_integrand = Pmusc/(zheta1^n * zheta2^n) + lambda1 * ddVE_dt^2;
    exp_integrand = ddVE_dt^2;

    

    %Solver based integration
    solver_based_integration = 0;
    if solver_based_integration
        insp_integrand_0 = all_global(index(init_keys, 'insp_integrand'), round(t0/dt) + 1);
        exp_integrand_0 = all_global(index(init_keys, 'exp_integrand'), round(t0/dt) + 1);

        dinsp_work_power = 1/Tresp * (insp_integrand - insp_integrand) * (t < TI);
        dexp_work_power = 1/Tresp * (exp_integrand - exp_integrand) * (t > TI);
        insp_work_power = y('insp_work_power');
        exp_work_power = y('exp_work_power');
        resp_work_power = insp_work_power + lambda2 * exp_work_power;
        Jfinal = resp_work_power;
    end
    %Testear después de correr los cambios implmentados

    

    all_global(index(init_keys, 'insp_integrand'), round(t/dt) + 1) = insp_integrand;
    all_global(index(init_keys, 'exp_integrand'), round(t/dt) + 1) = exp_integrand;
    
    try %this is to avoid 'index out of range', if that happens, it is forced to shrink Tresp
        insp_work_power = 1/Tresp * dt * sum(all_global(index(init_keys, 'insp_integrand'), round(t0/dt) + 1: round((t0 + TI)/dt) + 1)); 
        exp_work_power =  1/Tresp * dt * sum(all_global(index(init_keys, 'exp_integrand'), round((t0 + TI)/dt) + 1: round((t0 + Tresp)/dt) + 1)); %it should be until Tresp, but that happens at the end of the cycle
        resp_work_power = insp_work_power + lambda2 * exp_work_power;
        Jfinal = resp_work_power;

        
    catch
        Jfinal = all_global(index(init_keys, 'J'), round(t/dt));%10^5; %big M.
    end

    y('insp_integrand') = insp_integrand;
    y('exp_integrand') = exp_integrand;
    y('J') = Jfinal;
    all_global(index(init_keys, 'J'), round(t/dt) + 1) = Jfinal;

    y_ = y;

    

end

function [mean_dPaO2, mean_dPaCO2, mean_dPbCO2] = mean_values_computer(t, y, pars)
    
    %vars
    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');
    PbCO2 = y('PbCO2');

    mean_PaO2 = y('mean_PaO2');
    mean_PaCO2 = y('mean_PaCO2');
    mean_PbCO2 = y('mean_PbCO2');

    %means - we want a cut off filter that leaves signals at 0.1Hz range (to take the mean value)
    mean_dPaO2 = (PaO2 - mean_PaO2)/5;
    mean_dPaCO2 = (PaCO2 - mean_PaCO2)/5;
    mean_dPbCO2 = (PbCO2 - mean_PbCO2)/5;

    
end




end



