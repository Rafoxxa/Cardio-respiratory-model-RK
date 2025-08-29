%function xdot= model_basic_vascular(time, init_values, delays, pars, init_keys, taus_keys)
function xdot= model_basic_vascular(t, y, pars, init_keys)
    time = t;
    init_values = y;    
    global all_global;
    
    % Replace dictionary with containers.Map for MATLAB 2017 compatibility
    if isa(pars, 'containers.Map') 
        internal_variables = containers.Map(); %this array is for local variables that are used per iteration between different modules        
    else
        internal_variables = zeros(70);
    end
    
    %Basic parameters    
    dt = pars(248);
    %tau = pars(318);
    t0 = pars(317);
    t = time;
    

    if isa(pars, 'containers.Map') 
        y = containers.Map(init_keys, num2cell(init_values));
    else
        y = init_values;
    end

    %y = init_values

    %Functions to take variables values based on their names
    index = @(keys, key) find(strcmp(keys, key));
    index_fun = @(keys, key) find(strcmp(keys, key));

    % Delayed variables and necesary for integration variables storage
    [all_global, tiny_y_keys] = saving_in_globals(all_global, y);

 %% Testing variables: internal variables and other non state variables sometimes needs to be plotted, for that we create forced variables to stored them in the result's vector    
    test = containers.Map();
    test__vascular_checking__28_06_24 = containers.Map();

    

  %% Physiology Equations
    %''' General form:  [state_vector, internal_variables, dx1, dx2, ..., dxn, test] = physiological_fun(time, state_vector, parameters, internal_variables, test) '''
     
    %Adding the base value to control variables: as internal variables doesn't have a derivative, and because they are needed in multiple equations, we want to compute them at the beginning.  
    [internal_variables] = adding_base_values_in_control_variables(t, y, pars, internal_variables);
    
    %Respiratory (gases and lung's pressures)
    [internal_variables, test] =                               respiratory_pump(t, y, pars, internal_variables, test); % [5]
    [y, ddVua, dGaw, dVua, internal_variables, test] =         upper_airways(t, y, pars, internal_variables, test); % [2][6][7]
    [y, dPmusc, dPpl, dV, ddV, test] =                         pulmonary_mechanics(t, y, pars, internal_variables, test); % [1][2][6][7]
    [I, internal_variables, test] =                            metabolic_regulation(t, y, pars, internal_variables, test); % [0]   
    [y, dPAgas, ddPa, dP_1, dP_2, dP_3, dP_4, dP_5, a, test] = exchange_mixing(t, y, pars,index, test); % [1][2][9]
    [dv, dMRtgas, test] =                                      tissue(t, y, pars, test); % [0][2][8]
    [dPvbCO2, dPCSFCO2, internal_variables, test] =            brain(t, y, pars, internal_variables, test); %[2][3]
    [dMRR, dM_Rv, dMRv, test] =                                metabolism_dynamics(t, y, pars, test); % [8]

    %Mean values calculator
    [mean_dPaO2, mean_dPaCO2, mean_dPbCO2, mean_dP_sa, test] = mean_values_computer(t, y, pars, internal_variables, test); %[0]

    %Ventilatory Control
    [ddVE, test] =                                             ventilation_control(t, y, pars, internal_variables, test); %[3][8]    %VE computation
    [y, internal_variables, test] =                            neuromuscular_drive(t, y, pars, internal_variables, index, tiny_y_keys, test); %[?]     

    %Cardiovascular
    [internal_variables, test] =                                                                                  muscle_pump(t, y, pars, internal_variables, test); % [5]
    [internal_variables, test] =                                                                                  vena_cava(t, y, pars, internal_variables, test); % [2]       
    [internal_variables, dzheta_heart, dV_total_rv, dV_total_ra, dV_total_la, dV_total_lv, dQla, dTheart, test] = heart(t, y, pars, internal_variables, test); % [0][2][10]
    [internal_variables, dP_sa, dQ_sa, test] =                                                             systemic_arteries(t, y, pars, internal_variables, test); %[2][5] 
    [internal_variables, dQ_pa, dV_total_pp, dV_total_pv, dV_total_pa, dQpp, test] =                              pulmonary_circulation(t, y, pars, internal_variables, test); %[2][5]
    [internal_variables, dV_total_vc, dV_total_v, dP_sp, dQbp, test__vascular_checking__28_06_24] =                                            systemic_peripheric_and_venous_circulation(t, y, pars, internal_variables, test__vascular_checking__28_06_24); % [2][5]
    
    %Cardiovascular Control
    %Aferent
    [dP_mean, internal_variables] =                             afferent_barorreflex(t,y,pars, internal_variables, dP_sa); % [2][11]
    [df_ac] =                                                   afferent_chemoreceptor(t,y,internal_variables); % [2][4]
    [df_ap] =                                                   afferent_pulmonary_stretch(t,y,pars, internal_variables); % [1][2]
    %lLocal blood flow control
    [dxO2_b, dxCO2_b, dR_bp, internal_variables] =              cerebral_blood_flow(t, y, pars, internal_variables); % [2][4][11]
    [dxO2_e, dxO2_s, dxO2_p, internal_variables] = hipoxia_local_regulation(t, y, pars, internal_variables);  %propose by us
    [dxO2, dxCO2, dWh, internal_variables] =                    coronary_and_resting_muscle_blood_flow(t, y, pars, internal_variables); %[2][4][11]
    [dxO2_am, dx_met, dx_M, dphi_met, internal_variables] =     active_muscle_blood_flow(t,y,pars, internal_variables, index_fun); % [5]
    %Integration in SNC
    [dDThetaO2_s, dDThetaCO2_s, internal_variables] =           cns_ischemic_response(t, y, pars, internal_variables); %[2][4]
    %Eferent 
    [internal_variables] =                                      efferent_pathways(t, y, pars, internal_variables);    %[1][2][4][5][11]
    [dDTheta, dfh_s, dfp_s, dfv_s, internal_variables] =        reflex_control_R_Vu_E(t, y, pars, internal_variables, index_fun); %[2][11]
    [dDTsym, dDTvagal, dfv] =                                   reflex_control_HR(t, y, pars, internal_variables, index_fun); %[2][11]
    
  %% dxdt building   
    xdot_dict = containers.Map();
    if isa(y, 'containers.Map')
        xdot = zeros(1, length(y.keys));   
    else
        xdot = zeros(1, length(y));  
    end
    xdot(55) =                    dPmusc;
    xdot(56) =                      dPpl;
    xdot(82) =                        dV;
    xdot(106) =                       ddV;
    xdot(108) =                     ddVua;
    %xdot_dict('Vua') =                      dVua;
    xdot(21) =                      dGaw;   %fake derivative
    %xdot(57) =                     dPua;   
    %xdot(30) =                      y(30); %fake derivative
    xdot(107) =                      ddVE; %fake derivative
    xdot(23) =                        I;      %fake derivative
    xdot(132) =                      dv(1);
    xdot(131) =                     dv(2);
    xdot(27) =                    dMRtgas(1);
    xdot(26) =                   dMRtgas(2);
    xdot(58) =                   dPvbCO2;
    xdot(33) =                  dPCSFCO2;    
    xdot(32) =                     dPAgas(1);
    xdot(31) =                    dPAgas(2);
    xdot(105) =                    ddPa(1);
    xdot(104) =                   ddPa(2);
    xdot(53) =                     y(105);
    xdot(52) =                    y(104);
    xdot(103) =                      a(1);   %fake derivative
    xdot(102) =                     a(2);  %fake derivative
    xdot(35) =                    dP_1(1);
    xdot(34) =                   dP_1(2);
    xdot(37) =                    dP_2(1);
    xdot(36) =                   dP_2(2);
    xdot(39) =                    dP_3(1);
    xdot(38) =                   dP_3(2);
    xdot(41) =                    dP_4(1);
    xdot(40) =                   dP_4(2);
    xdot(43) =                    dP_5(1);
    xdot(42) =                   dP_5(2);
    xdot(25) =                      dMRR;
    xdot(29) =                     dM_Rv;
    xdot(28) =                      dMRv;
    xdot(123) =                mean_dPaO2; 
    xdot(122) =               mean_dPaCO2; 
    xdot(124) =               mean_dPbCO2; 
    xdot(121) =                mean_dP_sa; 
    xdot(113) =                  y(79); %fake derivative
    xdot(114) =               y(81); %fake derivative
    xdot(62) =                     dQ_pa; 
    xdot(91) =               dV_total_pp;
    xdot(92) =               dV_total_pv;
    xdot(90) =               dV_total_pa;
    xdot(97) =               dV_total_vc;
    xdot(65) =                     dQ_sa;
    xdot(45) =                     dP_sa;
    xdot(86) =              dV_total_v(1);
    xdot(96) =              dV_total_v(2);
    xdot(85) =              dV_total_v(3);
    xdot(87) =              dV_total_v(4);
    xdot(94) =             dV_total_v(5);
    xdot(84) =             dV_total_v(6);
    xdot(46) =                     dP_sp;
    xdot(95) =               dV_total_rv;
    xdot(93) =               dV_total_ra;
    xdot(88) =               dV_total_la;
    xdot(89) =               dV_total_lv;
    xdot(157) =              dzheta_heart;
    xdot(44) =                   dP_mean;
    xdot(111) =                     df_ac;
    xdot(112) =                     df_ap;
    xdot(149) =                    dxO2_b;
    xdot(142) =                   dxCO2_b;
    xdot(151) =                    dxO2(1);
    xdot(153) =                   dxO2(2);
    xdot(144) =                   dxCO2(1);
    xdot(146) =                  dxCO2(2);
    xdot(98) =                       dWh;
    xdot(148) =                   dxO2_am;
    xdot(156) =                    dx_met;
    xdot(155) =                      dx_M;
    xdot(4) =             dDThetaO2_s(1);
    xdot(5) =             dDThetaO2_s(2);
    xdot(6) =             dDThetaO2_s(3);
    xdot(1) =            dDThetaCO2_s(1);
    xdot(2) =            dDThetaCO2_s(2);
    xdot(3) =            dDThetaCO2_s(3);
    xdot(17) =                    dDTsym;
    xdot(18) =                  dDTvagal;
    xdot(10) =             dDTheta(1);
    xdot(12) =             dDTheta(2);
    xdot(11) =          dDTheta(3);
    xdot(9) =          dDTheta(4);
    xdot(14) =  dDTheta(5);
    xdot(16) =  dDTheta(6);
    xdot(15) = dDTheta(7);
    xdot(13) = dDTheta(8);
    xdot(7) =           dDTheta(9);
    xdot(8) =           dDTheta(10);
    xdot(125) =                  dphi_met;
    xdot(68) =                      dQpp;  %forced derivative
    xdot(66) =                      dQbp;  %forced derivative
    xdot(67) =                      dQla;  %forced derivative
    xdot(80) =                   dTheart;
    xdot(115) =                     dfh_s;    
    xdot(116) =                     dfp_s;    
    xdot(118) =                     dfv_s;
    xdot(117) =                       dfv;    
    xdot(72) =                     dR_bp; 
    xdot(150) =                     dxO2_e;
    xdot(154) =                     dxO2_s;
    xdot(152) =                     dxO2_p;   
    % xdot(60) =                      test__vascular_checking__28_06_24('dQ_e');
    % xdot(64) =                      test__vascular_checking__28_06_24('dQ_s');
    % xdot(61) =                      test__vascular_checking__28_06_24('dQ_h');
    % xdot(63) =                     test__vascular_checking__28_06_24('dQ_rm');
    % xdot(59) =                     test__vascular_checking__28_06_24('dQ_am');
    % xdot(48) =                    test__vascular_checking__28_06_24('dP_v_e');
    % xdot(51) =                    test__vascular_checking__28_06_24('dP_v_s');
    % xdot(49) =                    test__vascular_checking__28_06_24('dP_v_h');
    % xdot(50) =                   test__vascular_checking__28_06_24('dP_v_rm');
    % xdot(47) =                   test__vascular_checking__28_06_24('dP_v_am');
    % xdot(73) =                    test__vascular_checking__28_06_24('dR_e_p');
    % xdot(77) =                    test__vascular_checking__28_06_24('dR_s_p');
    % xdot(71) =                    test__vascular_checking__28_06_24('dR_b_p');
    % xdot(74) =                    test__vascular_checking__28_06_24('dR_h_p');
    % xdot(75) =                   test__vascular_checking__28_06_24('dR_rm_p');
    % xdot(69) =                   test__vascular_checking__28_06_24('dR_am_p'); 

    %Transforming dictionary to matlab vector, only for dict mode
    %%if length(xdot_dict) > 2
    %    
    %    for i = 1:length(init_keys)
    %        try
    %            xdot(i) = xdot_dict(init_keys(i));
    %            
    %        catch
    %            xdot(i) = 0;
    %        end
    %    end
    %    
    %    debugging_tools(xdot, y, init_keys);
    %%end
%
    xdot = xdot';
    if sum(imag(xdot)) > 0
        disp('error');
        klj = asdd(2);
        
    end

%
    %xdot = (xdot > 10^6)  * 10^6  + xdot * (xdot <= 10^6); %protection measure against high derivatives, mainly for the fitting process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55



%% Functions

%EQUATIONS (CHENG-SARMIENTO-SERNA)

%% Respiratory mechanics   Eq [x] --- [n]


function  [internal_variables_, test_] = respiratory_pump(t, y, pars, internal_variables, test)
        
    %pars definition
    Pthormax_n = pars(128);
    Pthormin_n = pars(129);
    Pabdmax_n = pars(122);
    Pabdmin_n = pars(123);
    VTn = pars(174);
    gthor = pars(288);
    gabd = pars(276);
    t0 = pars(317);

    %var definition
    TI = y(79);
    Tresp = y(81);
    VT = y(82);

    %Equations   
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
    elseif s >= (TI + TE)/Tresp 
        Ptor = Pthormax;
    end    

    %Pabd
    if s >= 0 && s < TI/2 * 1/Tresp
        Pabd = Pabdmax - (Pabdmax - Pabdmin) * Tresp/(TI/2) * s;
    elseif s >= TI/2 * 1/Tresp && s < TI/Tresp
        Pabd = Pabdmin;
    elseif s >= TI/Tresp && s < (TI + TE)/Tresp
        Pabd = Pabdmax - (Pabdmax - Pabdmin) * (TI + TE - Tresp*s)/TE;
    elseif s >= (TI + TE)/Tresp 
        Pabd = Pabdmax;
    end 
    
    internal_variables(44) = Pabd;
    internal_variables(35) = Ptor;
    internal_variables_ = internal_variables;    
    test_ = test;
end

function [y_, ddVua, dGaw, dVua, internal_variables_, test_] = upper_airways(t, y, pars, internal_variables, test)
    
    %pars
    R_trachea = pars(155);    
    Rl = pars(153);
    Rcw = pars(152);  
    Raw = pars(151);
    Cua = pars(30);
    bua = pars(231);
    Pcrit_min = pars(126);
    A0ua = pars(2);
    Kua = pars(78);
    
    %vars
    Ppl = y(56);
    dV = y(106);
    dVua = y(108);
    
    %Equations   
    Rrs = Raw + Rl + Rcw;
    dVla = dVua + dV;
    Pua = Ppl + dVla * Rrs;        
    dPua_dt = (Pua - y(57))/dt;
    ddVua = -1/R_trachea * (dPua_dt + dVua/Cua);    

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

    dGaw = (Gaw - y(21))/dt;
    internal_variables(4) = Gaw;
    internal_variables_ = internal_variables; 
    y_ = y;
    test_ = test;
end
function [y_, dPmusc, dPpl, dV, ddV, test_] = pulmonary_mechanics(t, y, pars, internal_variables, test)
    %pars 
    Ecw = pars(34);
    El = pars(35);
    kaw1  = pars(296);
    kaw2 = pars(297);
    Rcw = pars(152);
    Rrs = pars(154);
    Pao = pars(124);
    dt = pars(248);
    t0 = pars(317);

    %vars    
    V = y(82);    
    a0 = y(99);
    a1 = y(100);
    a2 = y(101);
    tau = y(129);
    TI = y(79);
    Gaw = internal_variables(4);  

    %Equations
    Ers = Ecw + El;    
    %Pmusc    
    t = t - t0;
    if 0 <= t && t <= TI
        Pmusc = a0 + a1*t + a2*t^2;        

    elseif  TI < t 
        PmuscTI = a0 + a1*TI + a2*TI^2;
        Pmusc = PmuscTI * exp(-(t - TI)/tau);    
    end
    %dV
    dV =  Gaw/Rrs * ((Pmusc - Pao) - Ers * V);   
    %dV =  1/Rrs * ((Pmusc - Pao) - Ers * V);   
    ddV = (dV - y(106))/dt;
    %Pcw
    if dV < 0
        Pcw = Ecw * V - 1;
        Pa_ = Pao;                               
    else
        Pcw = Ecw * V - 1 + Rcw * dV;
        Pa_ = Pao - kaw1 * dV - kaw2 * abs(dV)^2;
    end
    Pa = Pa_ * (Pa_ > 0);
    
    %Ppl
    Ppl = Pcw + Pa - Pmusc; 
 
    dPpl = (Ppl - y(56))/dt;
    dPmusc = (Pmusc - y(55))/dt;   

    y_ = y;
    test_ = test;
    
end


function [y_, internal_variables_, test_] = neuromuscular_drive(t, y, pars, internal_variables, index, tiny_y_keys, test)
    
    %pars
    dt = pars(248);
    t0 = pars(317);

    %vars
    TI = y(79);
    
    %Equations
    t_cycle = t - t0; 
    Nt = eps;
    if t_cycle < TI
        dVE_historic = all_global(1, round(t0/dt) + 1: round(t/dt) + 1);
        Nt = trapz(dVE_historic, 2) * dt;
        
    end 
    internal_variables(99) = Nt;  %dejar derivada en 0 y guardar todas las variables externas en una 
    internal_variables_ = internal_variables; 
    y_ = y;
    test_ = test;
    
end
function [I, internal_variables_, test_] = metabolic_regulation(t, y, pars, internal_variables, test)
    %pars
    MRtCO2_basal = pars(111);
    AT = pars(5);

    %vars
    MRtCO2 = y(26);

    %Equations
    %I = (MRtCO2 + 0.01 - MRtCO2_basal)/(AT - MRtCO2_basal);
    I = (MRtCO2 - MRtCO2_basal)/(AT - MRtCO2_basal);
    I = I * (I > 0) + 0 * (I <= 0); %this is a correction to avoid negative values in the metabolic rate, which is not possible
    %MRtCO2_ = MRtCO2 * (MRtCO2 > 0.3) + 0.3 * (MRtCO2 <= 0.3); %this is a correction to avoid negative values in the metabolic rate of CO2, which is not possible
    %I = (MRtCO2_ - 0.3)/(AT - 0.3);
    
    internal_variables(97) = I;
    internal_variables_ = internal_variables;
    test_ = test;

end

%%Gas exchange
function gas = dissociation(P, pars)

    %pars
    A1 = pars(3); % parameter in O2 dissociation equation
    A2 = pars(4); % parameter in CO3 dissociation equation
    alpha1 = pars(227); % parameter in O3 dissociation equation
    alpha2 = pars(228); % parameter in CO3 dissociation equation
    K1 = pars(59); % parameter in O3 dissociation equation
    K2 = pars(61); % parameter in CO3 dissociation equation
    beta1 = pars(229); % parameter in O3 dissociation equation
    beta2 = pars(230); % parameter in CO3 dissociation equation
    C1 = pars(10);
    C2 = pars(11);
    Z = pars(222);        
    CO2a_ = C2 * Z;
    O2a_ = C1 * Z;
    
    %vars
    PO2 = P(1);
    PCO2 = P(2);

    %Equations
    FCO2 = PCO2 * (1 + beta2 * PO2)/(K2 * (1 + alpha2 * PO2));
    CO2 = CO2a_ * FCO2^(1/A2)/(1 + FCO2^(1/A2));
    
    FO2 = PO2 * (1 + beta1 * PCO2)/(K1 * (1 + alpha1 * PCO2));
    O2 = O2a_ * FO2^(1/A1)/(1 + FO2^(1/A1));
    
    gas = [O2, CO2];    
    
end
function [y_, dPAgas, ddPa, dP_1, dP_2, dP_3, dP_4, dP_5, a, test_] = exchange_mixing(t, y, pars, index, test)

    %pars
    %fO2 = pars(250);
    fCO2 = pars(249);
    Patm = pars(125);
    Pws = pars(130);
    Vdead = pars(203);
    VLO2 = pars(173);
    VLCO2 = pars(172);
    T1 = pars(159);
    T2 = pars(160);
    LCTV = pars(79);
    settling_time = pars(316);
    type_of_input = pars(345);

    %vars    
    P_1O2 = y(35);
    P_1CO2 = y(34);
    P_2O2 = y(37);
    P_2CO2 = y(36);
    P_3O2 = y(39);
    P_3CO2 = y(38);
    P_4O2 = y(41);
    P_4CO2 = y(40);
    P_5O2 = y(43);
    P_5CO2 = y(42);
    PAO2 = y(32);
    PACO2 = y(31);
    dV = y(106);
    V = y(82);
    vO2 = y(132);
    vCO2 = y(131);    
    Qpp = y(68);
    Qla = y(67);
    PaO2 = y(53);
    PaCO2 = y(52);
    dPaO2 = y(105);
    dPaCO2 = y(104);

    %fgas = [fO2, fCO2];
    P_1 = [P_1O2, P_1CO2];
    P_2 = [P_2O2, P_2CO2];
    P_3 = [P_3O2, P_3CO2];
    P_4 = [P_4O2, P_4CO2];
    P_5 = [P_5O2, P_5CO2];

    Ta = 1000 * LCTV/(Qla + 250);    %this has a stop limit because Qla could drop, making the delay truly huge.  
    if t > Ta
        try
            PACO2_delayed = all_global(2, round((t-Ta)/dt) + 1);
            PAO2_delayed = all_global(3, round((t-Ta)/dt) + 1);

        catch
            PACO2_delayed = PACO2;
            PAO2_delayed = PAO2;
        end
    else
        PACO2_delayed = 40;
        PAO2_delayed = 88;
    end
    PAgas = [PAO2, PACO2];
    PAgas_delayed = [PAO2_delayed, PACO2_delayed];
    v = [vO2, vCO2];
    Pa = [PaO2, PaCO2];   
    dPa = [dPaO2, dPaCO2];   
    VL = [VLO2, VLCO2]; 

    %Equations
    
    %Inspired air
    %Computation from data
    fO2p_0 = pars(263);
    fO2p_1 = pars(264);
    fO2p_2 = pars(265);
    fO2p_3 = pars(266);
    fO2p_4 = pars(267);

    if type_of_input > 6
        if t >= settling_time
            tt = t - settling_time;
            
            fO2 = fO2p_0 + fO2p_1 * tt^1 + fO2p_2 * tt^2 + fO2p_3 * tt^3 + fO2p_4 * tt^4;
            if fO2 < 0.123
                fO2 = 0.123;
            end 
            fO2 = 100*fO2; %data is given in decimals
        else
            fO2 = 100 * 0.17;
        end
    else
        fO2 = pars(250);
    end
    
    fgas = [fO2, fCO2];

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
    % Dissociation
    a = dissociation(PAgas, pars); 
    % alveolar exchange 
    if dV >= 0
        dPAgas = (863 * Qpp/1000 * (v - a) + dV * (P_5 - PAgas)) * (eye(2) *  1./(VL + V)); %UNITS CORRECTION between flows and volumes mL to L
        
    else
        
        dPAgas = (863 * Qpp/1000 * (v - a)) * (eye(2) *  1./(VL + V));    %Here we have a 1/1000 to make units corrections (from ml to L)
    end

   
    % mixing
    ddPa = 1/(T1 * T2) * (PAgas_delayed - (T1 + T2)*dPa - Pa);    %momentary removal

    y(103) = a(1);
    y(102) = a(2);
    y_ = y;
    test_ = test;
end


%% 
function [dPvbCO2, dPCSFCO2, internal_variables_ , test_] = brain(t, y, pars, internal_variables, test)
    


    %pars
    KCSFCO2 = pars(66);
    KCCO2 = pars(63);
    dc = pars(234);
    h = pars(289);
    SCO2 = pars(156);
    SbCO2 = pars(157);
    MRbCO2 = pars(109);
    dt = pars(248);

    %vars
    PvbCO2 = y(58);
    PCSFCO2 = y(33);
    PaCO2 = y(52);    
    Qbp = y(66);

    %Equations
    dPvbCO2 = (MRbCO2 * 1 + Qbp * SCO2 * (PaCO2 - PvbCO2) - h)/SbCO2;  %UNIT CORRECTION, now is avoided 
    dPCSFCO2 = (PvbCO2 - PCSFCO2)/KCSFCO2;
    PbCO2 = PvbCO2 + (PCSFCO2 - PvbCO2) * exp(-dc * (Qbp * KCCO2)^0.5);
    %PbCO2 = 40;
    internal_variables(8) = PbCO2;
    internal_variables_ = internal_variables; 
    test_ = test;
end

function [dv, dMRtgas, test_] = tissue(t, y, pars, test)
    
    %Input variables- consumption-production rates
    pars = input_consumption(t, y, pars); %in this inside function we control the input

    %pars
    tauMR = pars(319);    
    Vtissue_CO2 = pars(204);
    Vtissue_O2 = pars(205);
    MRO2 = pars(99);   %This parameters correspond to the input values.  
    MRCO2 = pars(89); %This parameters correspond to the input values.

    %vars
    MRtO2 = y(27);
    MRtCO2 = y(26);    
    vO2 = y(132);
    vCO2 = y(131);
    aO2 = y(103);
    aCO2 = y(102);
    Qpp = y(68);
    Qbp = y(66);    
    
    MRtgas = [MRtO2, MRtCO2];  %l/min
    MRgas = [MRO2, MRCO2];
    a = [aO2, aCO2];
    v = [vO2, vCO2];
    Vtissue = [Vtissue_O2, Vtissue_CO2];

    Qt = Qpp - Qbp;    %ml/s 
    Qt_corrected = Qt/1000;   %UNIT CORRECTION, ml/s --> l/s
    MRtgas_corrected = MRtgas/60;   %UNIT CORRECTION, l/min --> l/s 
    dMRtgas = (MRgas - MRtgas)/tauMR;  %ATENTION to tauMR, if it's too big, MRtgas will not be able to follow MRgas in a short amount of time (so tauMR small for fast simulations (10-30s), and the actual value for large simulations (1 - 5 min))
    dv = (MRtgas_corrected * eye(2).*[-1,1] + Qt_corrected*(a-v)) * eye(2)./Vtissue; %we have to let things in lt/s, because the derivative has to be divided in seconds and the volumes must be in liters, due to Vtissue is in liters.

    function pars_ = input_consumption(t, y, pars)
        %Here you can define the input types        
        type_of_input = pars(345);

        %pars
        MRO2 = pars(99);
        MRCO2 = pars(89);
        MRtO2_basal = pars(112);
        MRtCO2_basal = pars(111);
        
        
        if type_of_input == 0
            pars(99) = MRO2;
            pars(89) = MRCO2;
            
        end

        if type_of_input == 1 %testing_simulation

            %if t < 30
                %barely basal
                MRO2 = 0.02 + MRtO2_basal; %10 * MRtO2_basal;
                MRCO2 = 0.02 + MRtCO2_basal; %10 * MRtCO2_basal;
            %elseif t > 30 && t < 60
            %%    %slightly over doubling it
            %    MRO2 = 0.1 + MRtO2_basal;
            %    MRCO2 = 0.1 + MRtCO2_basal;
            %elseif t > 60
            %
            %    MRO2 = 0.2 + MRtO2_basal;
            %    MRCO2 = 0.2 + MRtCO2_basal;
            %end
        
        end
        
        if type_of_input == 2 %paper_simulation: taken from the figures
            if t > 120 && t <= 360
                MRO2 = 0.7;
                MRCO2 = 0.6;            
            elseif t > 360 && t <= 660
                MRO2 =  0.875;
                MRCO2 = 0.828;
            elseif t > 660 && t <= 960
                MRO2 = 1.183;
                MRCO2 = 1.08;
            elseif t > 960 && t <= 1260
                MRO2 = 1.435;
                MRCO2 = 1.29;
            else
                MRO2 = MRtO2_basal;
                MRCO2 = MRtCO2_basal;
            end
        end

        if type_of_input == 3 %deploy_results
            
            pars(99) = MRO2;
            pars(89) = MRCO2;

        end

        if type_of_input == 4 %tiny paper_simulation: 10 times faster, but with 200 seconds of settling time
        tt = 10 * t - 200;
            if tt > 120 && tt <= 360
                MRO2 = 0.7;
                MRCO2 = 0.6;            
            elseif tt > 360 && tt <= 660
                MRO2 =  0.875;
                MRCO2 = 0.828;
            elseif tt > 660 && tt <= 960
                MRO2 = 1.183;
                MRCO2 = 1.08;
            elseif tt > 960 && tt <= 1260
                MRO2 = 1.435;
                MRCO2 = 1.29;
            else
                MRO2 = MRtO2_basal;
                MRCO2 = MRtCO2_basal;
            end
        end

        if type_of_input == 5 %fitting
            settling_time = pars(316);
            tt = t - settling_time;
            MRO2 = 0.02 + MRtO2_basal; %10 * MRtO2_basal;
            MRCO2 = 0.02 + MRtCO2_basal; %10 * MRtCO2_basal;
            if tt > 120 && t <= 360
                MRO2 = 0.7;
                MRCO2 = 0.6;            
            elseif tt > 360 && t <= 660
                MRO2 =  0.875;
                MRCO2 = 0.828;
            elseif tt > 660 && t <= 960
                MRO2 = 1.183;
                MRCO2 = 1.08;
            elseif tt > 960 && t <= 1260
                MRO2 = 1.435;
                MRCO2 = 1.29;
            else
                MRO2 = MRtO2_basal;
                MRCO2 = MRtCO2_basal;
            end
            

        end

        
        
        
        
        if type_of_input == 6 || type_of_input == 7 % vo2 and vco2 external 
            settling_time = pars(316);
            tt = t - settling_time;
            if t >= settling_time
                MRO2p_0 = pars(100);
                MRO2p_1 = pars(101);
                MRO2p_2 = pars(102);
                MRO2p_3 = pars(103);
                MRO2p_4 = pars(104);
                MRO2p_5 = pars(105);
                MRO2p_6 = pars(106);
                MRO2p_7 = pars(107);
                MRO2p_8 = pars(108);

                MRCO2p_0 = pars(90);
                MRCO2p_1 = pars(91);
                MRCO2p_2 = pars(92);
                MRCO2p_3 = pars(93);
                MRCO2p_4 = pars(94);
                MRCO2p_5 = pars(95);
                MRCO2p_6 = pars(96);
                MRCO2p_7 = pars(97);
                MRCO2p_8 = pars(98);
                

                MRO2 = MRO2p_0 + MRO2p_1*tt + MRO2p_2*tt^2 + MRO2p_3*tt^3 + MRO2p_4*tt^4 + MRO2p_5*tt^5 + MRO2p_6*tt^6 + MRO2p_7*tt^7 + MRO2p_8*tt^8;          
                MRCO2 = MRCO2p_0 + MRCO2p_1*tt + MRCO2p_2*tt^2 + MRCO2p_3*tt^3 + MRCO2p_4*tt^4 + MRCO2p_5*tt^5 + MRCO2p_6*tt^6 + MRCO2p_7*tt^7 + MRCO2p_8*tt^8;

                MRO2 = MRO2 * (MRO2 > MRtO2_basal) + MRtO2_basal * (MRO2 <= MRtO2_basal);
                MRCO2 = MRCO2 * (MRCO2 > MRtCO2_basal) + MRtCO2_basal * (MRCO2 <= MRtCO2_basal); 

                
            end               

        end


        
        pars(99) = MRO2;
        pars(89) = MRCO2;

        pars_ = pars;      
        
    end
    test_ = test;

end
function [dMRR, dM_Rv, dMRv, test_] = metabolism_dynamics(t, y, pars, test)
    %pars
    dt = pars(248);    
    MRbCO2 = pars(109);
    MRbO2 = pars(110);
    MRtCO2_basal = pars(111);
    MRtO2_basal = pars(112);
    tauMRv = pars(320);

    %vars 
    M_Rv = y(29);
    MRtO2 = y(27);
    MRtCO2 = y(26);
    %MRtO2_ = MRtO2 * (MRtO2 > 0.33) + 0.33 * (MRtO2 <= 0.33); %this is a correction to avoid negative values in the metabolic rate of O2, which is not possible
    %MRtCO2_ = MRtCO2 * (MRtCO2 > 0.3) + 0.3 * (MRtCO2 <= 0.3); %this is a correction to avoid negative values in the metabolic rate of CO2, which is not possible
    %Equations
    M_RR = (MRbCO2 + MRbO2 + MRtCO2 + MRtO2)/(MRbCO2 + MRbO2 + MRtCO2_basal + MRtO2_basal); %this variable then goes to dVE computing and it's an expression of the overall metabolic excersise
    %this is a correction to avoid negative values in the metabolic rate of O2, which is not possible
    %M_RR = (MRbCO2 + MRbO2 + MRtCO2_ + MRtO2_)/(MRbCO2 + MRbO2 + 0.3 + 0.33); %this variable then goes to dVE computing and it's an expression of the overall metabolic excersise

    if M_RR >= 1
        MRR = M_RR;
    else
        MRR = 1;
    end

    dMRR = (MRR - y(25))/dt;
    dM_Rv = ((MRR - 1) - M_Rv)/tauMRv;

    if M_Rv >= 0 && MRR > 1
        MRv = M_Rv;
    else
        MRv = 0;
    end

    dMRv = (MRv - y(28))/dt;
    test_ = test;
end

%Ventilatory controller
function  [ddVE, test_] = ventilation_control(t, y, pars, internal_variables, test)
    
    %pars    
    KpCO2 = pars(76);
    KpO2 = pars(77);
    KcMRv = pars(75);
    KcCO2 = pars(74);
    Kbg = pars(73);
    dV_rest = pars(233);
    V0dead = pars(203);
    GVdead = pars(40);
    dt = pars(248);

    %vars
    Tresp = y(81);
    MRv = y(28);
    mean_PaO2 = y(123);
    mean_PaCO2 = y(122);
    mean_PbCO2 = y(124);
    
    dVA_ = dV_rest * (KpCO2 * mean_PaCO2 + KcCO2 * mean_PbCO2 + (KpO2 * (104 - mean_PaO2)^4.9) * (mean_PaO2 < 104) + KcMRv * MRv - Kbg); %This should be alveolar minute volume, and it's part of the minute ventilation we want the lungs to adquire
    dVA = dVA_ * (dVA_ > 0); % as minute ventilation is a positive value, we must take the absolute value, we will never be able to remove air from the lungs, the minimum volume value will always be dead space volume 
    dVA = real(dVA);
    Vd = GVdead * dVA + V0dead; %the first part refers to dead space volume that changes because respiration (expansion of airway channels?)
    dVd = 1/Tresp * Vd;     %the amount of volume that is exchanged during a minute due to dead space, will be the breathing frequency times dead space volume
    
    dVE = dVA + dVd;    %dVE corresponds to minute ventilation (how much volume do I want to exchange over a minute, it's a indicator of respiration flow and its different than instant flow)
    %dVE = dVE * (dVE > 0.1) + 0.1 * (dVE <= 0.1); %this is a correction to avoid negative values in the minute ventilation, which is not possible, and also to avoid very low values that are not physiological (0.11 l/s, which are indeed 7 l/min, reported as common minute ventilation values)
    
    ddVE = (dVE - y(107))/dt;

    %all_global(1, round(t/dt) + 1) = dVE;
    
    
    
    test_ = test;
    %In the simulation is correct kind of the values we got (0.11 l/s, which are indeed 7 l/min, reported as common minute ventilation values)

end



function [mean_dPaO2, mean_dPaCO2, mean_dPbCO2, mean_dP_sa, test_] = mean_values_computer(t, y, pars, internal_variables, test)
    
    
    %vars
    PaO2 = y(53);
    PaCO2 = y(52);
    P_sa = y(45);    
    PbCO2 = internal_variables(8);
    mean_PaO2 = y(123);
    mean_PaCO2 = y(122);
    mean_P_sa = y(121);
    mean_PbCO2 = y(124);

    %means - we want a cut off filter that leaves signals at 0.1Hz range (to take the mean value)
    mean_dPaO2 = (PaO2 - mean_PaO2)/5;
    mean_dPaCO2 = (PaCO2 - mean_PaCO2)/5;
    mean_dPbCO2 = (PbCO2 - mean_PbCO2)/5;
    mean_dP_sa = (P_sa - mean_P_sa)/5;

    

    test_ = test;   
end

%Cardiovascular equations
%'''It's important to notice that different from previous modules, cardiovascular equations share a lot of variables between each other.
%   Even so, due to causal calculations some variables need to be computed outside their native modules, because they might be obtained from other variable which is calculated 3 modules after
% This could be avoided if modules weren't present, but for the sakes of order and easy understanding of the codes we keep it this way '''

function [internal_variables_, test_] = muscle_pump(t, y, pars, internal_variables, test)
    %-------------------------
    %pars
        Aim = pars(6);
        Tc = pars(164);
        Tim = pars(169);
    
    %-----------------------
    %vars
    
    %-----------------------
    %Equations
        alpha_ = mod(t, Tim);
        alpha = alpha_/Tim;
        phi = sin(pi * Tim/Tc * alpha) * (alpha <= Tc/Tim);
        P_im = Aim * phi;
    %---------------------------
    %Saving computations
        internal_variables(55) = P_im;
        internal_variables_ = internal_variables;
    test_ = test;
end

function [internal_variables_, test_] = vena_cava(t, y, pars, internal_variables, test)

    %----------------------
    %pars
    %vena_cava
        D1 = pars(32);
        D2 = pars(33);
        K1_vc = pars(60);
        K2_vc = pars(62);
        K_r_vc = pars(72);
        R_vc_n = pars(150);
        V_unstressed_vc = pars(200);
        V_vc_max = pars(201);
        V_vc_min = pars(202);
    
    %pulmonary
        V_unstressed_ra = pars(193);
        C_ra = pars(24);

    %-----------------------
    %vars
    %for other modules
        % Q_e_v = y('Q_e_v'); %extra-splanchnic venous
        % Q_s_v = y('Q_s_v'); %splanchnic venous
        % Q_b_v = y('Q_b_v'); %brain venous
        % Q_h_v = y('Q_h_v'); %coronary venous
        % Q_rm_v = y('Q_rm_v'); %resting muscular venous
        % Q_am_v = y('Q_am_v'); %active muscular venous    
        %Q_v = [Q_e_v, Q_s_v, Q_b_v, Q_h_v, Q_rm_v, Q_am_v];
    V_total_vc = y(97);
    V_total_ra = y(93);
    Ptor = internal_variables(35);  %this is in mmHg
    % Ptor = y('Ptor');

    %-----------------------
    %Equations
     
    %heart 
        V_ra = (V_total_ra - V_unstressed_ra) * ((V_total_ra - V_unstressed_ra) > 0);
        P_ra = V_ra/C_ra + Ptor;

    %vena_cava     
        V_vc = (V_total_vc - V_unstressed_vc) * ((V_total_vc - V_unstressed_vc) > 0);
        R_vc = K_r_vc * (V_vc_max/V_total_vc)^2 + R_vc_n;
        
        if V_total_vc >= V_unstressed_vc
            %P_vc_ = D1 + K1_vc * (V_vc - V_unstressed_vc);
            P_vc_ = D1 + K1_vc * (V_vc);
        else
            P_vc_ = D2 + K2_vc * exp(V_vc/V_vc_min);
        end
        P_vc = P_vc_ + Ptor;
        Q_ra = (P_vc - P_ra)/R_vc * (P_vc >= P_ra);
        
        %Computed in other modules (systemic_peripheric_and_venous_circulation)
            % Q_vc = sum(Q_v);             
            % dV_total_vc = Q_vc - Q_ra;   

    % y('P_vc') = P_vc;
    % y('V_vc') = V_vc;
    % y('Q_ra') = Q_ra;
    y_ = y;
    
    %-----------------------------
    %Saving computations

        internal_variables(56) = P_vc;
        internal_variables(53) = V_vc;
        internal_variables(45) = Q_ra;
        internal_variables(19) = P_ra;
        internal_variables(46) = V_ra;
        internal_variables_ = internal_variables;
    test_ = test;
end


function [internal_variables_,  dzheta_heart, dV_total_rv, dV_total_ra, dV_total_la, dV_total_lv, dQla, dTheart, test_] = heart(t, y, pars, internal_variables, test)

    %--------------------------------
    %pars
        %heart
            C_la = pars(20);
            K_E_lv = pars(70);
            K_E_rv = pars(71);
            KR_lv = pars(68);
            KR_rv = pars(69);
            ksys = pars(305);
            P0_lv = pars(114);
            P0_rv = pars(115);
            R_la = pars(139);
            R_ra = pars(144);
            Tsys_0 = pars(171);
            V_unstressed_la = pars(188);
            V_unstressed_lv = pars(189);
            V_unstressed_rv = pars(196);
            %for other modules
                %C_ra = pars(24);
                %V_unstressed_ra = pars(193);
        
        %pulmonary
            V_unstressed_pa = pars(190);
            C_pa = pars(21);
            V_unstressed_pv = pars(192);
            C_pv = pars(23);
            R_pv = pars(143);
        
        %general
        dt = pars(248);
    
    %------------------------------------- 
    %vars    
    %for other modules
        %V_total_ra = y(93);
     %this should come from internal variables
    V_total_lv = y(89);
    V_total_la = y(88);
    V_total_rv = y(95);
    V_total_pa = y(90);
    V_total_pv = y(92);
    P_sa = y(45);              
    zheta_heart = y(157);

    
    
    
    %------------------------------------
    %internal_variables
    Theart = internal_variables(132);
    Ptor = internal_variables(35);  
    Q_ra = internal_variables(45);
    P_ra = internal_variables(19);
    E_max_lv = internal_variables(130);
    E_max_rv = internal_variables(131);
    if t == 0
        E_max_lv = y(19);
        E_max_rv = y(20);
    end

    %--------------------------------------
    %Equations
    
    %Cardiac oscilator

        try  %to avoid the non existent first value
            t0_heart = all_global(5, round(t/dt));
            u_t0 = all_global(6, round(t/dt));
        catch
            t0_heart = 0;
            u_t0 = 0;
        end
        
        Tsys = Tsys_0 + ksys * 1/Theart; %systolic time

        %Activation function
            HR = 1/Theart;
            %all_global(7, round(t/dt) + 1) = HR;
            
            HR_integral = dt * sum(all_global(7, round(t0_heart/dt) + 1: round(t/dt) + 1));
            dzheta_heart = HR;
            U = HR_integral + u_t0 - floor(HR_integral + u_t0); %fractional part
            u = zheta_heart - floor(zheta_heart); %fractional part
            
            if U == 0        
                all_global(5, round(t/dt) + 1) = t;        
                all_global(6, round(t/dt) + 1) = u;
                
            else
                all_global(5, round(t/dt) + 1) = t0_heart;         
                all_global(6, round(t/dt) + 1) = u_t0;
                
            end 
                
            phi = sin(pi * Theart/Tsys * u)^2 * (u < Tsys/Theart);
            

    %Left ventricle    
        V_lv = (V_total_lv - V_unstressed_lv) * ((V_total_lv - V_unstressed_lv) > 0);
        P_max_lv = phi * E_max_lv * (V_total_lv - V_unstressed_lv) + (1 - phi) * P0_lv * (exp(K_E_lv * V_total_lv) - 1);  
        
        R_lv = KR_lv * P_max_lv;
        Q_lv = 1/R_lv * (P_max_lv - P_sa) * ((P_max_lv - P_sa) > 0);
        P_lv = P_max_lv - R_lv * Q_lv + Ptor;
        
        

    %Left atrium
        V_la = (V_total_la - V_unstressed_la) * ((V_total_la - V_unstressed_la) > 0);
        P_la = V_la/C_la + Ptor;
        Q_i_lv = 1/R_la * (P_la - P_lv) * ((P_la - P_lv) > 0);

        dV_lv_dt = (Q_i_lv - Q_lv) * (V_total_lv > V_unstressed_lv); %by chain's rule (dVlv_dt = piecewise'(V_total_v) * dV_total_v_dt), and piecewise'(x) is a step function
        Wh_lv = -(P_lv - Ptor) * dV_lv_dt;

    %pulmonary
        V_pa = (V_total_pa - V_unstressed_pa) * ((V_total_pa - V_unstressed_pa) > 0);
        P_pa = V_pa/C_pa + Ptor;

        V_pv = (V_total_pv - V_unstressed_pv) * ((V_total_pv - V_unstressed_pv) > 0);
        P_pv = V_pv/C_pv;

        Q_la = (P_pv + Ptor - P_la)/R_pv;
    
        

    %Right ventricle
        V_rv = (V_total_rv - V_unstressed_rv) * ((V_total_rv - V_unstressed_rv) > 0);
        P_max_rv = phi * E_max_rv * (V_total_rv - V_unstressed_rv) +(1 - phi) * P0_rv * (exp(K_E_rv * V_total_rv) - 1);
        R_rv = KR_rv * P_max_rv;
        Q_rv = 1/(R_rv) * (P_max_rv - P_pa) * ((P_max_rv - P_pa) > 0);
        P_rv = P_max_rv - R_rv * Q_rv + Ptor;
        
    
    %Right atrium
        %Computed in other modules (vena_cava)
            %V_ra = (V_total_ra - V_unstressed_ra) * ((V_total_ra - V_unstressed_ra) > 0);
            %P_ra = V_ra/C_ra + Ptor;
        Q_i_rv = 1/R_ra * (P_ra - P_rv) * ((P_ra - P_rv) > 0);
        dV_rv_dt = (Q_i_rv - Q_rv) * (V_total_rv > V_unstressed_rv); %again, chain's rule (is the same as before).
        Wh_rv = -(P_rv - Ptor) * dV_rv_dt;
        
    
    %Derivatives    
        dV_total_rv = Q_i_rv - Q_rv;
        dV_total_ra = Q_ra - Q_i_rv;
        dV_total_la = Q_la - Q_i_lv;
        dV_total_lv = Q_i_lv - Q_lv;
        
        %for correct internal_variables handling
        dQla = (Q_la - y(67))/dt;
        dTheart = (Theart - y(80))/dt;

    % y('Q_lv') = Q_lv;
    % y('Q_rv') = Q_rv;
    % y('P_la') = P_la;
    y_ = y;

    %------------------------------
    % Saving computations
        internal_variables(36) = Q_lv;
        internal_variables(39) = Q_la;
        internal_variables(38) = Q_rv;
        internal_variables(25) = P_la;
        internal_variables(40) = P_pv;
        internal_variables(41) = P_pa;
        internal_variables(50) = V_pa;
        internal_variables(52) = V_pv;
        internal_variables(48) = V_la;
        internal_variables(49) = V_lv;
        internal_variables(47) = V_rv;
        internal_variables(84) = Wh_lv;
        internal_variables(85) = Wh_rv;
        internal_variables_ = internal_variables;

    % -------------------------------
    %Auxiliar variables
        dvar1 = Wh_lv;
        dvar2 = Wh_rv;
    
    %Inner functions
    function frac = frac(x)    
        if x >= 1
            frac = 0;            
        else
            frac = x;
        end
    end
    test_ = test;
end

function [internal_variables_, dP_sa, dQ_sa, test_] = systemic_arteries(t, y, pars, internal_variables, test)
    
    %---------------------------
    %pars
    C_sa = pars(29);
    L_sa = pars(81);
    R_sa = pars(149);
    V_unstressed_sa = pars(199);
    
    %-----------------------------
    %vars    
    P_sa = y(45);
    Q_sa = y(65); 
    %Ptor = y('Ptor');
    P_sp = y(46);
    

    %------------------------------
    %internal_variables
    Ptor = internal_variables(35);
    Q_lv = internal_variables(36);
    
    %----------------------------
    % Equations
    %systemic_arteries
        V_sa = P_sa * C_sa;   
        V_total_sa = V_sa + V_unstressed_sa; 

        
    %Derivatives
        dvar3 = V_total_sa;
        dP_sa = 1/C_sa * (Q_lv - Q_sa); 
        dQ_sa = 1/L_sa * ((P_sa - Ptor) - R_sa * Q_sa - P_sp); %Careful with L_sa, it's too sensitive
        
        
    %-----------------------------
    %Saving computations
    internal_variables(54) = V_sa;
    internal_variables_ = internal_variables;
    
    y_ = y;
    test_ = test;
end

function [internal_variables_, dQ_pa, dV_total_pp, dV_total_pv, dV_total_pa, dQpp, test_] = pulmonary_circulation(t, y, pars, internal_variables, test)
    %------------------------------------
    %pars
    dt = pars(248);
    C_pp = pars(22);
    %R_pp = pars(142); now is calculated from internal variables (HIPOXIA)
    L_pa = pars(80);
    R_pa = pars(141);
    
    V_unstressed_pp = pars(191);
    % for other modules
        %C_pa = pars(21);
        %C_pv = pars(23);
        %R_pv = pars(143);
        %V_unstressed_pa = pars(190);
        %V_unstressed_pv = pars(192);

    %-------------------------------------
    %vars
    % for other modules
        %Ptor = y('Ptor');
        %V_total_pa = y(90);
        %V_total_pv = y(92);    
    Q_pa = y(62);
    V_total_pp = y(91);

    %------------------------------------
    %internal_variables   
    Q_rv = internal_variables(38);
    Q_la = internal_variables(39);  
    P_pv = internal_variables(40);
    P_pa = internal_variables(41);
    R_pp = internal_variables(139);
    
    %-----------------------------------
    %Equations
    
    %Pulmonary venous 
        %Computed in other modules (heart)
            %V_pv = (V_total_pv - V_unstressed_pv) * (V_total_pv - V_unstressed_pv > 0);
            %P_pv = V_pv/C_pv;
            %Q_la = (P_pv + Ptor - P_la)/R_pv;
    
    %Pulmonary peripheric
        V_pp = (V_total_pp - V_unstressed_pp) * (V_total_pp - V_unstressed_pp > 0);
        P_pp = V_pp/C_pp;
        Q_pp = (P_pp - P_pv)/R_pp;
    
    %Pulmonary arterial
        %Computed in other modules (heart)
            % V_pa = (V_total_pa - V_unstressed_pa) * (V_total_pa - V_unstressed_pa > 0);
            % P_pa = V_pa/C_pa + Ptor;
        dQ_pa = 1/L_pa * (P_pa - R_pa * Q_pa - P_pp);   
    
    %Derivatives
        dV_total_pp = Q_pa - Q_pp;
        dV_total_pv = Q_pp - Q_la;
        dV_total_pa = Q_rv - Q_pa;

        %for correct internal_variables handling
        dQpp = (Q_pp - y(68))/dt;
    
    %---------------------------------------
    %Saving computations
        internal_variables(51) = V_pp;
        internal_variables_ = internal_variables;

        y_ = y;
    test_ = test;

end

function [internal_variables_, dV_total_vc, dV_total_v, dP_sp, dQbp, test_] = systemic_peripheric_and_venous_circulation(t, y, pars, internal_variables, test)
    %---------------------------------
    %pars
    dt = pars(248);
    %compartment order: e,s,b,h,rm,am
        %systemic_peripheric_and_venous_circulation    
            k_r_am = pars(293);
            P0 = pars(113);
            C_e_p = pars(16);
            C_s_p = pars(27);
            C_b_p = pars(14);
            C_h_p = pars(18);
            C_rm_p = pars(25);
            C_am_p = pars(12);
            C_e_v = pars(17);
            C_s_v = pars(28);
            C_b_v = pars(15);
            C_h_v = pars(19);
            C_rm_v = pars(26);
            C_am_v = pars(13);
            R_e_n = pars(135);
            R_s_n = pars(147);
            R_b_n = pars(133);
            R_h_n = pars(137);
            R_rm_n = pars(145);
            R_am_n = pars(131);
            V_unstressed_e_p = pars(184);
            V_unstressed_s_p = pars(197); 
            V_unstressed_b_p = pars(182);
            V_unstressed_h_p = pars(186);
            V_unstressed_rm_p = pars(194);
            V_unstressed_am_p = pars(180);
            %V_unstressed_e_v = pars(185);    %this is actually a var from carfiac control
            %V_unstressed_s_v = pars(198);    %this is actually a var from carfiac control
            V_unstressed_b_v = pars(183);
            V_unstressed_h_v = pars(187);
            %V_unstressed_rm_v = pars(195);  %this is actually a var from carfiac control
            %V_unstressed_am_v = pars(181);  %this is actually a var from carfiac control
            V_tot = pars(175);

        %systemic arteries
            V_unstressed_sa = pars(199);
        %pulmonary
            V_unstressed_pa = pars(190);
            V_unstressed_pp = pars(191);
            V_unstressed_pv = pars(192);
        %heart
            V_unstressed_ra = pars(193);
            V_unstressed_la = pars(188);

    %-----------------------------------------        
    %vars
    Q_sa = y(65);
    I = y(23);
    V_total_e_v = y(86);
    V_total_s_v = y(96);
    V_total_b_v = y(85);
    V_total_h_v = y(87);
    V_total_rm_v = y(94);
    V_total_am_v = y(84);
    P_sp = y(46);
    
    %internal_variables
    Pabd = internal_variables(44);
    Q_ra = internal_variables(45);
    V_ra = internal_variables(46); 
    V_rv = internal_variables(47);
    V_la = internal_variables(48);
    V_lv = internal_variables(49);
    V_pa = internal_variables(50);
    V_pp = internal_variables(51); 
    V_pv = internal_variables(52); 
    V_vc = internal_variables(53); 
    V_sa = internal_variables(54); 
    P_im = internal_variables(55); 
    P_vc = internal_variables(56); 

    R_e_p = internal_variables(137);
    R_s_p = internal_variables(138);
    R_b_p = internal_variables(136);
    R_h_p = internal_variables(134);
    R_rm_p = internal_variables(135);
    R_am_p = internal_variables(133);

    V_unstressed_e_v = internal_variables(126);
    V_unstressed_s_v = internal_variables(127);
    V_unstressed_rm_v = internal_variables(128);
    V_unstressed_am_v = internal_variables(129);
    
    index = {'e','s','b','h','rm','am'};
    
    C_p = [C_e_p, C_s_p, C_b_p, C_h_p, C_rm_p, C_am_p];
    C_v = [C_e_v, C_s_v, C_b_v, C_h_v, C_rm_v, C_am_v];
    R_v_n = [R_e_n, R_s_n, R_b_n, R_h_n, R_rm_n, R_am_n];
    R_p = [R_e_p, R_s_p, R_b_p, R_h_p, R_rm_p, R_am_p];
    V_unstressed_p = [V_unstressed_e_p, V_unstressed_s_p, V_unstressed_b_p, V_unstressed_h_p, V_unstressed_rm_p, V_unstressed_am_p];
    V_unstressed_v = [V_unstressed_e_v, V_unstressed_s_v, V_unstressed_b_v, V_unstressed_h_v, V_unstressed_rm_v, V_unstressed_am_v];
    V_total_v = [V_total_e_v, V_total_s_v, V_total_b_v, V_total_h_v, V_total_rm_v, V_total_am_v];
   
    %Equations 
    %systemic_peripheric_and_venous_circulation
    P = Pabd * strcmp(index, 's');   %vectorial   
    
    V_v = (V_total_v - V_unstressed_v) .* (V_total_v >= V_unstressed_v);  %vectorial
    P_v = C_v.^(-1) .* V_v .* (V_total_v >= V_unstressed_v)  + strcmp(index, 'am') * P_im + (P0 * (1 - (V_total_v./V_unstressed_v).^(-3/2))).* (V_total_v < V_unstressed_v).* strcmp(index, 'am');  %vectorial
    
    V_u = V_unstressed_pp + V_unstressed_sa + V_unstressed_pa + V_unstressed_pv + V_unstressed_ra + V_unstressed_la + sum(V_unstressed_p) + sum(V_unstressed_v);
    
    V_s_v = sum( V_v.*strcmp(index, 's') );
    V_rm_v = sum( V_v.*strcmp(index, 'rm') );
    V_am_v = sum( V_v.*strcmp(index, 'am') );
    V_b_v = sum( V_v.*strcmp(index, 'b') );
    V_h_v = sum( V_v.*strcmp(index, 'h') );

    P_ev = 1/C_e_v * (V_tot - V_sa - V_ra - V_rv - V_la - V_lv - V_pa - V_pp - V_pv - V_s_v - V_rm_v - V_am_v - V_b_v - V_h_v - V_vc - V_u - P_sp * sum(C_p) ); % Check this with profe
    P_v = P_v .* ~strcmp(index, 'e') + P_ev .*  strcmp(index, 'e');  %we have to check this
    
    P_p = P_sp*ones(1,6);
    V_p = C_p .* P_p; %vectorial
    V_total_p = V_p + V_unstressed_p; %vectorial

    R__v = R_v_n .* (P_vc >= P) + R_v_n .* (P_v - P_vc)./(P_v - P) .* (P_vc < P); %vectorial
    R_v = R__v * (I <= 0) + (R__v .* ~strcmp(index, 'am') + k_r_am./V_total_v .* strcmp(index, 'am')) .* (I > 0); %vectorial        
    Q_p = (P_sp - P_v)./R_p; %vectorial.
    Q_v = (P_v - P_vc)./R_v .* (P_v >= P_vc); %vectorial    
    
    dV_total_v = Q_p - Q_v; %vectorial
    dP_sp = 1/(sum(C_p)) * (Q_sa - sum(Q_p)); %vectorial

    %vena_cava
        Q_vc = sum(Q_v);
        dV_total_vc = Q_vc - Q_ra;

    %for correct intenal_variables handling:
    Q_bp = Q_p(3);
    dQbp = (Q_bp - y(66))/dt;

    %['e','s','b','h','rm','am'];
    Q_e = Q_p(1);
    Q_s = Q_p(2);
    Q_h = Q_p(4);
    Q_rm = Q_p(5);
    Q_am = Q_p(6);
    dQ_e = (Q_e - y(60))/dt;
    dQ_s = (Q_s - y(64))/dt;
    dQ_h = (Q_h - y(61))/dt;
    dQ_rm = (Q_rm - y(63))/dt;
    dQ_am = (Q_am - y(59))/dt;
    P_v_e = P_v(1);
    P_v_s = P_v(2);
    P_v_h = P_v(4);
    P_v_rm = P_v(5);
    P_v_am = P_v(6);
    dP_v_e = (P_v_e - y(48))/dt;
    dP_v_s = (P_v_s - y(51))/dt;
    dP_v_h = (P_v_h - y(49))/dt;
    dP_v_rm = (P_v_rm - y(50))/dt;
    dP_v_am = (P_v_am - y(47))/dt;
    dR_e_p = (R_e_p - y(73))/dt;
    dR_s_p = (R_s_p - y(77))/dt;
    dR_b_p = (R_b_p - y(71))/dt;
    dR_h_p = (R_h_p - y(74))/dt;
    dR_rm_p = (R_rm_p - y(75))/dt;
    dR_am_p = (R_am_p - y(69))/dt;

    test('dQ_e') = dQ_e;
    test('dQ_s') = dQ_s;
    test('dQ_h') = dQ_h;
    test('dQ_rm') = dQ_rm;
    test('dQ_am') = dQ_am;    
    test('dP_v_e') = dP_v_e;
    test('dP_v_s') = dP_v_s;
    test('dP_v_h') = dP_v_h;
    test('dP_v_rm') = dP_v_rm;
    test('dP_v_am') = dP_v_am;
    test('dR_e_p') = dR_e_p;
    test('dR_s_p') = dR_s_p;
    test('dR_b_p') = dR_b_p;
    test('dR_h_p') = dR_h_p;
    test('dR_rm_p') = dR_rm_p;
    test('dR_am_p') = dR_am_p;

    %Saving computations    
    internal_variables(91) = Q_p(6);
    internal_variables(86) = Q_p(4);
    internal_variables(87) = Q_p(5);
    internal_variables(74) = Q_p(3);
    internal_variables(78) = Q_p(1);
    internal_variables(79) = Q_p(2);
    internal_variables_ = internal_variables;
    test_ = test;
end

function [dP_mean, internal_variables_] = afferent_barorreflex(t,y,pars, internal_variables, dP_sa)

    %pars
    f_ab_min = pars(252);
    f_ab_max = pars(251);
    kab = pars(294);
    P_n = pars(119);
    tau_p = pars(341);
    tau_z = pars(343);
    
    %vars
    P_mean = y(44);
    P_sa = y(45);


    %Equations

    %Derivatives
    dP_mean = 1/tau_p * (P_sa + tau_z * dP_sa - P_mean);    
    fab = 1/(1 + exp((P_mean - P_n)/kab)) * (f_ab_min + f_ab_max * exp((P_mean - P_n)/kab));

    %Saving computations
    internal_variables(98) = fab;
    internal_variables_ = internal_variables;
    
end

function  [df_ac] = afferent_chemoreceptor(t,y,internal_variables)

    %pars
    f_ac_CO2_n = pars(253);
    f_ac_max = pars(254);
    f_ac_min = pars(255);
    kac = pars(295);
    KH = pars(67);
    PaO2_ac_n = pars(121);
    PaCO2_n = pars(120);
    tau_ac = pars(336);


    %vars
    f_ac = y(111);
    PaO2 = y(53);
    PaCO2 = y(52);

    %Equations
    if PaO2 > 80
        K = KH;
    elseif PaO2 > 40 && PaO2 < 80
        K = KH - 1.2 * ((PaO2 - 80)/30);
    else
        K = KH - 1.6;
    end

    %Derivatives
    phi_ac = 1/(1 + exp((PaO2 - PaO2_ac_n)/kac)) * (f_ac_max + f_ac_min * exp((PaO2 - PaO2_ac_n)/kac)) * (K * log(PaCO2/PaCO2_n) + f_ac_CO2_n);
    df_ac = 1/tau_ac * (phi_ac - f_ac);

end

function [df_ap] = afferent_pulmonary_stretch(t,y,pars, internal_variables)

    %pars
    G_ap = pars(51);
    tau_ap = pars(337);

    %vars
    f_ap = y(112);
    V = y(82);

    %Equations
    phi_ap = G_ap * V;
    %Derivatives
    df_ap = 1/tau_ap * (phi_ap - f_ap);
end

function [dxO2_b, dxCO2_b, dR_bp, internal_variables_] = cerebral_blood_flow(t, y, pars, internal_variables)

    %pars
    A = pars(1);
    B = pars(7);
    C = pars(9);
    D = pars(31);
    dt = pars(248);
    vO2_b_n = pars(352);
    gO2_b = pars(269);
    MO2_bp = pars(83);
    R_bmp = pars(134);
    tau_CO2 = pars(321);
    tau_O2 = pars(325);
    PaCO2_n = pars(120);

    %vars
    xO2_b = y(149);
    xCO2_b = y(142);
    aO2 = y(103);
    PaCO2 = y(52);

    %internal variables
    %Q_b_p = internal_variables(74);
    Q_b_p = y(66);

    %Equations
    G_bp = 1/R_bmp * (1 + xO2_b + xCO2_b);    
    R_bp = 1/G_bp; 
    vO2_b = aO2 - MO2_bp/Q_b_p;
    phi_b = (A + B/(1 + C * exp(D * log(PaCO2))))/(A + B/(1 + C * exp(D * log(PaCO2_n) ))) - 1;
    %CC = PaCO2_n;
    %DD = 13.35; 
    %phi_b = ((A + B/(1 + (CC/PaCO2_n)^DD))/(A + B/(1 + (CC/PaCO2_n)^DD)) - 2)/2;
    %Derivatives
    dxO2_b = 1/tau_O2 * (-xO2_b - gO2_b * (vO2_b - vO2_b_n));
    dxCO2_b = 1/tau_CO2 * (-xCO2_b - phi_b);

    internal_variables(75) = R_bp;
    
    dR_bp =  (R_bp - y(72))/dt;
    internal_variables_ = internal_variables;
    
    
    

end

function [dxO2_e, dxO2_s, dxO2_p, internal_variables_] = hipoxia_local_regulation(t, y, pars, internal_variables) 

    % local control splacnic + extrasplacnic

    %pars    
    dt = pars(248);
    vO2_e_n = pars(353);
    vO2_s_n = pars(356);
    aO2_n = pars(226);

    gO2_e = pars(270);
    gO2_s = pars(274);
    gO2_p = pars(272);

    MO2_e = pars(84);
    MO2_s = pars(88);
    MO2_p = pars(86);

    R_e_p_n = internal_variables(122);
    R_s_p_n = internal_variables(123);
    R_p_p_n = pars(140);
    
    tau_CO2 = pars(321);
    tau_O2 = pars(325);

    %vars
    xO2_e = y(150);
    xCO2_e = y(143);
    xO2_s = y(154);
    xCO2_s = y(147);
    xO2_p = y(152);
    xCO2_p = y(145);
    aO2 = y(103);
    

    %internal variables
    Q_e_p = internal_variables(78);
    Q_s_p = internal_variables(79);
    %we are going asume that levels of oxygen that affect pulmonary resistances are not affected by pulmonary work   
      
    R_ep = R_e_p_n * 1/(1 + xO2_e);
    R_sp = R_s_p_n * 1/(1 + xO2_s);
    R_pp = R_p_p_n * (1 + xO2_p);
    
    vO2_e = aO2 - MO2_e/Q_e_p;
    vO2_s = aO2 - MO2_s/Q_s_p;
    %aO2 = dissociation(PAO2, pars);

    dxO2_e = 1/tau_O2 * (-xO2_e - gO2_e * (vO2_e - vO2_e_n) );
    dxO2_s = 1/tau_O2 * (-xO2_s - gO2_s * (vO2_s - vO2_s_n) );
    dxO2_p = 1/tau_O2 * (-xO2_p - gO2_p * (aO2 - aO2_n) );

    internal_variables(137) = R_ep;
    internal_variables(138) = R_sp;
    internal_variables(139) = R_pp;

    internal_variables_ = internal_variables;



      


end


function  [dxO2, dxCO2, dWh, internal_variables_] = coronary_and_resting_muscle_blood_flow(t, y, pars, internal_variables)

    %pars
    vO2_h_n = pars(354);   
    vO2_rm_n = pars(355);   
    gO2_h = pars(271);  
    gO2_rm = pars(273);  
    KCO2_h = pars(64);   
    KCO2_rm = pars(65);   
    MO2_h_p_n = pars(85);   
    MO2_rm_p = pars(87);   
    R_h_p_n = pars(138);   
    %R_rm_p_n = pars('R_rm_p_n');   
    R_rm_p_n = internal_variables(124);   
    tau_w = pars(342);   
    tau_O2 = pars(325);
    tau_CO2 = pars(321);
    Whn = pars(213); 
    PaCO2_n = pars(120); 
    

    vO2_n = [vO2_h_n, vO2_rm_n];
    gO2 = [gO2_h, gO2_rm];
    KCO2 = [KCO2_h, KCO2_rm];
    Rp_n = [R_h_p_n, R_rm_p_n];
    MO2_p = [MO2_h_p_n, MO2_rm_p]; %this changes afterwards in the code
    
    %vars
    xO2_h = y(151);
    xO2_rm = y(153);
    xCO2_h = y(144);
    xCO2_rm = y(146);
    Wh = y(98);
    PaCO2 = y(52);
    aO2 = y(103);

    xO2 = [xO2_h, xO2_rm];
    xCO2 = [xCO2_h, xCO2_rm];

    %internal variables
    Wh_lv = internal_variables(84);
    Wh_rv = internal_variables(85);
    Q_h_p = internal_variables(86);
    Q_rm_p = internal_variables(87);
    Q_p = [Q_h_p, Q_rm_p];

    %Equations
    Rp = Rp_n .* (1 + xCO2) ./ (1 + xO2); %vectorial
    MO2_h_p = MO2_h_p_n * Wh/Whn;
    
    MO2_p(1) = MO2_h_p;
    vO2 = aO2 - MO2_p ./ Q_p;   %vectorial   '''creo que está correcto porque por balance de flujos se tiene que expresar así: vo2 * Q - aO2 * Q = MO2p --> sale/t - entra/t = consumo/t. Si fuera la derivada del flujo además aparecería un segundo extra en el denominador
    phi = (1 - exp( (KCO2.^ -1) * (PaCO2 - PaCO2_n))) ./ (1 + exp( (KCO2.^ -1) * (PaCO2 - PaCO2_n)));  %vectorial
    wh = Wh_lv + Wh_rv;

    %Derivatives
    dxO2 = 1/tau_O2 * ( -xO2 - gO2 .* (vO2 - vO2_n)); %vectorial
    dxCO2 = 1/tau_CO2 * (-xCO2 + phi);   %vectorial
    dWh = 1/tau_w * (wh - Wh);
    
    internal_variables(134) = Rp(1);
    internal_variables(135) = Rp(2);
    internal_variables_ = internal_variables;


    
    
end

function [dxO2_am, dx_met, dx_M, dphi_met, internal_variables_] = active_muscle_blood_flow(t,y,pars, internal_variables, index_fun)
    
    %pars
    vO2_am_n = pars(351);
    delay_met = pars(247);
    gO2_am = pars(268);
    g_M = pars(275);
    I0_met = pars(54);
    kmet = pars(304);
    MO2_am_p_n = pars(82);
    phi_max = pars(314);
    phi_min = pars(315);
    tau_M = pars(324);
    tau_O2 = pars(325);
    tau_CO2 = pars(321);
    tau_met = pars(340);

    %vars
    aO2 = y(103);
    xO2_am  = y(148);
    x_M = y(155);
    x_met = y(156);
    %phi_met_delayed = all_global(:,round(t - delay_met/dt) + 1:end);  %we have to put the dictionary here
    phi_met_delayed = get_delayed_value(tiny_y_keys, t, delay_met, dt, all_global, index_fun, 0, 8);
    
    %internal variables
    I = internal_variables(97);
    Q_am_p = internal_variables(91);
    %Equations
    R_amp_n = internal_variables(125);
    R_am_p = R_amp_n/(1 + xO2_am + x_met);
    MO2_am_p = MO2_am_p_n * (1 + x_M);
    vO2_am = aO2 - MO2_am_p/Q_am_p;
    vO2_am = vO2_am * (vO2_am >= 0);
    
    phi_met = (phi_min + phi_max * exp((I - I0_met)/kmet))/(1 + exp((I - I0_met)/kmet));

    %Derivatives
    dxO2_am = 1/tau_O2 * (-xO2_am - gO2_am * (vO2_am - vO2_am_n));
    dx_M = 1/tau_M * (-x_M + g_M * I);
    dx_met = 1/tau_met * (-x_met + phi_met_delayed);
    internal_variables(133) = R_am_p;
    internal_variables_  = internal_variables;
    dphi_met = phi_met - y(125);

    function delayed_value = get_delayed_value(tiny_y_keys, t, delay, dt, all_global, index_fun, fj_init, fs)
        if delay > t
            % If delay is greater than current time, assign initial value
            delayed_value = fj_init;
        else
            % Otherwise, calculate the delayed value using the original code
            delayed_value = all_global(fs, round((t - delay)/dt) + 1);
        end
    


end

    

end

function [dDThetaO2_s, dDThetaCO2_s, internal_variables_] = cns_ischemic_response(t, y, pars, internal_variables)
     %pars
     gcc_h_s  = pars(285);
     gcc_p_s = pars(286);
     gcc_v_s = pars(287);
     k_isc_h_s = pars(290);
     k_isc_p_s = pars(291);
     k_isc_v_s = pars(292);
     PO2_ref_h_s = pars(116);
     PO2_ref_p_s = pars(117);
     PO2_ref_v_s = pars(118);
     tau_cc = pars(338);
     tau_isc = pars(339);
     Theta_h_s_n = pars(165);
     Theta_p_s_n = pars(166);
     Theta_v_s_n = pars(168);
     x_h_s = pars(357);
     x_p_s = pars(358);
     x_v_s = pars(359);
     PaCO2_n = pars(120);

     gcc_s = [gcc_h_s, gcc_p_s, gcc_v_s];
     k_isc_s = [k_isc_h_s, k_isc_p_s, k_isc_v_s];
     PO2_ref_s = [PO2_ref_h_s, PO2_ref_p_s, PO2_ref_v_s];
     Theta_s_n = [Theta_h_s_n, Theta_p_s_n, Theta_v_s_n];
     xs = [x_h_s, x_p_s, x_v_s];

     %vars
     PaO2 = y(53); 
     PaCO2 = y(52);
     DThetaO2_h_s = y(4);
     DThetaO2_p_s = y(5);
     DThetaO2_v_s = y(6);
     DThetaCO2_h_s = y(1);
     DThetaCO2_p_s = y(2);
     DThetaCO2_v_s = y(3);
     DThetaCO2_s = [DThetaCO2_h_s, DThetaCO2_p_s, DThetaCO2_v_s];
     DThetaO2_s = [DThetaO2_h_s, DThetaO2_p_s, DThetaO2_v_s];


     %internal_variables
     %Equations
     ws = xs ./ (1 + exp((PaO2 - PO2_ref_s) ./ k_isc_s)); %vectorial
     Theta_s = Theta_s_n - DThetaO2_s - DThetaCO2_s; %Theta_s = Theta_s_n - DThetaO2_s - DThetaCO2_s; %vectorial

     %Derivatives
     dDThetaCO2_s = 1/tau_cc * (-DThetaCO2_s + gcc_s * (PaCO2 - PaCO2_n)); %vectorial
     dDThetaO2_s = 1/tau_isc * (-DThetaO2_s + ws); %vectorial

     internal_variables(100) = Theta_s(1);
     internal_variables(101) = Theta_s(2);
     internal_variables(102) = Theta_s(3);

     internal_variables_ = internal_variables;

end

function [internal_variables_] = efferent_pathways(t, y, pars, internal_variables)
    
    %pars
    fab_0 = pars(256);
    fes_0 = pars(257);
    fes_inf = pars(258);
    fes_max = pars(259);
    fev_0 = pars(261);
    fev_inf = pars(262);
    kes = pars(302);
    kev = pars(303);
    I_0_h_s = pars(55);
    I_0_p_s = pars(56);
    I_0_v_s = pars(58);
    I_0_v = pars(57);
    kcc_h_s = pars(298);
    kcc_p_s = pars(299);
    kcc_v_s = pars(301);
    kcc_v = pars(300);
    gamma_h_s_max = pars(277);
    gamma_p_s_max = pars(279);
    gamma_v_s_max = pars(283);
    gamma_v_max = pars(281);
    gamma_h_s_min = pars(278);
    gamma_p_s_min = pars(280);
    gamma_v_s_min = pars(284);
    gamma_v_min = pars(282);
    Theta_v = pars(167);
    Wb_h_s = pars(206);
    Wb_p_s = pars(207);
    Wb_v_s = pars(208);
    Wc_h_s = pars(209);
    Wc_p_s = pars(210);
    Wc_v_s = pars(212);
    Wc_v = 0.9*pars(211);
    Wp_h_s = pars(214);
    Wp_p_s = pars(215);
    Wp_v_s = pars(217);
    Wp_v = pars(216);
    Wt_h_s = pars(218);
    Wt_p_s = pars(219);
    Wt_v_s = pars(221);
    Wt_v = pars(220);
    
    I_0 = [I_0_h_s, I_0_p_s, I_0_v_s];
    kcc = [kcc_h_s, kcc_p_s, kcc_v_s];
    gamma_max = [gamma_h_s_max, gamma_p_s_max, gamma_v_s_max];
    gamma_min = [gamma_h_s_min, gamma_p_s_min, gamma_v_s_min];
    Wb_s = [Wb_h_s, Wb_p_s, Wb_v_s];
    Wc_s = [Wc_h_s, Wc_p_s, Wc_v_s];
    Wp_s = [Wp_h_s, Wp_p_s, Wp_v_s];
    Wt_s = [Wt_h_s, Wt_p_s, Wt_v_s]; 

    %vars
    fac = y(111);
    fap = y(112);

    %internal variables
    I = internal_variables(97);
    fab = internal_variables(98);
    Nt = internal_variables(99);    
    
    Theta_h_s = internal_variables(100);
    Theta_p_s = internal_variables(101);
    Theta_v_s = internal_variables(102);

    Theta_s = [Theta_h_s, Theta_p_s, Theta_v_s];

    %Equations
    gamma_s = (gamma_min + gamma_max.* exp((I - I_0)./(kcc))) ./ (1 + exp((I - I_0)./(kcc))); %vectorial
    gamma_v = (gamma_v_min + gamma_v_max * exp((I - I_0_v)./(kcc_v))) ./ (1 + exp((I - I_0_v)./(kcc_v))); 
    fas = Wt_s * Nt + Wb_s * fab + Wc_s * fac + Wp_s * fap - Theta_s; %vectorial
    fs = fes_inf + (fes_0 + fes_inf) * exp(kes * fas) + gamma_s; %vectorial
    fs = fes_max.* (fs >= fes_max) + fs .* (fs < fes_max); %vectorial
    fv = (fev_0 + fev_inf * exp((fab - fab_0)/(kev)))/(1 + exp((fab - fab_0)/(kev))) - Wt_v * Nt + Wc_v * fac + Wp_v * fap - Theta_v + gamma_v;
    
    %internal variables
    internal_variables(120) = fs(1);
    internal_variables(108) = fs(2);
    internal_variables(109) = fs(3);
    internal_variables(121) = fv;

    internal_variables_ = internal_variables;
    
    %Derivatives
end

function [dDTheta, dfh_s, dfp_s, dfv_s, internal_variables_] = reflex_control_R_Vu_E(t, y, pars, internal_variables, index_fun)

    %pars
    delay_Emax_lv = pars(235);
    delay_Emax_rv = pars(236);
    delay_R_am_p = pars(237);
    delay_R_e_p = pars(238);
    delay_R_rm_p = pars(239);
    delay_R_s_p = pars(240);
    delay_V_u_am_v = pars(243);
    delay_V_u_e_v = pars(244);
    delay_V_u_rm_v = pars(245);
    delay_V_u_s_v = pars(246);
    Emax_lv_0 = pars(36);
    Emax_rv_0 = pars(37);
    R_am_p_0 = pars(132);
    R_e_p_0 = pars(136);
    R_rm_p_0 = pars(146);
    R_s_p_0 = pars(148);
    V_u_am_v_0 = pars(176);
    V_u_rm_v_0 = pars(178);
    V_u_e_v_0 = pars(177);
    V_u_s_v_0 = pars(179);
    G_Emax_lv = pars(41);
    G_Emax_rv = pars(42);
    G_R_am_p = pars(43);
    G_R_e_p = pars(44);
    G_R_rm_p = pars(45);
    G_R_s_p = pars(46);
    G_V_u_am_v = pars(47);
    G_V_u_e_v = pars(48);
    G_V_u_rm_v = pars(49);
    G_V_u_s_v = pars(50);
    tau_Emax_lv = pars(322);
    tau_Emax_rv = pars(323);
    tau_R_am_p = pars(326);
    tau_R_e_p = pars(327);
    tau_R_rm_p = pars(328);
    tau_R_s_p = pars(329);
    tau_V_u_am_v = pars(332);
    tau_V_u_e_v = pars(333);
    tau_V_u_rm_v = pars(334);
    tau_V_u_s_v = pars(335);
    fes_min = pars(260);
    
    tauTheta = [tau_R_e_p, tau_R_s_p, tau_R_rm_p, tau_R_am_p, tau_V_u_e_v, tau_V_u_s_v, tau_V_u_rm_v, tau_V_u_am_v, tau_Emax_lv, tau_Emax_rv];
    Theta0 = [R_e_p_0, R_s_p_0, R_rm_p_0, R_am_p_0, V_u_e_v_0, V_u_s_v_0, V_u_rm_v_0, V_u_am_v_0, Emax_lv_0, Emax_rv_0];
    GTheta = [G_R_e_p, G_R_s_p, G_R_rm_p, G_R_am_p, G_V_u_e_v, G_V_u_s_v, G_V_u_rm_v, G_V_u_am_v, G_Emax_lv, G_Emax_rv];

    %vars
    DTheta_R_e_p = y(10);
    DTheta_R_s_p = y(12);
    DTheta_R_rm_p_n = y(11);
    DTheta_R_am_p_n = y(9);
    DTheta_V_unstressed_e_v = y(14);
    DTheta_V_unstressed_s_v = y(16);
    DTheta_V_unstressed_rm_v = y(15);
    DTheta_V_unstressed_am_v = y(13);
    DTheta_Emax_lv = y(7);
    DTheta_Emax_rv = y(8);
    

    DTheta = [DTheta_R_e_p, DTheta_R_s_p, DTheta_R_rm_p_n, DTheta_R_am_p_n, DTheta_V_unstressed_e_v, DTheta_V_unstressed_s_v, DTheta_V_unstressed_rm_v, DTheta_V_unstressed_am_v, DTheta_Emax_lv, DTheta_Emax_rv];
    %internal variables
    f_h_s_actual = internal_variables(120);  
    f_p_s_actual = internal_variables(108);
    f_v_s_actual = internal_variables(109);

    f_h_s_delayed_Emaxlv = get_delayed_value(tiny_y_keys, t, delay_Emax_lv, dt, all_global, index_fun, f_h_s_actual, 9);
    f_h_s_delayed_Emaxrv = get_delayed_value(tiny_y_keys, t, delay_Emax_rv, dt, all_global, index_fun, f_h_s_actual, 9);
    f_p_s_delayed_Rep = get_delayed_value(tiny_y_keys, t, delay_R_e_p, dt, all_global, index_fun, f_p_s_actual, 10);
    f_p_s_delayed_Rsp = get_delayed_value(tiny_y_keys, t, delay_R_s_p, dt, all_global, index_fun, f_p_s_actual, 10);
    f_p_s_delayed_Rrmpn = get_delayed_value(tiny_y_keys, t, delay_R_rm_p, dt, all_global, index_fun, f_p_s_actual, 10);
    f_p_s_delayed_Rampn = get_delayed_value(tiny_y_keys, t, delay_R_am_p, dt, all_global, index_fun, f_p_s_actual, 10);
    f_v_s_delayed_Vuev = get_delayed_value(tiny_y_keys, t, delay_V_u_e_v, dt, all_global, index_fun, f_v_s_actual, 11);
    f_v_s_delayed_Vusv = get_delayed_value(tiny_y_keys, t, delay_V_u_s_v, dt, all_global, index_fun, f_v_s_actual, 11);
    f_v_s_delayed_Vurmv = get_delayed_value(tiny_y_keys, t, delay_V_u_rm_v, dt, all_global, index_fun, f_v_s_actual, 11);
    f_v_s_delayed_Vuamv = get_delayed_value(tiny_y_keys, t, delay_V_u_am_v, dt, all_global, index_fun, f_v_s_actual, 11);

    dfh_s = (f_h_s_actual - y(115))/dt;   
    dfp_s = (f_p_s_actual - y(116))/dt;    
    dfv_s = (f_v_s_actual - y(118))/dt;
    
    fs_delayed = [f_p_s_delayed_Rep, f_p_s_delayed_Rsp, f_p_s_delayed_Rrmpn, f_p_s_delayed_Rampn, f_v_s_delayed_Vuev, f_v_s_delayed_Vusv, f_v_s_delayed_Vurmv, f_v_s_delayed_Vuamv, f_h_s_delayed_Emaxlv, f_h_s_delayed_Emaxrv];
    fs_actual = [f_p_s_actual, f_p_s_actual, f_p_s_actual, f_p_s_actual, f_v_s_actual, f_v_s_actual, f_v_s_actual, f_v_s_actual, f_h_s_actual, f_h_s_actual];
    %testing
    %fs_delayed = fs_actual;

    %Equations
    
    sigmaTheta = (fs_actual >= fes_min) .* (GTheta.* log(fs_delayed - fes_min + 1)); 
    Theta = DTheta + Theta0;  %vectorial   (this are the actual effectors we have to change for the next iteration)
    
    %Derivatives 
    dDTheta = (tauTheta.^(-1)).* (-DTheta + sigmaTheta);  %vectorial

    internal_variables(137) = Theta(1);
    internal_variables(138) = Theta(2);
    internal_variables(124) = Theta(3);
    internal_variables(125) = Theta(4);
    internal_variables(126) = Theta(5);
    internal_variables(127) = Theta(6);
    internal_variables(128) = Theta(7);
    internal_variables(129) = Theta(8);
    internal_variables(130) = Theta(9);
    internal_variables(131) = Theta(10);

    internal_variables_ = internal_variables;


    function delayed_value = get_delayed_value(init_keys, t, delay, dt, all_global, index_fun, fj_init, fs)
        if delay > t
            % If delay is greater than current time, assign initial value
            delayed_value = fj_init;
        else
            % Otherwise, calculate the delayed value using the original code
            delayed_value = all_global(fs, round((t - delay)/dt) + 1);
        end
    end
end

function [dDTsym, dDTvagal, dfv, internal_variables_] = reflex_control_HR(t, y, pars, internal_variables, index_fun)
    %pars
    delay_Tsym = pars(241);
    delay_Tvagal = pars(242);
    GTsym = pars(38);
    GTvagal = pars(39);
    tau_Tsym = pars(330);
    tau_Tvagal = pars(331);
    fes_min = pars(260);
    %vars
    DTsym = y(17);
    DTvagal = y(18);
    %internal variables
    f_h_s_actual = internal_variables(120);
    f_v_actual = internal_variables(121);
    %if t < 10
    %    f_h_s_actual = 3.85;
    %    f_v_actual = 4.27;
    %end
    f_h_s_delayed = get_delayed_value(tiny_y_keys, t, delay_Tsym, dt, all_global, index_fun, f_h_s_actual, 9); %replace f_s_h_actual
    fv_delayed = get_delayed_value(tiny_y_keys, t, delay_Tvagal, dt, all_global, index_fun, f_v_actual, 12);    %replace f_v_actual

    dfv = (f_v_actual - y(117))/dt;

    %f_h_s_delayed = f_h_s_actual;
    %fv_delayed = f_v_actual;

    %Equations
    sigma_Tsym = GTsym * log(f_h_s_delayed - fes_min + 1) * (f_h_s_actual >= fes_min);
    sigma_Tvagal  = GTvagal * fv_delayed;
    %T = DTvagal + DTsym + T0;  %this is computed in the first module of the simulation
    %HR = 1/T;
    
    %Derivatives
    dDTsym = 1/tau_Tsym * (-DTsym + sigma_Tsym);
    dDTvagal = 1/tau_Tvagal * (-DTvagal + sigma_Tvagal);

    
    internal_variables_ = internal_variables;

    function delayed_value = get_delayed_value(init_keys, t, delay, dt, all_global, index_fun, fj_init, fs)
        if delay > t
            % If delay is greater than current time, assign initial value
            delayed_value = fj_init;
        else
            % Otherwise, calculate the delayed value using the original code
            delayed_value = all_global(fs, round((t - delay)/dt) + 1);
        end
    end

end

function [internal_variables_] = adding_base_values_in_control_variables(t, y, pars, internal_variables)


    Emax_lv_0 = pars(36);
    Emax_rv_0 = pars(37);
    R_am_p_0 = pars(132);
    R_e_p_0 = pars(136);
    R_rm_p_0 = pars(146);
    R_s_p_0 = pars(148);
    V_u_am_v_0 = pars(176);
    V_u_rm_v_0 = pars(178);
    V_u_e_v_0 = pars(177);
    V_u_s_v_0 = pars(179);
    T0 = pars(158);
    R_bmp = pars(134);
    R_h_p_n = pars(138);  
    
    %HIPOXIA    
    R_p_p_n =  pars(140); % pars(142);

    Theta0 = [R_e_p_0, R_s_p_0, R_rm_p_0, R_am_p_0, V_u_e_v_0, V_u_s_v_0, V_u_rm_v_0, V_u_am_v_0, Emax_lv_0, Emax_rv_0];
    
    %vars
    DTheta_R_e_p = y(10);
    DTheta_R_s_p = y(12);
    DTheta_R_rm_p_n = y(11);
    DTheta_R_am_p_n = y(9);
    DTheta_V_unstressed_e_v = y(14);
    DTheta_V_unstressed_s_v = y(16);
    DTheta_V_unstressed_rm_v = y(15);
    DTheta_V_unstressed_am_v = y(13);
    DTheta_Emax_lv = y(7);
    DTheta_Emax_rv = y(8);
    DTsym = y(17);
    DTvagal = y(18);
    xO2_am  = y(148);
    x_met = y(156);
    xO2_h = y(151);
    xO2_rm = y(153);
    xCO2_h = y(144);
    xCO2_rm = y(146);
    xO2_b = y(149);
    xCO2_b = y(142);
    xO2_e = y(150);
    xCO2_e = y(143);
    xO2_s = y(154);
    xCO2_s = y(147);
    xO2_p = y(152);
    xCO2_p = y(145);
    aO2 = y(103);

    xO2 = [xO2_h, xO2_rm];
    xCO2 = [xCO2_h, xCO2_rm];

    DTheta = [DTheta_R_e_p, DTheta_R_s_p, DTheta_R_rm_p_n, DTheta_R_am_p_n, DTheta_V_unstressed_e_v, DTheta_V_unstressed_s_v, DTheta_V_unstressed_rm_v, DTheta_V_unstressed_am_v, DTheta_Emax_lv, DTheta_Emax_rv];
    
    
    Theta = DTheta + Theta0;
    R_am_p_n = Theta(4);
    R_am_p = R_am_p_n/(1 + xO2_am + x_met);    
    Theart = DTvagal + DTsym + T0;    
    R_rm_p_n = Theta(3);
    Rp_n = [R_h_p_n, R_rm_p_n];   %%%THIS EQUATION MUST CHANGE
    Rp = Rp_n .* (1 + xCO2) ./ (1 + xO2);   
    G_bp = 1/R_bmp * (1 + xO2_b + xCO2_b);
    R_bp = 1/G_bp; 

    %HIPOXIA
    R_e_p_n = Theta(1);
    R_s_p_n = Theta(2);
    R_ep = R_e_p_n * 1/(1 + xO2_e);
    R_sp = R_s_p_n * 1/(1 + xO2_s);
    R_pp_ = R_p_p_n * (1 + xO2_p);  %This name is to avoid repetition




    internal_variables(122) = R_e_p_n;
    internal_variables(123) = R_s_p_n;
    internal_variables(124) = R_rm_p_n;
    internal_variables(125) = R_am_p_n;
    internal_variables(126) = Theta(5);
    internal_variables(127) = Theta(6);
    internal_variables(128) = Theta(7);
    internal_variables(129) = Theta(8);
    internal_variables(130) = Theta(9);
    internal_variables(131) = Theta(10);
    internal_variables(132) = Theart;    
    internal_variables(133) = R_am_p;
    internal_variables(134) = Rp(1);
    internal_variables(135) = Rp(2);
    internal_variables(136) = R_bp;
    %HIPOXIA
    internal_variables(137) = R_ep;
    internal_variables(138) = R_sp;
    internal_variables(139) = R_pp_;

    internal_variables_ = internal_variables;
end


function [all_global_, tiny_y_keys] = saving_in_globals(all_global, y)

    
    tiny_y_keys = ['dVE', 'PACO2', 'PAO2', 'Pmusc', 't0_heart', 'u_t0', 'HR', 'phi_met', 'fh_s', 'fp_s', 'fv_s', 'fv'];


    all_global_ = all_global;
   
end




end

% References:
% [1] Albanese A, Cheng L, Ursino M, Chbat NW. An integrated mathematical model of the human cardiopulmonary system: model development. Am J Physiol Circ Physiol 310: H899–H921, 2016.
% [2] Cheng L, Ivanova O, Fan H-H, Khoo MCK. An integrative model of respiratory and cardiovascular control in sleep-disordered breathing. Respir Physiol Neurobiol 174: 4–28, 2010.
% [3] Fincham WF, Tehrani FT. A mathematical model of the human respiratory system. J Biomed Eng 5: 125–33, 1983.
% [4] Magosso E, Ursino M. A mathematical model of CO2 effect on cardiovascular regulation. Am J Physiol Heart Circ Physiol 281: H2036-52, 2001.
% [5] Magosso E, Ursino M. Cardiovascular response to dynamic aerobic exercise: a mathematical model. Med Biol Eng Comput 40: 660–74, 2002.
% [6] Poon CS, Lin SL, Knudson OB. Optimization character of inspiratory neural drive. J Appl Physiol 72: 2005–17, 1992.
% [7] Serna Higuita LY, Mananas MA, Mauricio Hernandez A, Marina Sanchez J, Benito S. Novel Modeling of Work of Breathing for Its Optimization During Increased Respiratory Efforts. IEEE Syst J 10: 1003–1013, 2016.
% [8] Serna LY, Mañanas MA, Hernández AM, Rabinovich RA. An Improved Dynamic Model for the Respiratory Response to Exercise. Front Physiol 9: 69, 2018.
% [9] Spencer JL, Firouztale E, Mellins RB. Computational expressions for blood oxygen and carbon dioxide concentrations. Ann Biomed Eng 7: 59–66, 1979.
% [10] Ursino M. Interaction between carotid baroregulation and the pulsating heart: a mathematical model. Am J Physiol 275: H1733–H1747, 1998.
% [11] Ursino M, Magosso E. Acute cardiovascular response to isocapnic hypoxia. I. A mathematical model. Am J Physiol Heart Circ Physiol 279: H149-65, 2000.


