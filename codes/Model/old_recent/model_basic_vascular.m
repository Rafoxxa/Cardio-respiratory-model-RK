%function xdot= model_basic_vascular(time, init_values, delays, pars, init_keys, taus_keys)
function xdot= model_basic_vascular(t, y, pars, init_keys)
    time = t;
    init_values = y;    
    global all_global;
    if strcmp(class(pars), 'dictionary') 
        internal_variables = dictionary(); %this array is for local variables that are used per iteration between different modules
        
    else
        internal_variables = zeros(70);
    end
    
    %Basic parameters    
    dt = pars('dt');
    %tau = pars('tau');
    t0 = pars('t0');
    t = time;
    

    if strcmp(class(pars), 'dictionary') 
        y = dictionary(init_keys, init_values);
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
    test = dictionary();
    test__vascular_checking__28_06_24 = dictionary();

    

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
    xdot_dict = dictionary();
    if isa(y, 'dictionary')
        xdot = zeros(1, length(y.keys));   
    else
        xdot = zeros(1, length(y));  
    end
    xdot_dict("Pmusc") =                    dPmusc;
    xdot_dict("Ppl") =                      dPpl;
    xdot_dict("V") =                        dV;
    xdot_dict("dV") =                       ddV;
    xdot_dict("dVua") =                     ddVua;
    %xdot_dict("Vua") =                      dVua;
    xdot_dict("Gaw") =                      dGaw;   %fake derivative
    %xdot_dict("Pua") =                     dPua;   
    %xdot_dict("Nt") =                      y('Nt'); %fake derivative
    xdot_dict("dVE") =                      ddVE; %fake derivative
    xdot_dict("I") =                        I;      %fake derivative
    xdot_dict("vO2") =                      dv(1);
    xdot_dict("vCO2") =                     dv(2);
    xdot_dict("MRtO2") =                    dMRtgas(1);
    xdot_dict("MRtCO2") =                   dMRtgas(2);
    xdot_dict("PvbCO2") =                   dPvbCO2;
    xdot_dict("PCSFCO2") =                  dPCSFCO2;    
    xdot_dict("PAO2") =                     dPAgas(1);
    xdot_dict("PACO2") =                    dPAgas(2);
    xdot_dict("dPaO2") =                    ddPa(1);
    xdot_dict("dPaCO2") =                   ddPa(2);
    xdot_dict("PaO2") =                     y('dPaO2');
    xdot_dict("PaCO2") =                    y('dPaCO2');
    xdot_dict("aO2") =                      a(1);   %fake derivative
    xdot_dict("aCO2") =                     a(2);  %fake derivative
    xdot_dict("P_1O2") =                    dP_1(1);
    xdot_dict("P_1CO2") =                   dP_1(2);
    xdot_dict("P_2O2") =                    dP_2(1);
    xdot_dict("P_2CO2") =                   dP_2(2);
    xdot_dict("P_3O2") =                    dP_3(1);
    xdot_dict("P_3CO2") =                   dP_3(2);
    xdot_dict("P_4O2") =                    dP_4(1);
    xdot_dict("P_4CO2") =                   dP_4(2);
    xdot_dict("P_5O2") =                    dP_5(1);
    xdot_dict("P_5CO2") =                   dP_5(2);
    xdot_dict("MRR") =                      dMRR;
    xdot_dict("M_Rv") =                     dM_Rv;
    xdot_dict("MRv") =                      dMRv;
    xdot_dict('mean_PaO2') =                mean_dPaO2; 
    xdot_dict('mean_PaCO2') =               mean_dPaCO2; 
    xdot_dict('mean_PbCO2') =               mean_dPbCO2; 
    xdot_dict('mean_P_sa') =                mean_dP_sa; 
    xdot_dict('fake_TI') =                  y('TI'); %fake derivative
    xdot_dict('fake_Tresp') =               y('Tresp'); %fake derivative
    xdot_dict('Q_pa') =                     dQ_pa; 
    xdot_dict('V_total_pp') =               dV_total_pp;
    xdot_dict('V_total_pv') =               dV_total_pv;
    xdot_dict('V_total_pa') =               dV_total_pa;
    xdot_dict('V_total_vc') =               dV_total_vc;
    xdot_dict('Q_sa') =                     dQ_sa;
    xdot_dict('P_sa') =                     dP_sa;
    xdot_dict('V_total_e_v') =              dV_total_v(1);
    xdot_dict('V_total_s_v') =              dV_total_v(2);
    xdot_dict('V_total_b_v') =              dV_total_v(3);
    xdot_dict('V_total_h_v') =              dV_total_v(4);
    xdot_dict('V_total_rm_v') =             dV_total_v(5);
    xdot_dict('V_total_am_v') =             dV_total_v(6);
    xdot_dict('P_sp') =                     dP_sp;
    xdot_dict('V_total_rv') =               dV_total_rv;
    xdot_dict('V_total_ra') =               dV_total_ra;
    xdot_dict('V_total_la') =               dV_total_la;
    xdot_dict('V_total_lv') =               dV_total_lv;
    xdot_dict('zheta_heart') =              dzheta_heart;
    xdot_dict('P_mean') =                   dP_mean;
    xdot_dict('f_ac') =                     df_ac;
    xdot_dict('f_ap') =                     df_ap;
    xdot_dict('xO2_b') =                    dxO2_b;
    xdot_dict('xCO2_b') =                   dxCO2_b;
    xdot_dict('xO2_h') =                    dxO2(1);
    xdot_dict('xO2_rm') =                   dxO2(2);
    xdot_dict('xCO2_h') =                   dxCO2(1);
    xdot_dict('xCO2_rm') =                  dxCO2(2);
    xdot_dict('Wh') =                       dWh;
    xdot_dict('xO2_am') =                   dxO2_am;
    xdot_dict('x_met') =                    dx_met;
    xdot_dict('x_M') =                      dx_M;
    xdot_dict('DThetaO2_h_s') =             dDThetaO2_s(1);
    xdot_dict('DThetaO2_p_s') =             dDThetaO2_s(2);
    xdot_dict('DThetaO2_v_s') =             dDThetaO2_s(3);
    xdot_dict('DThetaCO2_h_s') =            dDThetaCO2_s(1);
    xdot_dict('DThetaCO2_p_s') =            dDThetaCO2_s(2);
    xdot_dict('DThetaCO2_v_s') =            dDThetaCO2_s(3);
    xdot_dict('DTsym') =                    dDTsym;
    xdot_dict('DTvagal') =                  dDTvagal;
    xdot_dict('DTheta_R_e_p') =             dDTheta(1);
    xdot_dict('DTheta_R_s_p') =             dDTheta(2);
    xdot_dict('DTheta_R_rm_p_n') =          dDTheta(3);
    xdot_dict('DTheta_R_am_p_n') =          dDTheta(4);
    xdot_dict('DTheta_V_unstressed_e_v') =  dDTheta(5);
    xdot_dict('DTheta_V_unstressed_s_v') =  dDTheta(6);
    xdot_dict('DTheta_V_unstressed_rm_v') = dDTheta(7);
    xdot_dict('DTheta_V_unstressed_am_v') = dDTheta(8);
    xdot_dict('DTheta_Emax_lv') =           dDTheta(9);
    xdot_dict('DTheta_Emax_rv') =           dDTheta(10);
    xdot_dict('phi_met') =                  dphi_met;
    xdot_dict('Qpp') =                      dQpp;  %forced derivative
    xdot_dict('Qbp') =                      dQbp;  %forced derivative
    xdot_dict('Qla') =                      dQla;  %forced derivative
    xdot_dict('Theart') =                   dTheart;
    xdot_dict('fh_s') =                     dfh_s;    
    xdot_dict('fp_s') =                     dfp_s;    
    xdot_dict('fv_s') =                     dfv_s;
    xdot_dict('fv') =                       dfv;    
    xdot_dict('R_bp') =                     dR_bp; 
    xdot_dict('xO2_e') =                     dxO2_e;
    xdot_dict('xO2_s') =                     dxO2_s;
    xdot_dict('xO2_p') =                     dxO2_p;   
    % xdot_dict("Q_e") =                      test__vascular_checking__28_06_24("dQ_e");
    % xdot_dict("Q_s") =                      test__vascular_checking__28_06_24("dQ_s");
    % xdot_dict("Q_h") =                      test__vascular_checking__28_06_24("dQ_h");
    % xdot_dict("Q_rm") =                     test__vascular_checking__28_06_24("dQ_rm");
    % xdot_dict("Q_am") =                     test__vascular_checking__28_06_24("dQ_am");
    % xdot_dict("P_v_e") =                    test__vascular_checking__28_06_24("dP_v_e");
    % xdot_dict("P_v_s") =                    test__vascular_checking__28_06_24("dP_v_s");
    % xdot_dict("P_v_h") =                    test__vascular_checking__28_06_24("dP_v_h");
    % xdot_dict("P_v_rm") =                   test__vascular_checking__28_06_24("dP_v_rm");
    % xdot_dict("P_v_am") =                   test__vascular_checking__28_06_24("dP_v_am");
    % xdot_dict("R_e_p") =                    test__vascular_checking__28_06_24("dR_e_p");
    % xdot_dict("R_s_p") =                    test__vascular_checking__28_06_24("dR_s_p");
    % xdot_dict("R_b_p") =                    test__vascular_checking__28_06_24("dR_b_p");
    % xdot_dict("R_h_p") =                    test__vascular_checking__28_06_24("dR_h_p");
    % xdot_dict("R_rm_p") =                   test__vascular_checking__28_06_24("dR_rm_p");
    % xdot_dict("R_am_p") =                   test__vascular_checking__28_06_24("dR_am_p"); 

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
%
    %xdot = (xdot > 10^6)  * 10^6  + xdot * (xdot <= 10^6); %protection measure against high derivatives, mainly for the fitting process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55



%% Functions

%EQUATIONS (CHENG-SARMIENTO-SERNA)

%% Respiratory mechanics   Eq [x] --- [n]


function  [internal_variables_, test_] = respiratory_pump(t, y, pars, internal_variables, test)
        
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
    VT = y('V');

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
    
    internal_variables('Pabd') = Pabd;
    internal_variables('Ptor') = Ptor;
    internal_variables_ = internal_variables;    
    test_ = test;
end

function [y_, ddVua, dGaw, dVua, internal_variables_, test_] = upper_airways(t, y, pars, internal_variables, test)
    
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
    
    %Equations   
    Rrs = Raw + Rl + Rcw;
    dVla = dVua + dV;
    Pua = Ppl + dVla * Rrs;        
    dPua_dt = (Pua - y('Pua'))/dt;
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

    dGaw = (Gaw - y('Gaw'))/dt;
    internal_variables('Gaw') = Gaw;
    internal_variables_ = internal_variables; 
    y_ = y;
    test_ = test;
end
function [y_, dPmusc, dPpl, dV, ddV, test_] = pulmonary_mechanics(t, y, pars, internal_variables, test)
    %pars 
    Ecw = pars('Ecw');
    El = pars('El');
    kaw1  = pars('kaw1');
    kaw2 = pars('kaw2');
    Rcw = pars('Rcw');
    Rrs = pars('Rrs');
    Pao = pars('Pao');
    dt = pars('dt');
    t0 = pars('t0');

    %vars    
    V = y('V');    
    a0 = y('a0');
    a1 = y('a1');
    a2 = y('a2');
    tau = y('tau');
    TI = y('TI');
    Gaw = internal_variables('Gaw');  

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
    ddV = (dV - y('dV'))/dt;
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
 
    dPpl = (Ppl - y('Ppl'))/dt;
    dPmusc = (Pmusc - y('Pmusc'))/dt;   

    y_ = y;
    test_ = test;
    
end


function [y_, internal_variables_, test_] = neuromuscular_drive(t, y, pars, internal_variables, index, tiny_y_keys, test)
    
    %pars
    dt = pars('dt');
    t0 = pars('t0');

    %vars
    TI = y('TI');
    
    %Equations
    t_cycle = t - t0; 
    Nt = eps;
    if t_cycle < TI
        dVE_historic = all_global(index(tiny_y_keys, 'dVE'), round(t0/dt) + 1: round(t/dt) + 1);
        Nt = trapz(dVE_historic, 2) * dt;
        
    end 
    internal_variables('Nt') = Nt;  %dejar derivada en 0 y guardar todas las variables externas en una 
    internal_variables_ = internal_variables; 
    y_ = y;
    test_ = test;
    
end
function [I, internal_variables_, test_] = metabolic_regulation(t, y, pars, internal_variables, test)
    %pars
    MRtCO2_basal = pars("MRtCO2_basal");
    AT = pars("AT");

    %vars
    MRtCO2 = y("MRtCO2");

    %Equations
    I = (MRtCO2 - MRtCO2_basal)/(AT - MRtCO2_basal);
    
    internal_variables('I') = I;
    internal_variables_ = internal_variables;
    test_ = test;

end

%%Gas exchange
function gas = dissociation(P, pars)

    %pars
    A1 = pars('A1'); % parameter in O2 dissociation equation
    A2 = pars('A2'); % parameter in CO3 dissociation equation
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
    %fO2 = pars('fO2');
    fCO2 = pars('fCO2');
    Patm = pars('Patm');
    Pws = pars('Pws');
    Vdead = pars('Vdead');
    VLO2 = pars('VLO2');
    VLCO2 = pars('VLCO2');
    T1 = pars('T1');
    T2 = pars('T2');
    LCTV = pars('LCTV');
    settling_time = pars("settling_time");
    type_of_input = pars("type_of_input");

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
    Qla = y('Qla');
    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');
    dPaO2 = y('dPaO2');
    dPaCO2 = y('dPaCO2');

    %fgas = [fO2, fCO2];
    P_1 = [P_1O2, P_1CO2];
    P_2 = [P_2O2, P_2CO2];
    P_3 = [P_3O2, P_3CO2];
    P_4 = [P_4O2, P_4CO2];
    P_5 = [P_5O2, P_5CO2];

    Ta = 1000 * LCTV/(Qla + 250);    %this has a stop limit because Qla could drop, making the delay truly huge.  
    if t > Ta
        try
            PACO2_delayed = all_global(index(tiny_y_keys, 'PACO2'), round((t-Ta)/dt) + 1);
            PAO2_delayed = all_global(index(tiny_y_keys, 'PAO2'), round((t-Ta)/dt) + 1);

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
    fO2p_0 = pars("fiO2_poly_0");
    fO2p_1 = pars("fiO2_poly_1");
    fO2p_2 = pars("fiO2_poly_2");
    fO2p_3 = pars("fiO2_poly_3");
    fO2p_4 = pars("fiO2_poly_4");

    if type_of_input > 6
        if t >= settling_time
            tt = t - settling_time;
            fO2 = fO2p_0 + fO2p_1 * tt^1 + fO2p_2 * tt^2 + fO2p_3 * tt^3 + fO2p_4 * tt^4;
            fO2 = 100*fO2; %data is given in decimals
        else
            fO2 = pars("fO2");
        end
    else
        fO2 = pars("fO2");
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

    y('aO2') = a(1);
    y('aCO2') = a(2);
    y_ = y;
    test_ = test;
end


%% 
function [dPvbCO2, dPCSFCO2, internal_variables_ , test_] = brain(t, y, pars, internal_variables, test)
    


    %pars
    KCSFCO2 = pars('KCSFCO2');
    KCCO2 = pars('KCCO2');
    dc = pars('dc');
    h = pars('h');
    SCO2 = pars('SCO2');
    SbCO2 = pars('SbCO2');
    MRbCO2 = pars('MRbCO2');
    dt = pars('dt');

    %vars
    PvbCO2 = y('PvbCO2');
    PCSFCO2 = y('PCSFCO2');
    PaCO2 = y('PaCO2');    
    Qbp = y('Qbp');

    %Equations
    dPvbCO2 = (MRbCO2 * 1 + Qbp * SCO2 * (PaCO2 - PvbCO2) - h)/SbCO2;  %UNIT CORRECTION, now is avoided 
    dPCSFCO2 = (PvbCO2 - PCSFCO2)/KCSFCO2;
    PbCO2 = PvbCO2 + (PCSFCO2 - PvbCO2) * exp(-dc * (Qbp * KCCO2)^0.5);
    %PbCO2 = 40;
    internal_variables('PbCO2') = PbCO2;
    internal_variables_ = internal_variables; 
    test_ = test;
end

function [dv, dMRtgas, test_] = tissue(t, y, pars, test)
    
    %Input variables- consumption-production rates
    pars = input_consumption(t, y, pars); %in this inside function we control the input

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
        type_of_input = pars('type_of_input');

        %pars
        MRO2 = pars('MRO2');
        MRCO2 = pars('MRCO2');
        MRtO2_basal = pars('MRtO2_basal');
        MRtCO2_basal = pars('MRtCO2_basal');
        
        
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
            
            pars('MRO2') = MRO2;
            pars('MRCO2') = MRCO2;

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
            settling_time = pars('settling_time');
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
            settling_time = pars('settling_time');
            tt = t - settling_time;
            if t >= settling_time
                MRO2p_0 = pars("MRO2_poly_0");
                MRO2p_1 = pars("MRO2_poly_1");
                MRO2p_2 = pars("MRO2_poly_2");
                MRO2p_3 = pars("MRO2_poly_3");
                MRO2p_4 = pars("MRO2_poly_4");
                MRO2p_5 = pars("MRO2_poly_5");
                MRO2p_6 = pars("MRO2_poly_6");
                MRO2p_7 = pars("MRO2_poly_7");
                MRO2p_8 = pars("MRO2_poly_8");

                MRCO2p_0 = pars("MRCO2_poly_0");
                MRCO2p_1 = pars("MRCO2_poly_1");
                MRCO2p_2 = pars("MRCO2_poly_2");
                MRCO2p_3 = pars("MRCO2_poly_3");
                MRCO2p_4 = pars("MRCO2_poly_4");
                

                MRO2 = MRO2p_0 + MRO2p_1*tt + MRO2p_2*tt^2 + MRO2p_3*tt^3 + MRO2p_4*tt^4 + MRO2p_5*tt^5 + MRO2p_6*tt^6 + MRO2p_7*tt^7 + MRO2p_8*tt^8;           
                MRCO2 = MRO2p_0 + MRO2p_1*tt + MRO2p_2*tt^2 + MRO2p_3*tt^3 + MRO2p_4*tt^4;
            end               

        end
        
        pars('MRO2') = MRO2;
        pars('MRCO2') = MRCO2;

        pars_ = pars;      
        
    end
    test_ = test;

end
function [dMRR, dM_Rv, dMRv, test_] = metabolism_dynamics(t, y, pars, test)
    %pars
    dt = pars('dt');    
    MRbCO2 = pars('MRbCO2');
    MRbO2 = pars('MRbO2');
    MRtCO2_basal = pars('MRtCO2_basal');
    MRtO2_basal = pars('MRtO2_basal');
    tauMRv = pars('tauMRv');

    %vars 
    M_Rv = y('M_Rv');
    MRtO2 = y('MRtO2');
    MRtCO2 = y('MRtCO2');
    
    %Equations
    M_RR = (MRbCO2 + MRbO2 + MRtCO2 + MRtO2)/(MRbCO2 + MRbO2 + MRtCO2_basal + MRtO2_basal); %this variable then goes to dVE computing and it's an expression of the overall metabolic excersise
    if M_RR >= 1
        MRR = M_RR;
    else
        MRR = 1;
    end

    dMRR = (MRR - y('MRR'))/dt;
    dM_Rv = ((MRR - 1) - M_Rv)/tauMRv;

    if M_Rv >= 0 && MRR > 1
        MRv = M_Rv;
    else
        MRv = 0;
    end

    dMRv = (MRv - y('MRv'))/dt;
    test_ = test;
end

%Ventilatory controller
function  [ddVE, test_] = ventilation_control(t, y, pars, internal_variables, test)
    
    %pars    
    KpCO2 = pars('KpCO2');
    KpO2 = pars('KpO2');
    KcMRv = pars('KcMRv');
    KcCO2 = pars('KcCO2');
    Kbg = pars('Kbg');
    dV_rest = pars('dV_rest');
    V0dead = pars('Vdead');
    GVdead = pars('GVdead');
    dt = pars('dt');

    %vars
    Tresp = y('Tresp');
    MRv = y('MRv');
    mean_PaO2 = y('mean_PaO2');
    mean_PaCO2 = y('mean_PaCO2');
    mean_PbCO2 = y('mean_PbCO2');
    
    dVA_ = dV_rest * (KpCO2 * mean_PaCO2 + KcCO2 * mean_PbCO2 + (KpO2 * (104 - mean_PaO2)^4.9) * (mean_PaO2 < 104) + KcMRv * MRv - Kbg); %This should be alveolar minute volume, and it's part of the minute ventilation we want the lungs to adquire
    dVA = dVA_ * (dVA_ > 0); % as minute ventilation is a positive value, we must take the absolute value, we will never be able to remove air from the lungs, the minimum volume value will always be dead space volume 
    dVA = real(dVA);
    Vd = GVdead * dVA + V0dead; %the first part refers to dead space volume that changes because respiration (expansion of airway channels?)
    dVd = 1/Tresp * Vd;     %the amount of volume that is exchanged during a minute due to dead space, will be the breathing frequency times dead space volume
    
    dVE = dVA + dVd;    %dVE corresponds to minute ventilation (how much volume do I want to exchange over a minute, it's a indicator of respiration flow and its different than instant flow)
    ddVE = (dVE - y("dVE"))/dt;
    all_global(index(tiny_y_keys, 'dVE'), round(t/dt) + 1) = dVE;
    
    test_ = test;
    %In the simulation is correct kind of the values we got (0.11 l/s, which are indeed 7 l/min, reported as common minute ventilation values)

end



function [mean_dPaO2, mean_dPaCO2, mean_dPbCO2, mean_dP_sa, test_] = mean_values_computer(t, y, pars, internal_variables, test)
    
    
    %vars
    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');
    P_sa = y('P_sa');    
    PbCO2 = internal_variables('PbCO2');
    mean_PaO2 = y('mean_PaO2');
    mean_PaCO2 = y('mean_PaCO2');
    mean_P_sa = y('mean_P_sa');
    mean_PbCO2 = y('mean_PbCO2');

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
        Aim = pars('Aim');
        Tc = pars('Tc');
        Tim = pars('Tim');
    
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
        internal_variables('P_im') = P_im;
        internal_variables_ = internal_variables;
    test_ = test;
end

function [internal_variables_, test_] = vena_cava(t, y, pars, internal_variables, test)

    %----------------------
    %pars
    %vena_cava
        D1 = pars('D1');
        D2 = pars('D2');
        K1_vc = pars('K1_vc');
        K2_vc = pars('K2_vc');
        K_r_vc = pars('K_r_vc');
        R_vc_n = pars('R_vc_n');
        V_unstressed_vc = pars('V_unstressed_vc');
        V_vc_max = pars('V_vc_max');
        V_vc_min = pars('V_vc_min');
    
    %pulmonary
        V_unstressed_ra = pars('V_unstressed_ra');
        C_ra = pars('C_ra');

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
    V_total_vc = y('V_total_vc');
    V_total_ra = y('V_total_ra');
    Ptor = internal_variables("Ptor");  %this is in mmHg
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
            P_vc_ = D1 + K1_vc * (V_vc - V_unstressed_vc);
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

        internal_variables('P_vc') = P_vc;
        internal_variables('V_vc') = V_vc;
        internal_variables('Q_ra') = Q_ra;
        internal_variables('P_ra') = P_ra;
        internal_variables('V_ra') = V_ra;
        internal_variables_ = internal_variables;
    test_ = test;
end


function [internal_variables_,  dzheta_heart, dV_total_rv, dV_total_ra, dV_total_la, dV_total_lv, dQla, dTheart, test_] = heart(t, y, pars, internal_variables, test)

    %--------------------------------
    %pars
        %heart
            C_la = pars('C_la');
            K_E_lv = pars('K_E_lv');
            K_E_rv = pars('K_E_rv');
            KR_lv = pars('KR_lv');
            KR_rv = pars('KR_rv');
            ksys = pars('ksys');
            P0_lv = pars('P0_lv');
            P0_rv = pars('P0_rv');
            R_la = pars('R_la');
            R_ra = pars('R_ra');
            Tsys_0 = pars('Tsys_0');
            V_unstressed_la = pars('V_unstressed_la');
            V_unstressed_lv = pars('V_unstressed_lv');
            V_unstressed_rv = pars('V_unstressed_rv');
            %for other modules
                %C_ra = pars('C_ra');
                %V_unstressed_ra = pars('V_unstressed_ra');
        
        %pulmonary
            V_unstressed_pa = pars('V_unstressed_pa');
            C_pa = pars('C_pa');
            V_unstressed_pv = pars('V_unstressed_pv');
            C_pv = pars('C_pv');
            R_pv = pars('R_pv');
        
        %general
        dt = pars('dt');
    
    %------------------------------------- 
    %vars    
    %for other modules
        %V_total_ra = y('V_total_ra');
     %this should come from internal variables
    V_total_lv = y('V_total_lv');
    V_total_la = y('V_total_la');
    V_total_rv = y('V_total_rv');
    V_total_pa = y('V_total_pa');
    V_total_pv = y('V_total_pv');
    P_sa = y('P_sa');              
    zheta_heart = y('zheta_heart');

    
    
    
    %------------------------------------
    %internal_variables
    Theart = internal_variables('Theart');
    Ptor = internal_variables('Ptor');  
    Q_ra = internal_variables('Q_ra');
    P_ra = internal_variables('P_ra');
    E_max_lv = internal_variables('Emax_lv');
    E_max_rv = internal_variables('Emax_rv');
    if t == 0
        E_max_lv = y('E_max_lv');
        E_max_rv = y('E_max_rv');
    end

    %--------------------------------------
    %Equations
    
    %Cardiac oscilator

        try  %to avoid the non existent first value
            t0_heart = all_global(index(tiny_y_keys, 't0_heart'), round(t/dt));
            u_t0 = all_global(index(tiny_y_keys, 'u_t0'), round(t/dt));
        catch
            t0_heart = 0;
            u_t0 = 0;
        end
        
        Tsys = Tsys_0 + ksys * 1/Theart; %systolic time

        %Activation function
            HR = 1/Theart;
            all_global(index(tiny_y_keys, 'HR'), round(t/dt) + 1) = HR;
            
            HR_integral = dt * sum(all_global(index(tiny_y_keys, 'HR'), round(t0_heart/dt) + 1: round(t/dt) + 1));
            dzheta_heart = HR;
            U = HR_integral + u_t0 - floor(HR_integral + u_t0); %fractional part
            u = zheta_heart - floor(zheta_heart); %fractional part
            
            if U == 0        
                all_global(index(tiny_y_keys, 't0_heart'), round(t/dt) + 1) = t;        
                all_global(index(tiny_y_keys, 'u_t0'), round(t/dt) + 1) = u;
                
            else
                all_global(index(tiny_y_keys, 't0_heart'), round(t/dt) + 1) = t0_heart;         
                all_global(index(tiny_y_keys, 'u_t0'), round(t/dt) + 1) = u_t0;
                
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
        dQla = (Q_la - y('Qla'))/dt;
        dTheart = (Theart - y('Theart'))/dt;

    % y('Q_lv') = Q_lv;
    % y('Q_rv') = Q_rv;
    % y('P_la') = P_la;
    y_ = y;

    %------------------------------
    % Saving computations
        internal_variables('Q_lv') = Q_lv;
        internal_variables('Q_la') = Q_la;
        internal_variables('Q_rv') = Q_rv;
        internal_variables('P_la') = P_la;
        internal_variables('P_pv') = P_pv;
        internal_variables('P_pa') = P_pa;
        internal_variables('V_pa') = V_pa;
        internal_variables('V_pv') = V_pv;
        internal_variables('V_la') = V_la;
        internal_variables('V_lv') = V_lv;
        internal_variables('V_rv') = V_rv;
        internal_variables('Wh_lv') = Wh_lv;
        internal_variables('Wh_rv') = Wh_rv;
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
    C_sa = pars('C_sa');
    L_sa = pars('L_sa');
    R_sa = pars('R_sa');
    V_unstressed_sa = pars('V_unstressed_sa');
    
    %-----------------------------
    %vars    
    P_sa = y('P_sa');
    Q_sa = y('Q_sa'); 
    %Ptor = y('Ptor');
    P_sp = y('P_sp');
    

    %------------------------------
    %internal_variables
    Ptor = internal_variables('Ptor');
    Q_lv = internal_variables('Q_lv');
    
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
    internal_variables('V_sa') = V_sa;
    internal_variables_ = internal_variables;
    
    y_ = y;
    test_ = test;
end

function [internal_variables_, dQ_pa, dV_total_pp, dV_total_pv, dV_total_pa, dQpp, test_] = pulmonary_circulation(t, y, pars, internal_variables, test)
    %------------------------------------
    %pars
    dt = pars('dt');
    C_pp = pars('C_pp');
    %R_pp = pars('R_pp'); now is calculated from internal variables (HIPOXIA)
    L_pa = pars('L_pa');
    R_pa = pars('R_pa');
    
    V_unstressed_pp = pars('V_unstressed_pp');
    % for other modules
        %C_pa = pars('C_pa');
        %C_pv = pars('C_pv');
        %R_pv = pars('R_pv');
        %V_unstressed_pa = pars('V_unstressed_pa');
        %V_unstressed_pv = pars('V_unstressed_pv');

    %-------------------------------------
    %vars
    % for other modules
        %Ptor = y('Ptor');
        %V_total_pa = y('V_total_pa');
        %V_total_pv = y('V_total_pv');    
    Q_pa = y('Q_pa');
    V_total_pp = y('V_total_pp');

    %------------------------------------
    %internal_variables   
    Q_rv = internal_variables('Q_rv');
    Q_la = internal_variables('Q_la');  
    P_pv = internal_variables('P_pv');
    P_pa = internal_variables('P_pa');
    R_pp = internal_variables('R_p_p');
    
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
        dQpp = (Q_pp - y('Qpp'))/dt;
    
    %---------------------------------------
    %Saving computations
        internal_variables('V_pp') = V_pp;
        internal_variables_ = internal_variables;

        y_ = y;
    test_ = test;

end

function [internal_variables_, dV_total_vc, dV_total_v, dP_sp, dQbp, test_] = systemic_peripheric_and_venous_circulation(t, y, pars, internal_variables, test)
    %---------------------------------
    %pars
    dt = pars('dt');
    %compartment order: e,s,b,h,rm,am
        %systemic_peripheric_and_venous_circulation    
            k_r_am = pars('k_r_am');
            P0 = pars('P0');
            C_e_p = pars('C_e_p');
            C_s_p = pars('C_s_p;');
            C_b_p = pars('C_b_p');
            C_h_p = pars('C_h_p');
            C_rm_p = pars('C_rm_p');
            C_am_p = pars('C_am_p');
            C_e_v = pars('C_e_v');
            C_s_v = pars('C_s_v');
            C_b_v = pars('C_b_v');
            C_h_v = pars('C_h_v');
            C_rm_v = pars('C_rm_v');
            C_am_v = pars('C_am_v');
            R_e_n = pars('R_e_n');
            R_s_n = pars('R_s_n');
            R_b_n = pars('R_b_n');
            R_h_n = pars('R_h_n');
            R_rm_n = pars('R_rm_n');
            R_am_n = pars('R_am_n');
            V_unstressed_e_p = pars('V_unstressed_e_p');
            V_unstressed_s_p = pars('V_unstressed_s_p'); 
            V_unstressed_b_p = pars('V_unstressed_b_p');
            V_unstressed_h_p = pars('V_unstressed_h_p');
            V_unstressed_rm_p = pars('V_unstressed_rm_p');
            V_unstressed_am_p = pars('V_unstressed_am_p');
            %V_unstressed_e_v = pars('V_unstressed_e_v');    %this is actually a var from carfiac control
            %V_unstressed_s_v = pars('V_unstressed_s_v');    %this is actually a var from carfiac control
            V_unstressed_b_v = pars('V_unstressed_b_v');
            V_unstressed_h_v = pars('V_unstressed_h_v');
            %V_unstressed_rm_v = pars('V_unstressed_rm_v');  %this is actually a var from carfiac control
            %V_unstressed_am_v = pars('V_unstressed_am_v');  %this is actually a var from carfiac control
            V_tot = pars('V_tot');

        %systemic arteries
            V_unstressed_sa = pars('V_unstressed_sa');
        %pulmonary
            V_unstressed_pa = pars('V_unstressed_pa');
            V_unstressed_pp = pars('V_unstressed_pp');
            V_unstressed_pv = pars('V_unstressed_pv');
        %heart
            V_unstressed_ra = pars('V_unstressed_ra');
            V_unstressed_la = pars('V_unstressed_la');

    %-----------------------------------------        
    %vars
    Q_sa = y('Q_sa');
    I = y('I');
    V_total_e_v = y('V_total_e_v');
    V_total_s_v = y('V_total_s_v');
    V_total_b_v = y('V_total_b_v');
    V_total_h_v = y('V_total_h_v');
    V_total_rm_v = y('V_total_rm_v');
    V_total_am_v = y('V_total_am_v');
    P_sp = y('P_sp');
    
    %internal_variables
    Pabd = internal_variables('Pabd');
    Q_ra = internal_variables('Q_ra');
    V_ra = internal_variables('V_ra'); 
    V_rv = internal_variables('V_rv');
    V_la = internal_variables('V_la');
    V_lv = internal_variables('V_lv');
    V_pa = internal_variables('V_pa');
    V_pp = internal_variables('V_pp'); 
    V_pv = internal_variables('V_pv'); 
    V_vc = internal_variables('V_vc'); 
    V_sa = internal_variables('V_sa'); 
    P_im = internal_variables('P_im'); 
    P_vc = internal_variables('P_vc'); 

    R_e_p = internal_variables('R_e_p');
    R_s_p = internal_variables('R_s_p');
    R_b_p = internal_variables('R_b_p');
    R_h_p = internal_variables('R_h_p');
    R_rm_p = internal_variables('R_rm_p');
    R_am_p = internal_variables('R_am_p');

    V_unstressed_e_v = internal_variables('V_unstressed_e_v');
    V_unstressed_s_v = internal_variables('V_unstressed_s_v');
    V_unstressed_rm_v = internal_variables('V_unstressed_rm_v');
    V_unstressed_am_v = internal_variables('V_unstressed_am_v');
    
    index = ["e","s","b","h","rm","am"];
    
    C_p = [C_e_p, C_s_p, C_b_p, C_h_p, C_rm_p, C_am_p];
    C_v = [C_e_v, C_s_v, C_b_v, C_h_v, C_rm_v, C_am_v];
    R_v_n = [R_e_n, R_s_n, R_b_n, R_h_n, R_rm_n, R_am_n];
    R_p = [R_e_p, R_s_p, R_b_p, R_h_p, R_rm_p, R_am_p];
    V_unstressed_p = [V_unstressed_e_p, V_unstressed_s_p, V_unstressed_b_p, V_unstressed_h_p, V_unstressed_rm_p, V_unstressed_am_p];
    V_unstressed_v = [V_unstressed_e_v, V_unstressed_s_v, V_unstressed_b_v, V_unstressed_h_v, V_unstressed_rm_v, V_unstressed_am_v];
    V_total_v = [V_total_e_v, V_total_s_v, V_total_b_v, V_total_h_v, V_total_rm_v, V_total_am_v];
   
    %Equations 
    %systemic_peripheric_and_venous_circulation
    P = Pabd * (index == "s");   %vectorial   
    
    V_v = (V_total_v - V_unstressed_v) .* (V_total_v >= V_unstressed_v);  %vectorial
    P_v = C_v.^(-1) .* V_v .* (V_total_v >= V_unstressed_v)  + (index == "am") * P_im + (P0 * (1 - (V_total_v./V_unstressed_v).^(-3/2))).* (V_total_v < V_unstressed_v).* (index == "am");  %vectorial
    
    V_u = V_unstressed_pp + V_unstressed_sa + V_unstressed_pa + V_unstressed_pv + V_unstressed_ra + V_unstressed_la + sum(V_unstressed_p) + sum(V_unstressed_v);
    
    V_s_v = sum( V_v.*(index == "s") );
    V_rm_v = sum( V_v.*(index == "rm") );
    V_am_v = sum( V_v.*(index == "am") );
    V_b_v = sum( V_v.*(index == "b") );
    V_h_v = sum( V_v.*(index == "h") );

    P_ev = 1/C_e_v * (V_tot - V_sa - V_ra - V_rv - V_la - V_lv - V_pa - V_pp - V_pv - V_s_v - V_rm_v - V_am_v - V_b_v - V_h_v - V_vc - V_u - P_sp * sum(C_p) ); % Check this with profe
    P_v = P_v .* (index ~= "e") + P_ev .*  (index == "e");  %we have to check this
    
    P_p = P_sp*ones(1,6);
    V_p = C_p .* P_p; %vectorial
    V_total_p = V_p + V_unstressed_p; %vectorial

    R__v = R_v_n .* (P_vc >= P) + R_v_n .* (P_v - P_vc)./(P_v - P) .* (P_vc < P); %vectorial
    R_v = R__v * (I == 0) + (R__v .* (index ~= "am") + k_r_am./V_total_v .* (index == "am")) .* (I > 0); %vectorial        
    Q_p = (P_sp - P_v)./R_p; %vectorial.
    Q_v = (P_v - P_vc)./R_v .* (P_v >= P_vc); %vectorial    
    
    dV_total_v = Q_p - Q_v; %vectorial
    dP_sp = 1/(sum(C_p)) * (Q_sa - sum(Q_p)); %vectorial

    %vena_cava
        Q_vc = sum(Q_v);
        dV_total_vc = Q_vc - Q_ra;

    %for correct intenal_variables handling:
    Q_bp = Q_p(3);
    dQbp = (Q_bp - y('Qbp'))/dt;

    %["e","s","b","h","rm","am"];
    Q_e = Q_p(1);
    Q_s = Q_p(2);
    Q_h = Q_p(4);
    Q_rm = Q_p(5);
    Q_am = Q_p(6);
    dQ_e = (Q_e - y('Q_e'))/dt;
    dQ_s = (Q_s - y('Q_s'))/dt;
    dQ_h = (Q_h - y('Q_h'))/dt;
    dQ_rm = (Q_rm - y('Q_rm'))/dt;
    dQ_am = (Q_am - y('Q_am'))/dt;
    P_v_e = P_v(1);
    P_v_s = P_v(2);
    P_v_h = P_v(4);
    P_v_rm = P_v(5);
    P_v_am = P_v(6);
    dP_v_e = (P_v_e - y('P_v_e'))/dt;
    dP_v_s = (P_v_s - y('P_v_s'))/dt;
    dP_v_h = (P_v_h - y('P_v_h'))/dt;
    dP_v_rm = (P_v_rm - y('P_v_rm'))/dt;
    dP_v_am = (P_v_am - y('P_v_am'))/dt;
    dR_e_p = (R_e_p - y('R_e_p'))/dt;
    dR_s_p = (R_s_p - y('R_s_p'))/dt;
    dR_b_p = (R_b_p - y('R_b_p'))/dt;
    dR_h_p = (R_h_p - y('R_h_p'))/dt;
    dR_rm_p = (R_rm_p - y('R_rm_p'))/dt;
    dR_am_p = (R_am_p - y('R_am_p'))/dt;

    test("dQ_e") = dQ_e;
    test("dQ_s") = dQ_s;
    test("dQ_h") = dQ_h;
    test("dQ_rm") = dQ_rm;
    test("dQ_am") = dQ_am;    
    test("dP_v_e") = dP_v_e;
    test("dP_v_s") = dP_v_s;
    test("dP_v_h") = dP_v_h;
    test("dP_v_rm") = dP_v_rm;
    test("dP_v_am") = dP_v_am;
    test("dR_e_p") = dR_e_p;
    test("dR_s_p") = dR_s_p;
    test("dR_b_p") = dR_b_p;
    test("dR_h_p") = dR_h_p;
    test("dR_rm_p") = dR_rm_p;
    test("dR_am_p") = dR_am_p;

    %Saving computations    
    internal_variables('Q_am_p') = Q_p(6);
    internal_variables('Q_h_p') = Q_p(4);
    internal_variables('Q_rm_p') = Q_p(5);
    internal_variables('Q_b_p') = Q_p(3);
    internal_variables('Q_e_p') = Q_p(1);
    internal_variables('Q_s_p') = Q_p(2);
    internal_variables_ = internal_variables;
    test_ = test;
end

function [dP_mean, internal_variables_] = afferent_barorreflex(t,y,pars, internal_variables, dP_sa)

    %pars
    f_ab_min = pars("f_ab_min");
    f_ab_max = pars("f_ab_max");
    kab = pars("kab");
    P_n = pars("P_n");
    tau_p = pars("tau_p");
    tau_z = pars("tau_z");
    
    %vars
    P_mean = y('P_mean');
    P_sa = y('P_sa');


    %Equations

    %Derivatives
    dP_mean = 1/tau_p * (P_sa + tau_z * dP_sa - P_mean);    
    fab = 1/(1 + exp((P_mean - P_n)/kab)) * (f_ab_min + f_ab_max * exp((P_mean - P_n)/kab));

    %Saving computations
    internal_variables('fab') = fab;
    internal_variables_ = internal_variables;
    
end

function  [df_ac] = afferent_chemoreceptor(t,y,internal_variables)

    %pars
    f_ac_CO2_n = pars('f_ac_CO2_n');
    f_ac_max = pars('f_ac_max');
    f_ac_min = pars('f_ac_min');
    kac = pars('kac');
    KH = pars('KH');
    PaO2_ac_n = pars('PaO2_ac_n');
    PaCO2_n = pars('PaCO2_n');
    tau_ac = pars('tau_ac');


    %vars
    f_ac = y('f_ac');
    PaO2 = y('PaO2');
    PaCO2 = y('PaCO2');

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
    G_ap = pars('G_ap');
    tau_ap = pars('tau_ap');

    %vars
    f_ap = y('f_ap');
    V = y('V');

    %Equations
    phi_ap = G_ap * V;
    %Derivatives
    df_ap = 1/tau_ap * (phi_ap - f_ap);
end

function [dxO2_b, dxCO2_b, dR_bp, internal_variables_] = cerebral_blood_flow(t, y, pars, internal_variables)

    %pars
    A = pars('A');
    B = pars('B');
    C = pars('C');
    D = pars('D');
    dt = pars('dt');
    vO2_b_n = pars('vO2_b_n');
    gO2_b = pars('gO2_b');
    MO2_bp = pars('MO2_bp');
    R_bmp = pars('R_bmp');
    tau_CO2 = pars('tau_CO2');
    tau_O2 = pars('tau_O2');
    PaCO2_n = pars('PaCO2_n');

    %vars
    xO2_b = y('xO2_b');
    xCO2_b = y('xCO2_b');
    aO2 = y('aO2');
    PaCO2 = y('PaCO2');

    %internal variables
    %Q_b_p = internal_variables('Q_b_p');
    Q_b_p = y("Qbp");

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

    internal_variables('R_bp') = R_bp;
    
    dR_bp =  (R_bp - y("R_bp"))/dt;
    internal_variables_ = internal_variables;
    
    
    

end

function [dxO2_e, dxO2_s, dxO2_p, internal_variables_] = hipoxia_local_regulation(t, y, pars, internal_variables) 

    % local control splacnic + extrasplacnic

    %pars    
    dt = pars('dt');
    vO2_e_n = pars('vO2_e_n');
    vO2_s_n = pars('vO2_s_n');
    aO2_n = pars('aO2_n');

    gO2_e = pars('gO2_e');
    gO2_s = pars('gO2_s');
    gO2_p = pars('gO2_p');

    MO2_e = pars('MO2_e');
    MO2_s = pars('MO2_s');
    MO2_p = pars('MO2_p');

    R_e_p_n = internal_variables('R_e_p_n');
    R_s_p_n = internal_variables('R_s_p_n');
    R_p_p_n = pars('R_p_p_n');
    
    tau_CO2 = pars('tau_CO2');
    tau_O2 = pars('tau_O2');

    %vars
    xO2_e = y('xO2_e');
    xCO2_e = y('xCO2_e');
    xO2_s = y('xO2_s');
    xCO2_s = y('xCO2_s');
    xO2_p = y('xO2_p');
    xCO2_p = y('xCO2_p');
    aO2 = y('aO2');
    

    %internal variables
    Q_e_p = internal_variables('Q_e_p');
    Q_s_p = internal_variables('Q_s_p');
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

    internal_variables('R_e_p') = R_ep;
    internal_variables('R_s_p') = R_sp;
    internal_variables('R_p_p') = R_pp;

    internal_variables_ = internal_variables;



      


end


function  [dxO2, dxCO2, dWh, internal_variables_] = coronary_and_resting_muscle_blood_flow(t, y, pars, internal_variables)

    %pars
    vO2_h_n = pars('vO2_h_n');   
    vO2_rm_n = pars('vO2_rm_n');   
    gO2_h = pars('gO2_h');  
    gO2_rm = pars('gO2_rm');  
    KCO2_h = pars('KCO2_h');   
    KCO2_rm = pars('KCO2_rm');   
    MO2_h_p_n = pars('MO2_h_p_n');   
    MO2_rm_p = pars('MO2_rm_p');   
    R_h_p_n = pars('R_h_p_n');   
    %R_rm_p_n = pars('R_rm_p_n');   
    R_rm_p_n = internal_variables('R_rm_p_n');   
    tau_w = pars('tau_w');   
    tau_O2 = pars('tau_O2');
    tau_CO2 = pars('tau_CO2');
    Whn = pars('Whn'); 
    PaCO2_n = pars('PaCO2_n'); 
    

    vO2_n = [vO2_h_n, vO2_rm_n];
    gO2 = [gO2_h, gO2_rm];
    KCO2 = [KCO2_h, KCO2_rm];
    Rp_n = [R_h_p_n, R_rm_p_n];
    MO2_p = [MO2_h_p_n, MO2_rm_p]; %this changes afterwards in the code
    
    %vars
    xO2_h = y('xO2_h');
    xO2_rm = y('xO2_rm');
    xCO2_h = y('xCO2_h');
    xCO2_rm = y('xCO2_rm');
    Wh = y('Wh');
    PaCO2 = y('PaCO2');
    aO2 = y('aO2');

    xO2 = [xO2_h, xO2_rm];
    xCO2 = [xCO2_h, xCO2_rm];

    %internal variables
    Wh_lv = internal_variables('Wh_lv');
    Wh_rv = internal_variables('Wh_rv');
    Q_h_p = internal_variables('Q_h_p');
    Q_rm_p = internal_variables('Q_rm_p');
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
    
    internal_variables('R_h_p') = Rp(1);
    internal_variables('R_rm_p') = Rp(2);
    internal_variables_ = internal_variables;


    
    
end

function [dxO2_am, dx_met, dx_M, dphi_met, internal_variables_] = active_muscle_blood_flow(t,y,pars, internal_variables, index_fun)
    
    %pars
    vO2_am_n = pars('vO2_am_n');
    delay_met = pars('delay_met');
    gO2_am = pars('gO2_am');
    g_M = pars('g_M');
    I0_met = pars('I0_met');
    kmet = pars('kmet');
    MO2_am_p_n = pars('MO2_am_p_n');
    phi_max = pars('phi_max');
    phi_min = pars('phi_min');
    tau_M = pars('tau_M');
    tau_O2 = pars('tau_O2');
    tau_CO2 = pars('tau_CO2');
    tau_met = pars('tau_met');

    %vars
    aO2 = y('aO2');
    xO2_am  = y('xO2_am');
    x_M = y('x_M');
    x_met = y('x_met');
    %phi_met_delayed = all_global(:,round(t - delay_met/dt) + 1:end);  %we have to put the dictionary here
    phi_met_delayed = get_delayed_value(tiny_y_keys, t, delay_met, dt, all_global, index_fun,0, 'phi_met');
    
    %internal variables
    I = internal_variables('I');
    Q_am_p = internal_variables('Q_am_p');
    %Equations
    R_amp_n = internal_variables('R_am_p_n');
    R_am_p = R_amp_n/(1 + xO2_am + x_met);
    MO2_am_p = MO2_am_p_n * (1 + x_M);
    vO2_am = aO2 - MO2_am_p/Q_am_p;
    
    phi_met = (phi_min + phi_max * exp((I - I0_met)/kmet))/(1 + exp((I - I0_met)/kmet));

    %Derivatives
    dxO2_am = 1/tau_O2 * (-xO2_am - gO2_am * (vO2_am - vO2_am_n));
    dx_M = 1/tau_M * (-x_M + g_M * I);
    dx_met = 1/tau_met * (-x_met + phi_met_delayed);
    internal_variables('R_am_p') = R_am_p;
    internal_variables_  = internal_variables;
    dphi_met = phi_met - y('phi_met');

    function delayed_value = get_delayed_value(tiny_y_keys, t, delay, dt, all_global, index_fun, fj_init, fs)
        if delay > t
            % If delay is greater than current time, assign initial value
            delayed_value = fj_init;
        else
            % Otherwise, calculate the delayed value using the original code
            delayed_value = all_global(index_fun(tiny_y_keys, fs), round((t - delay)/dt) + 1);
        end
    end

    

end

function [dDThetaO2_s, dDThetaCO2_s, internal_variables_] = cns_ischemic_response(t, y, pars, internal_variables)
     %pars
     gcc_h_s  = pars('gcc_h_s');
     gcc_p_s = pars('gcc_p_s');
     gcc_v_s = pars('gcc_v_s');
     k_isc_h_s = pars('k_isc_h_s');
     k_isc_p_s = pars('k_isc_p_s');
     k_isc_v_s = pars('k_isc_v_s');
     PO2_ref_h_s = pars('PO2_ref_h_s');
     PO2_ref_p_s = pars('PO2_ref_p_s');
     PO2_ref_v_s = pars('PO2_ref_v_s');
     tau_cc = pars('tau_cc');
     tau_isc = pars('tau_isc');
     Theta_h_s_n = pars('Theta_h_s_n');
     Theta_p_s_n = pars('Theta_p_s_n');
     Theta_v_s_n = pars('Theta_v_s_n');
     x_h_s = pars('x_h_s');
     x_p_s = pars('x_p_s');
     x_v_s = pars('x_v_s');
     PaCO2_n = pars('PaCO2_n');

     gcc_s = [gcc_h_s, gcc_p_s, gcc_v_s];
     k_isc_s = [k_isc_h_s, k_isc_p_s, k_isc_v_s];
     PO2_ref_s = [PO2_ref_h_s, PO2_ref_p_s, PO2_ref_v_s];
     Theta_s_n = [Theta_h_s_n, Theta_p_s_n, Theta_v_s_n];
     xs = [x_h_s, x_p_s, x_v_s];

     %vars
     PaO2 = y('PaO2'); 
     PaCO2 = y('PaCO2');
     DThetaO2_h_s = y('DThetaO2_h_s');
     DThetaO2_p_s = y('DThetaO2_p_s');
     DThetaO2_v_s = y('DThetaO2_v_s');
     DThetaCO2_h_s = y('DThetaCO2_h_s');
     DThetaCO2_p_s = y('DThetaCO2_p_s');
     DThetaCO2_v_s = y('DThetaCO2_v_s');
     DThetaCO2_s = [DThetaCO2_h_s, DThetaCO2_p_s, DThetaCO2_v_s];
     DThetaO2_s = [DThetaO2_h_s, DThetaO2_p_s, DThetaO2_v_s];


     %internal_variables
     %Equations
     ws = xs ./ (1 + exp((PaO2 - PO2_ref_s) ./ k_isc_s)); %vectorial
     Theta_s = Theta_s_n - DThetaO2_v_s - DThetaCO2_s; %vectorial

     %Derivatives
     dDThetaCO2_s = 1/tau_cc * (-DThetaCO2_s + gcc_s * (PaCO2 - PaCO2_n)); %vectorial
     dDThetaO2_s = 1/tau_isc * (-DThetaO2_s + ws); %vectorial

     internal_variables('Theta_h_s') = Theta_s(1);
     internal_variables('Theta_p_s') = Theta_s(2);
     internal_variables('Theta_v_s') = Theta_s(3);

     internal_variables_ = internal_variables;

end

function [internal_variables_] = efferent_pathways(t, y, pars, internal_variables)
    
    %pars
    fab_0 = pars('fab_0');
    fes_0 = pars('fes_0');
    fes_inf = pars('fes_inf');
    fes_max = pars('fes_max');
    fev_0 = pars('fev_0');
    fev_inf = pars('fev_inf');
    kes = pars('kes');
    kev = pars('kev');
    I_0_h_s = pars('I_0_h_s');
    I_0_p_s = pars('I_0_p_s');
    I_0_v_s = pars('I_0_v_s');
    I_0_v = pars('I_0_v');
    kcc_h_s = pars('kcc_h_s');
    kcc_p_s = pars('kcc_p_s');
    kcc_v_s = pars('kcc_v_s');
    kcc_v = pars('kcc_v');
    gamma_h_s_max = pars('gamma_h_s_max');
    gamma_p_s_max = pars('gamma_p_s_max');
    gamma_v_s_max = pars('gamma_v_s_max');
    gamma_v_max = pars('gamma_v_max');
    gamma_h_s_min = pars('gamma_h_s_min');
    gamma_p_s_min = pars('gamma_p_s_min');
    gamma_v_s_min = pars('gamma_v_s_min');
    gamma_v_min = pars('gamma_v_min');
    Theta_v = pars('Theta_v');
    Wb_h_s = pars('Wb_h_s');
    Wb_p_s = pars('Wb_p_s');
    Wb_v_s = pars('Wb_v_s');
    Wc_h_s = pars('Wc_h_s');
    Wc_p_s = pars('Wc_p_s');
    Wc_v_s = pars('Wc_v_s');
    Wc_v = pars('Wp_v');
    Wp_h_s = pars('Wp_h_s');
    Wp_p_s = pars('Wp_p_s');
    Wp_v_s = pars('Wp_v_s');
    Wp_v = pars('Wp_v');
    Wt_h_s = pars('Wt_h_s');
    Wt_p_s = pars('Wt_p_s');
    Wt_v_s = pars('Wt_v_s');
    Wt_v = pars('Wt_v');
    
    I_0 = [I_0_h_s, I_0_p_s, I_0_v_s];
    kcc = [kcc_h_s, kcc_p_s, kcc_v_s];
    gamma_max = [gamma_h_s_max, gamma_p_s_max, gamma_v_s_max];
    gamma_min = [gamma_h_s_min, gamma_p_s_min, gamma_v_s_min];
    Wb_s = [Wb_h_s, Wb_p_s, Wb_v_s];
    Wc_s = [Wc_h_s, Wc_p_s, Wc_v_s];
    Wp_s = [Wp_h_s, Wp_p_s, Wp_v_s];
    Wt_s = [Wt_h_s, Wt_p_s, Wt_v_s]; 

    %vars
    fac = y('f_ac');
    fap = y('f_ap');

    %internal variables
    I = internal_variables('I');
    fab = internal_variables('fab');
    Nt = internal_variables('Nt');    
    
    Theta_h_s = internal_variables('Theta_h_s');
    Theta_p_s = internal_variables('Theta_p_s');
    Theta_v_s = internal_variables('Theta_v_s');

    Theta_s = [Theta_h_s, Theta_p_s, Theta_v_s];

    %Equations
    gamma_s = (gamma_min + gamma_max.* exp((I - I_0)./(kcc))) ./ (1 + exp((I - I_0)./(kcc))); %vectorial
    gamma_v = (gamma_v_min + gamma_v_max * exp((I - I_0_v)./(kcc_v))) ./ (1 + exp((I - I_0_v)./(kcc_v))); 
    fas = Wt_s * Nt + Wb_s * fab + Wc_s * fac + Wp_s * fap - Theta_s; %vectorial
    fs = fes_inf + (fes_0 + fes_inf) * exp(kes * fas) + gamma_s; %vectorial
    fs = fes_max.* (fs >= fes_max) + fs .* (fs < fes_max); %vectorial
    fv = (fev_0 + fev_inf * exp((fab - fab_0)/(kev)))/(1 + exp((fab - fab_0)/(kev))) - Wt_v * Nt + Wc_v * fac + Wp_v * fap - Theta_v + gamma_v;
    
    %internal variables
    internal_variables('f_h_s') = fs(1);
    internal_variables('f_p_s') = fs(2);
    internal_variables('f_v_s') = fs(3);
    internal_variables('fv') = fv;

    internal_variables_ = internal_variables;
    
    %Derivatives
end

function [dDTheta, dfh_s, dfp_s, dfv_s, internal_variables_] = reflex_control_R_Vu_E(t, y, pars, internal_variables, index_fun)

    %pars
    delay_Emax_lv = pars('delay_Emax_lv');
    delay_Emax_rv = pars('delay_Emax_rv');
    delay_R_am_p = pars('delay_R_am_p');
    delay_R_e_p = pars('delay_R_e_p');
    delay_R_rm_p = pars('delay_R_rm_p');
    delay_R_s_p = pars('delay_R_s_p');
    delay_V_u_am_v = pars('delay_V_u_am_v');
    delay_V_u_e_v = pars('delay_V_u_e_v');
    delay_V_u_rm_v = pars('delay_V_u_rm_v');
    delay_V_u_s_v = pars('delay_V_u_s_v');
    Emax_lv_0 = pars('Emax_lv_0');
    Emax_rv_0 = pars('Emax_rv_0');
    R_am_p_0 = pars('R_am_p_0');
    R_e_p_0 = pars('R_e_p_0');
    R_rm_p_0 = pars('R_rm_p_0');
    R_s_p_0 = pars('R_s_p_0');
    V_u_am_v_0 = pars('V_u_am_v_0');
    V_u_rm_v_0 = pars('V_u_rm_v_0');
    V_u_e_v_0 = pars('V_u_e_v_0');
    V_u_s_v_0 = pars('V_u_s_v_0');
    G_Emax_lv = pars('G_Emax_lv');
    G_Emax_rv = pars('G_Emax_rv');
    G_R_am_p = pars('G_R_am_p');
    G_R_e_p = pars('G_R_e_p');
    G_R_rm_p = pars('G_R_rm_p');
    G_R_s_p = pars('G_R_s_p');
    G_V_u_am_v = pars('G_V_u_am_v');
    G_V_u_e_v = pars('G_V_u_e_v');
    G_V_u_rm_v = pars('G_V_u_rm_v');
    G_V_u_s_v = pars('G_V_u_s_v');
    tau_Emax_lv = pars('tau_Emax_lv');
    tau_Emax_rv = pars('tau_Emax_rv');
    tau_R_am_p = pars('tau_R_am_p');
    tau_R_e_p = pars('tau_R_e_p');
    tau_R_rm_p = pars('tau_R_rm_p');
    tau_R_s_p = pars('tau_R_s_p');
    tau_V_u_am_v = pars('tau_V_u_am_v');
    tau_V_u_e_v = pars('tau_V_u_e_v');
    tau_V_u_rm_v = pars('tau_V_u_rm_v');
    tau_V_u_s_v = pars('tau_V_u_s_v');
    fes_min = pars('fes_min');
    
    tauTheta = [tau_R_e_p, tau_R_s_p, tau_R_rm_p, tau_R_am_p, tau_V_u_e_v, tau_V_u_s_v, tau_V_u_rm_v, tau_V_u_am_v, tau_Emax_lv, tau_Emax_rv];
    Theta0 = [R_e_p_0, R_s_p_0, R_rm_p_0, R_am_p_0, V_u_e_v_0, V_u_s_v_0, V_u_rm_v_0, V_u_am_v_0, Emax_lv_0, Emax_rv_0];
    GTheta = [G_R_e_p, G_R_s_p, G_R_rm_p, G_R_am_p, G_V_u_e_v, G_V_u_s_v, G_V_u_rm_v, G_V_u_am_v, G_Emax_lv, G_Emax_rv];

    %vars
    DTheta_R_e_p = y('DTheta_R_e_p');
    DTheta_R_s_p = y('DTheta_R_s_p');
    DTheta_R_rm_p_n = y('DTheta_R_rm_p_n');
    DTheta_R_am_p_n = y('DTheta_R_am_p_n');
    DTheta_V_unstressed_e_v = y('DTheta_V_unstressed_e_v');
    DTheta_V_unstressed_s_v = y('DTheta_V_unstressed_s_v');
    DTheta_V_unstressed_rm_v = y('DTheta_V_unstressed_rm_v');
    DTheta_V_unstressed_am_v = y('DTheta_V_unstressed_am_v');
    DTheta_Emax_lv = y('DTheta_Emax_lv');
    DTheta_Emax_rv = y('DTheta_Emax_rv');
    

    DTheta = [DTheta_R_e_p, DTheta_R_s_p, DTheta_R_rm_p_n, DTheta_R_am_p_n, DTheta_V_unstressed_e_v, DTheta_V_unstressed_s_v, DTheta_V_unstressed_rm_v, DTheta_V_unstressed_am_v, DTheta_Emax_lv, DTheta_Emax_rv];
    %internal variables
    f_h_s_actual = internal_variables('f_h_s');  
    f_p_s_actual = internal_variables('f_p_s');
    f_v_s_actual = internal_variables('f_v_s');

    f_h_s_delayed_Emaxlv = get_delayed_value(tiny_y_keys, t, delay_Emax_lv, dt, all_global, index_fun, f_h_s_actual, 'fh_s');
    f_h_s_delayed_Emaxrv = get_delayed_value(tiny_y_keys, t, delay_Emax_rv, dt, all_global, index_fun, f_h_s_actual, 'fh_s');
    f_p_s_delayed_Rep = get_delayed_value(tiny_y_keys, t, delay_R_e_p, dt, all_global, index_fun, f_p_s_actual, 'fp_s');
    f_p_s_delayed_Rsp = get_delayed_value(tiny_y_keys, t, delay_R_s_p, dt, all_global, index_fun, f_p_s_actual, 'fp_s');
    f_p_s_delayed_Rrmpn = get_delayed_value(tiny_y_keys, t, delay_R_rm_p, dt, all_global, index_fun, f_p_s_actual, 'fp_s');
    f_p_s_delayed_Rampn = get_delayed_value(tiny_y_keys, t, delay_R_am_p, dt, all_global, index_fun, f_p_s_actual, 'fp_s');
    f_v_s_delayed_Vuev = get_delayed_value(tiny_y_keys, t, delay_V_u_e_v, dt, all_global, index_fun, f_v_s_actual, 'fv_s');
    f_v_s_delayed_Vusv = get_delayed_value(tiny_y_keys, t, delay_V_u_s_v, dt, all_global, index_fun, f_v_s_actual, 'fv_s');
    f_v_s_delayed_Vurmv = get_delayed_value(tiny_y_keys, t, delay_V_u_rm_v, dt, all_global, index_fun, f_v_s_actual, 'fv_s');
    f_v_s_delayed_Vuamv = get_delayed_value(tiny_y_keys, t, delay_V_u_am_v, dt, all_global, index_fun, f_v_s_actual, 'fv_s');

    dfh_s = (f_h_s_actual - y('fh_s'))/dt;   
    dfp_s = (f_p_s_actual - y('fp_s'))/dt;    
    dfv_s = (f_v_s_actual - y('fv_s'))/dt;
    
    fs_delayed = [f_p_s_delayed_Rep, f_p_s_delayed_Rsp, f_p_s_delayed_Rrmpn, f_p_s_delayed_Rampn, f_v_s_delayed_Vuev, f_v_s_delayed_Vusv, f_v_s_delayed_Vurmv, f_v_s_delayed_Vuamv, f_h_s_delayed_Emaxlv, f_h_s_delayed_Emaxrv];
    fs_actual = [f_p_s_actual, f_p_s_actual, f_p_s_actual, f_p_s_actual, f_v_s_actual, f_v_s_actual, f_v_s_actual, f_v_s_actual, f_h_s_actual, f_h_s_actual];
    %testing
    %fs_delayed = fs_actual;

    %Equations
    
    sigmaTheta = (fs_actual >= fes_min) .* (GTheta.* log(fs_delayed - fes_min + 1)); 
    Theta = DTheta + Theta0;  %vectorial   (this are the actual effectors we have to change for the next iteration)
    
    %Derivatives 
    dDTheta = (tauTheta.^(-1)).* (-DTheta + sigmaTheta);  %vectorial

    internal_variables('R_e_p') = Theta(1);
    internal_variables('R_s_p') = Theta(2);
    internal_variables('R_rm_p_n') = Theta(3);
    internal_variables('R_am_p_n') = Theta(4);
    internal_variables('V_unstressed_e_v') = Theta(5);
    internal_variables('V_unstressed_s_v') = Theta(6);
    internal_variables('V_unstressed_rm_v') = Theta(7);
    internal_variables('V_unstressed_am_v') = Theta(8);
    internal_variables('Emax_lv') = Theta(9);
    internal_variables('Emax_rv') = Theta(10);

    internal_variables_ = internal_variables;


    function delayed_value = get_delayed_value(init_keys, t, delay, dt, all_global, index_fun, fj_init, fs)
        if delay > t
            % If delay is greater than current time, assign initial value
            delayed_value = fj_init;
        else
            % Otherwise, calculate the delayed value using the original code
            delayed_value = all_global(index_fun(init_keys, fs), round((t - delay)/dt) + 1);
        end
    end
end

function [dDTsym, dDTvagal, dfv, internal_variables_] = reflex_control_HR(t, y, pars, internal_variables, index_fun)
    %pars
    delay_Tsym = pars('delay_Tsym');
    delay_Tvagal = pars('delay_Tvagal');
    GTsym = pars('GTsym');
    GTvagal = pars('GTvagal');
    tau_Tsym = pars('tau_Tsym');
    tau_Tvagal = pars('tau_Tvagal');
    fes_min = pars('fes_min');
    %vars
    DTsym = y('DTsym');
    DTvagal = y('DTvagal');
    %internal variables
    f_h_s_actual = internal_variables('f_h_s');
    f_v_actual = internal_variables('fv');
    f_h_s_delayed = get_delayed_value(tiny_y_keys, t, delay_Tsym, dt, all_global, index_fun, 3.8576, 'fh_s'); %replace f_s_h_actual
    fv_delayed = get_delayed_value(tiny_y_keys, t, delay_Tvagal, dt, all_global, index_fun, 4.2748, 'fv');    %replace f_v_actual

    dfv = (f_v_actual - y('fv'))/dt;

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
            delayed_value = all_global(index_fun(init_keys, fs), round((t - delay)/dt) + 1);
        end
    end

end

function [internal_variables_] = adding_base_values_in_control_variables(t, y, pars, internal_variables)


    Emax_lv_0 = pars('Emax_lv_0');
    Emax_rv_0 = pars('Emax_rv_0');
    R_am_p_0 = pars('R_am_p_0');
    R_e_p_0 = pars('R_e_p_0');
    R_rm_p_0 = pars('R_rm_p_0');
    R_s_p_0 = pars('R_s_p_0');
    V_u_am_v_0 = pars('V_u_am_v_0');
    V_u_rm_v_0 = pars('V_u_rm_v_0');
    V_u_e_v_0 = pars('V_u_e_v_0');
    V_u_s_v_0 = pars('V_u_s_v_0');
    T0 = pars('T0');
    R_bmp = pars('R_bmp');
    R_h_p_n = pars('R_h_p_n');  
    
    %HIPOXIA    
    R_p_p_n = pars('R_pp');

    Theta0 = [R_e_p_0, R_s_p_0, R_rm_p_0, R_am_p_0, V_u_e_v_0, V_u_s_v_0, V_u_rm_v_0, V_u_am_v_0, Emax_lv_0, Emax_rv_0];
    
    %vars
    DTheta_R_e_p = y('DTheta_R_e_p');
    DTheta_R_s_p = y('DTheta_R_s_p');
    DTheta_R_rm_p_n = y('DTheta_R_rm_p_n');
    DTheta_R_am_p_n = y('DTheta_R_am_p_n');
    DTheta_V_unstressed_e_v = y('DTheta_V_unstressed_e_v');
    DTheta_V_unstressed_s_v = y('DTheta_V_unstressed_s_v');
    DTheta_V_unstressed_rm_v = y('DTheta_V_unstressed_rm_v');
    DTheta_V_unstressed_am_v = y('DTheta_V_unstressed_am_v');
    DTheta_Emax_lv = y('DTheta_Emax_lv');
    DTheta_Emax_rv = y('DTheta_Emax_rv');
    DTsym = y('DTsym');
    DTvagal = y('DTvagal');
    xO2_am  = y('xO2_am');
    x_met = y('x_met');
    xO2_h = y('xO2_h');
    xO2_rm = y('xO2_rm');
    xCO2_h = y('xCO2_h');
    xCO2_rm = y('xCO2_rm');
    xO2_b = y('xO2_b');
    xCO2_b = y('xCO2_b');
    xO2_e = y('xO2_e');
    xCO2_e = y('xCO2_e');
    xO2_s = y('xO2_s');
    xCO2_s = y('xCO2_s');
    xO2_p = y('xO2_p');
    xCO2_p = y('xCO2_p');
    aO2 = y('aO2');

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




    internal_variables('R_e_p_n') = R_e_p_n;
    internal_variables('R_s_p_n') = R_s_p_n;
    internal_variables('R_rm_p_n') = R_rm_p_n;
    internal_variables('R_am_p_n') = R_am_p_n;
    internal_variables('V_unstressed_e_v') = Theta(5);
    internal_variables('V_unstressed_s_v') = Theta(6);
    internal_variables('V_unstressed_rm_v') = Theta(7);
    internal_variables('V_unstressed_am_v') = Theta(8);
    internal_variables('Emax_lv') = Theta(9);
    internal_variables('Emax_rv') = Theta(10);
    internal_variables('Theart') = Theart;    
    internal_variables('R_am_p') = R_am_p;
    internal_variables('R_h_p') = Rp(1);
    internal_variables('R_rm_p') = Rp(2);
    internal_variables('R_b_p') = R_bp;
    %HIPOXIA
    internal_variables('R_e_p') = R_ep;
    internal_variables('R_s_p') = R_sp;
    internal_variables('R_p_p') = R_pp_;

    internal_variables_ = internal_variables;
end


function [all_global_, tiny_y_keys] = saving_in_globals(all_global, y)

     %variables we want to save in all_global and external_global
    %arr = zeros(1,15);

    tiny_y_keys = ["dVE", "PACO2", "PAO2", "Pmusc", "t0_heart", "u_t0", "HR", "phi_met", "fh_s", "fp_s", "fv_s", "fv"];

    %round_time = round(time/dt);
    %round_time_plus_1 = round_time + 1;
    %disp(round_time_plus_1);
    %arr(:,1) = y(13);
    %arr(:,2) = y(52);
    %arr(:,3) = y(51);
    %arr(:,4) = y(15);
    %arr(:,5) = y(93);
    %arr(:,6) = y(95);
    %arr(:,7) = y(92);
    %arr(:,8) = y(127);
    %arr(:,9) = y(128);
    %arr(:,10) = y(129);
    %arr(:,11) = y(130);
    %arr(:,12) = y(131);
%
%
    %  
%
    %    
    %
    %all_global(:,  round_time_plus_1) = arr';
    %prev_all_global_ = all_global(:,1:round(time/dt) + 1);
    all_global_ = all_global;
    %all_global_(:,1:round(time/dt) + 1) = prev_all_global_ + (prev_all_global_ == 0).* arr';
    
    
    %%externals_global_(:,round(time/dt) + 1) = all_global(:,round(time/dt) + 1);
    %all_global_(:,round(time/dt) + 1) = init_values;
%
    %%prev_delays_global = delays_global(:,1:round(time/dt) + 1);
    %%delays_global(:,1:round(time/dt) + 1) = prev_delays_global + (prev_delays_global == 0).* delays;
    %
    %prev_all_global_ = all_global_(:,1:round(time/dt) + 1);
    %all_global_(:,1:round(time/dt) + 1) = prev_all_global_ + (prev_all_global_ == 0).* init_values;

end

function debugging_tools(xdot, y, init_keys)

    if(true)
    end


    %Small tests to check inf, imaginary, negative, zero values on variables

    %len_xdot = length(xdot);
%
    %if sum(imag(xdot)) ~= 0
    %    disp('alto');
    %end

    %% testing (it wont be possible until being ridd of dde23)
    % test_keys = test.keys(); 
    
    % for i=1:length(test_keys)
        
    %     if t == 0
    %         y(test_keys(i)) = 0;
    %         xdot(len_xdot + i) = test(test_keys(i));
    %     else
    %         xdot(len_xdot + i) = test(test_keys(i));
    %     end 
    % end
    

    

    
    %%%%%%%%%%%% Debugging tools %%%%%%%%%%%%%%%%%%%%

    % big_detector = xdot > 6000;
    % if sum(big_detector) 
    %     disp('variables ugly')
    % 
    %     for i = 1:length(big_detector)
    %         if big_detector(i) == 1
    %             disp(init_keys(i))
    %             disp(y(init_keys(i)));
    % 
    %         end
    %     end
    % 
    %     %disp(y("R_lv"))
    % 
    % 
    % 
    %     disp('tiempo');
    %     disp(t);
    % 
    % end
    % if sum(isnan(y.values))
    %         oo = 9;
    % end



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


