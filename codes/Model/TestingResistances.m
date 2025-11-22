%Testing Resistances
sv = struct_vars;

Rbp = 6.57./(1 + sv.xO2_b + sv.xCO2_b);
Ramp = (sv.DTheta_R_am_p_n + 4.36)./(1 + sv.xO2_am + sv.x_met);
Rrmp = (sv.DTheta_R_rm_p_n + 4.36).* (1 + sv.xCO2_rm) ./(1 + sv.xO2_rm);
Rhp = 22.27 * (1 + sv.xCO2_h)./(1 + sv.xO2_h);
Rep = (sv.DTheta_R_e_p + 2.87)./(1 + sv.xO2_e);
Rsp = (sv.DTheta_R_s_p + 2.31)./(1 + sv.xO2_s);
%t = sv.MRtO2;
figure; plot(t, Rbp); hold on;plot(t, Ramp); hold on;plot(t, Rrmp); hold on;plot(t, Rhp); hold on; plot(t, Rep); hold on; plot(t, Rsp); legend('b','am', 'rm', 'h', 'e', 's'); title('Resistances');

figure; plot(t, (1./Rbp + 1./Ramp + 1./Rrmp + 1./Rhp + 1./Rep + 1./Rsp).^-1); title('Total Parallel Resistance');
figure; hold on; plot(t, sv.DTheta_R_am_p_n + 4.36); title('Ram only simpathetic');

figure; hold on; plot(t, sv.xO2_am); title('xO2am');
hold on; yyaxis right; plot(t, sv.MRtCO2); legend('xO2am', 'CO2') ;title('MRtCO2 and xO2am');
hold on; yyaxis right; plot(t, sv.MRtO2); legend('xO2am', 'O2') ;title('MRtO2 and xO2am');

%Watching R_p_p
R_pp = s.pars('R_p_p_n') * 1./(1 + sv.xO2_p);
figure; plot(t, R_pp); title('R_pp');

%Plotting I and phimet variable

I = (sv.MRtCO2 - s.pars('MRtCO2_basal'))./(s.pars('AT') - s.pars('MRtCO2_basal'));
phimet = (s.pars('phi_min') + s.pars('phi_max') * exp((I - s.pars('I0_met'))./(s.pars('kmet'))))./(1 + exp((I - s.pars('I0_met'))./(s.pars('kmet'))));

figure; hold on; plot(t, sv.xO2_am); title('xO2am');
hold on; yyaxis right; plot(t, I); legend('xO2am', 'I') ;title('I and xO2am');




figure; hold on; plot(t, sv.x_met); title('xmet');
figure; hold on; plot(t, 1./(1 + sv.xO2_am + sv.x_met)); title('local effect on am');
figure; hold on; plot(t, 1./(1 + sv.xO2_am + 0.27)); title('local effect with xmet set to 0.27');

%Volumes
figure; plot(t, sv.V_total_b_v); hold on; plot(t, sv.V_total_am_v); hold on;plot(t, sv.V_total_rm_v); hold on;plot(t, sv.V_total_h_v); hold on;plot(t, sv.V_total_e_v); hold on; plot(t, sv.V_total_s_v); hold on; plot(t, Rsp); legend('b','am', 'rm', 'h', 'e', 's'); title('Venous Volumes');

%Cardiac Elastances
figure; plot(t, sv.DTheta_Emax_lv + s.pars('Emax_lv_0')); hold on; plot(t, sv.DTheta_Emax_rv + s.pars('Emax_rv_0')); legend('lv', 'rv'); title('lv and rv elastances');

