function  plot_porcentual()
long_transient_data = load('long_transient.mat');
sv = long_transient_data.struct_vars;
t = long_transient_data.t;
%plot(t, s_v.MRtCO2);



PACO2 = sv.PACO2;
PAO2 = sv.PAO2;
V = sv.V;
dVE = sv.dVE;
TE = sv.TE;
TI = sv.TI;
P_sa = sv.P_sa;
HR = sv.HR;

PACO2_percentual = take_percentage(PACO2, t);
PAO2_percentual = take_percentage(PAO2, t);
V_percentual = take_percentage(V, t);
dVE_percentual = take_percentage(dVE, t);
TE_percentual = take_percentage(TE, t);
TI_percentual = take_percentage(TI, t);
P_sa_percentual = take_percentage(P_sa, t);
HR_percentual = take_percentage(HR, t);

plot(t, PACO2_percentual);
title('PACO2');
figure;
plot(t, PAO2_percentual);
title('PAO2');
figure;
% plot(t, V_percentual);
% plot(t, dVE_percentual);
% plot(t, TE_percentual);
% plot(t, TI_percentual);
% plot(t, P_sa_percentual);
% plot(t, HR_percentual);


function var_percentual = take_percentage(var, t)
var_ = var;
var_400 = sum(var_.*(t > 390).*(t < 410))./sum((t > 390).*(t < 410));
%%disp(PAO2_400);
var_percentual = 100 * (var_ - var_400)/var_400;

end
end