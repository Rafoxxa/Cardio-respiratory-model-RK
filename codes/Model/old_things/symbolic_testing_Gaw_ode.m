% Define the parameters
Pao = 0;
Rrs = 3.02;
Ers = 20;
Ecw = 10.5;
a0 = 0;
a1 = 11;
a2 = -2;
tau = 0.2;
TI = 1.2;
TE = 1.3;
Pcrit = -40;
kaw1  = 1.85;
kaw2 = 0.43;
R_trachea = 10^6;    
Rl = 1.36;
Rcw = 0.83;  
Raw = 0.82;
Cua = 0.001;
A0ua = 1;
Kua = 1;

% Define symbolic variables
syms t V(t) dVua(t)

dV1 = diff(V, t);
dV2 = diff(V, t);

ddVua1 = diff(dVua, t);
ddVua2 = diff(dVua, t);
%%PIECEWISES
    % Define Pmusc
    Pmusc1 = a0 + a1*t + a2*t^2;
    Pmusc2 = (a0 + a1*TI + a2*TI^2) * exp(-(t - TI)/tau);

%Define Pleural 
    Pcw1 = Ecw * V - 1 + Rcw * dV1;
    Pa1 = Pao - kaw1 * dV1 - kaw2 * dV1^2;    
    Pcw2 = Ecw * dV2 - 1;
    Pa2 = Pao;   
    
    %Pa1 = max(Pa1_, 0);
    Ppl1 = Pcw1 + Pa1 - Pmusc1;
    Ppl2 = Pcw2 + Pa2 - Pmusc2; 

    Rrs = Raw + Rl + Rcw;
    Pua1 = Ppl1 + (dVua + dV1) * Rrs;        
    dPua1 = diff(Pua1, t);    
    odeVua1 = ddVua1 == -1/R_trachea * (dPua1 + dVua/Cua);    
    Gaw1 = A0ua * (1 - Pua1/Pcrit) * Kua;

    Pua2 = Ppl2 + (dVua + dV2) * Rrs;        
    dPua2 = diff(Pua2, t);    
    odeVua2 = ddVua2 == -1/R_trachea * (dPua2 + dVua/Cua);   
    Gaw2 = A0ua * (1 - Pua2/Pcrit) * Kua;

% Define the ODE for the first interval: 0 <= t <= TI
conds1 = [V(0) == 0, dVua(0) == 0];
odeV1 = dV1 == Gaw1/Rrs * ((Pmusc1 - Pao) - (Ers) * V(t));

% Solve for the first interval with initial condition V(0) == 0
sols1 = dsolve([odeV1, odeVua1], conds1);

% Define the ODE for the second interval: TI < t

ode2 = dV2 == Gaw2/Rrs * ((Pmusc2 - Pao) - (Ers) * V(t));

% Solve for the second interval, using sol1 evaluated at t=TI as initial condition
%solV = sols1.V;
%solVua = sols1.dVua;
%V_TI = subs(solV, t, TI);
%dVua_TI = subs(solVua, t, TI);
%conds2 = [V(TI) == V_TI, dVua(TI) == dVua_TI];
%sol2 = dsolve([odeV2, odeVua2], conds2);

% % Convert the solutions to numerical functions for plotting
% fV1 = matlabFunction(subs(sol1), 'Vars', t);
% fV2 = matlabFunction(subs(sol2), 'Vars', t);
% 
% % Define time intervals
% t1 = linspace(0, TI, 100);
% t2 = linspace(TI, TE + TI, 100);
% 
% % Evaluate the solutions
% V1 = fV1(t1);
% V2 = fV2(t2);
% 
% % Plot the solutions
% figure;
% hold on;
% plot(t1, V1, 'b-', 'LineWidth', 2); % Plot for 0 <= t <= TI
% plot(t2, V2, 'r-', 'LineWidth', 2); % Plot for TI < t
% xlabel('Time (s)');
% ylabel('Volume (V)');
% title('Solution of the ODE with Piecewise Input');
% legend('0 <= t <= TI', 'TI < t');
% grid on;
% hold off;
