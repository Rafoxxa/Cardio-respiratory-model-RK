% Define the parameters
Pao = 0;
Rrs = 3.02;
Ers = 20;
a0 = 0;
a1 = 11;
a2 = -2;
tau = 0.2;
TI = 1.2;
TE = 1.3;

% Define symbolic variables
syms t V(t)

% Define the piecewise function
Pmusc1 = a0 + a1*t + a2*t^2;
Pmusc2 = (a0 + a1*TI + a2*TI^2) * exp(-(t - TI)/tau);

% Define the ODE for the first interval: 0 <= t <= TI
dV1 = diff(V, t);
ode1 = dV1 == 1/Rrs * ((Pmusc1 - Pao) - (Ers) * V(t));

% Solve for the first interval with initial condition V(0) == 0
sol1 = dsolve(ode1, V(0) == 0);

% Define the ODE for the second interval: TI < t
dV2 = diff(V, t);
ode2 = dV2 == 1/Rrs * ((Pmusc2 - Pao) - (Ers) * V(t));

% Solve for the second interval, using sol1 evaluated at t=TI as initial condition
V_TI = subs(sol1, t, TI);
sol2 = dsolve(ode2, V(TI) == V_TI);

% Convert the solutions to numerical functions for plotting
fV1 = matlabFunction(subs(sol1), 'Vars', t);
fV2 = matlabFunction(subs(sol2), 'Vars', t);

% Define time intervals
t1 = linspace(0, TI, 100);
t2 = linspace(TI, TE + TI, 100);

% Evaluate the solutions
V1 = fV1(t1);
V2 = fV2(t2);

% Plot the solutions
figure;
hold on;
plot(t1, V1, 'b-', 'LineWidth', 2); % Plot for 0 <= t <= TI
plot(t2, V2, 'r-', 'LineWidth', 2); % Plot for TI < t
xlabel('Time (s)');
ylabel('Volume (V)');
title('Solution of the ODE with Piecewise Input');
legend('0 <= t <= TI', 'TI < t');
grid on;
hold off;


