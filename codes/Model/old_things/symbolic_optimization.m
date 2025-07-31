% Define symbolic variables
syms t V(t) dV_prev(t) dVua(t) Pua_prev(t) Pao Ecw El Rrs kaw1 kaw2 Rcw
syms R_trachea Rl Raw Cua bua Pcrit_min A0ua Kua
syms a0 a1 a2 tau TI TE

% Define state variables as symbolic
state_vars = [V; dV_prev; dVua; Pua_prev];

% Define symbolic parameters
params = [Pao, Ecw, El, Rrs, kaw1, kaw2, Rcw, ...
          R_trachea, Rl, Raw, Cua, bua, Pcrit_min, A0ua, Kua, ...
          a0, a1, a2, tau, TI];

% Define Pmusc as a symbolic expression
Pmusc = piecewise(0 <= t & t <= TI, a0 + a1*t + a2*t^2, ...
                  TI < t, (a0 + a1*TI + a2*TI^2) * exp(-(t - TI)/tau));

% Define Pcw and Pa_ as symbolic expressions
Pcw = piecewise(dV_prev(t) < 0, Ecw * V(t) - 1, ...
                dV_prev(t) >= 0, Ecw * V(t) - 1 + Rcw * dV_prev(t));

Pa_ = piecewise(dV_prev(t) < 0, Pao, ...
                dV_prev(t) >= 0, Pao - kaw1 * dV_prev(t) - kaw2 * abs(dV_prev(t))^2);

Pa = Pa_ * heaviside(Pa_);
Tresp = TI + TE;

% Define Ppl, Pua, and other symbolic expressions
Ppl = Pcw + Pa - Pmusc;
Rrs = Raw + Rl + Rcw;
dVla = dVua(t) + dV_prev(t);
Pua = Ppl + dVla * Rrs;
dPua = diff(Pua, t);
ddVua = -1/R_trachea * (dPua + dVua(t)/Cua);

% Define Pcrit and Gaw as symbolic expressions
Pcrit = piecewise(10 <= -1/(Cua * bua), 10, ...
                  Pcrit_min < -1/(Cua * bua) & -1/(Cua * bua) < 10, -1/(Cua * bua), ...
                  Pcrit_min >= -1/(Cua * bua), Pcrit_min);

Gaw = piecewise(Pua <= Pcrit, 0, ...
                Pcrit < Pua & Pua <= 0, A0ua * (1 - Pua/Pcrit) * Kua, ...
                Pua > 0, A0ua * Kua);

% Define ODE for volume wave as a symbolic expression
dV = Gaw/Rrs * ((Pmusc - Pao) - (Ecw + El) * V(t));
ddV = diff(dV, t);

% Define cost function J as a symbolic expression
syms Pmax dPmax lambda1 lambda2 n DT
Pmusc_sym = piecewise(t <= TI, a0 + a1*t + a2*t^2, ...
                      TI < t, (a0 + a1*TI + a2*TI^2) * exp(-(t - TI)/tau));
dPmusc_dt_sym = diff(Pmusc_sym, t);

% Define zheta1 and zheta2
zheta1 = (1 - Pmusc_sym/Pmax) * heaviside(Pmax - Pmusc_sym) + 0.01 * heaviside(Pmusc_sym - Pmax);
zheta2 = (1 - dPmusc_dt_sym/dPmax) * heaviside(dPmax - dPmusc_dt_sym) + 0.01 * heaviside(dPmusc_dt_sym - dPmax);

% Define integrals for inspiratory and expiratory work power
insp_integrand = Pmusc_sym./(zheta1.^n .* zheta2.^n) + lambda1 * ddV.^2;
exp_integrand = ddV.^2;

% Define J using integrals
J = (1/(TI + TE)) * int(insp_integrand, t, 0, TI) + ...
    lambda2 * (1/(TI + TE)) * int(exp_integrand, t, TI, Tresp);

% Compute symbolic gradients
grad_J = gradient(J, state_vars);

