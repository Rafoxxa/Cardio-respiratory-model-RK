function [t, x_dot, x_vars, X_KEYS, INDEX] = run_ode(model, pars, init, simulation_time, dt)
    
    
    global all_global;
    all_global = zeros(15, round(10 * simulation_time/dt) + 1) + 0;
    global iter_step_prev;
    iter_step_prev = 0;
    control_on = 1;
    
    % Loop invariant variables
    %disp("simulation_time");
    %disp(simulation_time);
    if isa(init, "dictionary")
        X_KEYS = init.keys;
        INDEX = dictionary(X_KEYS, [1:length(X_KEYS)]');    
    else
        X_KEYS = 0;
        INDEX = 0;
    end

    %TINY_X_KEYS = ["dVE", "PACO2", "PAO2", "Pmusc", "insp_integrand", "exp_integrand", "J", "t0_heart", "u_t0", "HR", "phi_met", "fh_s", "fp_s", "fv_s", "fv"]; 
    %TINY_INDEX = dictionary(TINY_X_KEYS, [1:length(TINY_X_KEYS)]);
    %OPTIONS = ddeset('AbsTol',1e-8,'RelTol',1e-5); 
    if isa(pars, "dictionary")
        OPTIONS_CENTRAL = odeset('AbsTol',1e-3,'RelTol',1e-2, 'OutputFcn', @(t, y, flag, varargin) simple_output_fcn(t, y, flag, pars.values, varargin) );
    else
        OPTIONS_CENTRAL = odeset('AbsTol',1e-3,'RelTol',1e-2, 'OutputFcn', @(t, y, flag, varargin) simple_output_fcn(t, y, flag, pars, varargin) );
    end
    OPTIONS = odeset('AbsTol',1e-3,'RelTol',1e-2);
    SIMULATION_TIME = simulation_time;
    DT = dt;
    
    
if isa(pars, 'dictionary')
pars = pars.values;
end
if isa(init, 'dictionary')
init = init.values;
end
PARS = pars;   
    
    
    % cycle periods  
    Tresp = init(5);

    % define cycle-elements
    cycle = struct();
    cycle.t_ini = 0;
    cycle.t_end = Tresp;
    cycle.time_interval = cycle.t_ini:DT:cycle.t_end;
    cycle.init_vars = init;
    %cycle.init_taus = taus;
    cycle.x_vars = [];
    cycle.progress = 0;
    cycle.progress_dec = 0;
    cycle.real_time = 0;
    
    %define stacked/historical-elements
    stacked = struct();
    stacked.t = [];
    stacked.vars = [];
    stacked.t0 = [];
 

    while cycle.t_end < SIMULATION_TIME            
                tic;
                sol = ode23(model, [cycle.t_ini, cycle.t_end], cycle.init_vars, OPTIONS_CENTRAL, PARS, X_KEYS);
                cycle.x_vars = deval(sol, cycle.time_interval);         
        
        [PARS, stacked, cycle] = cycle_driver(stacked, cycle, DT, INDEX, PARS);
        
        cycle.real_time = cycle.real_time +  toc;
        %disp('tiempo real:');
        %disp(cycle.real_time);
        

        
        
    end

    x_vars = stacked.vars;
    t = stacked.t;
    
    stacked.t0 = [stacked.t0 t(end)]; 
    t0 = stacked.t0(1);

    x_dot = zeros(size(x_vars));  
    cnt = 1;      
      
    for i = 1:length(t)
        if t(i) >= stacked.t0(cnt+1);
            cnt = cnt + 1;
            t0 = stacked.t0(cnt);
        end
        PARS(77) = t0;
        x_dot(:,i) = model(t(i), x_vars(:,i), PARS, X_KEYS);
        %x_dot(:,i) = model(t(i), x_vars(:,i), delays_global(:, :, round(t(i)/DT) + 1), PARS, init.keys, taus.keys);
        
    end    
    
    
    
    function [PARS_, stacked, cycle] = cycle_driver(stacked, cycle, DT, INDEX, PARS)       
        
        
        Tresp = cycle.x_vars(5, end); 
        dVE = cycle.x_vars(13, end);
        stacked.t = [stacked.t cycle.time_interval(1:end-1)];
        stacked.vars = [stacked.vars cycle.x_vars(:,1:end-1)];
        
        stacked.t0 = [stacked.t0 cycle.t_ini];   
        cycle.t_ini = cycle.t_end;
        PARS(77) = cycle.t_ini;
        
        cycle.init_vars = cycle.x_vars(:,end);
        %cycle.init_vars_values = [cycle.x_vars(:, end)];% ; squeeze(delays_global(:,:,end))];        
        %cycle.init_taus = cycle.init_taus;  %this is unnecesary, because tau is not variable
        
        PARS_ = PARS;

        if control_on
        
            %disp('dVE');
            %disp(dVE);            
            args_optimal = control_optimization(cycle, OPTIONS, PARS_, dVE);
            %disp('fin optim')
            %disp('args_optimal')
            %disp(args_optimal);

            cycle.init_vars(4) = args_optimal(1);
            cycle.init_vars(6) = args_optimal(2);
            %cycle.init_vars(19) = args_optimal(3);
            cycle.init_vars(20) = args_optimal(3);
            cycle.init_vars(21) = args_optimal(4);
            cycle.init_vars(24) = args_optimal(5);
        end
        cycle.init_vars(5) = cycle.init_vars(4) + cycle.init_vars(6); 
        Tresp = cycle.init_vars(5);
        cycle.t_end = cycle.t_end + Tresp;
        cycle.time_interval = cycle.t_ini: DT: cycle.t_end;
        %cycle.progress = PARS(77)/SIMULATION_TIME * 100;
        %if cycle.progress >= cycle.progress_dec 
        %    fprintf('The percentage is: %.2f%%\n', cycle.progress);
        %    cycle.progress_dec = cycle.progress_dec + 10;
        %end

        
        
        
    end


    function J = cost_fun(cycle, OPTIONS, PARS, args)
        
        %DT = pars(2);
        %optim_init = cycle.init_vars;
        optim_pars = PARS;
        optim_pars(352) = args(1);
        optim_pars(353) = args(2);
        optim_pars(354) = 0;  %fixed by authors
        optim_pars(355) = args(3);
        optim_pars(356) = args(4);
        optim_pars(357) = args(5);
        
        optim_pars(358) = optim_pars(352) + optim_pars(353);
        state_var = [cycle.init_vars(12), cycle.init_vars(10), cycle.init_vars(11), cycle.init_vars(16)];

        J = respiratory_work(state_var, OPTIONS, optim_pars);
        J = log(J);
        
        % Cost function equations

        function ydot = tiny_system_ua(t, init_value, pars)
            %this should be the corect equation, that diminshes the speed of the ode incrporating upper airways dynamics, giving a lowe volume for the same muscular presion.

            %parameters
            DT = pars(3);
            Pao = pars(18);
            Ecw = pars(31);
            El = pars(32);
            Rrs = pars(28);
            Ers = Ecw + El;          
            
            kaw1  = pars(38);
            kaw2 = pars(39);
            Rcw = pars(25);           

            R_trachea = pars(29);    
            Rl = pars(26);
            Rcw = pars(25);  
            Raw = pars(27);
            Cua = pars(30);
            bua = pars(35);
            Pcrit_min = pars(15);
            A0ua = pars(36);
            Kua = pars(37);
            
            %define state vars dict           
            
            %obtain state vars
            V = init_value(1);
            dV_prev = init_value(2);
            dVua = init_value(3);
            Pua_prev = init_value(4);

            a0 = 0;
            a1 = pars(355);
            a2 = pars(356);
            tau = pars(357);
            TI = pars(352);

            %Pmusc function
            if 0 <= t && t <= TI
                Pmusc = a0 + a1*t + a2*t^2;        
        
            elseif  TI < t 
                PmuscTI = a0 + a1*TI + a2*TI^2;
                Pmusc = PmuscTI * exp(-(t - TI)/tau);    %cmH2O            
            end

            %Pleural pressure function
            if dV_prev < 0
                Pcw = Ecw * V - 1;
                Pa_ = Pao;                               
            else
                Pcw = Ecw * V - 1 + Rcw * dV_prev;
                Pa_ = Pao - kaw1 * dV_prev - kaw2 * abs(dV_prev)^2;
            end
            Pa = Pa_ * (Pa_ > 0);
            Ppl = Pcw + Pa - Pmusc; 

            %Upper airways computation
            Rrs = Raw + Rl + Rcw;
            dVla = dVua + dV_prev;
            Pua = Ppl + dVla * Rrs;        
            dPua = (Pua - Pua_prev)/DT;
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

            %ODE for volume wave
            dV =  Gaw/Rrs * ((Pmusc - Pao) - Ers * V);
            %dV =  1/Rrs * ((Pmusc - Pao) - Ers * V);
            ddV = (dV - dV_prev)/DT;

            ydot = [dV, ddV, ddVua, dPua]';
            
        end
        
        
        function dV = tiny_system(t, init_value, pars)
            % For the moment, based on the graphs of Gaw, we are taking an average value, trying to not compromise copmutational speed.
            

            %parameters
            Pao = pars(18);
            Ecw = pars(31);
            El = pars(32);
            Rrs = pars(28);
            Ers = Ecw + El;
            
            %define state vars dict
            
            
            %obtain state vars
            V = init_value;
            a0 = 0;
            a1 = pars(355);
            a2 = pars(356);
            tau = pars(357);
            TI = pars(352);

            %pressure function
            if 0 <= t && t <= TI
                Pmusc = a0 + a1*t + a2*t^2;        
        
            elseif  TI < t 
                PmuscTI = a0 + a1*TI + a2*TI^2;
                Pmusc = PmuscTI * exp(-(t - TI)/tau);    %cmH2O            
            end

            %ode for volume wave
            dV =  1/Rrs * ((Pmusc - Pao) - Ers * V);
            
        end

        function J = respiratory_work(state_vars, OPTIONS, optim_pars)
            
            %parameters for cost function
            Pmax = optim_pars(75);
            dPmax = optim_pars(76);
            lambda1 = optim_pars(73);
            lambda2 = optim_pars(74);
            n = optim_pars(72);
            DT = optim_pars(3);

            TI = optim_pars(352);
            TE  = optim_pars(353);
            a0 = optim_pars(354);
            a1 = optim_pars(355);
            a2  = optim_pars(356);
            tau = optim_pars(357);

            TI_check = var_inside_boundries(TI, optim_pars(347), optim_pars(342));
            TE_check = var_inside_boundries(TE, optim_pars(348), optim_pars(343)); 
            a1_check = var_inside_boundries(a1, optim_pars(349), optim_pars(344));
            a2_check = var_inside_boundries(a2, optim_pars(350), optim_pars(345));
            tau_check = var_inside_boundries(tau, optim_pars(351), optim_pars(346));
            vars_inside_boundries = TI_check && TE_check && a1_check && a2_check && tau_check;
            if vars_inside_boundries
                Tresp = TI + TE;
                time_interval = 0: DT: Tresp;

                %run the volume wave ode, given that this tiny system doesn't require previous time of the simulation, we don't need information of the cycle.
                tiny_model = @(t, init) tiny_system_ua(t, init, optim_pars); 
                %sol = ode23(tiny_model, [0, optim_init(5)], optim_init.values, OPTIONS, PARS, optim_init.keys);
                try
                    sol = ode23(tiny_model, [0, Tresp], state_vars, OPTIONS);
                    y = deval(sol, time_interval);  %the only variable is volume
                    %V = y(1);
                    dV_dt = y(2,:);


                    %obtain time vector
                    t = time_interval; 
                    %recompute pressure wave             
                    Pmusc = (t <= TI).* (a0 + a1*t + a2*t.^2) + (TI < t).*((a0 + a1*TI + a2*TI^2) * exp(-(t - TI)/tau));   
                
                    %obtain finite difference aproximation for Pmusc
                    dPmusc_dt = diff(Pmusc) ./ DT;
                    
                    %obtain finite difference aproximaiton for ddV/ddt
                    %dV_dt = diff(V) ./ diff(t);
                    ddV_ddt = diff(dV_dt) ./ DT;        

                    %establish same sizes for all arrays:
                    %dPmusc_dt = dPmusc_dt(1:end-1);     
                    Pmusc = Pmusc(1:end-1);     
                    t = t(1:end-1);
                    
                    %compute zheta 1 and zheta 2 parameters for optimization
                    zheta1 = (1 - Pmusc/Pmax) .* (Pmax > Pmusc) + 0.01 * (Pmax <= Pmusc); %this is for handling off limits cases
                    zheta2 = (1 - dPmusc_dt/dPmax) .* (dPmax > dPmusc_dt) + 0.01 * (dPmax <= dPmusc_dt); %this is for handling off limits cases

                    % obtain the integrals and the respiratory work: all vectorial
                    insp_integrand = Pmusc./(zheta1.^n .* zheta2.^n) + lambda1 * ddV_ddt.^2;
                    exp_integrand = ddV_ddt.^2;
                    insp_work_power = 1/Tresp * DT * sum(insp_integrand .* (t <= TI)); 
                    exp_work_power =  1/Tresp * DT * sum(exp_integrand .* (t > TI)); %it should be until Tresp, but that happens at the end of the cycle
                    J = insp_work_power + lambda2 * exp_work_power; 
                catch
                    J = 10^10;
                end 
            else
                J = 10^10;
            end          
        end



    end

    function args_optimal = control_optimization(cycle, OPTIONS, PARS, dVE)
        % Define initial values for optimization parameters
        initial_TI_value = cycle.init_vars(4)  ;
        initial_TE_value = cycle.init_vars(6);
        initial_a1_value = cycle.init_vars(20) ;
        initial_a2_value = cycle.init_vars(21)  ;
        initial_tau_value = cycle.init_vars(24) ;
        args_init = [initial_TI_value, initial_TE_value, initial_a1_value, initial_a2_value, initial_tau_value];
    
        % Set up optimization options
        %options = optimoptions("fmincon", "Display", "iter", "Algorithm", "sqp");
        options = optimset('fmincon');
        options.Algorithm = 'sqp';
        options.Display = 'off';
        options.Diagnostics = 'off';
        %%options.TolX = 10^-3;
        %%options.TolFun = 10^-3; 
        options.MaxIter = 3;
    
        % Set lower and upper bounds for the optimization parameters
        lb = [pars(342), pars(343), pars(344), pars(345), pars(346)];
        ub = [pars(347), pars(348), pars(349), pars(350), pars(351)];
    
        %Contraints
        % Tresp > TI
        % A = [1 -1 0 0 0 0];
        % b = 0;

        % Perform optimization
        nonlcon = @(x) ventilation_constraint(x, dVE, PARS, OPTIONS);
        args_optimal = fmincon(@(args) cost_fun(cycle, OPTIONS, PARS, args), args_init, [], [], [], [], lb, ub, nonlcon, options);

        %Non Linear Constraint

        function ydot = fast_tiny_system(t, init_value, pars, a1, a2)
            %parameters
            %Pao = pars(18);
            %Ecw = pars(31);
            %El = pars(32);
            %Rrs = pars(28);
            %Ers = Ecw + El;
            %%obtain state vars
            %V = init_value;
            %a0 = 0;           
            %%pressure function
            %Pmusc = a0 + a1*t + a2*t^2; 
            %%ode for volume wave
            %dV =  Gaw/Rrs * ((Pmusc - Pao) - Ers * V);  %I think Gaw gots to be well computed instead of putting this static value o 0.3
            DT = PARS(3);
            Pao = pars(18);
            Ecw = pars(31);
            El = pars(32);
            Ers = Ecw + El;
            kaw1  = pars(38);
            kaw2 = pars(39);
            R_trachea = pars(29);    
            Rl = pars(26);
            Rcw = pars(25);  
            Raw = pars(27);
            Cua = pars(30);
            Pcrit_min = pars(15);
            A0ua = pars(36);
            Kua = pars(37);
            V = init_value(1);
            dV_prev = init_value(2);
            dVua = init_value(3);
            Pua_prev = init_value(4);

            a0 = 0;
            
            %Pmusc function
            Pmusc = a0 + a1*t + a2*t^2;        
            %Pleural pressure function
            Pcw = Ecw * V - 1 + Rcw * dV_prev;
            Pa_ = Pao - kaw1 * dV_prev - kaw2 * abs(dV_prev)^2;
            Pa = Pa_ * (Pa_ > 0);
            Ppl = Pcw + Pa - Pmusc; 
            %Upper airways computation
            Rrs = Raw + Rl + Rcw;
            dVla = dVua + dV_prev;
            Pua = Ppl + dVla * Rrs;        
            dPua = (Pua - Pua_prev)/dt;
            ddVua = -1/R_trachea * (dPua + dVua/Cua);    
            Pcrit = Pcrit_min;
            Gaw = A0ua * (1 - Pua/Pcrit) * Kua;
            %ODE for volume wave
            dV =  Gaw/Rrs * ((Pmusc - Pao) - Ers * V);
            
            ddV = (dV - dV_prev)/dt;
            ydot = [dV, ddV, ddVua, dPua]';
        end

        function [c, ceq] = ventilation_constraint(args, dVE, PARS, OPTIONS)
            DT = PARS(3);
            Pao = PARS(18);
            Ecw = PARS(31);
            El = PARS(32);            
            Ers = Ecw + El;

            TI = args(1);
            TE = args(2);
            a1 = args(3);
            a2 = args(4);

            TI_check = var_inside_boundries(TI, PARS(347), PARS(342));
            TE_check = var_inside_boundries(TE, PARS(348), PARS(343));
            a1_check = var_inside_boundries(a1, PARS(349), PARS(344));
            a2_check = var_inside_boundries(a2, PARS(350), PARS(345));

            vars_inside_boundries = TI_check && TE_check && a1_check && a2_check;
            if vars_inside_boundries

                fast_tiny_model = @(t, init) fast_tiny_system(t, init, PARS, a1, a2); 
                time_interval = 0:DT:TI;
                sol = ode23(fast_tiny_model, [0, TI], [0,0,0,0], OPTIONS);
                try
                    y = deval(sol, time_interval);
                    %disp(V(end));
                    c = [];
                    VT = y(1,end);  %Volume at TIs
                catch 
                    VT = max(sol.y(1,:));
                    c = [];
                    
                    
                end
                %VT = 0.9*(a1*TI + a2*TI^2 - Pao)/Ers; %we are multiplicating by a proxy of Gaw to avoid running the whole ode
                %ceq = VT - 0.1;
                %ceq = VT/(TE + TI) - dVE; %dVE is already in lt/s
                
                ceq = VT/(TE + TI) - dVE;
                
                %ceq = (TE + TI)*dVE * Ers + Pao - a1 * TI + a2 * TI^2; 
            else
                c = [];
                ceq = 10^10;
            end

        end
     
    end


    function status = simple_output_fcn(t, y, flag, pars, varargin)
        
    
        % Initialize the persistent variable during initialization
        init = y;
        if strcmp(flag, 'init')
            dt = pars(2); % Set the time step         
            
        elseif isempty(flag)
            % Regular call for each solver step
            round_time = round(t/dt);
            iter_step = round_time;
            
            delta_step = iter_step - iter_step_prev;
            iter_step_plus_1 = iter_step + 1;
            iter_step_prev_plus_1 = iter_step_prev + 1;
            
            if delta_step > 0
                all_global(1,iter_step_prev_plus_1:iter_step_plus_1) = init(13);
                all_global(2,iter_step_prev_plus_1:iter_step_plus_1) = init(52);
                all_global(3,iter_step_prev_plus_1:iter_step_plus_1) = init(51);
                all_global(4,iter_step_prev_plus_1:iter_step_plus_1) = init(15);
                all_global(5,iter_step_prev_plus_1:iter_step_plus_1) = init(93);
                all_global(6,iter_step_prev_plus_1:iter_step_plus_1) = init(95);
                all_global(7,iter_step_prev_plus_1:iter_step_plus_1) = init(92);
                all_global(8,iter_step_prev_plus_1:iter_step_plus_1) = init(127);
                all_global(9,iter_step_prev_plus_1:iter_step_plus_1) = init(128);
                all_global(10,iter_step_prev_plus_1:iter_step_plus_1) = init(129);
                all_global(11,iter_step_prev_plus_1:iter_step_plus_1) = init(130);
                all_global(12,iter_step_prev_plus_1:iter_step_plus_1) = init(131);
            
            elseif delta_step <= 0

                all_global(1,iter_step_plus_1) = init(13);
                all_global(2,iter_step_plus_1) = init(52);
                all_global(3,iter_step_plus_1) = init(51);
                all_global(4,iter_step_plus_1) = init(15);
                all_global(5,iter_step_plus_1) = init(93);
                all_global(6,iter_step_plus_1) = init(95);
                all_global(7,iter_step_plus_1) = init(92);
                all_global(8,iter_step_plus_1) = init(127);
                all_global(9,iter_step_plus_1) = init(128);
                all_global(10,iter_step_plus_1) = init(129);
                all_global(11,iter_step_plus_1) = init(130);
                all_global(12,iter_step_plus_1) = init(131);
            
            end

            
            
            

            iter_step_prev = iter_step;
            
       
        end
    
        status = 0; % Continue the solver
    end

    function logical_output = var_inside_boundries(value, upper, lower)
        d = sign(value) * 0.001;
        logical_output = (value > upper * (1 + d)) || (value < lower * (1 - d));         
        logical_output = ~logical_output;
    end

    

 
end