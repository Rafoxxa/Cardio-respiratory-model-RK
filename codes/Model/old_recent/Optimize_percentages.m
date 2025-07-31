[optimized_params] = optimize_parameters();
%save_excel(optimized_params);
function [optimized_params] = optimize_parameters()
    % Load original parameter values from file
    param_data = readtable('excel_para_ori.xlsx');
    
    %original_params = param_data{:, 2};
    
    % Define percentages that need to be optimized (from volumes and pressures)
    
    
    
    % Lower and upper bounds (allow variation between 50% and 150%)
    lb = [0.8 1 0.25 0.6 0.7 0.9      0.15 0.06 0.05 0.13 0.1  0.12 0.03           0.7  0.45 0.5 0.75 0.7 0.65        0.65 0.2 0.5 0.75 0.6 0.4       0.05 0.1 0.05 0.05 0.05 0.1     0.15 0.04 1 0.25 0.2  0.12 1];
    ub = [0.9 1 0.3  0.7 0.8 1        0.3  0.09 0.12  0.3  0.3  0.22  0.04         0.85  0.75 0.7 0.9  0.8  0.7       0.7  0.3 0.6 0.9 0.9 0.5        0.5  0.5 0.5  0.5   0.5 0.5    0.2  0.05 1 0.3  0.25 0.15 1];
    
    %ori, pressure perc, second chunk (x(7:))0.35 0.07 0.05 0.25 0.1  0.12 0.03
    %0.35 0.1  0.1  0.3  0.12 0.15 0.04

    %ori, venous perc, third chunk0.7  0.45 0.6 0.85 0.75 0.65
%                                      0.8  0.65 0.7 0.9  0.8  0.7

% ori, unstressed venous, fourth 0.65 0.2 0.5 0.75 0.6 0.4
%                                0.7  0.3 0.6 0.9 0.9 0.5 

%ori, unestressde perif, fifht  0.05 0.1 0.05 0.05 0.05 0.1
%                               0.15 0.2 0.1  0.1  0.1  0.15
    %
%    P_m_p
%P_h_p
%P_p_p
%P_s_p
%P_e_p
%P_b_p
    x0 = (lb + ub)/2;
    

    % Define Aeq and beq for the sum constraint
    n = length(lb);
    Aeq = zeros(1, n);
    Aeq([7 8 9 10 11 12 13]) = 1; % Set 1s at selected variable positions
    beq = 1; % The sum should equal 1    

    % Run optimization
    %options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    %optimized_params = fmincon(@(x) costFunction(x, param_data), x0, [], [], Aeq, beq, lb, ub, [], options);



       % Set options for pattern search
       numSeeds = 10;  % Number of starting points
       results = cell(numSeeds, 1);
       fvals = zeros(numSeeds, 1);
       options = optimoptions('patternsearch', ...
    'UseCompletePoll', true, ...       % Ensures full polling at each iteration
    'UseCompleteSearch', true, ...     % Searches the full pattern
    'Display', 'iter', ...             % Show iteration progress
    'MaxIterations', 100, ...         % Increase iterations
    'MeshTolerance', 1e-6, ...         % Fine-tune search precision
    'StepTolerance', 1e-8, ...         % Stop criterion
    'InitialMeshSize', 1.0);           % Controls exploration range
      
      parfor i = 1:numSeeds
          x0 = lb + (ub - lb) .* rand(size(lb));  % Random initial point in bounds
          [x_opt, fval] = patternsearch(@(x) costFunction(x, param_data), ...
                                        x0, [], [], Aeq, beq, lb, ub, [], options);
          results{i} = x_opt;
          fvals(i) = fval;
      end
      
      % Choose the best result
      [bestFval, bestIdx] = min(fvals);
      bestSolution = results{bestIdx};
      optimized_params = bestSolution;
    
    % Save optimized parameters
    
    %param_data.EstimationResults = optimized_params;
    %writetable(param_data, 'optimized_parameters.xlsx');
    %disp('Optimized parameters saved to optimized_parameters.xlsx');
end

function cost = costFunction(x, param_data)
    % Compute estimated parameters based on x multipliers
    [cost, ~] = estimate_newton_ohm(x, param_data);
end
function save_excel(x)
    param_data = readtable('excel_para_ori.xlsx');
    [~, data_dict] = estimate_newton_ohm(x, param_data);
    % Save optimized parameters
    param_names = param_data{:, 1};   
    
    % Match parameter names and retrieve values
    n_params = length(param_names);
    estimation_results = cell(n_params, 1); % Placeholder for estimated results
    
    for i = 1:n_params
        param_name = param_names{i};
        if isfield(data_dict, param_name)
            estimation_results{i} = data_dict.(param_name); % Fetch value if exists
        else
            estimation_results{i} = NaN; % Mark as missing if not found
        end
    end
    
    % Add a new column to the table
    param_data.EstimationResults = estimation_results;
    
    % Save the updated table to a new Excel file
    output_file_name = 'updated_parameters_opt.xlsx';
    writetable(param_data, output_file_name);
    
    disp(['Results have been saved to ', output_file_name]);
end
function [cost, data_dict] = estimate_newton_ohm(x, param_data)
    % Compute estimated parameters based on x multipliers


    BW        = 70; %data.BW;                         % Body weight (kg)                
    Hgt       = 172;%data.Hgt;                        % Height (cm)                     
    Gender    = 2;%data.Gender;                     % Gender (1=female, 2=male)       

    %% Total blood volume (mL)
    if (Gender == 2)  
           ToTBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * BW) + 0.6041) * 1000; % Total blood volume (mL)
           BSA   = 0.000579479 * BW^0.38 * Hgt^1.24;                          % Body surface area (m^2)
       else
           ToTBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * BW) + 0.1833) * 1000; % Total blood volume (mL)
           BSA   = 0.000975482 * BW^0.46 * Hgt^1.08;                          % Body surface area (m^2)
    end
    
    %% Cardiac Output (mL/s)

    Ave_HR        = 60;%data.Ave_HR;                 % Average heart rate (beats/min)    
    %TotFlow       = TotBV/60;                    % Total flow
    %CO            = TotFlow;                     % Cardiac output 

    HI            = Ave_HR/60;                   % Heart rate (s)
    %HI              =100;                   % Heart rate

    PM = 100; %1/3 * 120 + 2/3*75;
    ToTBV = 5027.6;  
    TotFlow = 5027.6/60;  
    CO = TotFlow;
    %% Blood Pressure
    P_m_p = x(1) * PM;
    P_h_p = x(2) * PM;
    P_p_p = x(3) * PM;
    P_s_p = x(4) * PM;
    P_e_p = x(5) * PM;
    P_b_p = x(6) * PM - 10; %Intracraneal pressure has to be substracted
    

    P_m_v = 3;
    P_h_v = 4.5;
    P_p_v = 5;
    P_s_v = 6;
    P_e_v = 7.5;
    P_b_v = 5;
    P_vc_v = 2.5;   %This value has to be correctly acquired,.

    %% Volume (Delta Volume = VTotal - Vunestressed)
    V_m_Total =  x(7)  * ToTBV;  % Muscle Circulation (Active + Resting, Arterial + Venous): 15%
    V_h_Total =  x(8)  * ToTBV;   % Heart: 7%
    V_p_Total =  x(9)  * ToTBV;   % Pulmonary Circulation: 9%
    V_s_Total =  x(10) * ToTBV;   % Splanchnic Circulation (Arterial + Venous): 20%
    V_e_Total =  x(11) * ToTBV;   % Extrasplanchnic Circulation (Arterial + Venous): 15%
    V_b_Total =  x(12) * ToTBV;   % Brain Circulation (Arterial + Venous): 8%
    V_vc_Total = x(13) * ToTBV;  % Vena Cava (Superior & Inferior): 5%   
   
    V_m_v = x(14) * V_m_Total;  
    V_h_v =  x(15) * V_h_Total;  
    V_p_v =  x(16) * V_p_Total;  
    V_s_v =  x(17) * V_s_Total;  
    V_e_v =  x(18) * V_e_Total;  
    V_b_v =  x(19) * V_b_Total;  
    V_vc_v =  V_vc_Total;  
    
    V_m_p =  V_m_Total  - V_m_v   ;
    V_h_p =   V_h_Total  - V_h_v  ; 
    V_p_p =   V_p_Total  - V_p_v  ; 
    V_s_p =   V_s_Total  - V_s_v  ; 
    V_e_p =   V_e_Total  - V_e_v  ; 
    V_b_p =   V_b_Total  - V_b_v  ; 
    V_vc_p = 0; 

    V_m_v_un = x(20) * V_m_v  ;
    V_h_v_un = x(21) * V_h_v  ;
    V_p_v_un = x(22) * V_p_v  ;
    V_s_v_un = x(23) * V_s_v  ;
    V_e_v_un = x(24) * V_e_v  ;
    V_b_v_un = x(25) * V_b_v  ;
    V_vc_v_un =  V_vc_Total ;

    V_m_p_un = x(26) * V_m_p  ;
    V_h_p_un = x(27) * V_h_p  ;
    V_p_p_un = x(28) * V_p_p  ;
    V_s_p_un = x(29) * V_s_p  ;
    V_e_p_un = x(30) * V_e_p  ;
    V_b_p_un = x(31) * V_b_p  ;
    %V_vc_p_un =  V_vc_Total ;

    V_m_v_delta = V_m_v - V_m_v_un;
    V_h_v_delta =   V_h_v -  V_h_v_un;
    V_p_v_delta =   V_p_v -  V_p_v_un;
    V_s_v_delta =  V_s_v - V_s_v_un;
    V_e_v_delta =   V_e_v -  V_e_v_un;
    V_b_v_delta =   V_b_v -  V_b_v_un;


    V_m_p_delta = V_m_p - V_m_p_un;
    V_h_p_delta =   V_h_p -  V_h_p_un;
    V_p_p_delta =   V_p_p -  V_p_p_un;
    V_s_p_delta =  V_s_p - V_s_p_un;
    V_e_p_delta =   V_e_p -  V_e_p_un;
    V_b_p_delta =   V_b_p -  V_b_p_un;
    
    
    % V_rm_p = 0.45 * ToTBV * 0.01;
    % V_h_p = 0.103 * ToTBV * 0.01;
    % V_p_p = 0.459 * ToTBV * 0.01;
    % V_s_p = 1.19 * ToTBV * 0.01;
    % V_e_p = 0.57 * ToTBV * 0.01;
    % V_b_p = 0.3 * ToTBV * 0.01;
    % 
    % 
    % V_rm_v = 2.18 * ToTBV * 0.01;
    % V_h_v = 0.42 * ToTBV * 0.01;
    % V_p_v = 0.4556 * ToTBV * 0.01;
    % V_s_v = 6.2 * ToTBV * 0.01;
    % V_e_v = 2.77 * ToTBV * 0.01;
    % V_b_v = 1.27 * ToTBV * 0.01;
    % V_vc_v = 0.534 * ToTBV * 0.01;

    %% Flows
    Q_m =  x(32) * TotFlow;
    Q_h =  x(33) * TotFlow;
    Q_p =  x(34) * TotFlow;
    Q_s =  x(35) * TotFlow;
    Q_e =  x(36) * TotFlow;
    Q_b =  x(37) * TotFlow;
    Q_vc = x(38) * TotFlow;

    %% Resistances
    
    R_m_p = (P_m_p - P_m_v )/Q_m;
    R_h_p = (P_h_p - P_h_v )/Q_h;
    R_p_p = (P_p_p - P_p_v )/Q_p;
    R_s_p = (P_s_p - P_s_v )/Q_s;
    R_e_p = (P_e_p - P_e_v )/Q_e;
    R_b_p = (P_b_p - P_b_v )/Q_b;
    %R_vc = (P_vc_p - P_vc_v )/Q_vc;
    R_m_v = (P_m_v - P_vc_v )/Q_m;
    R_h_v = (P_h_v - P_vc_v )/Q_h;
    R_p_v = (P_p_v - P_vc_v )/Q_p;
    R_s_v = (P_s_v - P_vc_v )/Q_s;
    R_e_v = (P_e_v - P_vc_v )/Q_e;
    R_b_v = (P_b_v - P_vc_v )/Q_b;


    %% Capacitances

    %C_m_p = V_m_p_delta / P_m_p ;
    %C_h_p =  V_h_p_delta / P_h_p;
    %C_p_p =  V_p_p_delta / P_p_p;
    %C_s_p =  V_s_p_delta / P_s_p;
    %C_e_p =  V_e_p_delta / P_e_p;
    %C_b_p =  V_b_p_delta / P_b_p;
    %
    %C_m_v = V_m_v_delta / P_m_v;
    %C_h_v =  V_h_v_delta / P_h_v;
    %C_p_v =  V_p_v_delta / P_p_v;
    %C_s_v =  V_s_v_delta / P_s_v;
    %C_e_v =  V_e_v_delta / P_e_v;
    %C_b_v =  V_b_v_delta / P_b_v;
    %C_vc_v = V_vc_v / P_vc_v;
%
    %Renaming
    V_unstressed_am_p = V_m_p_un * 0.6;
    %V_unstressed_am_v = V_m_v_un;
    V_unstressed_b_p = V_b_p_un;
    V_unstressed_b_v = V_b_v_un;
    V_unstressed_e_p = V_e_p_un;
    %V_unstressed_e_v = V_e_v_un;
    V_unstressed_h_p = V_h_p_un;
    V_unstressed_h_v = V_h_v_un;
    %V_unstressed_la = 
    %V_unstressed_lv
    %V_unstressed_pa
    V_unstressed_pp = V_p_p_un;
    V_unstressed_pv = V_p_v_un;
    %V_unstressed_ra
    V_unstressed_rm_p = V_m_p_un * 0.4;
    %V_unstressed_rm_v = V_m_v_un;
    %V_unstressed_rv 
    V_unstressed_s_p = V_s_p_un;
    %V_unstressed_s_v = V_s_v_un;
    %V_unstressed_sa = 
    V_unstressed_vc = V_vc_v_un;

    V_u_am_v_0 = V_m_v_un * 0.6;
    V_u_e_v_0 = V_e_v_un;
    V_u_rm_v_0 = V_m_v_un * 0.4;
    V_u_s_v_0 = V_s_v_un;

    R_am_n = R_m_v;
    R_am_p_0 = R_m_p;
    R_b_n = R_b_v;
    R_e_n = R_e_v;
    R_e_p_0 = R_e_p;
    R_h_n = R_h_v;
    R_h_p_n = R_h_p;
    %R_la = ;
    %R_pa = ;
    R_pp = R_p_p;
    R_pv = R_p_v;
    %R_ra = ;
    R_rm_n = R_m_v;
    R_rm_p_0 = R_m_p;
    R_s_n = R_s_v;
    R_s_p_0 = R_s_p;

    %C_am_p = C_m_p;
    %C_rm_p = C_m_p;
    %C_am_v = C_m_v;
    %C_rm_v = C_m_v;
    %C_pp = C_p_p;
    %C_pv = C_p_v;


    vars = whos;  % Get all variables in the workspace
    data_dict = struct();
    
    for i = 1:length(vars)
        var_name = vars(i).name;
        data_dict.(var_name) = eval(var_name); % Save all variables into struct
    end

    param_names = param_data{:, 1};
    param_values = param_data{:, 2};
    
    % Match parameter names and retrieve values
    n_params = length(param_names);
    estimation_results = cell(n_params, 1); % Placeholder for estimated results

    error = 0;
    for i = 1:n_params
        param_name = param_names{i};
        if isfield(data_dict, param_name)
            if isnan(data_dict.(param_name)) || isnan(param_values(i))
                disp("error");
            else
                error = error + ((data_dict.(param_name) -  param_values(i))/param_values(i))^2; % Compute difference
            end
            
        end
    end
    
    
    % Compute Euclidean distance between estimated and original parameters
    cost = sqrt(error/n_params);
end
