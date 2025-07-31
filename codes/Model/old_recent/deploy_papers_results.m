function [mrco2, vars] = deploy_papers_results(label, solver_fun, data_structurer_fun, varargin)
    global delays_global;
    global all_global;
    global externals_global;
    global dV_out;

if label == 'Sarmiento-2023-fig-3'
    [mrco2, vars] = steady_state_simulation(solver_fun, data_structurer_fun, varargin{:})    

end




function [mrco2, vars] = steady_state_simulation(solver_fun, data_structurer_fun, varargin)
    
    p = inputParser;
    p.KeepUnmatched = true;
    
    addParameter(p, 'model', '');  
    addParameter(p, 'pars', '');   
    addParameter(p, 'init', '');   
    addParameter(p, 'taus', '');   
    addParameter(p, 'simulation_time', '');   
    addParameter(p, 'dt', '');   
    addParameter(p, 'control_on', '');   
    parse(p, varargin{:});
    params = p.Results;

    structure_data = data_structurer_fun;
    solver = solver_fun;
    
    %let's put the initial conditions of the variables in the same range as the first steady state 
    %params.init("dV") = 11;
    params.init("PACO2") = 35;
    params.init("PAO2") = 88;
    %from the steady state plots
    TI_series = [1.4,1.4,1.35,1.3, 1.3, 1.3, 1.3, 1.3];
    Tresp_series = [3.5, 3, 2.85, 2.85, 2.85, 2.85, 2.72, 2.72];
    
    steps = 8;
    params.simulation_time = 10;   
    VARS =     zeros(10, steps);
    VCO2_DOT = zeros(1, steps);

    for step = 1:steps
        %let's take higher MO2 and MCO2 value, maybe there's is some sensitivity issue around.
        

        %params.init('TI') = TI_series(step);
        %params.init('Tresp') = Tresp_series(step);
        %params.init('TE') = params.init('Tresp') - params.init('TI');
        [t, x_dot, x_vars, x_keys, index] = solver(params.model, params.pars, params.init, params.taus, params.simulation_time, params.dt, params.control_on);
        
        struct_vars = structure_data(x_vars, x_keys);
        struct_dots = structure_data(x_dot, x_keys);
        [P_sa_systolic_end, P_sa_diastolic_end] = find_hill_valley_last_T(struct_vars.P_sa, t, 4); %we want the last 2 seconds signal 
        [VT_end, ~] = find_hill_valley_last_T(struct_vars.V, t, 4);  %we want the last 4 seconds
        mean_P_sa_end = struct_vars.mean_P_sa(end);
        dVE_end = VT_end * 1/(struct_dots.fake_Tresp(end)) * 60;    %remember this must be lt/min. VT_end is in lt (because volume is in lt) and Tresp is in seconds-
        %Extras       
        %struct_vars.dVE = struct_dots.dVE;        
        %disp(struct_vars.dVE);
        struct_vars.fake_TI = struct_dots.fake_TI; 
        struct_vars.fake_Tresp = struct_dots.fake_Tresp; 
        
        %el dV(end) hay que resolver

        %cambiar HR, PM, dVE!!!

        %[P_sa_systolic_end, P_sa_diastolic_end] = find_hill_valley_last_T(struct_vars.P_sa, t, 3); 
        try
        VARS(:,step) = [dVE_end, VT_end, 60/struct_vars.fake_Tresp(end), struct_vars.fake_TI(end),   struct_vars.PACO2(end), 60/struct_vars.Theart(end),  P_sa_systolic_end , P_sa_diastolic_end, (0.66*P_sa_diastolic_end + 0.33*P_sa_systolic_end), struct_vars.PAO2(end)];
        catch
            disp('error of min peak');
        end
        VCO2_DOT(:,step) = params.pars('MRCO2');

        params.pars('MRO2') =  params.pars('MRO2') + 0.1;
        params.pars('MRCO2') =  params.pars('MRCO2') + 0.1; 
        params.init('MRtO2') =  params.pars('MRO2');
        params.init('MRtCO2') =  params.pars('MRCO2');
    end    
        
       
    mrco2 = VCO2_DOT;
    vars = VARS;
    % for index = 1:8 
    %     var =   VARS(index,:);
    %     MRCO2 = VCO2_DOT(1,:);
    %     figure;
    % 
    %     plot(MRCO2, var);
    % end


    function [hill, valley] = find_hill_valley_last_T(wave, t_vector, T)
        
        
        from_T = t_vector(end) - T;     %if we want the last 3 seconds, we should take the t_end - 3 seconds index
        if from_T <= 0
             error("T must not be lower or equal than total time simulation");
        end

        index_T = find_index_t_seconds(from_T,t_vector);

        %pressure_wave_final_50_percent = pressure_wave(length(pressure_wave) - uint16(length(pressure_wave) * 0.5) : length(pressure_wave));
        wave_from_T = wave(index_T:end);
        [peaks, peak_locs] = findpeaks(wave_from_T, 'MinPeakHeight', max(wave_from_T) * 0.5);
        [valleys, valley_locs] = findpeaks(-wave_from_T, 'MinPeakHeight', min(-wave_from_T) * 0.5);
        hill = max(peaks);
        valley = -min(valleys);
    end
        
    function index_T = find_index_t_seconds(T, t_vector)
        diff = abs(t_vector - T);
        min_diff = min(diff);
        index_T = find(diff == min_diff);

    end
 end


end

