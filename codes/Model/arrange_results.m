function struct_vars = arrange_results(x_dot, x_vars, x_keys, t) 

    struct_vars = structure_data(x_vars, x_keys);
    struct_dots = structure_data(x_dot, x_keys);

    struct_vars.dV = struct_dots.V;
    struct_vars.dV_total_rv = struct_dots.V_total_rv;
    struct_vars.dV_total_lv = struct_dots.V_total_lv;
    struct_vars.dV_total_ra = struct_dots.V_total_ra;
    struct_vars.dV_total_la = struct_dots.V_total_la;
    %Fake derivatives (this are created in oder to have the history acumulated of important internal variables that in esence are not state variables)
    %struct_vars.Gaw = struct_dots.Gaw;
    struct_vars.aO2 = struct_dots.aO2;
    struct_vars.aCO2 = struct_dots.aCO2;
    struct_vars.Nt = struct_dots.Nt;    
    struct_vars.ddPaO2 = struct_dots.dPaO2;  
    
    struct_vars.TI = struct_dots.fake_TI;
    struct_vars.Tresp = struct_dots.fake_Tresp;  
    
    %Derived from operations
    struct_vars.HR = struct_vars.Theart.^-1;      

    struct_vars.VT = data_processing("volume", struct_vars.V, t); %volume
    out_pressure = data_processing("pressure", struct_vars.P_sa, t);
    struct_vars.ps = out_pressure{2};
    struct_vars.pd = out_pressure{3};
    struct_vars.pm = struct_vars.mean_P_sa;
    
    
    %Experimental
    % struct_vars.dff = (struct_vars.Qpp - struct_vars.Qbp).*(struct_vars.aO2 - struct_vars.vO2);
    % struct_vars.perfusionpart = 863 .* struct_vars.Qpp./1000 .* (struct_vars.vO2 - struct_vars.aO2);
    % struct_vars.airpart = struct_vars.dV .* (struct_vars.P_5O2 - struct_vars.PAO2);
    % struct_vars.press_diff = struct_vars.P_sa - struct_vars.P_sp;
    struct_vars.Vla = struct_vars.dVua + struct_vars.dV;
    struct_vars.Pua = struct_vars.Ppl + struct_vars.Vla * 3.02;
    struct_vars.gaw = 1 - struct_vars.Pua*(-1/40);
    struct_vars.pseudo_dVE = struct_vars.V ./ struct_vars.Tresp;

    struct_vars.Vtotal = struct_vars.V_total_e_v + struct_vars.V_total_s_v + struct_vars.V_total_b_v + struct_vars.V_total_h_v + struct_vars.V_total_rm_v + struct_vars.V_total_am_v + struct_vars.V_total_vc + struct_vars.V_total_lv + struct_vars.V_total_la + struct_vars.V_total_rv + struct_vars.V_total_ra + struct_vars.V_total_pa + struct_vars.V_total_pp + struct_vars.V_total_pv;

    %Data structure
    function results_table =  structure_data(x_values, xkeys)
        
        results_table = struct();

        for i = 1:length(xkeys)
            results_table.(xkeys{i}) = x_values(i, :);
        end
    end

end