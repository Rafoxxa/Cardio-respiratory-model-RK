function plot_global_vars_in_time(global_var, keys, var_names_array, dt, final_time)
    %try not to show more than 3 vars, to avoid time collapsing of the solver
    for index = 1:length(var_names_array)

        var_name = var_names_array(index);
        var_i_array = select_global_var(global_var, var_name, keys);
        var_i_InTime = var_in_time(var_i_array, dt, final_time);
        figure;
        plot(0:dt:final_time, var_i_InTime);
        title(var_name);
    end
end