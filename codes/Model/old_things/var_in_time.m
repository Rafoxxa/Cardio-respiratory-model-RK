function varInTime = var_in_time(var, dt, final_time)
    time_vector = 0:dt:finla_time;
    if length(time_vector) < length(var)
        varInTime = var(1:length(time_vector));
    else
        disp('time vector length is equal or greater than the var variable');
    end

end