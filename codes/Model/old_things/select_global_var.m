function var = select_global_var(global_var, var_name, keys)

    index = find(keys == var_name);
    var = global_var(index, :);
end