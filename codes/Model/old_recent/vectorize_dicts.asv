function vectorize_dicts(in_run_ode, in_model, out_run_ode, out_model)
%function transform_files(in_run_ode, in_model, out_run_ode, out_model)
transform_model(in_model, out_model);
transform_run_ode(in_run_ode, out_run_ode);
%end

function transform_run_ode(in_run_ode, out_run_ode)

    [pars, init] = load_global_easy();    
    keys_pars = keys(pars); 
    keys_init = keys(init); 
    % Crear un diccionario que mapea llaves a sus índices
    index_map_pars = containers.Map(keys_pars, 1:length(keys_pars));
    
    % Hacer lo mismo para el diccionario init si es necesario    
    index_map_init = containers.Map(keys_init, 1:length(keys_init));

    fid = fopen(in_run_ode, 'rt');
    original_code = fread(fid, '*char')';
    fclose(fid);

    expr_pars = 'pars\(["'']([^"''"]+)["'']\)';
    expr_PARS = 'PARS\(["'']([^"''"]+)["'']\)';
    expr_optim = 'optim_pars\(["'']([^"''"]+)["'']\)';

    expr_init = 'init\(["'']([^"''"]+)["'']\)';
    expr_cycle = 'cycle.init_vars\(["'']([^"''"]+)["'']\)';
    expr_INDEX = 'INDEX\(["'']([^"''"]+)["'']\)';

    [matched_parts_pars, tokens_pars] = regexp(original_code, expr_pars, 'match', 'tokens');
    [matched_parts_PARS, tokens_PARS] = regexp(original_code, expr_PARS, 'match', 'tokens');
    [matched_parts_optim, tokens_optim] = regexp(original_code, expr_optim, 'match', 'tokens');

    [matched_parts_init, tokens_init] = regexp(original_code, expr_init, 'match', 'tokens');
    [matched_parts_cycle, tokens_cycle] = regexp(original_code, expr_cycle, 'match', 'tokens');
    [matched_parts_INDEX, tokens_INDEX] = regexp(original_code, expr_INDEX, 'match', 'tokens');

    new_code = original_code;
    for i = 1:length(matched_parts_pars)
        key = tokens_pars{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_pars, key)
            % Create the replacement string
            
            replacement = sprintf('pars(%d)', index_map_pars(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            warning('Key "%s" not present in index_map_pars.', key);
            replacement = matched_parts_pars{i};
        end
        new_code = strrep(new_code, matched_parts_pars{i}, replacement);
    end 

    for i = 1:length(matched_parts_PARS)
        key = tokens_PARS{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_pars, key)
            % Create the replacement string
            
            replacement = sprintf('PARS(%d)', index_map_pars(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            warning('Key "%s" not present in index_map_pars.', key);
            replacement = matched_parts_PARS{i};
        end
        new_code = strrep(new_code, matched_parts_PARS{i}, replacement);
    end 

    for i = 1:length(matched_parts_optim)
        key = tokens_optim{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_pars, key)
            % Create the replacement string
            
            replacement = sprintf('optim_pars(%d)', index_map_pars(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            warning('Key "%s" not present in index_map_pars.', key);
            replacement = matched_parts_optim{i};
        end
        new_code = strrep(new_code, matched_parts_optim{i}, replacement);
    end 

    for i = 1:length(matched_parts_init)
        key = tokens_init{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_init, key)
            % Create the replacement string
            
            replacement = sprintf('init(%d)', index_map_init(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts_init{i};
        end
        new_code = strrep(new_code, matched_parts_init{i}, replacement);
    end 

    for i = 1:length(matched_parts_cycle)
        key = tokens_cycle{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_init, key)
            % Create the replacement string
            
            replacement = sprintf('cycle.init_vars(%d)', index_map_init(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts_cycle{i};
        end
        new_code = strrep(new_code, matched_parts_cycle{i}, replacement);
    end 

    for i = 1:length(matched_parts_INDEX)
        key = tokens_INDEX{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_init, key)
            % Create the replacement string
            
            replacement = sprintf('%d', index_map_init(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts_INDEX{i};
        end
        new_code = strrep(new_code, matched_parts_INDEX{i}, replacement);
    end 

        % 1. Add lines after "DT = dt"
    pattern1 = '(DT\s*=\s*dt\s*;\s*)';
    replacement1 = '$1\nif isa(pars, ''dictionary'')\npars = pars.values;\nend\nif isa(init, ''dictionary'')\ninit = init.values;\nend\n';
    new_code = regexprep(new_code, pattern1, replacement1);

    % 2. Replace "cycle.init_vars.values" with "cycle.init_vars"
    pattern2 = 'cycle\.init_vars\.values';
    replacement2 = 'cycle.init_vars';
    new_code = regexprep(new_code, pattern2, replacement2);

  % 3. Change "cycle.init_vars = dictionary(...)" to "cycle.init_vars = cycle.x_vars"
    pattern3 = 'cycle\.init_vars\s*=\s*dictionary\([^)]*\)';
    replacement3 = 'cycle.init_vars = cycle.x_vars(:,end';
    new_code = regexprep(new_code, pattern3, replacement3);

    fid = fopen(out_run_ode, 'wt');
    fwrite(fid, new_code);
    fclose(fid);
    

      
      
end

function transform_model(in_model, out_model)
    
    [pars, init] = load_global_easy();    
    keys_pars = keys(pars); 
    keys_init = keys(init);    
    
    % Crear un diccionario que mapea llaves a sus índices
    index_map_pars = containers.Map(keys_pars, 1:length(keys_pars));
    
    % Hacer lo mismo para el diccionario init si es necesario    
    index_map_init = containers.Map(keys_init, 1:length(keys_init));

    fid = fopen(in_model, 'rt');
    original_code = fread(fid, '*char')';
    fclose(fid);

    expr_pars = 'pars\(["'']([^"''"]+)["'']\)';
   
    % Find all matches in the original code
    [matched_parts, tokens] = regexp(original_code, expr_pars, 'match', 'tokens');

    % Loop through all matched parts and tokens
    new_code = original_code;
    for i = 1:length(matched_parts)
        key = tokens{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_pars, key)
            % Create the replacement string
            
            replacement = sprintf('pars(%d)', index_map_pars(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            %warning('Key "%s" not present in index_map_pars.', key);
            replacement = matched_parts{i};
        end
        new_code = strrep(new_code, matched_parts{i}, replacement);
    end 

    expr_init = 'y\(["'']([^"''"]+)["'']\)';
    
   
    % Find all matches in the original code
    [matched_parts, tokens] = regexp(new_code, expr_init, 'match', 'tokens');
    


    % Loop through all matched parts and tokens
    
    for i = 1:length(matched_parts)
        key = tokens{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_init, key)
            % Create the replacement string
            
            replacement = sprintf('y(%d)', index_map_init(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            %warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts{i};
        end
        new_code = strrep(new_code, matched_parts{i}, replacement);
    end 

    expr_xdot = 'xdot_dict\(["'']([^"''"]+)["'']\)';
   
    % Find all matches in the original code
    [matched_parts, tokens] = regexp(new_code, expr_xdot, 'match', 'tokens');

    for i = 1:length(matched_parts)
        key = tokens{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_init, key)
            % Create the replacement string
            
            replacement = sprintf('xdot(%d)', index_map_init(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            %warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts{i};
        end
        new_code = strrep(new_code, matched_parts{i}, replacement);
    end 

    expr_internal = 'internal_variables\(["'']([^"''"]+)["'']\)';
   
    % Find all matches in the original code
    [matched_parts, tokens] = regexp(original_code, expr_internal, 'match', 'tokens');
    strArray = strings(0);
    for i=1:length(matched_parts)
        
        strArray(i) = tokens{i};
        
    end

    index_map_internal = containers.Map(strArray, 1:length(strArray));

    for i = 1:length(matched_parts)
        key = tokens{i}{1};  % Extract the key (e.g., "R" or 'R')
        if isKey(index_map_internal, key)
            % Create the replacement string
            
            replacement = sprintf('internal_variables(%d)', index_map_internal(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            %warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts{i};
        end
        new_code = strrep(new_code, matched_parts{i}, replacement);
    end 



    %%%%%%%%%%%%%%%%% para los index

    %map es el tiny
    
    % Regular expression to find the line where tiny_y_keys is defined
    pattern = 'tiny_y_keys\s*=\s*\[([^\]]+)\]';

    % Find the tiny_y_keys line and extract the keys
    tokens = regexp(new_code, pattern, 'tokens');

    % Check if the array was found
    if ~isempty(tokens)
        % Extract the keys as a single string
        keys_str = tokens{1}{1};

        % Use another regex to extract individual keys from the string
        key_pattern = '"([^"]+)"';
        key_matches = regexp(keys_str, key_pattern, 'tokens');
        
        % Flatten the cell array of tokens
        tiny_y_keys = [key_matches{:}];

        % Create a mapping from keys to indices
        index_map_tiny = containers.Map(tiny_y_keys, 1:numel(tiny_y_keys));  
    end

    %para cambiar el index_fun
    expr_indexfun = 'index_fun\(\s*\w+,\s*(\w+)\s*\)';
    new_code = regexprep(new_code, expr_indexfun, '$1');

    % idearme cómo reemplazar las llaves que entran en índices, el patron que parte get_delayed_valued(abc,...,z), que z me lo cambie
    expr_delay = 'get_delayed_value\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^,]+)\s*,\s*(''[^'']+?'')\s*\)\s*;';
    [matched_parts, tokens] = regexp(new_code, expr_delay, 'match', 'tokens');
    for i = 1:length(matched_parts)
        final_arg = strrep(tokens{i}{end}, '''', '');  
        if isKey(index_map_tiny, final_arg)
            % Create the replacement string
            format_str = 'get_delayed_value(%s, %s, %s, %s, %s, %s, %s, %d);';
            index_value = index_map_tiny(final_arg);
            replacement = sprintf(format_str, tokens{i}{1}, tokens{i}{2}, tokens{i}{3}, tokens{i}{4}, tokens{i}{5}, tokens{i}{6}, tokens{i}{7}, index_value);          
           
        else
            % If the key is not present, return the original match and issue a warning
            %warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts{i};
        end
        new_code = strrep(new_code, matched_parts{i}, replacement);
    end 
    

    %ahora vamos a cambiar el index general    

    expr_index = 'index\(\s*([^,]+)\s*,\s*''([^'']+)''\s*\)';
    [matched_parts, tokens] = regexp(new_code, expr_index, 'match', 'tokens');

    for i = 1:length(matched_parts)
        key = tokens{i}{end};  
        if isKey(index_map_tiny, key)
            % Create the replacement string
            
            replacement = sprintf('%d', index_map_tiny(key));
           
        else
            % If the key is not present, return the original match and issue a warning
            %warning('Key "%s" not present in index_map_init.', key);
            replacement = matched_parts{i};
        end
        new_code = strrep(new_code, matched_parts{i}, replacement);
    end 


  

    

  

    % Apply the replacement using regexprep
    %new_code = regexprep(original_code, expr_pars,  @(match, key) replacement_pars(match, key, index_map_pars));

    %expr_init = 'init\(["'']([^"''"]+)["'']\)';
    %new_code = regexprep(new_code, expr_init, '${sprintf("y(%d)", index_map_init($1))}');

    fid = fopen(out_model, 'wt');
    fwrite(fid, new_code);
    fclose(fid);
    
   
    
    

end
end