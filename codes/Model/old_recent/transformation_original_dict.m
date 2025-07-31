
transform_files("load_global_easy.m", "model_basic_vascular.m", "load_global_vec.m", "model_basic_vec.m");

function transform_files(input_params_file, input_model_file, output_params_file, output_model_file)
    % Transform the parameters file
    %transform_parameters(input_params_file, output_params_file);
    
    % Transform the model file
    transform_model(input_model_file, output_model_file);
end

function transform_parameters(input_file, output_file)
    % Read the content of the parameters file
    fid = fopen(input_file, 'rt');
    original_code = fread(fid, '*char')';
    fclose(fid);
    %original_code = strrep(original_code, '\r\n', '\n');
    %filtered_text = regexprep(original_code, '^\s*%.*\n', '', 'lineanchors');

    % Regular expression to find dictionary parameter assignments
    
    dict_pattern = '(?<dict_name>\w+)\(("|'')(?<key>\w+)\2\)\s*=\s*(?<value>[^;]+)';


    % Find all dictionary assignments
    dict_matches = regexp(original_code, dict_pattern, 'names');

    % Create vectors for the parameters
    vector_names = unique({dict_matches.dict_name});  % Dict names
    vector_lines = '';

    % Loop through each dictionary name
    for i = 1:length(vector_names)
        vec_name = vector_names{i};  % Dict name size vector
        vec_values = {dict_matches(strcmp({dict_matches.dict_name}, vec_name)).value}; % Vector values

        % Resolve references
        for j = 1:length(vec_values)
            vec_values{j} = resolve_reference(vec_values{j}, dict_matches);
        end

        % Create the vector line
        vec_line = sprintf('\n%s = [%s];\n', vec_name, strjoin(vec_values, ', '));  % This creates the complete vector in one row
        vector_lines = strcat(vector_lines, vec_line);  % Create the code
    end

    % Write the transformed code to the output file
    fid = fopen(output_file, 'wt');
    fwrite(fid, vector_lines);
    fclose(fid);
end

function resolved_value = resolve_reference(value, dict_matches)
    % Pattern to match any dictionary reference
    ref_pattern = '(?<dict_name>\w+)\((["''])(?<key>\w+)\2\)';
    
    % Recursively resolve references
    while contains(value, '(')
        %disp(value)
        key_match = regexp(value, ref_pattern, 'names');
        if isempty(key_match)
            break; % No more references to resolve
        end
        
        for k = 1:length(key_match)
            dict_name = key_match(k).dict_name;
            referenced_key = key_match(k).key;

            % Find the value of the referenced key in the specified dictionary
            referenced_value = {dict_matches(strcmp({dict_matches.dict_name}, dict_name) & ...
                                             strcmp({dict_matches.key}, referenced_key)).value};
            if ~isempty(referenced_value)
                % Replace the reference with its resolved value
                new_value_simple = strrep(value, ...
                    sprintf('%s(''%s'')', dict_name, referenced_key), ...
                    resolve_reference(referenced_value{1}, dict_matches));
                if strcmp(value, new_value_simple)
                
                    new_value_double = strrep(value, ...
                    sprintf('%s("%s")', dict_name, referenced_key), ...
                    resolve_reference(referenced_value{1}, dict_matches));
                    value = new_value_double;
                else
                    value = new_value_simple;
                    
                end
                
                
            end
        end
    end
    
    resolved_value = value;  % Return the fully resolved value
end


 
function transform_model(input_file, output_file)
    % Read the content of the model file
    fid = fopen(input_file, 'rt');
    original_code = fread(fid, '*char')';
    fclose(fid);

    % Regular expression to find dictionary parameter usage
    dict_usage_pattern = '(?<dict_name>\w+)\(("|'')(?<key>\w+)\2\)';

    % Find all dictionary usages
    dict_usages = regexp(original_code, dict_usage_pattern, 'names');

    % Create a unique list of dictionary names and their fields
    vector_names = unique({dict_usages.dict_name});
    new_code = original_code;

    for i = 1:length(vector_names)
        vec_name = vector_names{i};
        vec_fields = unique({dict_usages(strcmp({dict_usages.dict_name}, vec_name)).field});

        % Replace dictionary references with vector references
        for j = 1:length(vec_fields)
            field = vec_fields{j};
            % Replace all occurrences of dict_name.field with vec_name(j)
            field_pattern = sprintf('%s\\.%s', vec_name, field);
            replacement = sprintf('%s(%d)', vec_name, j);
            new_code = regexprep(new_code, field_pattern, replacement);
        end
    end

    % Write the transformed code to the output file
    fid = fopen(output_file, 'wt');
    fwrite(fid, new_code);
    fclose(fid);
end
