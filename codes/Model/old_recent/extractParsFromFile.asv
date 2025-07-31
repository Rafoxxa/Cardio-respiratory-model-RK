function extractParsFromFile(filepath)
    % Check if file exists
    if ~isfile(filepath)
        error('File not found: %s', filepath);
    end
    
    % Read file content using importdata
    data = importdata(filepath);
    
    % Ensure data was loaded correctly (cell array with strings)
    if isempty(data)
        warning('No data found in file: %s', filepath);
        return;
    end
    
    % Extract the desired string (assuming data is a cell array)
    text = char(data{1});  % Assuming the first cell contains the text
    
    % Find matches of pars("...") or pars('...')
     matches = regexp(text, 'pars\(([^)]+)\)', 'tokens');


    
    % Ensure matches are correctly extracted
    if isempty(matches)
        warning('No matches found in file: %s', filepath);
        return;
    end
    
    % Convert cell array of cells into a single string array
    matches = string([matches{:}]);
    
    % Remove duplicates from the matches
    matches = unique(matches);

    % Open the file for writing
    fid = fopen('output.txt', 'w');
    
    % Check if the file was opened successfully
    if fid == -1
        error('Unable to open the file for writing.');
    end
    
    % Write each unique match to the file on a new line
    for i = 1:length(matches)
        fprintf(fid, '%s\n', matches(i));
    end
    
    % Close the file
    fclose(fid);
    
    fprintf('Matches saved to: output.txt\n');
end
