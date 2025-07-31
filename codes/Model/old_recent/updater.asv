% Read the main code text file
main_code_text = fileread('model_basic.txt');

% Find occurrences of 'pars' and 'vars' in the text
matches_pars = regexp(main_code_text, 'pars\.[a-zA-Z0-9_]+', 'match');
matches_vars = regexp(main_code_text, 'vars\.[a-zA-Z0-9_]+', 'match');

% Generate updated load_global.m content
updated_load_global_content = [
    'function [pars, init] = load_global_easy(pars, vars)\n', ...
    '    % Your updated code here\n\n', ...
    '    disp(''Global parameters:'');\n', ...
    '    disp(pars);\n', ...
    '    disp(''Global variables:'');\n', ...
    '    disp(vars);\n', ...
    'end\n'
];

% Replace placeholders with actual code
for i = 1:length(matches_pars)
    updated_load_global_content = strrep(updated_load_global_content, ...
        ['pars.' matches_pars{i}(6:end)], matches_pars{i});
end

for i = 1:length(matches_vars)
    updated_load_global_content = strrep(updated_load_global_content, ...
        ['vars.' matches_vars{i}(6:end)], matches_vars{i});
end

% Write updated load_global.m file
fid = fopen('load_global.m', 'w');
fprintf(fid, updated_load_global_content);
fclose(fid);

disp('Updated load_global.m file created.');
