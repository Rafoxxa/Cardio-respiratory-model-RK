%plotting_data_for_manuscript
%DATA
% 
% [t1_N, y1_N, poly1VO2_N, poly1VCO2_N, poly1fiO2_N, poly1TI_N, poly1Tresp_N, ~] = data_preprocessing2(1, "normoxia", ".", 0, "ori");
% [t4_N, y4_N, poly4VO2_N, poly4VCO2_N, poly4fiO2_N, poly4TI_N, poly4Tresp_N, ~] = data_preprocessing2(4, "normoxia", ".", 0, "ori");
% [t5_N, y5_N, poly5VO2_N, poly5VCO2_N, poly5fiO2_N, poly5TI_N, poly5Tresp_N, ~] = data_preprocessing2(5, "normoxia", ".", 0, "ori");
% [t6_N, y6_N, poly6VO2_N, poly6VCO2_N, poly6fiO2_N, poly6TI_N, poly6Tresp_N, ~] = data_preprocessing2(6, "normoxia", ".", 0, "ori");
% 
% [t1_H, y1_H, poly1VO2_H, poly1VCO2_H, poly1fiO2_H, poly1TI_H, poly1Tresp_H, ~] = data_preprocessing2(1, "hipoxia", "mix", 0, "ori");
% [t4_H, y4_H, poly4VO2_H, poly4VCO2_H, poly4fiO2_H, poly4TI_H, poly4Tresp_H, ~] = data_preprocessing2(4, "hipoxia", "mix", 0, "ori");
% [t5_H, y5_H, poly5VO2_H, poly5VCO2_H, poly5fiO2_H, poly5TI_H, poly5Tresp_H, ~] = data_preprocessing2(5, "hipoxia", "mix", 0, "ori");
% [t6_H, y6_H, poly6VO2_H, poly6VCO2_H, poly6fiO2_H, poly6TI_H, poly6Tresp_H, ~] = data_preprocessing2(6, "hipoxia", "mix", 0, "ori");




try
hipoxia = 1;
normoxia = 1 - hipoxia;

if normoxia
    fast_mask1 = t1_N < 908;
    y1_N(fast_mask1, 7:10) = nan;
    fast_mask6 = t6_N > 1748;
    y6_N(fast_mask6, 7:10) = nan;
elseif hipoxia
    fast_mask4 = t4_H > (35.05 * 60);
    y4_H(fast_mask4, 8:10) = nan;
    fast_mask1 = t1_H > (39.5 * 60);
    y1_H(fast_mask1, 8:10) = nan;

end



variables_to_show = [7,8,9,10];%[5,6];%[1, 2,3,4];[11, 12]
subject = 4;
width = 2;
height = 2;
with_mask = 0;

y_cells_H = {y1_H', y4_H', y5_H', y6_H'};
t_cells_H = {t1_H', t4_H', t5_H', t6_H'};

y_cells_N = {y1_N', y4_N', y5_N', y6_N'};
t_cells_N = {t1_N, t4_N, t5_N, t6_N};

VO2_cells_H = {poly1VO2_H, poly4VO2_H, poly5VO2_H, poly6VO2_H};
VCO2_cells_H = {poly1VCO2_H, poly4VCO2_H, poly5VCO2_H, poly6VCO2_H};


VO2_cells_N = {poly1VO2_N, poly4VO2_N, poly5VO2_N, poly6VO2_N};
VCO2_cells_N = {poly1VCO2_N, poly4VCO2_N, poly5VCO2_N, poly6VCO2_N};

TI_cells_N = {poly1TI_N, poly4TI_N, poly5TI_N, poly6TI_N};
Tresp_cells_N = {poly1Tresp_N, poly4Tresp_N, poly5Tresp_N, poly6Tresp_N};

TI_cells_H = {poly1TI_H, poly4TI_H, poly5TI_H, poly6TI_H};
Tresp_cells_H = {poly1Tresp_H, poly4Tresp_H, poly5Tresp_H, poly6Tresp_H};

VT1_N = {791.00, 705.00, 668.00, 766.00}; %t@VT1
VT1_H = {2864.00, 2721.00, 2799.00, 3184.00}; %t@VT1

tmax_N = 2300;
tmax_H = 3400;

variable_names = ["dVE", "VT", "TI", "Tresp","PAO2", "PACO2", "HR", "PS", "PD", "PM", "VO2", "VCO2", "Pow"];
subject_ids = [1,2,3,4];%[1,4,5,6];
units_table = readtable("variables_units.xlsx");

% if hipoxia
%    plot_data(y_cells_H, variables_to_show ,t_cells_H, variable_names , {}, units_table,  VO2_cells_H, VCO2_cells_H, VT1_H, subject_ids, tmax_H, with_mask);
% else
%    plot_data(y_cells_N, variables_to_show ,t_cells_N, variable_names , {}, units_table,  VO2_cells_N, VCO2_cells_N, VT1_N, subject_ids, tmax_N, with_mask);
% end

%average = 0;

%plot_masked_comparison(y_cells_N, t_cells_N, y_cells_H, t_cells_H, variables_to_show, variable_names, units_table, width, height, subject);

%plot_poly(t_cells_N, TI_cells_N, subject_ids);
old_hipoxia_paths = {'3300_sec_hipoxia-28-07-2025.mat','3300_sec_hipoxia-28-07-2025.mat','3300_sec_hipoxia-23-07-2025.mat','3300_sec_hipoxia-28-07-2025.mat'};
old_normoxia_paths = {'1200_sec_normoxia-28-07-2025.mat','1200_sec_normoxia-28-07-2025.mat','1200_sec_normoxia-23-07-2025.mat','1200_sec_normoxia-28-07-2025.mat'};

catch
    disp('no data loaded');
end

Csa1 = 0.28;
Csa2 = 0.18;
Csa3 = 0.08;
Lsa1 = 0.2;
Lsa2 = 0.02;
Lsa3 = 0.0022;

%Rsa = 4;
%bode_damping_analysis([Csa1, Lsa1, Csa2, Lsa2, Csa3, Lsa3]);
%plot_resistances('..\Simulations\simulation_after_fitting\5\1200_sec_normoxia-28-10-2025.mat', '..\Simulations\simulation_after_fitting\5\3300_sec_hipoxia-28-10-2025.mat');
%plot_topCV_params('parameters_for_manuscript.csv', {'gO2_e', 'gO2_p', 'gO2_s', 'vO2_e_n', 'vO2_s_n', 'aO2_n'});
%plot_topCV_params_box('parameters_for_manuscript.csv');

%paramGroups = {{'gO2_p'},{'gO2_s', 'gO2_e'}, {'vO2_b_n', 'vO2_e_n', 'vO2_s_n'}, {'aO2_n'}, {'GTsym', 'Wb_h_s', 'Wb_p_s', 'Wb_v_s', 'Wp_p_s', 'Wp_v'}, {'Wt_v_s'}, {'kcc_p_s'}, {'ub_a1'}, {'P_n', 'PaO2_ac_n'}, {'L_sa'}, {'Tsys_0'}, {'tau_M'}, {'f_ab_max', 'kev'}, {'x_h_s', 'x_p_s'}, {'Kbg'}, {'phi_max'}};
paramGroups = {{'gO2_p','gO2_s', 'gO2_e'}, {'vO2_b_n', 'vO2_e_n', 'vO2_s_n', 'aO2_n'},{'GTsym', 'Wb_h_s', 'Wb_p_s', 'Wb_v_s', 'Wp_p_s', 'Wp_v'}, {'P_n', 'PaO2_ac_n'}, {'f_ab_max', 'kev'}, {'x_h_s', 'x_p_s'}};
%paramGroups = {{'aO2_n'}}%, {'GTsym', 'Wb_h_s', 'Wb_p_s', 'Wb_v_s', 'Wp_p_s', 'Wp_v'}, {'Wt_v_s'}, {'kcc_p_s'}, {'ub_a1'}, {'P_n', 'PaO2_ac_n'}, {'L_sa'}, {'Tsys_0'}, {'tau_M'}, {'f_ab_max', 'kev'}, {'x_h_s', 'x_p_s'}, {'Kbg'}, {'phi_max'}};
plot_param_groups_box('parameters_for_manuscript.csv', paramGroups)

%data = validation_section();

%plot_J_iterations();

%old_hipoxia_paths = {'3300_sec_hipoxia-28-07-2025.mat','3300_sec_hipoxia-28-07-2025.mat','3300_sec_hipoxia-23-07-2025.mat','3300_sec_hipoxia-28-07-2025.mat'};
%old_normoxia_paths = {'1200_sec_normoxia-28-07-2025.mat','1200_sec_normoxia-28-07-2025.mat','1200_sec_normoxia-23-07-2025.mat','1200_sec_normoxia-28-07-2025.mat'};

%show_old_simulations_control([1,4,5,6], old_hipoxia_paths, old_normoxia_paths)

%SIMULATIONS
%subject_ids = [1,4,5,6];

%autonomic_vars = ["fh_s", "fp_s", "fv_s", "fv"];
% local_vars = ["xO2_b", "xO2_h", "xO2_rm", "xO2_am"];
% resp_vars = ["HR"];
% gas_vars = ["PACO2"];
% 
% autonomic_colors = [
%     [0.00, 0.45, 0.74];   % Blue
%     [0.53, 0.81, 0.98];   % Sky blue
%     [0.49, 0.18, 0.56];   % Purple
%     [0.40, 0.26, 0.13]    % Brown
% ];

% local_colors = [
%     0.1216, 0.4667, 0.7059;  % Azul (similar al default de MATLAB)
%     1.0000, 0.4980, 0.0549;  % Naranja
%     0.1725, 0.6275, 0.1725;  % Verde
%     0.8392, 0.1529, 0.1569;  % Rojo
% ];
% common_vars = ["MRtO2", "MRtCO2"];
% title_name = "Autonomic response";
% plot_simulations(subject_ids, autonomic_vars, autonomic_colors, title_name, common_vars);

% thresh_arr = [0.6, 0.8, 0.6, 0.5];
% plot_LSA(subject_ids, "all", thresh_arr);
% %The Sall was computed in the shared_params code


function plot_LSA(subject_ids, plot_type, thresh_arr)
    figure;
    for idx = 1:length(subject_ids)
        subplot(2,2,idx)
        parameter_analysis_fun_plot(subject_ids(idx), plot_type, thresh_arr(idx));
        hold on;
    end

end

function plot_simulations(subjects, vars, colors, title_name, common)

    units_table = readtable("variables_units.xlsx");    
    for idx = 1:length(subjects)

        path_n = sprintf("../Simulations/only_simulation/%d/1200_sec_normoxia-28-07-2025.mat", subjects(idx));
        path_h = sprintf("../Simulations/only_simulation/%d/3300_sec_hipoxia-28-07-2025.mat", subjects(idx));
        
        sim_n = load(path_n);
        sim_h = load(path_h);

        struct_vars_n = sim_n.struct_vars;
        time_n = sim_n.t;

        struct_vars_h = sim_h.struct_vars;
        time_h = sim_h.t;
        
        
        custom_plot("multiple_to_show", {vars, struct_vars_n, time_n, units_table, colors, common}); 
        str_title_n = sprintf("Subject %d: %s during normoxia", subjects(idx), title_name); 
        %title("holi");
        title(char(str_title_n));
        


        custom_plot("multiple_to_show", {vars, struct_vars_h, time_h, units_table, colors, common}); 
        str_title_h = sprintf("Subject %d: %s during hypoxia", subjects(idx), title_name);
        title(char(str_title_h));
        


    end
end

function plot_poly(time_cell, vo_cell,subject_ids)
    for index = 1:length(vo_cell)
        subplot(2,2,index);
        bestPolynomialFit(time_cell{index}, vo_cell{index}, 8, "VO2(l/min)", 2);
        name_title = sprintf("subject %d", subject_ids(index));
        title(name_title);
    end

end









%set(findall(gca, 'Type', 'Line'), 'LineWidth', 0.1);
function plot_data(y_cells, V, t_cells, nombres, together, units_table, vo2_cells, vco2_cells, vt1_times, subject_ids, tmax, with_mask)
     for ci = 1:length(vt1_times)    
         vt1_times{ci} = vt1_times{ci}/60;
         t_cells{ci} = t_cells{ci}/60;
     end

    tmax = tmax/60;
% y_cells: celda de matrices VxT, una por sujeto
% V: vector de índices de variables a graficar
% t_cells: celda de vectores de tiempo, uno por sujeto
% nombres: nombres de las variables (cell array de strings)
% together: cell array agrupando variables por subplot
% units_table: tabla con columnas 'Variable' y 'MeasureUnit'
% vo2_cells, vco2_cells: no se usan, mantenidos por compatibilidad
% vt1_times: celda con tiempos VT1 por sujeto (en segundos)
% subject_ids: vector con los identificadores reales de cada sujeto (ej: [1 4 5 6])

    if ~iscell(y_cells)
        error("y_cells debe ser una celda de matrices VxT.");
    end
    if ~iscell(t_cells)
        error("t_cells debe ser una celda de vectores de tiempo.");
    end
    if length(y_cells) ~= length(t_cells)
        error("y_cells y t_cells deben tener el mismo largo.");
    end
    if isempty(together)
        together = num2cell(V(:)');
    end
    if nargin < 10 || isempty(subject_ids)
        subject_ids = 1:length(y_cells);  % valores por defecto si no se entrega
    end
    if length(subject_ids) ~= length(y_cells)
        error("subject_ids debe tener la misma longitud que y_cells.");
    end

    T_units = units_table;
    num_subjects = length(y_cells);
    num_plots = numel(together);

    base_colores = lines(num_subjects);
    %base_colores(3,1) = min(base_colores(3,1) * 1.8, 1);

    rgb = base_colores(3,:);
    rgb(1) = min(rgb(1) * 1.1, 1);   % un poco más rojo
    rgb(2) = rgb(2) * 0.75;          % menos verde
    rgb(3) = rgb(3) * 0.5;           % un toque más oscuro (menos azul)

    base_colores(3,:) = rgb;

    line_styles = {'-', '--', ':', '-.'};
    y1 = y_cells{1};
    t1 = t_cells{1};

    power_series = y1(13, :);
    time_for_power = t1;
    figure;
    for i = 1:num_plots
        subplot(num_plots, 1, i);
        hold on;

        vars = together{i};
        num_vars = length(vars);

        if num_vars > numel(line_styles)
            warning('Más variables que estilos de línea disponibles; algunos estilos se repetirán.');
        end

        for s = 1:num_subjects
            y = y_cells{s};
            t = t_cells{s};
            subj_id = subject_ids(s);

            % Obtener máscara (última columna)
            mask = y(end - 1, :);  % vector lógico por variable (una fila por var)
            mask = logical(mask);  % asegúrate de que sea lógico

            for v = 1:num_vars
                var_idx = vars(v);
                style_idx = mod(v-1, length(line_styles)) + 1;
                style = line_styles{style_idx};

                legend_name = sprintf('Subject %d', subj_id);

                % Datos originales
                data = y(var_idx, 1:end);
                time = t;

                % Aplicar máscara (pon NaN donde mask = 0)
                %valid_mask = logical(mask(var_idx));
                if with_mask
                    data(~mask) = NaN;   
                end

                plot(time, data', ...
                     'Color', base_colores(s,:), ...
                     'LineStyle', style, ...
                     'LineWidth', 1.2, ...
                     'DisplayName', legend_name);
            end

            

            % Línea vertical en tiempo de VT1 (sin leyenda)
            if ~isempty(vt1_times) && ~isempty(vt1_times{s}) && ~isnan(vt1_times{s})
                t_vt1 = vt1_times{s};

                xline(t_vt1, '--', num2str(subj_id), ...
                      'Color', base_colores(s,:), ...
                      'LabelOrientation', 'horizontal', ...
                      'Alpha', 0.8, ...
                      'LineWidth', 1.5, ...
                      'HandleVisibility', 'off');
            end

            legend('show', 'Location', 'northwest');
        end

        yyaxis right
        ax = gca;
        ax.YColor = 'k';  % Cambia el color del eje derecho a negro
        plot(time_for_power, power_series, 'Color', [0.3 0.3 0.3], 'DisplayName', "power-ladder");
        ylabel('W')  % o la unidad que corresponda
        yyaxis left  % volver al eje original para futuros plots


        hold off;
        legend('show');
        grid on;

        if isfinite(tmax)
            xlim([0 tmax]);
        end

        nombre_vars = nombres(vars);
        unidades = strings(size(vars));
        for k = 1:length(vars)
            try
                unidades(k) = find_unit(T_units, nombre_vars{k});
            catch
                unidades(k) = "";
            end
        end

        if all(unidades == unidades(1)) && unidades(1) ~= ""
            ylabel(unidades(1));
        else
            ylabel('Valor');
        end

        title(strjoin(nombre_vars, ', '));
        xlabel('Tiempo [min]');
    end
    
end


function show_old_simulations_control(subject_ids, hypoxia_files, normoxia_files, normalize, second_axis_autonomic, second_axis_local)
% SHOW_OLD_SIMULATIONS_CONTROL procesa simulaciones en hipoxia y normoxia,
% aplica un filtro pasabajo y genera subplots comparativos (original y normalizado opcionalmente).
%
% Uso:
%   show_old_simulations_control(subject_ids, hypoxia_files, normoxia_files, normalize, second_axis_autonomic, second_axis_local)
%
% Ejemplo:
%   show_old_simulations_control([5], ...
%       {'3300_sec_hipoxia-23-07-2025.mat'}, ...
%       {'3300_sec_normoxia-23-07-2025.mat'}, ...
%       true, {'vO2'}, {'VCO2', 'otra_var'})

    if nargin < 4, normalize = false; end
    if nargin < 5, second_axis_autonomic = {'fv_s'}; end
    %if nargin < 5, second_axis_autonomic = {'aO2', 'vO2'}; end
    if nargin < 6, second_axis_local = {''}; end
    
    % Convertir a cell si es char
    if ischar(second_axis_autonomic)
        second_axis_autonomic = {second_axis_autonomic};
    end
    if ischar(second_axis_local)
        second_axis_local = {second_axis_local};
    end

    base_path = "../Simulations/only_simulation/";
    x_limits = [0.55, 0.8]; % límites comunes para VO2

    for i = 1:numel(subject_ids)
        subj = subject_ids(i);
        fprintf('Procesando sujeto %d...\n', subj);

        % === Construir rutas completas ===
        hypoxia_path = fullfile(base_path, num2str(subj), hypoxia_files{i});
        normoxia_path = fullfile(base_path, num2str(subj), normoxia_files{i});

        % === Cargar datos ===
        load_hypoxia = load(hypoxia_path);
        load_normoxia = load(normoxia_path);

        struct_hypoxia = load_hypoxia.struct_vars;
        struct_normoxia = load_normoxia.struct_vars;
        t_hypoxia = load_hypoxia.t;
        t_normoxia = load_normoxia.t;

        % === Filtro pasabajo ===
        fields_autonomic = {'fh_s', 'fp_s', 'fv'};%{'PAO2'};
        fields_local = {'xO2_b', 'xO2_h', 'xO2_rm', 'xO2_am'};%{'PACO2'};
        cutoff = 0.2;  order = 4;

        struct_filt_hypoxia = lowpass_struct(struct_hypoxia, t_hypoxia, cutoff, order, [fields_autonomic, fields_local]);
        struct_filt_normoxia = lowpass_struct(struct_normoxia, t_normoxia, cutoff, order, [fields_autonomic, fields_local]);

        % === Calcular min y max conjuntos para normalización ===
        [struct_normoxia_norm, struct_hypoxia_norm] = normalize_pair(struct_filt_normoxia, struct_filt_hypoxia, ...
            [fields_autonomic, fields_local]);

        % === Calcular límites verticales conjuntos (hipoxia y normoxia)
        [auto_ylim, local_ylim, second_ylim_autonomic, second_ylim_local] = calculate_common_y_limits(...
            struct_filt_normoxia, struct_filt_hypoxia, fields_autonomic, fields_local, ...
            second_axis_autonomic, second_axis_local);

        % === Crear figura para el sujeto ===
        figure('Name', sprintf('Sujeto %d', subj), 'Color', 'w');
        if normalize
            tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
        else
            tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
        end
        

        % === Fila 1: respuestas originales ===
        nexttile(1);
        plot_autonomic(struct_filt_normoxia, false, second_axis_autonomic, auto_ylim, second_ylim_autonomic); 
        title('A. Normoxia - Autonomic');
        xlim(x_limits); 

        nexttile(2);
        plot_autonomic(struct_filt_hypoxia, false, second_axis_autonomic, auto_ylim, second_ylim_autonomic); 
        title('B. Hypoxia - Autonomic');
        xlim(x_limits); 

        nexttile(3);
        plot_local(struct_filt_normoxia, false, second_axis_local, local_ylim, second_ylim_local); 
        title('C. Normoxia - Local');
        xlim(x_limits);

        nexttile(4);
        plot_local(struct_filt_hypoxia, false, second_axis_local, local_ylim, second_ylim_local); 
        title('D. Hypoxia - Local');
        xlim(x_limits);

        % === Si se pidió normalización, agregar segunda fila ===
        if normalize
            nexttile(5);
            plot_autonomic(struct_normoxia_norm, true, second_axis_autonomic);
            title('E. Normoxia - Autonomic (normalized)');
            xlim(x_limits); ylim([0 1]);

            nexttile(6);
            plot_autonomic(struct_hypoxia_norm, true, second_axis_autonomic);
            title('F. Hypoxia - Autonomic (normalized)');
            xlim(x_limits); ylim([0 1]);

            nexttile(7);
            plot_local(struct_normoxia_norm, true, second_axis_local);
            title('G. Normoxia - Local (normalized)');
            xlim(x_limits); ylim([0 1]);

            nexttile(8);
            plot_local(struct_hypoxia_norm, true, second_axis_local);
            title('H. Hypoxia - Local (normalized)');
            xlim(x_limits); ylim([0 1]);
        end
    end
end

% -------------------------------------------------------------------------
function [structA_norm, structB_norm] = normalize_pair(structA, structB, field_list)
    structA_norm = structA; structB_norm = structB;
    for i = 1:numel(field_list)
        f = field_list{i};
        if isfield(structA, f) && isfield(structB, f)
            valsA = structA.(f); valsB = structB.(f);
            global_min = min([valsA(:); valsB(:)]);
            global_max = max([valsA(:); valsB(:)]);
            if global_max > global_min
                structA_norm.(f) = (valsA - global_min) / (global_max - global_min);
                structB_norm.(f) = (valsB - global_min) / (global_max - global_min);
            end
        end
    end
end

% -------------------------------------------------------------------------
function [auto_ylim, local_ylim, second_ylim_autonomic, second_ylim_local] = calculate_common_y_limits(...
    structA, structB, fields_auto, fields_local, second_axis_autonomic, second_axis_local)
    
    all_auto = []; 
    all_local = [];
    all_second_autonomic = [];
    all_second_local = [];

    maskA = logical((structA.MRtO2 > 0.58) .* (structA.MRtO2 < 0.8));
    maskB = logical((structB.MRtO2 > 0.58) .* (structB.MRtO2 < 0.8));

    % --- Campos autonómicos
    for f = fields_auto
        if isfield(structA, f{1}), all_auto = [all_auto; structA.(f{1})(maskA)]; end
        if isfield(structB, f{1}), all_auto = [all_auto; structB.(f{1})(maskB)]; end
    end

    % --- Campos locales
    for f = fields_local
        if isfield(structA, f{1}), all_local = [all_local; structA.(f{1})(maskA)]; end
        if isfield(structB, f{1}), all_local = [all_local; structB.(f{1})(maskB)]; end
    end

    % --- Variables del segundo eje autonómico
    if nargin > 4 && ~isempty(second_axis_autonomic)
        for i = 1:numel(second_axis_autonomic)
            var = second_axis_autonomic{i};
            if isfield(structA, var)
                all_second_autonomic = [all_second_autonomic; structA.(var)(maskA)];
            end
            if isfield(structB, var)
                all_second_autonomic = [all_second_autonomic; structB.(var)(maskB)];
            end
        end
    end

    % --- Variables del segundo eje local
    if nargin > 5 && ~isempty(second_axis_local)
        for i = 1:numel(second_axis_local)
            var = second_axis_local{i};
            if isfield(structA, var)
                all_second_local = [all_second_local; structA.(var)(maskA)];
            end
            if isfield(structB, var)
                all_second_local = [all_second_local; structB.(var)(maskB)];
            end
        end
    end

    % --- Cálculo de límites
    auto_ylim = [min(all_auto)*0.9, max(all_auto)*1.1];
    local_ylim = [min(all_local)*0.9, max(all_local)*1.1];

    if ~isempty(all_second_autonomic)
        second_ylim_autonomic = [min(all_second_autonomic)*0.9, max(all_second_autonomic)*1.1];
    else
        second_ylim_autonomic = [];
    end

    if ~isempty(all_second_local)
        second_ylim_local = [min(all_second_local)*0.9, max(all_second_local)*1.1];
    else
        second_ylim_local = [];
    end
end


function plot_autonomic(struct_filt, normalize, second_axis_vars, auto_ylim, second_ylim)
    if nargin < 2, normalize = false; end
    if nargin < 3, second_axis_vars = {}; end
    if nargin < 4, auto_ylim = []; end
    if nargin < 5, second_ylim = []; end
    
    % Convertir a cell si es char
    if ischar(second_axis_vars)
        second_axis_vars = {second_axis_vars};
    end

    
    

    hold on;

    vars = {'fh_s', 'fp_s', 'fv'};
    
    colors = [
        0.9290 0.6940 0.1250;   % amarillo
        0.0000 0.4470 0.7410;   % azul
        0.4940 0.1840 0.5560    % morado
    ];
    
    plot_handles = [];
    legend_labels = {};
    
    % --- Eje izquierdo (autonómico)
    yyaxis left
    for i = 1:numel(vars)
        if isfield(struct_filt, vars{i})
            y = struct_filt.(vars{i});
            
            % Selección de estilo
            if strcmp(vars{i}, 'fv')
                line_style = '--';   % fv punteado
            else
                line_style = '-';    % el resto sólido
            end
            
            p = plot(struct_filt.MRtO2, y, line_style, ...
                     'LineWidth', 1.5, ...
                     'Color', colors(i,:));
            
            plot_handles(end+1) = p;
            legend_labels{end+1} = vars{i};
        end
    end

    if ~isempty(auto_ylim)
        ylim(auto_ylim);
    end
    if normalize
        ylabel('Normalized response');
    else
        ylabel('Autonomic Response [spike/s]');
    end

    % --- Eje derecho (múltiples variables secundarias)
    if ~isempty(second_axis_vars)
        yyaxis right
        color_orange = [0.85, 0.33, 0.10]; % naranjo
        %color_purple = [0.4940, 0.1840, 0.5560];
        line_styles = {'-', '--', ':', '-.'}; % estilos de línea para distinguir
        
        ylabel_strs = {};
        for i = 1:numel(second_axis_vars)
            var = second_axis_vars{i};
            if isfield(struct_filt, var)
                y2 = struct_filt.(var);
                style_idx = mod(i-1, length(line_styles)) + 1;
                p2 = plot(struct_filt.MRtO2, y2, line_styles{style_idx}, 'LineWidth', 1.5, 'Color', color_orange);
                plot_handles(end+1) = p2;
                legend_labels{end+1} = var;
                ylabel_strs{end+1} = var;
            end
        end
        
        if ~isempty(second_ylim)
            ylim(second_ylim);
        end
        
        if ~isempty(ylabel_strs)
            ylabel_combined = strjoin(cellfun(@(x) strrep(x, '_', '\_'), ylabel_strs, 'UniformOutput', false), ', ');
            ylabel(sprintf('%s [spike/s]', ylabel_combined));
        end
    end

    % --- Eje X
    xlabel('VO2 (l/min)');

    % --- Leyenda combinada
    legend(plot_handles, legend_labels, 'Location', 'best');

    hold off;
end


function plot_local(struct_filt, normalize, second_axis_vars, local_ylim, second_ylim)
    % Dibuja las variables locales, con opción a normalización y segundo eje
    if nargin < 2, normalize = false; end
    if nargin < 3, second_axis_vars = {}; end
    if nargin < 4, local_ylim = []; end
    if nargin < 5, second_ylim = []; end
    
    % Convertir a cell si es char
    %if ischar(second_axis_vars)
    %    second_axis_vars = {second_axis_vars};
    %end
    
    hold on;
    %yyaxis left
        
    vars = {'xO2_b', 'xO2_h', 'xO2_rm', 'xO2_am'};
    colors = lines(numel(vars));
    plot_handles = [];
    legend_labels = {};

    % --- Eje izquierdo (variables locales)
    %yyaxis left
    for i = 1:numel(vars)
        v = vars{i};
        if isfield(struct_filt, v)
            y = struct_filt.(v);
            p = plot(struct_filt.MRtO2, y, '-', 'LineWidth', 1.5, 'Color', colors(i,:));  
            
            plot_handles(end+1) = p;
            legend_labels{end+1} = v;
        end
    end
    
    if ~isempty(local_ylim)
        ylim(local_ylim);
    end

    if normalize
        ylabel('Normalized response');
    else
        ylabel('Local Response');
    end

    % --- Eje derecho (múltiples variables secundarias)
    % if ~isempty(second_axis_vars)
    %     yyaxis right
    %     color_orange = [0.85, 0.33, 0.10]; % naranjo
    %     line_styles = {'-', '--', ':', '-.'}; % estilos de línea para distinguir
    % 
    %     ylabel_strs = {};
    %     % for i = 1:numel(second_axis_vars)
    %     %     var = second_axis_vars{i};
    %     %     if isfield(struct_filt, var)
    %     %         y2 = struct_filt.(var);
    %     %         style_idx = mod(i-1, length(line_styles)) + 1;
    %     %         p2 = plot(struct_filt.MRtO2, y2, line_styles{style_idx}, 'LineWidth', 1.5, 'Color', color_orange);
    %     %         plot_handles(end+1) = p2;
    %     %         legend_labels{end+1} = var;
    %     %         ylabel_strs{end+1} = var;
    %     %     end
    %     % end
    % 
    %     if ~isempty(second_ylim)
    %         ylim(second_ylim);
    %     end
    % 
    %     if ~isempty(ylabel_strs)
    %         ylabel_combined = strjoin(cellfun(@(x) strrep(x, '_', '\_'), ylabel_strs, 'UniformOutput', false), ', ');
    %         ylabel(sprintf('%s', ylabel_combined));
    %     end
    % end

    % --- Eje X
    xlabel('VO2 (l/min)');

    % --- Leyenda combinada
    legend(plot_handles, legend_labels, 'Location', 'best');

    hold off;
end

function plot_J_iterations()

    mainFig = figure; % figura principal
    for i = 1:4
        % Ejecuta la función que genera la figura
        ForwardFittingModel('hard:only-loading', i);

        % Captura la última figura creada por la función
        h = gcf;

        % Busca los ejes (axes) generados
        ax = findall(h, 'type', 'axes');

        if isempty(ax)
            warning('No se encontraron ejes en la figura %d.', i);
            close(h);
            continue;
        end

        % Crea el subplot destino en la figura principal
        subplot(4,1,i);
        hold on;

        % Copia todos los hijos de los ejes (líneas, textos, etc.)
        for k = 1:numel(ax)
            newAx = gca;
            % Copiamos todos los objetos gráficos del eje original al subplot
            copyobj(allchild(ax(k)), newAx);
        end

        % Copia etiquetas y título (opcional)
        title(newAx, get(get(ax(1), 'Title'), 'String'));
        xlabel(newAx, get(get(ax(1), 'XLabel'), 'String'));
        ylabel(newAx, get(get(ax(1), 'YLabel'), 'String'));

        hold off;

        % Cierra la figura temporal
        close(h);
    end
end

function plot_topCV_params_box(csvFile, selectedParams)
% plot_topCV_params_box - Boxplots for the 10 parameters with highest CV
%
% USAGE:
%   plot_topCV_params_box('estimated_parameters.csv')
%   plot_topCV_params_box('estimated_parameters.csv', {'go2e','go2s'})

    % -----------------------------
    % 1. Load CSV
    % -----------------------------
    if ~isfile(csvFile)
        error('File not found: %s', csvFile);
    end
    T = readtable(csvFile);

    expectedCols = {'Parameters','Units','S1','S2','S3','S4','Mean','StdDv','StdDv_Mean','OriginalValue'};
    if ~all(ismember(expectedCols, T.Properties.VariableNames))
        missing = setdiff(expectedCols, T.Properties.VariableNames);
        error('Missing required columns in CSV: %s', strjoin(missing, ', '));
    end

    % -----------------------------
    % 2. If user selected specific params
    % -----------------------------
    if nargin >= 2 && ~isempty(selectedParams)
        mask = ismember(T.Parameters, selectedParams);
        T = T(mask, :);
        if isempty(T)
            error('None of the specified parameters were found in the table.');
        end
    end

    % -----------------------------
    % 3. Sort by CV and take top 10
    % -----------------------------
    T = sortrows(T, 'StdDv_Mean', 'descend'); % highest CV first
    maxParams = min(10, height(T));
    T = T(1:maxParams, :); % keep only top 10 (or fewer if <10)

    nParams = height(T);
    if nParams == 0
        error('No parameters to plot.');
    end

    % -----------------------------
    % 4. Prepare data
    % -----------------------------
    paramNames = string(T.Parameters);
    units      = T.Units;
    dataMat    = [T.S1, T.S2, T.S3, T.S4];
    meanVals   = T.Mean;
    origVals   = T.OriginalValue;
    CV         = T.StdDv_Mean;

    subjectLabels = {'S1','S2','S3','S4'};
    nSubjects = numel(subjectLabels);
    subjColors = lines(nSubjects);

    % -----------------------------
    % 5. Layout
    % -----------------------------
    %nCols = ceil(sqrt(nParams * 1.5));
    %nRows = ceil(nParams / nCols);

    nRows = 2;
    nCols = 5;

    figW = min(1600, 300*nCols);
    figH = min(1200, 240*nRows);
    figure('Color','w','Position',[100 100 figW figH]);

    t = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % -----------------------------
    % 6. Plot each parameter
    % -----------------------------
    for i = 1:nParams
        ax = nexttile(t);
        vals = dataMat(i, :);

        % Boxplot
        try
            boxplot(ax, vals, 'Widths', 0.45, 'Symbol', '');
        catch
            patch(ax, [0.7 1.3 1.3 0.7], [min(vals) min(vals) max(vals) max(vals)], ...
                0.9*[1 1 1], 'EdgeColor','none');
        end
        hold(ax, 'on');

        % Scatter for each subject
        jitterWidth = 0.18;
        for s = 1:nSubjects
            x = 1 + (rand-0.5)*jitterWidth;
            scatter(ax, x, vals(s), 60, subjColors(s,:), 'filled', 'MarkerEdgeColor','k');
        end

        % Y limits + padding
        yMin = min(vals);
        yMax = max(vals);
        yPad = max(0.05*(yMax - yMin), 0.01*abs(yMax)+eps);
        ylim(ax, [yMin-yPad, yMax+yPad]);

        % Mean line
        plot(ax, [0.6 1.4], [meanVals(i) meanVals(i)], ':k', 'LineWidth', 1.2);

        % ---- Nominal reference (clipped if outside) ----
        nom = origVals(i);
        clipped_nom = nom;
        outside = false;

        if nom > yMax
            clipped_nom = yMax;
            outside = true;
        elseif nom < yMin
            clipped_nom = yMin;
            outside = true;
        end

        if outside
            plot(ax, [0.6 1.4], [clipped_nom clipped_nom], '-m', 'LineWidth', 2.5);
            text(ax, 1.42, clipped_nom, sprintf('%.2g', nom), 'Color','m', ...
                'FontSize',11, 'FontWeight','bold', 'VerticalAlignment','middle');
        else
            plot(ax, [0.6 1.4], [clipped_nom clipped_nom], '-r', 'LineWidth', 1.5);
        end

        % Title includes CV
        title(ax, sprintf('%s (CV = %.2f)', paramNames(i), CV(i)), ...
            'Interpreter','none','FontWeight','bold','FontSize',10);

        ylabel(ax, units{i}, 'Interpreter','none');

        % Remove x-axis labels
        set(ax, 'XTick', 1, 'XTickLabel', {''});

        hold(ax, 'off');
    end

    % -----------------------------
    % 7. Global legend
    % -----------------------------
    % ax_legend = axes('Position',[0.35 0.02 0.3 0.05],'Visible','off');
    % hold(ax_legend, 'on');
    % for s = 1:nSubjects
    %     scatter(ax_legend, s, 1, 80, subjColors(s,:), 'filled','MarkerEdgeColor','k');
    % end
    % hold(ax_legend, 'off');

    lg = legend(subjectLabels, ...
        'Orientation','horizontal', 'Location','southoutside', 'Box','off');
    lg.FontSize = 14;
end

function plot_param_groups_box(csvFile, paramGroups)
% Plot boxplots for groups of parameters, each group in its own subplot (tiledlayout).
% Ensures Y-limits include both data and nominal values and shows nominal markers.

    if ~isfile(csvFile)
        error('File not found: %s', csvFile);
    end

    T = readtable(csvFile);

    expectedCols = {'Parameters','Units','S1','S2','S3','S4','Mean','StdDv','StdDv_Mean','OriginalValue'};
    if ~all(ismember(expectedCols, T.Properties.VariableNames))
        missing = setdiff(expectedCols, T.Properties.VariableNames);
        error('Missing required columns: %s', strjoin(missing, ', '));
    end

    subjectLabels = {'S1','S2','S3','S4'};
    subjColors = lines(numel(subjectLabels));

    % --- tiledlayout ---
    figure('Color','w');
    nRows = 2;
    nCols = 3;
    t = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    % Handles para leyenda global
    globalSubjectHandles = gobjects(1,4);
    globalNominalHandle = gobjects(1,1);
    legendCaptured = false;

    for g = 1:numel(paramGroups)
        group = paramGroups{g};

        mask = ismember(T.Parameters, group);
        Tg = T(mask, :);

        if isempty(Tg)
            warning('Group %d has no matching parameters.', g);
            continue;
        end

        % Keep order
        [~, idx] = ismember(group, Tg.Parameters);
        idx(idx==0) = [];
        Tg = Tg(idx, :);

        paramNames = string(Tg.Parameters);
        units = Tg.Units;
        dataMat = [Tg.S1 Tg.S2 Tg.S3 Tg.S4];
        meanVals = Tg.Mean;
        origVals = Tg.OriginalValue;

        nParams = height(Tg);

        % --- global Y limits ---
        finiteVals = dataMat(isfinite(dataMat));
        finiteNom = origVals(isfinite(origVals));

        globalMin = min([finiteVals(:); finiteNom(:)]);
        globalMax = max([finiteVals(:); finiteNom(:)]);
        delta = 0.04 * (globalMax - globalMin + eps);
        yLimits = [globalMin - delta, globalMax + delta];

        ax = nexttile(t);
        hold(ax,'on');

        % --- prepare data for boxplot ---
        dataMat_t = dataMat';
        y = dataMat_t(:);
        paramIndex = reshape(repmat(1:nParams, 4, 1), [], 1);

        % --- boxplot ---
        boxplot(ax, y, paramIndex, 'Widths', 0.55, 'Symbol', '');

        % --- scatter subject points ---
        jitterWidth = 0.18;
        localSubjectHandles = gobjects(1,4);

        for p = 1:nParams
            for s = 1:4
                jitterX = p + (rand-0.5)*jitterWidth;
                h = scatter(ax, jitterX, dataMat(p,s), 55, subjColors(s,:), ...
                        'filled', 'MarkerEdgeColor','k');
                if ~legendCaptured
                    localSubjectHandles(s) = h;
                end
            end
        end

        % --- mean lines ---
        for p = 1:nParams
            plot(ax, [p-0.25 p+0.25], [meanVals(p) meanVals(p)], ':k', 'LineWidth', 1.2);
        end

        % --- nominal values ---
        localNomHandle = gobjects(1,1);
        for p = 1:nParams
            nom = origVals(p);
            if ~isfinite(nom), continue; end

            if nom >= yLimits(1) && nom <= yLimits(2)
                hL = plot(ax, [p-0.25 p+0.25], [nom nom], '-', 'Color',[0.8 0 0.8], 'LineWidth',1.8);
                plot(ax, p, nom, 'd', 'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0 0.8], 'MarkerSize',7);
            else
                clipped = max(min(nom, yLimits(2)), yLimits(1));
                hL = plot(ax, [p-0.12 p+0.12], [clipped clipped], '-', 'Color',[0.8 0 0.8], 'LineWidth',1.6);
                text(ax, p+0.08, clipped, sprintf('%.3g', nom), 'Color',[0.8 0 0.8], ...
                    'FontSize',9,'FontWeight','bold','Interpreter','none');
            end

            if ~legendCaptured
                localNomHandle = hL;
            end
        end

        % --- Store handles for global legend (first subplot only) ---
        if ~legendCaptured
            globalSubjectHandles = localSubjectHandles;
            globalNominalHandle = localNomHandle;
            legendCaptured = true;
        end

        % --- axes formatting ---
        xlim(ax, [0.5, nParams+0.5]);
        ylim(ax, yLimits);

        xticks(ax, 1:nParams);
        xticklabels(ax, paramNames);
        xtickangle(ax, 35);

        if all(strcmp(units, units{1}))
            ylabel(ax, units{1}, 'Interpreter','none');
        else
            ylabel(ax, 'Value', 'Interpreter','none');
        end

        title(ax, sprintf('Parameter Group %d', g), 'FontSize', 14, 'FontWeight','bold');

        % NO legend here
        grid(ax,'on');
        hold(ax,'off');
    end

    % --- GLOBAL LEGEND BELOW ---
    lg = legend([globalSubjectHandles globalNominalHandle], ...
        [subjectLabels {'Nominal'}], ...
        'Orientation','horizontal', 'Box','off');
    
    lg.Layout.Tile = 'south';
end















function plot_resistances(file_normoxia, file_hipoxia)
    %% ====== CARGA DE DATOS ======
    lfN = load(file_normoxia);
    lfH = load(file_hipoxia);

    svN = lowpass_struct(lfN.struct_vars_normoxia, lfN.t_normoxia, 0.2, 4, ...
        {'MRtO2','xO2_b','xCO2_b','DTheta_R_am_p_n','xO2_am','x_met','DTheta_R_rm_p_n','xCO2_rm','xO2_rm','xCO2_h','xO2_h','DTheta_R_e_p','xO2_e','DTheta_R_s_p','xO2_s','pd'});
    svH = lowpass_struct(lfH.struct_vars_hipoxia, lfH.t_hipoxia, 0.2, 4, ...
        {'MRtO2','xO2_b','xCO2_b','DTheta_R_am_p_n','xO2_am','x_met','DTheta_R_rm_p_n','xCO2_rm','xO2_rm','xCO2_h','xO2_h','DTheta_R_e_p','xO2_e','DTheta_R_s_p','xO2_s','pd'});

    % Eje X basado en MRtO2
    xN = svN.MRtO2;
    xH = svH.MRtO2;

    %% ====== NORMOXIA ======
    RbpN  = 6.57 ./ (1 + svN.xO2_b + svN.xCO2_b);
    RampN = (svN.DTheta_R_am_p_n + 4.36) ./ (1 + svN.xO2_am + svN.x_met);
    RrmpN = (svN.DTheta_R_rm_p_n + 4.36) .* (1 + svN.xCO2_rm) ./ (1 + svN.xO2_rm);
    RhpN  = 22.27 * (1 + svN.xCO2_h) ./ (1 + svN.xO2_h);
    RepN  = (svN.DTheta_R_e_p + 2.87) ./ (1 + svN.xO2_e);
    RspN  = (svN.DTheta_R_s_p + 2.31) ./ (1 + svN.xO2_s);

    RtotalN = (1./RbpN + 1./RampN + 1./RrmpN + 1./RhpN + 1./RepN + 1./RspN).^(-1);
    Ram_sympN = svN.DTheta_R_am_p_n + 4.36;
    local_effectN = 1 ./ (1 + svN.xO2_am + svN.x_met);

    %% ====== HIPOXIA ======
    RbpH  = 6.57 ./ (1 + svH.xO2_b + svH.xCO2_b);
    RampH = (svH.DTheta_R_am_p_n + 4.36) ./ (1 + svH.xO2_am + svH.x_met);
    RrmpH = (svH.DTheta_R_rm_p_n + 4.36) .* (1 + svH.xCO2_rm) ./ (1 + svH.xO2_rm);
    RhpH  = 22.27 * (1 + svH.xCO2_h) ./ (1 + svH.xO2_h);
    RepH  = (svH.DTheta_R_e_p + 2.87) ./ (1 + svH.xO2_e);
    RspH  = (svH.DTheta_R_s_p + 2.31) ./ (1 + svH.xO2_s);

    RtotalH = (1./RbpH + 1./RampH + 1./RrmpH + 1./RhpH + 1./RepH + 1./RspH).^(-1);
    Ram_sympH = svH.DTheta_R_am_p_n + 4.36;
    local_effectH = 1 ./ (1 + svH.xO2_am + svH.x_met);

    %% ====== FIGURA 2x4 ======
    figure('Name', 'Comparison: Normoxia vs Hypoxia', 'Color', 'w');
    xlim_all = [0.5, 1.2]; % rango común para el eje X
    lw = 2; % grosor común de líneas

    %% ---------- 1–2: Resistencias individuales ----------
    subplot(2,4,1);
    plot(xN, RbpN, 'LineWidth', lw); hold on;
    plot(xN, RampN, 'LineWidth', lw);
    plot(xN, RrmpN, 'LineWidth', lw);
    plot(xN, RhpN, 'LineWidth', lw);
    plot(xN, RepN, 'LineWidth', lw);
    plot(xN, RspN, 'LineWidth', lw);
    title('A. Individual Resistances (Normoxia)');
    ylabel('R (mmHgs/ml)');
    legend('b','am','rm','h','e','s', 'Location', 'best');
    grid on; xlim(xlim_all);
    xlabel('VO2 (l/min)');

    subplot(2,4,2);
    plot(xH, RbpH, 'LineWidth', lw); hold on;
    plot(xH, RampH, 'LineWidth', lw);
    plot(xH, RrmpH, 'LineWidth', lw);
    plot(xH, RhpH, 'LineWidth', lw);
    plot(xH, RepH, 'LineWidth', lw);
    plot(xH, RspH, 'LineWidth', lw);
    title('B. Individual Resistances (Hypoxia)');
    ylabel('R (mmHgs/ml)');
    grid on; xlim(xlim_all);
    xlabel('VO2 (l/min)');

    yl = [min([ylim(subplot(2,4,1)), ylim(subplot(2,4,2))]), ...
          max([ylim(subplot(2,4,1)), ylim(subplot(2,4,2))])];
    subplot(2,4,1); ylim(yl);
    subplot(2,4,2); ylim(yl);

    %% ---------- 3–4: Resistencia total ----------
    subplot(2,4,3);
    plot(xN, RtotalN, 'k', 'LineWidth', lw);
    yyaxis right; p1 = plot(xN, svN.pd, 'LineWidth', lw); ylabel('pd');
    yyaxis left; ylabel('R (mmHgs/ml)');
    title('C. Total R and PD (Normoxia)');
    grid on; xlim(xlim_all);
    xlabel('VO2 (l/min)');

    subplot(2,4,4);
    plot(xH, RtotalH, 'k', 'LineWidth', lw);
    yyaxis right; p2 = plot(xH, svH.pd, 'LineWidth', lw); ylabel('Diastolic Pressure');
    yyaxis left; ylabel('R (mmHgs/ml)');
    title('D. Total R and PD (Hypoxia)');
    grid on; xlim(xlim_all);
    xlabel('VO2 (l/min)');

    % Igualar límites de R y de pd
    ylR = [min([ylim(subplot(2,4,3)), ylim(subplot(2,4,4))]), ...
           max([ylim(subplot(2,4,3)), ylim(subplot(2,4,4))])];
    subplot(2,4,3); ylim(ylR);
    subplot(2,4,4); ylim(ylR);

    % Igualar límites de pd (yyaxis right)
    yyaxis(subplot(2,4,3), 'right'); yl_pdN = ylim;
    yyaxis(subplot(2,4,4), 'right'); yl_pdH = ylim;
    yl_pd = [min([yl_pdN, yl_pdH]), max([yl_pdN, yl_pdH])];
    yyaxis(subplot(2,4,3), 'right'); ylim(yl_pd);
    yyaxis(subplot(2,4,4), 'right'); ylim(yl_pd);

    %% ---------- 5–6: Simpático ----------
    subplot(2,4,5);
    plot(xN, Ram_sympN, 'r', 'LineWidth', lw);
    title('E. Sympathetic Ram (Normoxia)');
    ylabel('R (mmHgs/ml)');
    xlabel('VO2 (l/min)');
    grid on; xlim(xlim_all);

    subplot(2,4,6);
    plot(xH, Ram_sympH, 'r', 'LineWidth', lw);
    title('F. Sympathetic Ram (Hypoxia)');
    ylabel('R (mmHgs/ml)');
    xlabel('VO2 (l/min)');
    grid on; xlim(xlim_all);

    yl = [min([ylim(subplot(2,4,5)), ylim(subplot(2,4,6))]), ...
          max([ylim(subplot(2,4,5)), ylim(subplot(2,4,6))])];
    subplot(2,4,5); ylim(yl);
    subplot(2,4,6); ylim(yl);

    %% ---------- 7–8: Efecto local ----------
    subplot(2,4,7);
    plot(xN, local_effectN, 'b', 'LineWidth', lw);
    title('G. Local Effect (Normoxia) Ram');
    ylabel('Scaling factor');
    xlabel('VO2 (l/min)');
    grid on; xlim(xlim_all);

    subplot(2,4,8);
    plot(xH, local_effectH, 'b', 'LineWidth', lw);
    title('H. Local Effect Ram (Hypoxia)');
    ylabel('Scaling factor');
    xlabel('VO2 (l/min)');
    grid on; xlim(xlim_all);

    yl = [min([ylim(subplot(2,4,7)), ylim(subplot(2,4,8))]), ...
          max([ylim(subplot(2,4,7)), ylim(subplot(2,4,8))])];
    subplot(2,4,7); ylim(yl);
    subplot(2,4,8); ylim(yl);

    %% ---------- Título general ----------
    sgtitle('Comparison of Vascular Resistances under Normoxia vs Hypoxia');
end

function bode_damping_analysis(params)
    % ==========================================================
    % Bode magnitude plots for 3 damping conditions and 3 (C,L) sets
    % Each row = damping condition; each column = (C,L) combination
    % Shows f_n, f_z, and HR band (1–2 Hz)
    % Each subplot labeled A–I
    % ==========================================================

    if numel(params) ~= 6
        error('Input vector must have 6 elements: [C1 L1 C2 L2 C3 L3]');
    end

    % Extract parameter sets
    C_values = params(1:2:end);
    L_values = params(2:2:end);

    % Labels for damping types
    damping_labels = {'Underdamped', 'Critically damped', 'Overdamped'};

    % Frequency setup
    f = logspace(-2, 2, 2000);      % 0.01–100 Hz
    w = 2*pi*f;
    f_hr_band = [1, 2];             % HR band (1–2 Hz)

    % Letters for subplot labels
    letters = 'ABCDEFGHI';

    figure('Name','Bode Damping Analysis (3x3)','Color','w','Position',[100 100 1200 800]);

    for j = 1:3  % columns = parameter sets
        Csa = C_values(j);
        Lsa = L_values(j);

        % Compute critical resistance and its variants
        Rcrit = 2 * sqrt(Lsa / Csa);
        R_values = [Rcrit/2, Rcrit, Rcrit*1.7]; % under, critical, over

        % Natural frequency
        wn = 1 / sqrt(Lsa * Csa);
        f_n = wn / (2*pi);

        for i = 1:3  % rows = damping conditions
            idx = (i-1)*3 + j; % subplot index
            Rsa = R_values(i);
            f_z = Rsa / ((2*pi) * Lsa);

            % Transfer function P(s)/Qlv(s)
            num = [Lsa Rsa];
            den = [Csa*Lsa, Rsa*Csa, 1];
            sys = tf(num, den);

            [mag, ~] = bode(sys, w);
            mag = squeeze(20*log10(mag));

            % Subplot position
            subplot(3,3,idx)
            semilogx(f, mag, 'b', 'LineWidth', 1.4); hold on;

            % Vertical lines for key frequencies
            xline(f_n, '--r', 'f_n', 'LabelVerticalAlignment','bottom');
            xline(f_z, '--k', 'f_z', 'LabelVerticalAlignment','bottom');
            xline(f_hr_band(1), '--g', '1 Hz', 'LabelVerticalAlignment','bottom');
            xline(f_hr_band(2), '--g', '2 Hz', 'LabelVerticalAlignment','bottom');

            grid on; xlim([0.01 10]);
            xlabel('Frequency [Hz]');
            ylabel('|P/Q_{lv}| [dB]');

            % Title with current parameters
            title(sprintf('%s\nR = %.3f, C = %.3f, L = %.3f', ...
                damping_labels{i}, Rsa, Csa, Lsa), 'FontSize', 10);

            % Add subplot letter (A–I) in top-left corner (visible over plot)
            text(0.02, 0.92, sprintf('%s.', letters(idx)), ...
                 'Units','normalized', 'FontWeight','bold', ...
                 'FontSize',12, 'Color','k', ...
                 'BackgroundColor','w', 'Margin',2, ...
                 'EdgeColor','k', 'VerticalAlignment','top', ...
                 'HorizontalAlignment','left');
        end
    end

    % --- Display summary in console ---
    fprintf('\nParameter summary:\n');
    for j = 1:3
        Csa = C_values(j);
        Lsa = L_values(j);
        Rcrit = 2 * sqrt(Lsa / Csa);
        fprintf('\nCase %d: C = %.3f, L = %.3f\n', j, Csa, Lsa);
        fprintf('  Underdamped:    R = %.3f\n', Rcrit/2);
        fprintf('  Critical:       R = %.3f\n', Rcrit);
        fprintf('  Overdamped:     R = %.3f\n', Rcrit*1.7);
    end
end

function data = validation_section()

    subjects = [1, 4, 5, 6];

    base_path = 'C:\Users\Rafael\Desktop\universidad\magister\Cardio-respiratory-model-RK\codes\Simulations\simulation_after_fitting\';

    data = struct();

    for s = subjects
        folder = fullfile(base_path, num2str(s));

        % Archivos a cargar
        file_hipoxia  = fullfile(folder, '3300_sec_hipoxia-04-11-2025.mat');
        file_normoxia = fullfile(folder, '1200_sec_normoxia-04-11-2025.mat');

        % Cargar
        H = load(file_hipoxia);
        N = load(file_normoxia);

        subj = ['subj' num2str(s)];

        % ----------- HIPOXIA -----------
        if isfield(H, 'struct_vars_hipoxia')
            H.struct_vars_hipoxia = compute_resistances(H.struct_vars_hipoxia);
            data.(subj).hipoxia.struct_vars = H.struct_vars_hipoxia;            
        elseif isfield(H, 'struct_vars')
            H.struct_vars= compute_resistances(H.struct_vars);
            data.(subj).hipoxia.struct_vars = H.struct_vars;
        else
            error('No se encontró struct_vars_hipoxia ni struct_vars en %s', file_hipoxia)
        end

        if isfield(H, 't_hipoxia')
            data.(subj).hipoxia.t = H.t_hipoxia;
        elseif isfield(H, 't')
            data.(subj).hipoxia.t = H.t;
        else
            error('No se encontró t_hipoxia ni t en %s', file_hipoxia)
        end

        % ----------- NORMOXIA -----------
        if isfield(N, 'struct_vars_normoxia')
            N.struct_vars_normoxia = compute_resistances(N.struct_vars_normoxia);
            data.(subj).normoxia.struct_vars = N.struct_vars_normoxia;
        elseif isfield(N, 'struct_vars')
            N.struct_vars = compute_resistances(N.struct_vars);
            data.(subj).normoxia.struct_vars = N.struct_vars;
        else
            error('No se encontró struct_vars_normoxia ni struct_vars en %s', file_normoxia)
        end

        if isfield(N, 't_normoxia')
            data.(subj).normoxia.t = N.t_normoxia;
        elseif isfield(N, 't')
            data.(subj).normoxia.t = N.t;
        else
            error('No se encontró t_normoxia ni t en %s', file_normoxia)
        end

    end

    disp('Carga completada. Estructura "data" creada.');

    plot_variable_vs_variable(data, 'Qh', 'MRtO2', 0.55, 0.9, 'Coronary flow (ml/s)');
    
    
    %plot_3vars_hypoxia_normoxia(data, 0.55, 0.9)

end

function plot_variable_vs_variable(data, varY, varX, minX, maxX, units)
% plot_variable_vs_variable - Grafica cualquier variable Y vs X
% comparando hipoxia y normoxia para cada sujeto.
%
% Sintaxis:
%   plot_variable_vs_variable(data, 'dVE', 'MRtO2', 0.5)
%
% Inputs:
%   data  - estructura generada por validation_section()
%   varY  - string con el nombre de la variable en struct_vars para eje Y
%   varX  - string con el nombre de la variable en struct_vars para eje X
%   minX  - valor mínimo de X para filtrar los datos

subjects = fieldnames(data);  % Extrae nombres de sujetos
figure;
numSubjects = length(subjects);
letters = {'A', 'B', 'C', 'D'};
% Colores clásicos de MATLAB
color_hipoxia  = [0 0.4470 0.7410];   % azul
color_normoxia = [0.8500 0.3250 0.0980]; % naranjo




for i = 1:numSubjects
    subj = subjects{i};

    %data.(subj).hipoxia.struct_vars = lowpass_struct(data.(subj).hipoxia.struct_vars, data.(subj).hipoxia.t, 0.5, 4, [{varY}, {}]);
    %data.(subj).normoxia.struct_vars = lowpass_struct(data.(subj).normoxia.struct_vars, data.(subj).normoxia.t, 0.5, 4, [{varY}, {}]);
    
    % -------- Hipoxia --------
    if isfield(data.(subj).hipoxia.struct_vars, varY) && isfield(data.(subj).hipoxia.struct_vars, varX)
        Y_h = data.(subj).hipoxia.struct_vars.(varY);
        X_h = data.(subj).hipoxia.struct_vars.(varX);
        %PAO2_h = data.(subj).hipoxia.struct_vars.('PAO2');
        
        % Filtrar por minX
        idx_h = logical( (X_h >= minX) .* (X_h <= maxX));
        X_h = X_h(idx_h);
        Y_h = Y_h(idx_h);
        %PAO2_h = PAO2_h(idx_h);
    else
        warning('No se encontró %s o %s en hipoxia de %s', varY, varX, subj);
        continue
    end
    
    % -------- Normoxia --------
    if isfield(data.(subj).normoxia.struct_vars, varY) && isfield(data.(subj).normoxia.struct_vars, varX)
        Y_n = data.(subj).normoxia.struct_vars.(varY);
        X_n = data.(subj).normoxia.struct_vars.(varX);
        %PAO2_n = data.(subj).normoxia.struct_vars.('PAO2');
        
        idx_n = logical( (X_n >= minX) .* (X_n <= maxX));
        X_n = X_n(idx_n);
        Y_n = Y_n(idx_n);
        %PAO2_n = PAO2_n(idx_n);
    else
        warning('No se encontró %s o %s en normoxia de %s', varY, varX, subj);
        continue
    end
    
    % -------- Crear subplot --------
    subplot(2,2,i)
    plot(X_h, Y_h, 'Color', color_hipoxia, 'LineWidth', 1.5); hold on
    %plot(PAO2_h, Y_h, 'Color', color_hipoxia, 'LineWidth', 1.5); hold on
    plot(X_n, Y_n, 'Color', color_normoxia, 'LineWidth', 1.5);
    %plot(PAO2_n, Y_n, 'Color', color_normoxia, 'LineWidth', 1.5);
    
    ylabel(units); xlabel('VO2 (l/min)');
    name = sprintf('%s. S %d', letters{i}, i);
    title([name]);
    legend({'Hypoxia', 'Normoxia'}, 'Location', 'best');
    grid on;
end

%sgtitle(sprintf('Comparación %s vs %s por sujeto (Hipoxia vs Normoxia)', varY, varX));

end

function plot_3vars_hypoxia_normoxia(data, minX, maxX)

%vars  = {'fh_s', 'fp_s', 'fv_s'};
vars  = {'V_total_rv', 'V_total_lv'};
%units = {'fh_s (spike/s)', 'fp_s (spike/s)', 'fv_s (spike/s)'}; % edítalas si quieres
units = {'V_total_rv (ml)', 'V_total_lv (ml)'}; % edítalas si quieres
%letters = {'A', 'D', 'G', 'J', 'B', 'E', 'H', 'K', 'C', 'F', 'I', 'L'};
letters = {'A', 'E', 'B', 'F', 'C', 'G', 'D', 'H'};
subjects = fieldnames(data);
numSubjects = length(subjects);

if numSubjects ~= 4
    warning('Se esperaban 4 sujetos, se encontraron %d', numSubjects);
end

color_hipoxia  = [0 0.4470 0.7410];       % azul
color_normoxia = [0.8500 0.3250 0.0980];  % naranjo

figure;

for i = 1:numSubjects
    subj = subjects{i};

    for v = 1:length(vars)
        varY = vars{v};
        varX = 'MRtO2';  % como antes

        % ---- Hipoxia ----
        if isfield(data.(subj).hipoxia.struct_vars, varY) && isfield(data.(subj).hipoxia.struct_vars, varX)
            Y_h = data.(subj).hipoxia.struct_vars.(varY);
            X_h = data.(subj).hipoxia.struct_vars.(varX);
            idx_h = (X_h >= minX) & (X_h <= maxX);
            X_h = X_h(idx_h);
            Y_h = Y_h(idx_h);
        else
            warning('Falta %s o %s en hipoxia de %s', varY, varX, subj);
            continue
        end

        % ---- Normoxia ----
        if isfield(data.(subj).normoxia.struct_vars, varY) && isfield(data.(subj).normoxia.struct_vars, varX)
            Y_n = data.(subj).normoxia.struct_vars.(varY);
            X_n = data.(subj).normoxia.struct_vars.(varX);
            idx_n = (X_n >= minX) & (X_n <= maxX);
            X_n = X_n(idx_n);
            Y_n = Y_n(idx_n);
        else
            warning('Falta %s o %s en normoxia de %s', varY, varX, subj);
            continue
        end

        % ---- Ubicación en grilla 6x2 ----
        row = (i-1)*length(vars) + v;  % 1 → 12
        subplot(2,6,row)

        % ---- Plot ----
        plot(X_h, Y_h, 'Color', color_hipoxia,  'LineWidth', 1.5); hold on
        plot(X_n, Y_n, 'Color', color_normoxia, 'LineWidth', 1.5);

        ylabel(units{v});
        xlabel('VO₂ (L/min)');
        title(sprintf('%s. S %d — %s', letters{4*(v-1) + i}, i, varY));
        legend({'Hypoxia', 'Normoxia'}, 'Location', 'best');
        grid on;
    end
end

%sgtitle('Hypoxia vs Normoxia — 3 Variables per Subject')

end



function sv = compute_resistances(sv)
    sv.Rbp  = 6.57 ./ (1 + sv.xO2_b + sv.xCO2_b);
    sv.local_b = 1./ (1 + sv.xO2_b + sv.xCO2_b);   
    sv.Ramp = (sv.DTheta_R_am_p_n + 4.36) ./ (1 + sv.xO2_am + sv.x_met);
    sv.local_am = 1./ (1 + sv.xO2_am + 0*sv.x_met);    
    sv.Rrmp = (sv.DTheta_R_rm_p_n + 4.36) .* (1 + sv.xCO2_rm) ./ (1 + sv.xO2_rm);
    sv.local_rm = 1./ (1 + sv.xO2_rm);
    sv.Rhp  = 22.27 * (1 + sv.xCO2_h) ./ (1 + sv.xO2_h);
    sv.local_h = 1./ (1 + sv.xO2_h);
    sv.Rep  = (sv.DTheta_R_e_p + 2.87) ./ (1 + sv.xO2_e);
    sv.local_e = 1./ (1 + sv.xO2_e);
    sv.Rsp  = (sv.DTheta_R_s_p + 2.31) ./ (1 + sv.xO2_s);
    sv.local_s = 1./ (1 + sv.xO2_s);

    sv.local_total = sv.local_b.^-1 + sv.local_am.^-1 + sv.local_rm.^-1 + sv.local_h.^-1 + sv.local_e.^-1 + sv.local_s.^-1;
    

    sv.Rtotal = (1./sv.Rbp + 1./sv.Ramp + 1./sv.Rrmp + 1./sv.Rhp + 1./sv.Rep + 1./sv.Rsp).^(-1);
    sv.Ram_symp = sv.DTheta_R_am_p_n + 4.36;
    sv.Vtotal_bigOnes = sv.V_total_e_v +  sv.V_total_s_v;

    sv.HR = sv.HR * 60;
    sv.R_pp = 0.0894 * (1 + sv.xO2_p); 

    sv.Qh = sv.P_sp ./ sv.Rhp;
    
    %sv.local_effect = 1 ./ (1 + sv.xO2_am + sv.x_met);
end



































function unit = find_unit(table, var)
    row_index = strcmp(table.Variable, var);
    unit = table.MeasureUnit{row_index};    
end