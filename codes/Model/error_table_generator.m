%% === Configuración ===
subjectIDs = 1:4;
date1 = '28-10-2025';
date2 = '04-11-2025';

% Labels de las variables en orden
varLabels = {
    'PAO2','PACO2','dVE','VT','TI','Tresp','HR','PM','PS','PD'
};

% Matrices para almacenar resultados
Jnorm_all = zeros(length(subjectIDs), length(varLabels));
Jhypo_all = zeros(length(subjectIDs), length(varLabels));

%% === Correr modelo por sujeto ===
for i = 1:length(subjectIDs)
    s = subjectIDs(i);
    fprintf('Procesando sujeto %d...\n', s);

    out = ForwardFittingModel('only-compute-J', s, date1, date2);

    Jnorm_all(i, :) = out.J_normoxia(:)';
    Jhypo_all(i, :) = out.J_hipoxia(:)';
end

%% === Construir tablas ===
SubjectCol = strcat("S", string(subjectIDs')) ;

T_norm = array2table(Jnorm_all, 'VariableNames', varLabels);
T_norm = addvars(T_norm, SubjectCol, 'Before', 1);
T_norm.Properties.VariableNames{1} = 'Subject';

T_hypo = array2table(Jhypo_all, 'VariableNames', varLabels);
T_hypo = addvars(T_hypo, SubjectCol, 'Before', 1);
T_hypo.Properties.VariableNames{1} = 'Subject';

%% === Guardar a Excel ===
excelFile = 'J_subjects_matrix.xlsx';
writetable(T_norm, excelFile, 'Sheet', 'Normoxia', 'WriteMode','overwrite');
writetable(T_hypo, excelFile, 'Sheet', 'Hipoxia', 'WriteMode','overwrite');

fprintf('\n✅ Excel generado correctamente: %s\n', excelFile);
