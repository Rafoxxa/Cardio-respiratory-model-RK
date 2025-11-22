%% Nombre del PDF final
pdf_filename = 'Sujeto_1_analisis_presiones.pdf';

%% Hacer todas las figuras invisibles
%set(0, 'DefaultFigureVisible', 'off');

%% Ejecutar función que genera múltiples figuras
%for person = 1:3
%    ForwardFittingModel('only-plot', person, '04-10-2025', '06-10-2025');
%end

%% Guardar todas las figuras en un PDF
figs = findall(0,'Type','figure');
for k = 1:numel(figs)
    exportgraphics(figs(k), pdf_filename, 'Append', true);
    close(figs(k));
end

%% Restaurar visibilidad normal
%set(0, 'DefaultFigureVisible', 'on');
