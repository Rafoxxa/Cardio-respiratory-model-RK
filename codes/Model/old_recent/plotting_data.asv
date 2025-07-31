 data_path_1 = "../../data/CPETS_septiembre_2024_nuevoProtocolo/CPET-1_AG_cleaned.xlsx";
 data_path_2 = "../../data/CPETS_septiembre_2024_nuevoProtocolo/CPET-2_RK_cleaned.xlsx";
 data_path_3 = "../../data/CPETS_septiembre_2024_nuevoProtocolo/CPET-3_AM_cleaned.xlsx";

 table1 = readtable(data_path_1);
 table2 = readtable(data_path_2);
 table3 = readtable(data_path_3);

 plot_tables(table1, "Sujeto 1");
 figure;
 plot_tables(table2, "Sujeto 2");
 figure;
 plot_tables(table3, "Sujeto 3");

function plot_tables(table, sujeto)
 subplot(4,2,1);
 plot(table.t, table.VCO2);
 ylabel("VCO2");
 subplot(4,2,2);
 plot(table.t, table.BR);
 ylabel("BF");
 subplot(4,2,3);
 plot(table.t, table.VE);
 ylabel("VE");
 subplot(4,2,4);
 plot(table.t, table.VT);
 ylabel("VT");
 subplot(4,2,5);
 plot(table.t, table.Ti);
 ylabel("TI");
 subplot(4,2,6);
 plot(table.t, table.PetCO2);
 ylabel("PACO2");
 subplot(4,2,7);
 plot(table.t, table.PetO2);
 ylabel("PAO2");
 subplot(4,2,8);
 plot(table.t, table.Potencia);
 ylabel("Potencia");
 sgtitle(sujeto);
end