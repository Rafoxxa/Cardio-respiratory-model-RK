%%Plotting for document

%Plot for protocol visualization
tt = 1:44*60;
power = 0 * (tt < 2*60) + 1 * (tt >= 2*60) .* (tt < 6*60) + 25 * (tt >= 6*60) .* (tt < 11*60) + 50 * (tt >= 11*60) .* (tt < 16*60) + 75 * (tt >= 16*60) .* (tt < 21*60) + 100 * (tt >= 21*60) .* (tt < 26*60) + 125 * (tt >= 26*60) .* (tt < 31*60) + 100 * (tt >= 31*60) .* (tt < 32*60) + 75 * (tt >= 32*60) .* (tt < 33*60) + 50 * (tt >= 33*60) .* (tt < 34*60) + 25 * (tt >= 34*60) .* (tt < 35*60) + 1 * (tt >= 35*60) .* (tt < 39*60) + 0 * (tt >= 39*60) .* (tt < 44*60);

figure;
hold on;
area(tt, power);
xline_at = xline(1200, ':k', 'LineWidth', 1.5); % Dotted vertical line

grid on;
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

xlabel('Time (s)');
ylabel('Power (W)');

legend(xline_at, 'AT'); % Add legend for the dotted line

hold off;

%% Plot for AT detection
AT = 1.2;
VO2 = 0:0.05:3;
VCO2 = 0.8 * VO2 .* (VO2 < AT) + (1.5 * VO2 + 0.8*AT - 1.5*AT) .* (VO2 >= AT);

figure;
hold on;
plot(VO2, VCO2);
xline_at = xline(AT, ':k', 'LineWidth', 1.5); % Dotted vertical line

grid on;
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

xlabel("VO2'");
ylabel("VCO2'");

legend(xline_at, 'AT'); % Add legend for the dotted line

hold off;

%% Computations for AT
VO2 = zeros(4, 1000);
VCO2 = zeros(4, 1000);
times = zeros(4, 1000);
for patient = [1,4,5,6]
    %try
        [t, ~, VO2_, VCO2_, ~] = data_preprocessing(patient, "hipoxia", "exercise", 0, "ori");
        
        VO2(patient, 1:size(VO2_, 1)) = VO2_;
        VCO2(patient, 1:size(VO2_, 1)) = VCO2_;
        times(patient, 1:size(t, 1)) = t;
    %catch
    %    disp("error");







        %[~, ~, VO2_, VCO2_, ~] = data_preprocessing(patient, "hipoxia", "exercise", 0, "ori");
        
        %VO2(patient, 1:size(VO2_, 1)) = VO2_;
        %VCO2(patient, 1:size(VO2_, 1)) = VCO2_;
    %end
end

%The plot will consist in 4 curves, all in the VCO2 vs VO2 plot:
% 1. data points
% 2. Tendency curve from the first 30 points
% NOT THIS, TO COMPLICATE 3. Tendency curve from the last 30 points
% RER = 1 line, dividing plane


manually_AT = [1.07, 1.5, 1, 1.18, 1.36];
manually_AT_for_hipoxia = [1.6, 0, 0, 2.6, 2.44, 1.11];

figure; % Create a new figure
n = 1;
for i = [1, 4, 5, 6]%1:4
    subplot(2, 2, n); % Create a 3x2 grid and place plots in first 5 positions
    hold on; % Keep all plots on the same axes

    %Scatter plot of VO2 vs VCO2
    h1 = plot(VO2(i,:), VCO2(i,:), '.'); 

    % Line plot of VO2 vs itself
    h2 = plot(VO2(i,:), VO2(i,:)); 

    % Line plot with initial slope
    slope_ini = mean(VCO2(i,1:30) ./ VO2(i,1:30));
    h3 = plot(VO2(i,:), slope_ini * VO2(i,:), '-'); 

    %Vertical line at manually_AT(i)
    %xline_AT = xline(manually_AT(i), ':k', 'LineWidth', 1.5);

    %Add all legend entries at once
    %legend([h1, h2, h3, xline_AT], "data", "RER = 1", "first slope", "AT");
    legend([h1, h2], "data", "RER = 1");


    % Labels and title
    xlabel("VO2' l/min ");
    ylabel("VCO2' l/min ");
    title(['Sujeto ', num2str(i)]);
    grid on;

    
   
     % vo2 = VO2(i, VO2(i,:) > 0);
     % timess = times(i,VO2(i,:) > 0);
     % yline_AT = yline(manually_AT_for_hipoxia(i), ':k', 'LineWidth', 1.5);
     % plot(timess, vo2);
    n = n + 1;
    
end



% Adjust layout to avoid overlapping titles
sgtitle('VO2 vs VCO2 Subplots (NORMOXIA)');



