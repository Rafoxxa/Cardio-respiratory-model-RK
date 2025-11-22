
clear

%data_file = ['C:\Work\uc_lab\samples\050923LFmarcha01.csv']; % 20 subframes%
%data_file = ['C:\Work\uc_lab\samples\080923VBmarcha01.csv']; % 20 subframes
%data_file = ['C:\Users\Gonzalo Lab\OneDrive - Universidad Católica de Chile\Documentos\UC\Lab\Fondecyt\Resultados Preliminares\Datos preliminares poster chicago\Marcha Vicon\POST\290923PZmarcha02.csv']; % 20 subframes


%data_file = ['C:\Work\uc_lab\samples\290923PZmarcha04.csv']; % 20 subframes %%%%

%data_file = ['C:\Work\uc_lab\samples\150923RSmarcha02.csv']; % 20 subframes
%data_file = ['C:\Work\uc_lab\samples\150923RSmarcha09.csv']; % 20 subframes
%data_file = ['C:\Work\uc_lab\samples\150923RSmarcha19.csv']; % 20 subframes
%data_file = ['C:\Work\uc_lab\samples\250723MM.csv']; %% This Data file have 10 Subframe per frame, remember Change it in parameters!
%data_file = '../../../../../090924RM MARCHA 08.csv';
data_file = '../../../../../050724MF marcha 04.csv';
%% Parameters
% Sampling data, some data files have 10 subframe per frame
frame_freq_hz = 100; % Frame rate, 
subframe_per_frame = 20; % Numbers of subframe per frame

% Data start (Sample 0 0 Row number)
sample_0_row = 36;

% Ground Forces in axis Z Column for each foot:
gf_1_col = 5;
gf_2_col = 14;
% Center of Mass in Axis Z column (This is the one used in the original code): 
comz_col = 20; 

%% Getting only the interval of interest
scan_interval = 1000;
max_scan_interval = 20;
force_treshold = 15; % Newtons threshold
for i=1:max_scan_interval
    index_1 = sample_0_row+scan_interval*(i-1);
    index_2 = sample_0_row+scan_interval*i;
    dgf_1 = readmatrix(data_file,'Delimiter',{',',';'},'Range',[index_1 gf_1_col index_2 gf_1_col]); % 
    dgf_2 = readmatrix(data_file,'Delimiter',{',',';'},'Range',[index_1 gf_2_col index_2 gf_2_col]); % 
    gf_1_start_index = find(dgf_1 < -force_treshold, 1); % Trigger Threshold setted on 20 [N]
    gf_2_start_index = find(dgf_2 < -force_treshold, 1);
    if (isempty(gf_1_start_index)) && (isempty(gf_2_start_index))
        if i == max_scan_interval %Mean we reach the end of the scan without finding anything
            ME = MException('Stride interval not found');
            throw(ME)
        end
    else
        % Get interval range
        if isempty(gf_1_start_index)
           front_gf_start_index = gf_2_start_index + index_1;
           back_foot_col = gf_1_col; %gf_2 is front foot
           front_foot_col = gf_2_col;

        else
           front_gf_start_index = gf_1_start_index + index_1;
           back_foot_col = gf_2_col; %gf_1 is front foot
           front_foot_col = gf_1_col;
        end
        % Scan for the end of back_gf, 4 seconds of interval seems enough
        scan_interval_seconds = 4;
        stride_start_frame = readmatrix(data_file,'Delimiter',{',',';'},'Range',[front_gf_start_index 1 front_gf_start_index 1]);

        frame_offset = readmatrix(data_file,'Delimiter',{',',';'},'Range',[1 1 100 1]);
        frame_offset_index = find(frame_offset >= 1, 3);
        sample_0_row = frame_offset_index(3)+1;
        frame_offset = frame_offset(sample_0_row);
        stride_start_index = (stride_start_frame-frame_offset)*subframe_per_frame+sample_0_row;

        back_foot_gf = readmatrix(data_file,'Delimiter',{',',';'},'Range',[front_gf_start_index back_foot_col front_gf_start_index+frame_freq_hz*subframe_per_frame*scan_interval_seconds back_foot_col]);
        back_gf_start_index = find(back_foot_gf < -force_treshold, 1); % start first
        back_gf_end_index = find(back_foot_gf(back_gf_start_index+subframe_per_frame:end,1) > -force_treshold, 1) + back_gf_start_index+subframe_per_frame*2 + front_gf_start_index; % Search end
        stride_end_frame = readmatrix(data_file,'Delimiter',{',',';'},'Range',[back_gf_end_index 1 back_gf_end_index 1]);
        stride_end_index = (stride_end_frame-frame_offset)*subframe_per_frame+sample_0_row;
        gf_1 = -readmatrix(data_file,'Delimiter',{',',';'},'Range',[stride_start_index  front_foot_col   stride_end_index    front_foot_col]);
        gf_2 = -readmatrix(data_file,'Delimiter',{',',';'},'Range',[stride_start_index  back_foot_col    stride_end_index    back_foot_col]);
        % stride_end_index = back_gf_end_index + front_gf_start_index;
        % stride_start_index = front_gf_start_index-subframe_per_frame;
        

        
        clear back_foot_gf % Clear memory
        break
    end
end

%Repeat with CoMz

for i=1:max_scan_interval
    index_1 = stride_end_index+scan_interval*(i-1);
    index_2 = stride_end_index+scan_interval*i;
    comz = readmatrix(data_file,'Delimiter',{',',';'},'Range',[index_1 1 index_2 1]);
    comz_start_index = find(comz == stride_start_frame, 1); % Sync frame
     if isempty(comz_start_index)
        if i == max_scan_interval %Mean we reach the end of the scan without finding anything
            ME = MException('ComZ not found');
            throw(ME)
        end
     else
        try
         comz_start_index = comz_start_index + index_1 - 1;

        comz = readmatrix(data_file,'Delimiter',{',',';'},'Range',[comz_start_index comz_col comz_start_index+(stride_end_frame-stride_start_frame) comz_col]);
        comvz = [0;diff(comz)]; % 
        comvz(1,1) = comvz(2,1); % Fill first data
        catch
            disp('jiji');
        end
        break
     end

end


%% Generate Time Vectors & ComZ derivate
freq_gf = frame_freq_hz*subframe_per_frame;
T_gf = 1/freq_gf;
t_gf = 0:T_gf:(size(gf_1)/freq_gf - T_gf);
% Get diff GF vector for analisys 


freq_com = frame_freq_hz;
T_com = 1/freq_com;
t_com = 0:T_com:(size(comz)/freq_com - T_com);
comvz = comvz./T_com; % Transform to mm/s. 
abs_comv = abs(comvz);

% Low pass filter (Butterworth 8th order) for Ground forces
fs = freq_gf; % Sample Frequency (Need validation)
fc = 100; % Cutoff freqcuency of the filter, should investigate why this fc is chosen

[b,a] = butter(8,fc/(fs/2)); 
%[b,a] = besself(8,fc/(fs/2)); 

gf_1 = filtfilt(b,a,gf_1);
gf_2 = filtfilt(b,a,gf_2);

gf1_dt = filtfilt(b,a,[0;diff(gf_1)]); % 
gf2_dt = filtfilt(b,a,[0;diff(gf_2)]); % 


% First find max of both ground forces to have a reference
max_gf1 = max(gf_1);
max_gf2 = max(gf_2);


p1 = find(gf_1 > gf_2 - 10 & gf_1 < gf_2 + 10 & gf_1 > max_gf1*0.3, 1); % This is the middle point of the transition, we should find interval between the minimun and max comz around this dot
t_mid = t_gf(p1);

% FootContact of the front leg and ToesOff of the backleg:
FC_index = find(gf_2 > max_gf2/50 ,1) ;
FC_time = t_gf(FC_index);
TO_index = find(gf_1(p1:end) < max_gf1/50, 1)  + p1;
TO_time = t_gf(TO_index);

p2 = find(gf_2(TO_index:end) < max_gf2/50, 1) + TO_index;


% 
% 
% gf2_dt_aux = []
% for i=1:size(gf2_dt, 1)
%     if gf2_dt(i) > 0 
%         gf2_dt_aux(i) = 0;
%     else
%         gf2_dt_aux(i) = gf2_dt(i);
%     end
% end
% figure
% plot(gf2_dt_aux)
% 
% p2 = find(gf2_dt < -1.9, 1);

for i=1:p2
    if gf2_dt(p2-i) >= 0
        p2 = p2 - i;
        break
    end
end

%% Find ComV peaks and filter the consecutive peaks
[comv_peaks,comv_pk_index] = findpeaks(abs_comv); 
n = 0;
while n == 0
    [comv_peaks2,comv_pk_index2] = erase_extrapeak(comvz(comv_pk_index),comv_pk_index);
    size1 = size(comv_pk_index,1);
    size2 = size(comv_pk_index2,1);
    if(size1 == size2)
        break
    else
        comv_pk_index = comv_pk_index2;
        continue
    end
end


%% Get Step to step interval

% The Step-to-Step transition is the area where the CoM reaches its minimal vertical velocity before FC and its maximal vertical velocity after TO
index_2 = 0;
index_1 = 0;
for i=1:size(comv_pk_index,1)
    if t_com(comv_pk_index(i)) > t_mid
        index_2 = comv_pk_index(i);
        index_1 = comv_pk_index(i-1);
        break
    end

    if index_2 == size(comv_pk_index,1)
        ME = MException('LINE 177: Error finding step-to-step interval');
        throw(ME)
    end
end


sts_t1 = 0;
sts_t2 = 0;

if(t_com(index_1) < FC_time)
    sts_t1 = t_com(index_1);
else
    sts_t1 = FC_time;
end

if(t_com(index_2) > TO_time)
    sts_t2 = t_com(index_2);
else
    sts_t2 = TO_time;
end

i1 = find(t_gf >= sts_t1, 1);
i2 = find(t_gf >= sts_t2, 1);


%% Plot

figure
subplot(2,1,1)
hold on
plot(t_gf,gf_1,'b')
plot(t_gf,gf_2,'r')
xline(FC_time,'--r');
xline(TO_time,'--b');
xline(t_gf(p2),'black');
xregion(sts_t1, sts_t2, 'FaceColor', 'g')


legend('Back leg force','Front leg force')
xlabel('Time [s]')
ylabel('Ground Force [N]')
subplot(2,1,2)
hold on
plot(t_com, comvz,'b')
xline(FC_time,'--r');
xline(TO_time,'--b');
xregion(sts_t1, sts_t2, 'FaceColor', 'g')
plot(t_com(comv_pk_index), comvz(comv_pk_index), 'blacko')
legend('CoMv')
xlabel('Time [s]')
ylabel('Center of Mass velocity [mm/s]')


 figure
 hold on 
 plot(t_gf(i1:i2),gf_1(i1:i2),'b')
 plot(t_gf(i1:i2),gf_2(i1:i2),'r')


ForceRatio = mean(gf_2(i1:i2))/mean(gf_1(i1:i2))
sts_vvmin_and_vvmax_percentage = [100*sts_t1/(t_gf(p2)), 100*sts_t2/(t_gf(p2))]
sts_duration = sts_t2 - sts_t1
Stride_time = t_gf(p2)





% Function to clear a peaks arround the "real peak" 
function [ Y,I ] = erase_extrapeak(x, z) % 
    for i=1:(size(x,1)-1)
        value = x(i);
        next_value = x(i+1);
        if isempty(value)
            break
        else
            if value < 0
                if next_value < 0
                    if next_value <= value
                        z(i) = [];
                        x(i) = [];
                        Y = x;
                        I = z;
                        break
                    else
                        x(i+1) = [];
                        z(i+1) = [];
                        Y = x;
                        I = z;
                        break
                    end
    
                else
                    continue
                end
    
            else
                if next_value > 0
                    if next_value >= value
                        z(i) = [];
                        x(i) = [];
                        Y = x;
                        I = z;
                        break;
                    else
                        x(i+1) = [];
                        z(i+1) = [];
                        Y = x;
                        I = z;
                        break
                    end
    
                else
                    continue
                end
    
            end
    
        end
    end
    if i == (size(x,1)-1)
        Y = x;
        I = z;
        return
    end
end


function [ I, y ] = second_max( x ) % Used to get the 2nd largest value of the array, used to calculate the end of the stride with the right force
   [I, y] = max(x(x<max(x)));
end



%force_ratio_2 (3).m
%Mostrando force_ratio_2 (3).m.