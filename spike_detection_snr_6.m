%%
% it calculates spikes upon ROI manual selection
% last uptate: 11-04-2021, Victor Hugo Cornejo
close all; clear; clc

%% files load

basedir = '\...\';                                                         %% enter folder
denoised = loadtiff([basedir 'mov_denoisedCh2_mc_crop.tif']);              %% denoised tiff file
% load('collected_data.mat','M1');                                         %% M1 non-filtered matrix
load('collected_data.mat','M1');                                           %% M filtered matrix
% xml_file_name = 'RFP-spASAP-anesth-stim50ms-V1-013';                     %% .xml raw file

%% Initial processing
[h,w,ss] = size(denoised);
denoised=double(denoised);

data_raw = zeros(h,w,ss);                                                  %% non filtered Matrix
for i=1:h*w*ss
	data_raw(i)=M1(i);
end

% data_filt = zeros(h,w,ss);                                               %% filtered Matrix
% for i=1:h*w*ss
% 	data_filt(i)=M(i);
% end

%% select ROI

beep
figure(1);
subplot(1,3,1)
imagesc(mean(denoised,3));
title('Select spine')
axis image; 
colormap(jet);
roi_spine = choose_polygon2(w,h); 

subplot(1,3,2)
imagesc(mean(denoised,3));
title('Select dendrite')
axis image; 
colormap(jet);
roi_dendrite = choose_polygon2(w,h); 

subplot(1,3,3)
imagesc(mean(denoised,3));
title('Select background')
axis image; 
colormap(jet);
roi_bg = choose_polygon2(w,h); 

mean_spine_filt = zeros(1,ss);
mean_dendrite_filt = zeros(1,ss);
mean_bg_filt = zeros(1,ss);
for frame = 1:size(data_raw,3)
    temp = data_raw(:,:,frame);
    mean_spine_filt(frame) = mean(temp(roi_spine));
    mean_dendrite_filt(frame) = mean(temp(roi_dendrite));
    mean_bg_filt(frame) = mean(temp(roi_bg));
end

%% Detect spikes in spines
v_thr = find(mean_spine_filt>=std(mean_bg_filt(20:end))*3);                %% select SNR threshold: SNR=3
f = find(diff(v_thr)>1);

periods = zeros(length(f),2);
periods(1,1) = v_thr(1);
for i = 1:size(periods,1)
    periods(i,2) = v_thr(f(i));
    periods(i+1,1) = v_thr(f(i)+1);
end
periods(end,2) = v_thr(end);
for i = 1:size(periods,1)
    [m,temp] = max(mean_spine_filt(periods(i,1):periods(i,2)));
    spike_times_s(i) = temp+periods(i,1)-1;
end

ex_first_s=find(spike_times_s<=20);
spike_times_s(ex_first_s)=[];
ex_end_s=find(spike_times_s>=(length(mean_spine_filt)-20));
spike_times_s(ex_end_s)=[];

spike_s_diff=diff(spike_times_s);
spike_s_diff=find(spike_s_diff<3)+1;                                       %% frames two adjacent spikes
spike_times_s(spike_s_diff)=[];

s = length(spike_times_s);

%% Detect spikes in ROI dendrite
v_thrd = find(mean_dendrite_filt>=std(mean_bg_filt)*3);                    %% select SNR threshold: SNR=3
fd = find(diff(v_thrd)>1);

periodsd = zeros(length(fd),2);
periodsd(1,1) = v_thrd(1);
for i = 1:size(periodsd,1)
    periodsd(i,2) = v_thrd(fd(i));
    periodsd(i+1,1) = v_thrd(fd(i)+1);
end
periodsd(end,2) = v_thrd(end);
for i = 1:size(periodsd,1)
    [md,tempd] = max(mean_dendrite_filt(periodsd(i,1):periodsd(i,2)));
    spike_times_d(i) = tempd+periodsd(i,1)-1;
end

ex_first_d=find(spike_times_d<=20);
spike_times_d(ex_first_d)=[];
ex_end_d=find(spike_times_d>=(length(mean_dendrite_filt)-20));
spike_times_d(ex_end_d)=[];

spike_d_diff=diff(spike_times_d);
spike_d_diff=find(spike_d_diff<3)+1;  %% frames two adjacent spikes
spike_times_d(spike_d_diff)=[];

sd = length(spike_times_d);

%% Spike separation Spine vs dendrite

frame_window= 1;                                                           % frames between spikes detection in sp vs dend

frame_search_d=zeros((frame_window*2+1),length(spike_times_d));
for r=1:length(spike_times_d)
    r2=(spike_times_d(r)-frame_window):1:(spike_times_d(r)+frame_window);
    frame_search_d(:,r)=r2;
end

r4=zeros(1,length(spike_times_s));
for r=1:length(spike_times_s)
    r3=find(frame_search_d==spike_times_s(r));
        if isempty(r3)==1
           r3=1;
        end
    r4(r)=frame_search_d(r3);
end
era_dis=find(r4 > frame_search_d(1));
spike_times_EPSP_AP=r4(era_dis);

spike_times_EPSP_only=spike_times_s;
for r5=1:length(spike_times_EPSP_AP)
    r6=find(spike_times_EPSP_only==spike_times_EPSP_AP(r5));
    spike_times_EPSP_only(r6)=0;
end
r7=find(spike_times_EPSP_only==0);
spike_times_EPSP_only(r7)=[];

%
frame_search_s=zeros((frame_window*2+1),length(spike_times_s));
for r=1:length(spike_times_s)
    r8=(spike_times_s(r)-frame_window):1:(spike_times_s(r)+frame_window);
    frame_search_s(:,r)=r8;
end

r9=zeros(1,length(spike_times_d));
for r=1:length(spike_times_d)
    r10=find(frame_search_s==spike_times_d(r));
        if isempty(r10)==1
           r10=1;
        end
    r9(r)=frame_search_s(r10);
end
era_dis1=find(r9 > frame_search_s(1));
spike_times_AP_EPSP=r9(era_dis1);

spike_times_AP_only=spike_times_d;
for r=1:length(spike_times_AP_EPSP)
    r11=find(spike_times_AP_only==spike_times_AP_EPSP(r));
    spike_times_AP_only(r11)=0;
end
r12=find(spike_times_AP_only==0);
spike_times_AP_only(r12)=[];

r13=spike_times_AP_EPSP-spike_times_EPSP_AP;
den_to_sp=spike_times_AP_EPSP(r13<=0);
sp_to_den=spike_times_EPSP_AP(r13>0);

%% Figure Representation

figure(2);
yyaxis left;
plot(mean_spine_filt, 'b');
hold on;
plot(spike_times_EPSP_only,mean_spine_filt(spike_times_EPSP_only),...
    'o','MarkerSize',7,'MarkerEdgeColor','b'); hold on;
plot(spike_times_EPSP_AP,mean_spine_filt(spike_times_EPSP_AP),'*',...
    'MarkerSize',8,'MarkerEdgeColor','k');

yyaxis right;
plot(mean_dendrite_filt, 'r');
hold on;
plot(spike_times_AP_only,mean_dendrite_filt(spike_times_AP_only),'o',...
    'MarkerSize',7,'MarkerEdgeColor','r'); hold on;
plot(spike_times_AP_EPSP,mean_dendrite_filt(spike_times_AP_EPSP),'*',...
    'MarkerSize',8,'MarkerEdgeColor','k');

legend('Spine', 'Spine Only', 'Both', 'Dendrite', 'Dendrite Only');

%
figure(3);
yyaxis left;
plot(mean_spine_filt, 'b');
hold on;
plot(sp_to_den,mean_spine_filt(sp_to_den),'*','MarkerSize',7,'MarkerEdgeColor','b');

yyaxis right;
plot(mean_dendrite_filt, 'r');
hold on;
plot(den_to_sp,mean_dendrite_filt(den_to_sp),'*','MarkerSize',7,'MarkerEdgeColor','r');

legend('Spine', 'Spine->Dendrite', 'Dendrite', 'Dendrite->Spine' );

%% Calculate peak dF/F and SNR

st_spine_filt = calculate_sta_ia(mean_spine_filt',spike_times_s,15,15,1);
st_dendrite_filt = calculate_sta_ia(mean_dendrite_filt',spike_times_d,15,15,1);

%spine

peak_dfof_sp=zeros(1,length(spike_times_s));
fo_sp=zeros(1,length(spike_times_s));
fluo_change_sp=zeros(1,length(spike_times_s));
SNR_sp=zeros(1,length(spike_times_s));
area_sp=zeros(1,length(spike_times_s));

for i=1:length(spike_times_s)
    peak_dfof_sp(1,i)=max(st_spine_filt(1:30,i));
    fo_sp(1,i)=min(st_spine_filt(1:15,i));
    fluo_change_sp(1,i)= peak_dfof_sp(1,i)-fo_sp(1,i);
    SNR_sp(1,i)=peak_dfof_sp(1,i)/std(mean_bg_filt);
    area_sp(1,i)=trapz(st_spine_filt(1:30,i));
end

% dendrite

peak_dfof_den=zeros(1,length(spike_times_d));
fo_den=zeros(1,length(spike_times_d));
fluo_change_den=zeros(1,length(spike_times_d));
SNR_den=zeros(1,length(spike_times_d));
area_den=zeros(1,length(spike_times_d));

for i=1:length(spike_times_d)
    peak_dfof_den(1,i)=max(st_dendrite_filt(1:30,i));
    fo_den(1,i)=min(st_dendrite_filt(1:15,i));
    fluo_change_den(1,i)= peak_dfof_den(1,i)-fo_den(1,i);
    SNR_den(1,i)=peak_dfof_den(1,i)/std(mean_bg_filt);
    area_den(1,i)=trapz(st_dendrite_filt(1:30,i));
end

%%
table_SP_only=zeros(length(spike_times_s),1);
for i=1:length(spike_times_s)
    r14=(spike_times_s(i)==spike_times_EPSP_only);
    if sum(r14)==0
       table_SP_only(i)=0;
    else 
       table_SP_only(i)=1;
    end
end

table_SP_DEN=zeros(length(spike_times_s),1);
for i=1:length(spike_times_s)
    r15=(spike_times_s(i)==spike_times_EPSP_AP);
    if sum(r15)==0
       table_SP_DEN(i)=0;
    else 
       table_SP_DEN(i)=1;
    end
end

table_SP_first_DEN=zeros(length(spike_times_s),1);
for i=1:length(spike_times_s)
    r16=(spike_times_s(i)==sp_to_den);
    if sum(r16)==0
       table_SP_first_DEN(i)=0;
    else 
       table_SP_first_DEN(i)=1;
    end
end

%%
table_DEN_only=zeros(length(spike_times_d),1);
for i=1:length(spike_times_d)
    r14=(spike_times_d(i)==spike_times_AP_only);
    if sum(r14)==0
       table_DEN_only(i)=0;
    else 
       table_DEN_only(i)=1;
    end
end

table_DEN_SP=zeros(length(spike_times_d),1);
for i=1:length(spike_times_d)
    r15=(spike_times_d(i)==spike_times_AP_EPSP);
    if sum(r15)==0
       table_DEN_SP(i)=0;
    else 
       table_DEN_SP(i)=1;
    end
end

table_DEN_first_SP=zeros(length(spike_times_d),1);
for i=1:length(spike_times_d)
    r16=(spike_times_d(i)==den_to_sp);
    if sum(r16)==0
       table_DEN_first_SP(i)=0;
    else 
       table_DEN_first_SP(i)=1;
    end
end

%% Compiled info

COMPILED_TABLE_SPINE=table(spike_times_s', peak_dfof_sp', fo_sp', fluo_change_sp',...
    SNR_sp', area_sp', table_SP_only, table_SP_DEN, table_SP_first_DEN,...
    'VariableNames',{'Frame','Peak_dfof','F0','Fluo_change','SNR','Area','SP_only',... 
    'SP_DEN', 'SP_first'});

COMPILED_TABLE_DENDRITE=table(spike_times_d', peak_dfof_den', fo_den', fluo_change_den',...
    SNR_den', area_den', table_DEN_only, table_DEN_SP, table_DEN_first_SP,...
    'VariableNames',{'Frame','Peak_dfof','F0','Fluo_change','SNR','Area','DEN_only',... 
    'DEN_SP', 'DEN_first'});

%% Save files

% files
clear data_raw
clear denoised
clear M1

mkdir spike_analysis_spine;
cd spike_analysis_spine;
save('collected_data_spikes.mat');
writetable(COMPILED_TABLE_SPINE,'COMPILED_TABLE_SPINE.xlsx')
writetable(COMPILED_TABLE_DENDRITE,'COMPILED_TABLE_DENDRITE.xlsx')
savefig(1,'figure1_ROI.fig');
savefig(2,'figure2_traces1.fig');
savefig(3,'figure3_traces2.fig');
savefig(4,'figure4_trial_average_sp.fig');
savefig(5,'figure5_trial_average_den.fig');