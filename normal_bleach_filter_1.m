%%
% it transforms raw denoised series into beached corrected, normalized and
% filtered tiff series.
% last uptate: 11-04-2021, Victor Hugo Cornejo
close all; clear; clc

%%

file_name='mov_denoisedCh2_mc_crop';
Freq_filt = 10;                                                            %% filtering frequency cuttoff in Hz
frame_rate=60;

%% files load
tic
data = loadtiff([file_name '.tif']);                                       %% denoised tiff file

%% Create new image normalized, filtered and bleach corrected

% figure(1);
% subplot(1,2,1)
% imagesc(mean(data,3));
% title('Original Image', 'FontSize',12);
% colormap (jet);
% axis image;
% set(gca,'YTickLabel',[],'XTickLabel',[]);

[h,w,ss] = size(data);
data=double(data);

M=zeros(h*w,ss);
parfor i=1:h*w*ss
M(i)=data(i);
end

%wb = waitbar(0, 'Calculating bleach correction...');
parfor ii=1:h*w
%   f_bleach=fit( (1:ss)', M(ii,:)','exp2');                               % exp2 for decay
    f_bleach=fit( (1:ss)', M(ii,:)','smoothingspline','SmoothingParam',0.000001);  
    M(ii,:) = ((  f_bleach(1:ss) - M(ii,:)' ) ./ f_bleach(1:ss)  )' ;
    %waitbar(ii/(h*w), wb);
end
%close(wb);

M1 = M;                                                                    %% M1 non-filtered matrix, M filtered matrix

[b,a] = butter(5,Freq_filt/(frame_rate/2));                                            
% freqz(b,a)
parfor ii=1:h*w
   M(ii,:)=filtfilt(b,a,M(ii,:));
end

parfor i=1:h*w*ss
	data(i)=M(i);
end

% subplot(1,2,2)
% imagesc(mean(data,3));
% title('Processed image', 'FontSize',12);
% colormap (jet);
% axis image;
% set(gca,'YTickLabel',[],'XTickLabel',[]);

%% Image parameters for tiff construction

zero_real_filt = min(data(:));                                             %% converted to 0 tiff file: 0 intensity
max_real_filt = max(data(:));                                              %% converted to 1 tiff file: 255 intensity

%% Save files

% mkdir data_mov_denoisedCh2_mc_df_f0_bleach_filt;
% cd data_mov_denoisedCh2_mc_df_f0_bleach_filt;
% savefig(1,'Original_processed_image.fig');
save('collected_data.mat','zero_real_filt','max_real_filt','M','M1','Freq_filt','frame_rate');
toc
save_tif_stack_YS(data, [file_name '_df_f0_bleach_filt_' num2str(Freq_filt) 'hz.tif']);
% sound(sin(1:3000)); 
disp('work done');