function [st] = calculate_sta_ia(signal,spikes_frames,frames_before_spike,frames_after_spike,time_per_frame)

input_mat = signal; % the fluorescence matrix
ff = 1:length(input_mat); % change this if you want to choose a shorter trace than the orignial data

sta = zeros(frames_before_spike+frames_after_spike,1);
vec_ref = zeros(frames_before_spike+frames_after_spike,1);

st = zeros(frames_before_spike+frames_after_spike,length(spikes_frames));
for s = 1:length(spikes_frames) % loop over spikes in a trial
    if spikes_frames(s) > ff(end) || spikes_frames(s) < ff(1)
        continue
    elseif spikes_frames(s) < frames_before_spike
        frames = 1:spikes_frames(s)+frames_after_spike;
        st(end-length(frames)+1:end,s) = input_mat(frames);
        sta = sta + st(:,s);
        vec_ref(end-length(frames)+1:end) = vec_ref(end-length(frames)+1:end)+1;
    elseif spikes_frames(s) >= size(input_mat,1) - frames_after_spike
        frames = spikes_frames(s)-(frames_before_spike-1):size(input_mat,1);
        st(1:length(frames),s) = input_mat(frames);
        sta = sta + st(:,s);
        vec_ref(1:length(frames)) = vec_ref(1:length(frames))+1;
    else
        frames = spikes_frames(s)-(frames_before_spike-1):spikes_frames(s)+frames_after_spike;
        st(:,s) = input_mat(frames);
        sta = sta + st(:,s);
        vec_ref = vec_ref+1;
    end
end
sta = sta./(vec_ref);

% figure;
% x = -1*(frames_before_spike-1):frames_after_spike;
% % %errorbar(x*time_per_frame,sta,nanstd(st,[],2)/sqrt(length(spikes_frames)),'b'); hold on
% shadedErrorBar(x*time_per_frame,sta,nanstd(st,[],2)/sqrt(length(spikes_frames))); hold on
% % xlabel('time')


%% do the same for a flourescence shuffled matrix
% all_sta = zeros(1000,frames_before_spike+frames_after_spike);
% for i = 1:1000
%     shuf_ind = randperm(length(input_mat)); % shuffle the data
%     shuf_mat = input_mat(shuf_ind);
%     
%     sta = zeros(frames_before_spike+frames_after_spike,1);
%     vec_ref = zeros(frames_before_spike+frames_after_spike,1);
%     
%     for s = 1:length(spikes_frames) % loop over spikes in a trial
%         st_shuf = zeros(frames_before_spike+frames_after_spike,1);
%         if spikes_frames(s) > ff(end) || spikes_frames(s) < ff(1)
%             continue
%         elseif spikes_frames(s) < frames_before_spike
%             frames = 1:spikes_frames(s)+frames_after_spike;
%             st_shuf(end-length(frames)+1:end) = shuf_mat(frames);
%             sta = sta + st_shuf;
%             vec_ref(end-length(frames)+1:end) = vec_ref(end-length(frames)+1:end)+1;
%         elseif spikes_frames(s) >= size(input_mat,1) - frames_after_spike
%             frames = spikes_frames(s)-(frames_before_spike-1):size(shuf_mat,1);
%             st_shuf(1:length(frames)) = shuf_mat(frames);
%             sta = sta + st_shuf;
%             vec_ref(1:length(frames)) = vec_ref(1:length(frames))+1;
%         else
%             frames = spikes_frames(s)-(frames_before_spike-1):spikes_frames(s)+frames_after_spike;
%             st_shuf = shuf_mat(frames);
%             sta = sta + st_shuf;
%             vec_ref = vec_ref+1;
%         end
%     end
%     sta = sta./(vec_ref);
%     all_sta(i,:) = sta;
% end

% errorbar(x*time_per_frame,mean(all_sta),std(all_sta),'r')
% shadedErrorBar(x*time_per_frame,mean(all_sta),std(all_sta),'lineprops','r')
% plot(x,sta,'r');
% legend('real','shuffle')    
