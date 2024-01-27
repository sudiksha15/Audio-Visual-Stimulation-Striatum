%% Compute spectral power changes in movement speed in baseline vs stim 

close all
clear all

%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';
stim_freq=1;
aSPs=[];aDx=[];
if stim_freq==1
% 10Hz
allSES{1}='\608448\01062020\608448';  % no LFP
allSES{2}='\608451\01062020\608451';
allSES{3}='\608452\01062020\608452';
allSES{4}='\608450\03092020\608450';
allSES{5}='\602101\03092020\602101';
allSES{6}='\611311\03092020\611311';
allSES{7}='\608448\01172020\608448';
allSES{8}='\608451\07082020\608451';
allSES{9}='\608452\07082020\608452';
allSES{10}='\611111\07082020\611111';
allSES{11}='\602101\07082020\602101';
allSES{12}='\612535\07082020\612535';
allSES{13}='\615883\02052021\10Hz\615883';
allSES{14}='\615883\03122021\10Hz\615883';
else
% 140Hz
allSES{1}='\608448\01102020\608448';
allSES{2}='\608451\01102020\608451';
allSES{3}='\608452\01102020\608452';
allSES{4}='\602101\03112020\602101';
allSES{5}='\611311\03112020\611311';
allSES{6}='\611111\03112020\611111';
allSES{7}='\608450\01172020\608450';
allSES{8}='\608451\07012020\608451';
allSES{9}='\608452\07012020\608452';
allSES{10}='\611111\07012020\611111';
allSES{11}='\612535\07022020\612535';
allSES{12}='\602101\07022020\602101'; %exclude LFP
allSES{13}='\615883\02052021\145Hz\615883';
allSES{14}='\615883\03122021\145Hz\615883';
end
p_all_sess=[];
for ses=1:length(allSES)
if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual\motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual\motion_corrected'])   
end

motion=h5read('processed_motion.h5','/raw_speed_trace');
moving_period=h5read('processed_motion.h5','/moving_period');
stn_period=h5read('processed_motion.h5','/stationary_period');

load('LFP_ts.mat')
stim_onsets=LFP_data.Stim_onset;
stim_offsets=LFP_data.Stim_offset;

v=motion;
idx=find(isnan(v)==0);
%Remove NANs
v = v(~isnan(v));
moving_period=moving_period(idx);
stn_period=stn_period(idx);
%% Taking all periods regardless of baseline and stim
hi_idx=find(moving_period==1);
low_idx=find(stn_period==1);

%vamp=abs(hilbert(v));
vamp=fastsmooth(v,50,1,1); % motion smoothed

%% Baseline vs stim - taking all baseline periods
stim_vec=zeros(1,size(motion,1));
for i=1:length(stim_onsets)
timsel=(stim_onsets(i)-0:stim_offsets(i));
stim_vec(timsel)=1;
end
stim_vec=stim_vec(~isnan(v));
stim_vec=stim_vec';

bs_vec=zeros(1,size(motion,1));
timsel=(1:stim_onsets(1)-1);
bs_vec(timsel)=1;
for i=1:length(stim_onsets)-1
timsel=(stim_offsets(i)+1200:stim_onsets(i+1)-1);
bs_vec(timsel)=1;
end
timsel=(stim_offsets(5)+1200:size(motion,1));
bs_vec(timsel)=1;
bs_vec=bs_vec(~isnan(v));
bs_vec=bs_vec';
moving_bs= moving_period.*bs_vec;
moving_stim= moving_period.*stim_vec;
hi_bs_idx=find(moving_bs==1);
hi_stim_idx=find(moving_stim==1);

stn_bs= stn_period.*bs_vec;
stn_stim= stn_period.*stim_vec;
low_bs_idx=find(stn_bs==1);
low_stim_idx=find(stn_stim==1);
%% Comparing a minute before to a minute during stim exactly 

stim_vec2=zeros(1,size(motion,1));
for i=1:length(stim_onsets)
timsel=(stim_onsets(i):stim_onsets(i)+1200-1);
stim_vec2(timsel)=1;
end
stim_vec2=stim_vec2(~isnan(v));
stim_vec2=stim_vec2';


bs_vec2=zeros(1,size(motion,1));
for i=1:length(stim_onsets)
timsel=(stim_onsets(i)-1200:stim_onsets(i)-1);
bs_vec2(timsel)=1;
end
bs_vec2=bs_vec2(~isnan(v));
bs_vec2=bs_vec2';
moving_bs2= moving_period.*bs_vec2;
moving_stim2= moving_period.*stim_vec2;
hi_bs_idx2=find(moving_bs2==1);
hi_stim_idx2=find(moving_stim2==1);

stn_bs2= stn_period.*bs_vec2;
stn_stim2= stn_period.*stim_vec2;
low_bs_idx2=find(stn_bs2==1);
low_stim_idx2=find(stn_stim2==1);

%% Compute spectrogram
lfp=[];
lfp.trial{1}= v';
lfp.time{1}= (1:length(v))/20;
lfp.label= {'A'}

cfg = []; %block_type == cfg.blk
cfg.method ='mtmconvol'; %'mvar';
cfg.output ='fourier';
cfg.taper='hanning';
cfg.keeptapers ='yes';
cfg.keeptrials ='yes';
cfg.trials='all';cfg.tapsmofrq =6;%
cfg.channel= ['all']; %chans=cfg.channel;
cfg.foi= [0.2:0.2:10];
cfg.t_ftimwin =[ones(1,length(cfg.foi))*5];  
cfg.toi=lfp.time{1}(1:1:end) ;
freq2 = ft_freqanalysis(cfg, lfp);

s= abs(squeeze(freq2.fourierspctrm));
z=nanmean(abs(s(:,hi_idx)),2)./nanmean(abs(s(:,low_idx)),2);
z_hi=nanmean(abs(s(:,hi_idx)),2);
z_low=nanmean(abs(s(:,low_idx)),2);
p_all_sess(:,ses)=z;
p_hi_all_sess(:,ses)=z_hi;
p_low_all_sess(:,ses)=z_low;

z2=nanmean(abs(s(:,hi_bs_idx)),2)./nanmean(abs(s(:,low_bs_idx)),2);
z3=nanmean(abs(s(:,hi_stim_idx)),2)./nanmean(abs(s(:,low_stim_idx)),2);
z_hi_bs=nanmean(abs(s(:,hi_bs_idx)),2);
z_hi_stim=nanmean(abs(s(:,hi_stim_idx)),2);
p_hi_bs_all_sess(:,ses)=z_hi_bs;
p_hi_stim_all_sess(:,ses)=z_hi_stim;
p_all_sess2(:,ses)=z2;
p_all_sess3(:,ses)=z3;

z4=nanmean(abs(s(:,hi_bs_idx2)),2)./nanmean(abs(s(:,low_bs_idx2)),2);
z5=nanmean(abs(s(:,hi_stim_idx2)),2)./nanmean(abs(s(:,low_stim_idx2)),2);
z_hi_bs2=nanmean(abs(s(:,hi_bs_idx2)),2);
z_hi_stim2=nanmean(abs(s(:,hi_stim_idx2)),2);
p_hi_bs_all_sess2(:,ses)=z_hi_bs2;
p_hi_stim_all_sess2(:,ses)=z_hi_stim2;
p_all_sess4(:,ses)=z4;
p_all_sess5(:,ses)=z5;
end

%% Plotting all periods
figure('Color','w'),
plot(freq2.freq,nanmean(p_hi_all_sess,2),'g')
hold on,plot(freq2.freq,nanmean(p_low_all_sess,2),'r--')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
set(gca,'Xtick',[  0.5 1 2 3 4 5 6 7 8 9 10], 'Xticklabel', [ 0.5 1 2 3 4 5 6 7 8 9 10])
axis tight
legend run rest 
xlabel('Frequency (log)') , ylabel('Amplitude (log)')

figure('COlor','w','Position', [ 300 200 130 pheight],'Renderer', 'painters') 
plot(freq2.freq,nanmean(p_all_sess,2),'k','LineWidth',2.0)
fill_error_area2(freq2.freq,nanmean(p_all_sess,2), nanstd(p_all_sess,[],2)./sqrt(size(p_all_sess,2)),[ 0.5 0.5 0.5])
axis tight
xlim([0.6,10])
%ylim([3.5,8])
yticks([4:8])
xticks([2:2:10])
xlabel('Frequency') ;ylabel('Relative Amplitude')
title('Relative amplitude (mov/rest)')

savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\'
%% Compare all baseline vs all stim
% figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters') 
% plot(freq2.freq,nanmean(p_hi_bs_all_sess,2),'b')
% fill_error_area2(freq2.freq,nanmean(p_hi_bs_all_sess,2), nanstd(p_hi_bs_all_sess,[],2)./sqrt(size(p_hi_bs_all_sess,2)),'b')
% hold on,plot(freq2.freq,nanmean(p_hi_stim_all_sess,2),'r--')
% fill_error_area2(freq2.freq,nanmean(p_hi_stim_all_sess,2), nanstd(p_hi_stim_all_sess,[],2)./sqrt(size(p_hi_stim_all_sess,2)),'r')
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
% set(gca,'Xtick',[  0.5 1 2 3 4 5 6 7 8 9 10], 'Xticklabel', [ 0.5 1 2 3 4 5 6 7 8 9 10])
% axis tight
% legend baseline stim
% xlabel('Frequency (log)') , ylabel('Amplitude (log)')
% 
% figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters') 
% plot(freq2.freq,nanmean(p_hi_bs_all_sess,2),'b')
% fill_error_area2(freq2.freq,nanmean(p_hi_bs_all_sess,2), nanstd(p_hi_bs_all_sess,[],2)./sqrt(size(p_hi_bs_all_sess,2)),'b')
% hold on,plot(freq2.freq,nanmean(p_hi_stim_all_sess,2),'r--')
% fill_error_area2(freq2.freq,nanmean(p_hi_stim_all_sess,2), nanstd(p_hi_stim_all_sess,[],2)./sqrt(size(p_hi_stim_all_sess,2)),'r')
% % set(gca,'Xscale','log')
% % set(gca,'Yscale','log')
% set(gca,'Xtick',[  0.5 1 2 3 4 5 6 7 8 9 10], 'Xticklabel', [ 0.5 1 2 3 4 5 6 7 8 9 10])
% axis tight
% legend baseline stim
% xlabel('Frequency') , ylabel('Amplitude')
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'mov_vs_rest_amp_all_bs_stim_145Hz.pdf'])
% 
idx=find(freq2.freq>=3);
idx2=find(freq2.freq<=4);
bs_stat=nanmean(p_hi_bs_all_sess(idx(1):idx2(end),:),1);
stim_stat=nanmean(p_hi_stim_all_sess(idx(1):idx2(end),:),1);
[h,p] = ttest(bs_stat,stim_stat);
% significant 
p1 = signrank(bs_stat,stim_stat);
pheight=160;
figure('COlor','w','Position', [ 300 200 130 pheight],'Renderer', 'painters')
line([ 0.8 2.3], [ -0.01 -0.01],'COlor', [ 1 1 1 ],'Linewidth',0.5)
violinplot2(bs_stat',[1.2 ],'ViolinColor', [ 0.3 0.5 0.9; 0 0 0.9])
violinplot2(stim_stat',[1.9 ],'ViolinColor', [ 0.9 0.3 0.7; 0 0 0.9]);axis tight
line([ 0.9 1.5], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
line([ 1.6 2.2], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
%ylim([-0.8 25]);
set(gca,'Xtick',[1.2 1.9],'Xticklabel',{'Base','10 Hz'})
%set(gca,'Xtick',[],'Xticklabel',{})
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'VIOL_mov_peak_spect_10Hz_all_periods.pdf'])

% x=[ones(14,1);ones(14,1).*2];
% y=[bs_stat';stim_stat'];
% figure
% boxplot(y,x)
% 
figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters')  
plot(freq2.freq,nanmean(p_all_sess2,2),'b')
fill_error_area2(freq2.freq,nanmean(p_all_sess2,2), nanstd(p_all_sess2,[],2)./sqrt(size(p_all_sess2,2)),'b')
plot(freq2.freq,nanmean(p_all_sess3,2),'r')
fill_error_area2(freq2.freq,nanmean(p_all_sess3,2), nanstd(p_all_sess3,[],2)./sqrt(size(p_all_sess3,2)),'r')
% Add standard deviation
axis tight
xlabel('Frequency') ;ylabel('Relative Amplitude')
title('Relative amplitude (mov/rest)- Baseline vs stim')
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'mov_rest_amp_all_bs_stim_145Hz.pdf'])
% 
% idx=find(freq2.freq>=3);
% idx2=find(freq2.freq<=4);
% bs_stat=nanmean(p_all_sess2(idx(1):idx2(end),:),1);
% stim_stat=nanmean(p_all_sess3(idx(1):idx2(end),:),1);
% [h,p] = ttest(bs_stat,stim_stat);
% p1 = signrank(bs_stat,stim_stat);

%% Compare a minute before stim to a minute during stim 
figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters') 
plot(freq2.freq,nanmean(p_hi_bs_all_sess2,2),'g')
fill_error_area2(freq2.freq,nanmean(p_hi_bs_all_sess2,2), nanstd(p_hi_bs_all_sess2,[],2)./sqrt(size(p_hi_bs_all_sess2,2)),'g')
hold on,plot(freq2.freq,nanmean(p_hi_stim_all_sess2,2),'r--')
fill_error_area2(freq2.freq,nanmean(p_hi_stim_all_sess2,2), nanstd(p_hi_stim_all_sess2,[],2)./sqrt(size(p_hi_stim_all_sess2,2)),'r')
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')
set(gca,'Xtick',[  0.5 1 2 3 4 5 6 7 8 9 10], 'Xticklabel', [ 0.5 1 2 3 4 5 6 7 8 9 10])
axis tight
legend baseline stim
xlabel('Frequency') , ylabel('Amplitude')
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'mov_vs_rest_amp_1min_bs_stim_145Hz.pdf'])

idx=find(freq2.freq>=3);
idx2=find(freq2.freq<=4);
bs_stat=nanmean(p_hi_bs_all_sess2(idx(1):idx2(end),:),1);
stim_stat=nanmean(p_hi_stim_all_sess2(idx(1):idx2(end),:),1);
[h,p] = ttest(bs_stat,stim_stat);
% significant 
p1 = signrank(bs_stat,stim_stat);
pheight=160;
figure('COlor','w','Position', [ 300 150 130 pheight],'Renderer', 'painters')
line([ 0.8 2.3], [ 0.2 0.2],'COlor', [ 1 1 1 ],'Linewidth',0.5)
violinplot2(bs_stat',[1.2 ],'ViolinColor', [ 0.3 0.5 0.9; 0 0 0.9])
violinplot2(stim_stat',[1.9 ],'ViolinColor', [ 0.9 0.3 0.7; 0 0 0.9]);axis tight
%line([ 0.9 1.5], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
%line([ 1.6 2.2], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
ylim([0.2 1.2]);
set(gca,'Xtick',[1.2 1.9],'Xticklabel',{'Base','145 Hz'})
%set(gca,'Xtick',[],'Xticklabel',{})
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'VIOL_mov_peak_spect_145Hz_final.pdf'])


% cd('/home/hanlabadmins/Codes_git/Stimulation-Project/Trace_analysis/Final_codes')
% filename='bs_stim_delta_power_movement_145Hz_2_4.h5';
% dataset = bs_stat; 
% h5create(filename, ['/1_min_bs' ], size(dataset));
% h5write(filename, ['/1_min_bs'], dataset);
% dataset = stim_stat; 
% h5create(filename, ['/1_min_stim' ], size(dataset));
% h5write(filename, ['/1_min_stim'], dataset);
% dataset = bs_stat2; 
% h5create(filename, ['/all_bs' ], size(dataset));
% h5write(filename, ['/all_bs'], dataset);
% dataset = stim_stat2; 
% h5create(filename, ['/all_stim' ], size(dataset));
% h5write(filename, ['/all_stim'], dataset);
%     
figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters') 
fill_error_area2(freq2.freq,nanmean(p_all_sess4,2), nanstd(p_all_sess4,[],2)./sqrt(size(p_all_sess4,2)),'b')
fill_error_area2(freq2.freq,nanmean(p_all_sess5,2), nanstd(p_all_sess5,[],2)./sqrt(size(p_all_sess5,2)),'r')
plot(freq2.freq,nanmean(p_all_sess4,2),'b')
plot(freq2.freq,nanmean(p_all_sess5,2),'r')
% Add standard deviation
axis tight
xlabel('Frequency') ;ylabel('Relative Amplitude')
title('Relative amplitude (mov/rest)- Baseline vs stim')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'mov_rest_amp_1min_bs_stim_10Hz.pdf'])

% idx=find(freq2.freq>=3);
% idx2=find(freq2.freq<=4);
% bs_stat=nanmean(p_all_sess4(idx(1):idx2(end),:),1);
% stim_stat=nanmean(p_all_sess5(idx(1):idx2(end),:),1);
% [h,p] = ttest(bs_stat,stim_stat);
% p1 = signrank(bs_stat,stim_stat);
% %save('motion_spect_10Hz.mat','p_all_sess')
% 


%% All sessions across 10 Hz and 145 Hz
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar')
% load('motion_spect_10Hz.mat')
% p_10Hz=p_all_sess;
% load('motion_spect_145Hz.mat')
% p_145Hz=p_all_sess;
% p_tot=[p_10Hz,p_145Hz];
% figure('Color','w') 
% plot(freq2.freq,nanmean(p_tot,2),'k')
% fill_error_area2(freq2.freq,nanmean(p_tot,2), nanstd(p_tot,[],2)./sqrt(size(p_tot,2)),[ 0.5 0.5 0.5])
% axis tight
% % Add standard deviation
% xlabel('Frequency') ;ylabel('Relative Amplitude')
% title('Relative amplitude (mov/rest)')