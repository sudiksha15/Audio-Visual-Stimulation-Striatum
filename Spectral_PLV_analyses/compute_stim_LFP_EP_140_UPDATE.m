clear all

%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';
stim_freq=2;
aEP=[];aEP2=[];aEP3=[];
 clear dataT2
mm2=0
mm=0;
if stim_freq==1
% 10Hz
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
removal_stim{1}= [];
removal_stim{2}= [];
removal_stim{3}= [];
removal_stim{4}= [1 2];
removal_stim{5}= [3 4 5];
removal_stim{6}= [1 2 3 4 5];
removal_stim{7}= [1 2 3 4 5];
removal_stim{8}= [1 2 3 4 5];
removal_stim{9}= [1 2];
removal_stim{10}= [1 2 3 4 5];
removal_stim{11}= [1 2 3 4 5];
removal_stim{12}= [1 2 3 4 5];
removal_stim{13}= [1 2 3 4 5];
removal_stim{14}= [1 2 4 5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% 6 problematic
for ses=[1:5 9 14]%length(allSES)
    ses
if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual/motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual/motion_corrected'])   
end
rem_stim= removal_stim{ses};
%cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Sudi_Sridhar\608452\07012020\608452_145Hz_AudioVisual\motion_corrected')




motion=h5read('processed_motion.h5','/raw_speed_trace');
tracesF=h5read('processed_trace.h5','/trace');
traces=h5read('processed_trace.h5','/onset_binary_trace');
load('LFP_ts.mat')
stim_onsets=LFP_data.Stim_onset;
stim_offsets=LFP_data.Stim_offset;
v=motion;
idx=find(isnan(v)==0);
%Remove NANs
v = v(~isnan(v));  % motion signal
% Align traces and motiontraces=traces(idx,:);
traces=traces(idx,:);tracesF=tracesF(idx,:);


Sampling_freq=20;
% Mean across all neuron
mean_traces=zscore(nanmean(traces,2));  %% Here is use median instead of mean as a better population average
% High-speed periods
moving_period=h5read('processed_motion.h5','/moving_period');
moving_period=moving_period(idx);  % moving period in 0 and 1

Sampling_freq=20;

%2099.95 should be used for some sessions 
time_vect1  = 0:1/Sampling_freq:2099.9;% Original time vector (20Hz)
signal1=v ; % signal
time_vect2 = 0:1/1000:2099.9; % time vector for interpolation (1000Hz)
%signal1_Intp=interp1(time_vect1,signal1,time_vect2);%  interpolated signal
stim_vec=zeros(1,size(tracesF,1));
for i=1:length(stim_onsets)
timsel=(stim_onsets(i)-0:stim_onsets(i)+1200-1) -idx(1);
stim_vec(timsel)=1;
end

% Below is from start of imaging...
stim_onsets=LFP_data.Stim_onset;
stim_offsets=LFP_data.Stim_offset;
stim_onsets= stim_onsets/20*1000;
stim_offsets=stim_offsets/20*1000;
stim_vecLFP=zeros(1,size(LFP_data.LFP,1));
for i=1:length(LFP_data.LFP_stim_onset)
timsel=(stim_onsets(i)-0:stim_offsets(i)) ;
stim_vecLFP(timsel)=1;
end

LFP=LFP_data.LFP;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 vv=LFP(1:50:end);
 vv=vv(idx);
%% Align LFP and motion
start_frame=LFP_data.Start_Imaging;
stim_onsets=LFP_data.Stim_onset;
stim_offsets=LFP_data.Stim_offset;
delay_frames=start_frame+ idx(1);
delay_frame_LFP = ceil(delay_frames*1000/20);
%aligned_LFP = LFP(delay_frame_LFP:delay_frame_LFP+length(signal1_Intp)-1);
shifted_stim_onsets=(stim_onsets-idx(1))/20*1000;
shifted_stim_offsets=(stim_offsets-idx(1))/20*1000;
aligned_LFP=LFP(delay_frame_LFP:1:end); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim_vecLFP=stim_vecLFP(delay_frame_LFP:1:end);
%aligned_LFP=aligned_LFP(1:length(moving_period));

v_Intp=interp1(time_vect1,v(1:length(time_vect1)),time_vect2);%  interpolated signal
moving_period=interp1(time_vect1,moving_period(1:length(time_vect1)),time_vect2);%  interpolated signal
aligned_LFP=aligned_LFP(1:length(v_Intp));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim_vec=interp1(time_vect1,stim_vec(1:length(time_vect1)),time_vect2);%  interpolated signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matFiles = dir('*.mat') ; 
matFiles_name = {matFiles.name} ; 
idx = find(~cellfun(@isempty,strfind(matFiles_name,'plex'))) ;
file=matFiles_name(idx);
load(file{1})
% Check if the timestamps are from imaging start or starting LFP recording 
stim_TS= plx.Timestamp_stim;
% Seconds to timestamps 
stim_TS=stim_TS*1000;
% frame diff = 0.1 s = 100
stim_TS=stim_TS-delay_frame_LFP;

stim_vecTS=zeros(1,size(aligned_LFP,1));
stim_vecTS(round(stim_TS))=1;


aligned_LFP=zscore(aligned_LFP);
%%%%
Fn = 1000/2;FB=[ 58 63 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
LFP_A=filtfilt(B,A,   aligned_LFP);
aligned_LFP=aligned_LFP-LFP_A;

aligned_LFP=aligned_LFP-fastsmooth(aligned_LFP,10,1,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REMOVAL STIM PERIODS%%%%%%%%%
onsets2=find(stim_vecTS);
stim_trans=find(diff(onsets2)>100);
stim_trans=[ 0 onsets2(stim_trans)-1 onsets2(end)+1];
for iter1=rem_stim
    onsets2(onsets2>stim_trans(iter1) & onsets2< stim_trans(iter1+1))=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clear dataT2
mm2=0
sig=v_Intp';
sig=aligned_LFP;
;wind=40;
for id=1:length(onsets2)
    if onsets2(id)>wind  & onsets2(id)+wind <size(sig,1) 
        mm2=mm2+1;
   dataT2(:,mm2)= zscore(sig(onsets2(id)-wind:  onsets2(id)+wind  ));  
      dataT2(:,mm2)=   dataT2(:,mm2)-nanmean(   dataT2(:,mm2),1);
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   FS=1000;
        lfp_all=[];
    lfp_all.trial{1}(1,:)= (aligned_LFP);%zscore(v-fastsmooth(v,400,1,1));
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP'};
    
 cfg = []; %block_type == cfg.blk
cfg.method ='wavelet';%mtmfft'; %'mvar';
cfg.output ='fourier';
cfg.taper='hanning';
cfg.keeptapers ='yes';
cfg.keeptrials ='yes';
cfg.trials='all';cfg.tapsmofrq =6;%
cfg.channel=1;% ['all']; %chans=cfg.channel;
cfg.foi= [20:5:170];%[61:2:160]; %[2:2:40];%[61:2:160];    
cfg.t_ftimwin =[ones(1,length(cfg.foi))*1];
cfg.toi= lfp_all.time{1};   cfg.width =5;
freq2 = ft_freqanalysis(cfg, lfp_all);
cwt_out=squeeze(abs(freq2.fourierspctrm));
%size(cwt_out)


 clear dataP
mm2=0
%onsets2=find(stim_vecTS);
sig=aligned_LFP;
;wind2=40;
for id=1:length(onsets2)
    if onsets2(id)>wind2  & onsets2(id)+wind2 <size(sig,1) 
        mm2=mm2+1;
  dataP(:,:,mm2)=cwt_out(:,onsets2(id)-wind2:  onsets2(id)+wind2  );      
    end;    
end
%figure,imagesc(-wind:wind,freq2.freq,zscore(nanmean(dataP,3),[],2));axis xy


aEP= [aEP, nanmean(dataT2,2)];
%aEP2= [aEP2, nanmean(dataTM,2)];%
%aEP3= [aEP3, nanmean(dataTNM,2)];
mm=mm+1;
aEP_P(:,:,mm)=nanmean(dataP,3);
end
 %savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\SUDI\filt_delt_stats\'
pheight=160

% figure('COlor','w','Position',[300 300 250 180])
% plot((-wind:wind)./1000,nanmean(dataT2,2),'k')
% fill_error_area2((-wind:wind)./1000,nanmean(dataT2,2), nanstd(dataT2,[],2)./sqrt(size(dataT2,2)),[0.4 0.4 0.4])
% axis tight
% line([ 0 0], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
% line([ 0.1 0.1], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
% line([ -.1 -.1], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
% xlim([-0.11 0.13])


figure('COlor','w','Position',[300 300 250 180])
for id=1:size(aEP,2)
plot((-wind:wind)./1000,aEP(:,id)+id*0.3);hold on;
end
% line([ 0 0], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
% line([ 0.1 0.1], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
% line([ -.1 -.1], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
axis tight
xlim([-0.04 0.04])

figure('COlor','w','Position',[300 300 250 180])
plot((-wind:wind)./1000,nanmean(aEP,2),'k')
fill_error_area2((-wind:wind)./1000,nanmean(aEP,2), nanstd(aEP,[],2)./sqrt(size(aEP,2)),[0.4 0.4 0.4])
axis tight
%line([ 0 0], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
%line([ 0.1 0.1], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
%line([ -.1 -.1], [ -.07 .07],'Color', [ 0.1 0.1 0.1])
xlim([-0.02 0.02])
ylim([-0.06 0.06])
%%%

% figure('COlor','w','Position',[300 300 250 180])
% plot((-wind:wind)./1000,nanmean(aEP2,2),'r')
% fill_error_area2((-wind:wind)./1000,nanmean(aEP2,2), nanstd(aEP2,[],2)./sqrt(size(aEP2,2)),[0.4 0.4 0.4])
% plot((-wind:wind)./1000,nanmean(aEP3,2),'b')
% fill_error_area2((-wind:wind)./1000,nanmean(aEP3,2), nanstd(aEP3,[],2)./sqrt(size(aEP3,2)),[0.4 0.4 0.4])
% axis tight
% xlim([-0.11 0.13])
% 
%%%%%%%%%%%

% col=colormap(jet(15));
% figure('COlor','w','Position',[300 300 250 180])
% for id=1:size(aEP,2)
% plot((-wind:wind)./1000,aEP2(:,id)+id*0.7,'Color', col(id,:));hold on;
% plot((-wind:wind)./1000,aEP3(:,id)+id*0.7,'--','Color', col(id,:));hold on;
% end
% axis tight
% xlim([-0.11 0.13])
 % 
 % figure('COlor','w','Position',[300 300 250 180])
 % imagesc((-wind2:wind2)./1000,freq2.freq,bsxfun(@rdivide, nanmean(aEP_P,3), nanmean(nanmean(aEP_P,3),2)  ))
 % colormap(jet)
 % axis xy
 % set(gca,'Clim',[ 0.98 1.02])
 
 