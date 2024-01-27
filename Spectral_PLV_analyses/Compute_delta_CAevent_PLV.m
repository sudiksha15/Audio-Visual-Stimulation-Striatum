clear all
close all

%% Ca-speed and Ca-LFP PLV - why fail to deconvolve trace??
%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
addpath(genpath('U:\eng_research_handata\EricLowet'))
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';

aCOHs=[];aSPs=[];aCOHs2=[];aPHs=[];aNT=[];corrLFP=[];
cha=1;% 2= movement, 1=LFP
for stim_freq=2  % 1= 10Hz, 2=145Hz
    
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

if stim_freq==1
ses_sel=[ 2:13];
else
  ses_sel=[ 1:11 13 14];  
end


for ses=  ses_sel


if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual\motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual\motion_corrected'])   
end




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

% resp_cells=h5read('Motion_resp_cells_sustained_10Hz.h5',allSES{4});
% resp_cells=squeeze(resp_cells);
% % Python to MATLAB
% resp_cells=resp_cells+1;
% 
% zscore(lfp_all.trial{1}(1,:))

Sampling_freq=20;
% Mean across all neuron
mean_traces=zscore(nanmean(traces,2));  %% Here is use median instead of mean as a better population average
% High-speed periods
moving_period=h5read('processed_motion.h5','/moving_period');
moving_period=moving_period(idx);  % moving period in 0 and 1

%2099.95 should be used for some sessions 
%time_vect1  = 0:1/Sampling_freq:2099.9;% Original time vector (20Hz)
%signal1=v ; % signal
%time_vect2 = 0:1/1000:2099.9; % time vector for interpolation (1000Hz)
%signal1_Intp=interp1(time_vect1,signal1,time_vect2);%  interpolated signal


LFP=LFP_data.LFP;
%asel=fastsmooth(zscore(fastsmooth(abs(hilbert(LFP)),5,1,1))>4,300,1,1);
%LFP(asel>0)=median(LFP)+randn(1, length(find(asel>0))).*std(LFP);
%figure(1)
%plot(LFP)
 vv=LFP(1:50:end); %down-sampling
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
aligned_LFP=LFP(delay_frame_LFP:50:end);
aligned_LFP=aligned_LFP(1:length(moving_period));
%%%%%%%%%%%%%%%%%%%%%%%
if 1
 
    lfp_all=[];FS=20
    lfp_all.trial{1}(1,:)= (aligned_LFP);%zLFP
        lfp_all.trial{1}(2,:)= zscore(v-fastsmooth(v,400,1,1)); %Motion signal
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP', 'motion'};
      
%%%%%%%%%%%%%%%
deltP=lfp_all.trial{1}(2,:);
Fn = FS/2;FB=[ 3 4 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
deltP= abs(hilbert(filtfilt(B,A,    deltP)));
thresD= prctile(deltP(moving_period==1),30);   
% thresD is a threshold on movement delta power as additional criterion. The idea is
%that to compute PLV only in time periods when there is suffcient movement delta
%power. Currently it is set to 0

corrLFP=[corrLFP,xcorr(zscore(v-fastsmooth(v,400,1,1)), aligned_LFP,2500,'Coeff')]
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  cfg = []; %block_type == cfg.blk
% cfg.method ='wavelet';%mtmfft'; %'mvar';
% cfg.output ='fourier';
% cfg.taper='hanning';
% cfg.keeptapers ='yes';
% cfg.keeptrials ='yes';
% cfg.trials='all';cfg.tapsmofrq =5;%
% cfg.channel= ['all']; %chans=cfg.channel;
% cfg.foi= [1:0.5:9.5];    cfg.t_ftimwin =[ones(1,length(cfg.foi))*1];
% cfg.toi= lfp_all.time{1};   cfg.width =6;
% freq2 = ft_freqanalysis(cfg, lfp_all);

[cwt_out,frs]=runcwt(lfp_all, [1 20],FS);

freq2.freq=frs;
wavA = abs(squeeze(cwt_out(1,cha,:,:)));
wavD = angle(squeeze(cwt_out(1,cha,:,:)));

%wavD=squeeze(angle(freq2.fourierspctrm(1,cha,:,:)));
shuf=randperm(size(wavD,2));
%wavD=wavD(:,end:-1:1);
%wavD2=wavD;wavD2(:,asel')=NaN;
wavD2=wavD(:,shuf);
wavD(:,~(deltP>thresD  & (moving_period==1)') )=NaN;
wavD2(:,~(deltP>thresD & (moving_period==1)'))=NaN;

% POW=nanmean(squeeze(abs(freq2.fourierspctrm(:,:,:,deltP>thresD))),2);
% POWB=nanmean(squeeze(abs(freq2.fourierspctrm(:,:,:,moving_period==1))),2);
% figure,plot(freq2.freq, POW)
% hold on,,plot(freq2.freq, POWB)

mm=0;clear allCOHS allC allCOHS2 allSP allPH1 allNT
for neuron=1:size(tracesF,2)
 [ s_oasis, onsets2 ,onsets]= deconvSudi(tracesF(:,neuron));
minimum_spikes=20;
    if sum(~isnan(wavD(1,onsets2))) >minimum_spikes
        mm=mm+1
        
        %%%%%%%% PLV %%%%%%%%%%%%%%
        T=    abs(nanmean(exp(1i.*wavD(:,onsets2)),2)).^2;
        NT=sum(~isnan(wavD(1,onsets2)));
        allCOHS(:,mm)=    (((1/(NT-1))*((T.*NT-1)))) ;% abs(nanmean(exp(1i.*wavD(:,onsets2)),2));
        T=    abs(nanmean(exp(1i.*wavD2(:,onsets2)),2)).^2;
        allCOHS2(:,mm)=    (((1/(NT-1))*((T.*NT-1)))); 
        allNT(mm)=NT;
        %%%%%%%%% PHASE %%%%%%%%%%%%%%%%%%%%%%%
        allPH1(:,mm)=    angle(nanmean(exp(1i.*wavD(:,onsets2)),2));
      
 clear dataT2
mm2=0;wind=40;
for id=1:length(onsets2)
    if onsets2(id)>wind  & onsets2(id)+wind <size(traces,1) 
        mm2=mm2+1;
   dataT2(:,mm2)= zscore(tracesF(onsets2(id)-wind:  onsets2(id)+wind,neuron  ));  
    end
    
end
allSP(:,mm)=nanmean(dataT2,2);
    end
end
end

if mm>0
    Dx=max(diff(allSP));

calc_transient_perc=0;
aCOHs= [aCOHs, allCOHS(:,Dx>prctile(Dx,calc_transient_perc))];
aCOHs2= [aCOHs2, allCOHS2(:,Dx>prctile(Dx,calc_transient_perc))];
aPHs= [aPHs, allPH1(:,Dx>prctile(Dx,calc_transient_perc))];
aSPs= [aSPs, allSP(:,Dx>prctile(Dx,calc_transient_perc))];
aNT=[aNT, allNT(Dx>prctile(Dx,calc_transient_perc))];
end
end
end
 
savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\'
% Threshold for number of spikes 
nt=20;
   figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
 fill_error_area2(freq2.freq,nanmean(aCOHs(:,aNT>nt),2), nanstd(aCOHs(:,aNT>nt),[],2)./sqrt(size(aCOHs(:,aNT>nt),2)),[ 0.5 0.5 0.5])
  fill_error_area2(freq2.freq,nanmean(aCOHs2(:,aNT>nt),2), nanstd(aCOHs2(:,aNT>nt),[],2)./sqrt(size(aCOHs2(:,aNT>nt),2)),[ 0.5 0.5 0.5])
     plot(freq2.freq,nanmean(aCOHs(:,aNT>nt),2),'r')
plot(freq2.freq,nanmean(aCOHs2(:,aNT>nt),2),'b','Color', [ 0.2 0.2 .2])
axis tight
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LFP_PLV_lock' num2str(stim_freq) '.pdf'])


%    figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
%  fill_error_area2(freq2.freq,nanmean(aCOHs(:,aNT>nt)-(aCOHs2(:,aNT>nt)),2), nanstd(aCOHs(:,aNT>nt)-(aCOHs2(:,aNT>nt)),[],2)./sqrt(size(aCOHs(:,aNT>nt),2)),[ 0.5 0.5 0.5])
%   %fill_error_area2(freq2.freq,nanmean(aCOHs2(:,aNT>nt),2), nanstd(aCOHs2(:,aNT>nt),[],2)./sqrt(size(aCOHs2(:,aNT>nt),2)),[ 0.5 0.5 0.5])
%      plot(freq2.freq,nanmean(aCOHs(:,aNT>nt)-(aCOHs2(:,aNT>nt)),2),'r','Linewidth',1)
% %plot(freq2.freq,nanmean(aCOHs2(:,aNT>nt),2),'b','Color', [ 0.2 0.2 .2])
% hold on,line([1 10], [ 0 0],'Color', [ 0 0 0])
% axis tight
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LFP_PLV_lock1' num2str(stim_freq) '.pdf'])

 
frsel=find(freq2.freq>3 & freq2.freq<4) 
%frsel=find(freq2.freq>0 & freq2.freq<1.5)  
delt_angs= circ_mean(aPHs(frsel,:));
delt_PLV= nanmean(aCOHs(frsel,:));
delt_PLV_SH= nanmean(aCOHs2(frsel,:));
%figure,rose(delt_angs)

[h,p,ci,stats] =ttest(delt_PLV,delt_PLV_SH)
p=signrank(delt_PLV,delt_PLV_SH)

   figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters') 
polarhistogram(delt_angs(delt_PLV>0.0),[-pi:0.33:pi],'FaceColor','red','Normalization','probability','Facealpha',0.3);hold on,
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LFP_PLV_PHASE' num2str(stim_freq) '.pdf'])


