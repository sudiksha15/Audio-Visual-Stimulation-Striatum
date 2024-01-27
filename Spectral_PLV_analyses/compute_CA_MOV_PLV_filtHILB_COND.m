clear all

%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
addpath(genpath('U:\eng_research_handata\EricLowet'))
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';

stim_freq=2; % 1=10Hz, 2=145Hz

aCOHs=[];aSPs=[];aCOHs2=[];aPHs=[];
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



for ses=1:length(allSES)
if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual\motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual\motion_corrected'])   
end

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










%%%%%%%%%%%%%%%%%%%%%%%
if 1
    FS=1000;
        lfp_all=[];
    lfp_all.trial{1}(1,:)= (aligned_LFP);%zscore(v-fastsmooth(v,400,1,1));
        lfp_all.trial{1}(2,:)= zscore(v_Intp);%-fastsmooth(v_Intp,10000,1,1));
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP', 'motion'};
    
  
      
%%%%%%%%%%%%%%%
deltP=lfp_all.trial{1}(2,:);
Fn = FS/2;FB=[ 2 4 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
%deltP= abs(hilbert(filtfilt(B,A,    deltP)));
deltphase= angle(hilbert(filtfilt(B,A,    deltP)));
thresD= prctile(deltP(moving_period==1),0);   %deltP(moving_period==1)
%deltphase(    moving_period==1)=NaN;                                  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavD=deltphase;wavDN=wavD;
shuf=randperm(size(wavD,2));
wavD2=wavD;%(:,shuf);wavDN2=wavD(:,shuf);
wavD(~(moving_period==1 &stim_vec==1 ))=NaN;
wavD2(~(moving_period==1&stim_vec==0 ))=NaN;


num_event_th=20;
mm=0;mm2=0;clear allCOHS allC allCOHS2 allSP allPH1
for neuron=1:size(tracesF,2)
 [ s_oasis, onsets2 ,onsets]= deconvSudi(tracesF(:,neuron));
 onsets2=onsets2.*50; onsets2(onsets2> length(wavD))=[];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(~isnan(wavD(1,onsets2))) >num_event_th
        mm=mm+1
     T=    abs(nanmean(exp(1i.*wavD(onsets2)),2)).^2;
     NT=sum(~isnan(wavD(onsets2)));
    allCOHS(mm)=    (((1/(NT-1))*((T.*NT-1)))) ;% 
    allPH1(mm)=    angle(nanmean(exp(1i.*wavD(onsets2)),2));
   % end
    %    if sum(~isnan(wavD2(1,onsets2))) >num_event_th
    mm2=mm2+1;  T=    abs(nanmean(exp(1i.*wavD2(onsets2)),2)).^2;
    NT2=sum(~isnan(wavD2(onsets2))); allCOHS2(mm2)=    (((1/(NT2-1))*((T.*NT2-1))));    
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%            
end

end

aCOHs= [aCOHs, allCOHS];
aCOHs2= [aCOHs2, allCOHS2];
aPHs= [aPHs, allPH1];
end
savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\'
pheight=160

 [h,p,ci,stats] = ttest(aCOHs,aCOHs2)
 
 %1620, >20, tstat=-1.5530, p=  0.1206; % stim_frque=2;
 % 683, >20, tstat=-1.07, p=0.28%stim_freq=2
% 904, >20, tstat=0.9, p=  0.3651 %stim_freq=1
aCOHs2=aCOHs2(~isnan(aCOHs2));
aCOHs=aCOHs(~isnan(aCOHs2));
aCOHs2=aCOHs2(~isinf(aCOHs2));
aCOHs=aCOHs(~isinf(aCOHs2));
V1=aCOHs';V2=aCOHs2';
p1=signrank(aCOHs,aCOHs2)
[h,p,ci,stats] = ttest(aCOHs,aCOHs2)
figure('COlor','w','Position', [ 300 150 130 pheight],'Renderer', 'painters')
%violinplot2(V1- V2,[1.3 ],'ViolinColor', [ 0.4 0.7 0.6; 0 0 0.9])
violinplot2(V1,[1.2 ],'ViolinColor', [ 0.3 0.5 0.9; 0 0 0.9])
violinplot2( V2,[1.9 ],'ViolinColor', [ 0.9 0.3 0.7; 0 0 0.9])
line([ 0.8 2.3], [ 0.2 0.2],'COlor', [ 1 1 1 ],'Linewidth',0.5)
axis tight;
set(gca,'Xtick',[1.2 1.9],'Xticklabel',{'Base','145 Hz'})
ylim([-0.2 0.9])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'CA_MOV_COND_STIM_PLV' num2str(stim_freq) '.pdf'])
%  
% figure('COlor','w','Position',[300 300 180 160])
% polarhistogram(aPHs,20,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',.3,'Normalization' , 'probability')
% axis tight
%   print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'CA__LFP_COND_STIM_PLV' num2str(stim_freq) '.pdf'])
%   
%   figure('COlor','w','Position',[300 300 180 160])
% polarhistogram(aPHs( aCOHs>prctile(aCOHs,80)),20,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',.3,'Normalization' , 'probability')
% axis tight
%   print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'CA__LFP_COND_STIM_PLV80perc' num2str(stim_freq) '.pdf'])
% % %
% %
% figure('COlor','w','Position',[300 300 180 160])
% hist(aPHs,20,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',.3,'Normalization' , 'probability')
% axis tight
%   print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'CA__LFP_STIM_PLVHIST' num2str(stim_freq) '.pdf'])

%   figure,plot(nanmean(aSPs,2))

