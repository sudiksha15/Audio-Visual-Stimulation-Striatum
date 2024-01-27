close all
clear all
%% Calculate phase locking between stim triggers and LFP to check entrainment.

%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';
stim_freq=2;
aCOHs=[];aSPs=[];aCOHs2=[];aPHs=[];aNT=[];
   aAs=[];   aAs2=[];
   mm=0;
if stim_freq==1
% 10Hz
allSES{1}='\608448\01062020\608448'; allRES{1}='\608448_01062020'; % no LFP
allSES{2}='\608451\01062020\608451';allRES{2}='\608451_01062020';
allSES{3}='\608452\01062020\608452';allRES{3}='\608451_01062020';
allSES{4}='\608450\03092020\608450';allRES{4}='\608450_03092020';
allSES{5}='\602101\03092020\602101';allRES{5}='\602101_03092020';
allSES{6}='\611311\03092020\611311';allRES{6}='\611311_03092020';
allSES{7}='\608448\01172020\608448';allRES{7}='\608448_01172020';
allSES{8}='\608451\07082020\608451';allRES{8}='\608451_07082020';
allSES{9}='\608452\07082020\608452';allRES{9}='\608452_07082020';
allSES{10}='\611111\07082020\611111';allRES{10}='\611111_07082020';
allSES{11}='\602101\07082020\602101';allRES{11}='\602101_07082020';
allSES{12}='\612535\07082020\612535';allRES{12}='\612535_07082020';
allSES{13}='\615883\02052021\10Hz\615883';allRES{13}='\615883_02052021';
allSES{14}='\615883\03122021\10Hz\615883';allRES{14}='\615883_03122021';
removal_stim{1}= [];
removal_stim{2}= [];
removal_stim{3}= [];
removal_stim{4}= [];
removal_stim{5}= [];
removal_stim{6}= [];
removal_stim{7}= [];
removal_stim{8}= [];
removal_stim{9}= [];
removal_stim{10}= [];
removal_stim{11}= [];
removal_stim{12}= [];
removal_stim{13}= [];
removal_stim{14}= [];
else
% 140Hz
allSES{1}='/608448/01102020/608448';allRES{1}='/608448_01102020'; 
allSES{2}='/608451/01102020/608451';allRES{2}='/608451_01102020';
allSES{3}='/608452/01102020/608452';allRES{3}='/608452_01102020';
allSES{4}='/602101/03112020/602101';allRES{4}='/602101_03112020';
allSES{5}='/611311/03112020/611311';allRES{5}='/611311_03112020';
allSES{6}='/611111/03112020/611111';allRES{6}='/611111_03112020';
allSES{7}='/608450/01172020/608450';allRES{7}='/608450_01172020';
allSES{8}='/608451/07012020/608451';allRES{8}='/608451_07012020';
allSES{9}='/608452/07012020/608452';allRES{9}='/608452_07012020';
allSES{10}='/611111/07012020/611111';allRES{10}='/611111_07012020';
allSES{11}='/612535/07022020/612535';allRES{11}='/612535_07022020';
allSES{12}='/602101/07022020/602101'; allRES{12}='/602101_07022020';%exclude LFP
allSES{13}='/615883/02052021/145Hz/615883'; allRES{13}='/615883_02052021';
allSES{14}='/615883/03122021/145Hz/615883';allRES{14}='/615883_03122021';
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
end
if stim_freq==1
ses_sel=[ 2:14]%14];
else
  ses_sel=[1:5 9 14];
end

for ses= ses_sel
    ses
if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual\motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual\motion_corrected'])   
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



%  figure,plot(aligned_LFP); hold on,plot(stim_vecTS,'.')
%  hold on,plot(moving_period.*.5)

 %figure,plot(stim_vecLFP); hold on,plot(stim_vecTS,'.')
 %hold on,plot(moving_period.*.5)
%%%%%%%%%%%%%%%%%%%%%


% v_Intp=v_Intp(1:5:end);
% moving_period=moving_period(1:5:end);
% aligned_LFP=aligned_LFP(1:5:end);stim_vec=stim_vec(1:5:end);
% stim_vecLFP=stim_vecLFP(1:5:end);
% 



%%%%%%%%%%%%%%%%%%%%%%%
if 1

        FS=1000;
 Fn = FS/2;FB=[ 7 8 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
thetaS= ((filtfilt(B,A,   aligned_LFP)));
 %aligned_LFP= aligned_LFP-thetaS;
   aligned_LFP=zscore(aligned_LFP);
    devianz=( abs(   aligned_LFP)>4);
    lfp_all=[];
    lfp_all.trial{1}(1,:)= (aligned_LFP);%zscore(v-fastsmooth(v,400,1,1));
        lfp_all.trial{1}(2,:)= zscore(v_Intp);%-fastsmooth(v_Intp,10000,1,1));
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP', 'motion'};
      
%%%%%%%%%%%%%%%
deltP=lfp_all.trial{1}(2,:);
Fn = FS/2;FB=[ 2.5 5 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
deltP= abs(hilbert(filtfilt(B,A,    deltP)));
deltPh= angle(hilbert(filtfilt(B,A,    deltP)));
thresD= prctile(deltP(moving_period==1),0);   %deltP(moving_period==1)
%    
%%%%%%%%%%%%
dtrigs=zeros(1, length(deltPh) );
 timer1=0;
for xi=1:length(deltPh)
     timer1=timer1+1;
    if deltPh(xi)>0  & deltPh(xi)<0.2  &  timer1>40
      dtrigs(xi)=1;
      timer1=0;
    else
        
    end
end

%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cfg = []; %block_type == cfg.blk
cfg.method ='mtmconvol';%mtmfft'; %'mvar';
cfg.output ='fourier';
cfg.taper='hanning';
cfg.keeptapers ='yes';
cfg.keeptrials ='yes';
cfg.trials='all';cfg.tapsmofrq =6;%
cfg.channel=1;% ['all']; %chans=cfg.channel;
if stim_freq==1
cfg.foi= [2:2:40];%[61:2:160]; %[2:2:40];%[61:2:160]; 
else
cfg.foi=[61:2:160];end
cfg.t_ftimwin =[ones(1,length(cfg.foi))*1];
cfg.toi= lfp_all.time{1};   cfg.width =6;
freq2 = ft_freqanalysis(cfg, lfp_all);
cwt_out=freq2.fourierspctrm;
%[cwt_out,frs]=runcwt(lfp_all, [1 20],FS);
cha=1;
%cwt_out(:,:,:,fastsmooth(devianz,50,1,1)>0)=NaN;
%freq2.freq=frs;
wavD = angle(squeeze(cwt_out(1,cha,:,:)));
%wavA = abs(squeeze(cwt_out(1,1,:,:)));%%%%
%wavD2 = angle(squeeze(cwt_out(1,1,:,end:-1:1)));
% wavD(:,~(  deltP>thresD  & (stim_vec==0) & (moving_period==1)))=NaN;
% wavD2(:,~( deltP>thresD  & (stim_vec==0) &  (moving_period==1)))=NaN;
% wavA(:,~(  deltP>thresD  & (stim_vec==0) & (moving_period==1)))=NaN;


%% REMOVAL STIM PERIODS%%%%%%%%%
onsets2=find(stim_vecTS);
stim_trans=find(diff(onsets2)>100);
stim_trans=[ 0 onsets2(stim_trans)-1 onsets2(end)+1];
for iter1=rem_stim
    onsets2(onsets2>stim_trans(iter1) & onsets2< stim_trans(iter1+1))=[];
end

 stimS=onsets2  ;

  sP=length( stimS);
 rsel= randperm(length(stim_vecTS));

%  figure,imagesc(bsxfun(@rdivide,nanmean(data,3),nanmean(nanmean(data,3),2)))
%  axis xy
%%%%%%%%%%%%%
%PLV=abs(nanmean(exp(1i.*wavD(:,stim_vecTS==1)),2));

PLV=abs(nanmean(exp(1i.*wavD(:, stimS)),2));
PLVN= abs(nanmean(exp(1i.*wavD(:,rsel(1:sP))),2));
%NT= sum(~isnan(wavD(4,1:100:end)));
%NT2= sum(~isnan(wavDN(4,1:100:end)));
%PLVu=(((1/(NT-1))*(((PLV.^2).*NT-1))));
%PLVu2=(((1/(NT2-1))*(((PLVN.^2).*NT2-1))));
%figure,plot(freq2.freq,PLV),hold on,plot(freq2.freq,PLVN)
%calc_transient_perc=25;
 aCOHs= [aCOHs, PLV];  % stim off
  aCOHs2= [aCOHs2, PLVN]; % stim on
   %aAs= [aAs, nanmean(wavA,2)];  % stim off
  %aAs2= [aAs2, nanmean(wavAN,2)]; % stim on
  mm=mm+1;
end
end
 
savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\';
% nt=20;
   figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
 fill_error_area2(freq2.freq,nanmean(aCOHs,2), nanstd(aCOHs,[],2)./sqrt(size(aCOHs,2)),[ 0.5 0.5 0.5])
  fill_error_area2(freq2.freq,nanmean(aCOHs2,2), nanstd(aCOHs2,[],2)./sqrt(size(aCOHs2,2)),[ 0.5 0.5 0.5])
     plot(freq2.freq,nanmean(aCOHs,2),'r','Linewidth',1)
plot(freq2.freq,nanmean(aCOHs2,2),'k','Color', [ 0.1 0.1 .1],'Linewidth',1)
axis tight
ylim([0,0.2])
   print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LFP__STIM_PLV' num2str(stim_freq) '.pdf'])
% 

if stim_freq==1
fsel=find(freq2.freq>9 & freq2.freq <11); % 10Hz
else
    fsel=find(freq2.freq>144 & freq2.freq <146); % 10Hz
end

V1=nanmean(aCOHs(fsel,:),1);
V2=nanmean(aCOHs2(fsel,:),1);

[h,p,ci,stats] =ttest(V1,V2)
p2=signrank(V1,V2);
p2
% 10Hz, df=12,tstat: 4.2176,p= 0.0012;
% 145Hz, df=5,tstat: 3.6783,  p= 0.0143