clear all
close all
% Getting CFC stats and violin plots

% Why is gamma not working??

%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';
aCOHs=[];aSPs=[];aCOHs2=[];aPHs=[];aNT=[];PLVtrue=[];PLVshuf=[];PLVtrue2=[];PLVshuf2=[];
   aAs=[];   aAs2=[];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   freqband =2; %1=gamma, 2=beta
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mm=0;
   for stim_freq=1:2;
if stim_freq==1
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
else
% 140Hz
allSES{1}='\608448\01102020\608448';allRES{1}='\608448_01102020'; 
allSES{2}='\608451\01102020\608451';allRES{2}='\608451_01102020';
allSES{3}='\608452\01102020\608452';allRES{3}='\608452_01102020';
allSES{4}='\602101\03112020\602101';allRES{4}='\602101_03112020';
allSES{5}='\611311\03112020\611311';allRES{5}='\611311_03112020';
allSES{6}='\611111\03112020\611111';allRES{6}='\611111_03112020';
allSES{7}='\608450\01172020\608450';allRES{7}='\608450_01172020';
allSES{8}='\608451\07012020\608451';allRES{8}='\608451_07012020';
allSES{9}='\608452\07012020\608452';allRES{9}='\608452_07012020';
allSES{10}='\611111\07012020\611111';allRES{10}='\611111_07012020';
allSES{11}='\612535\07022020\612535';allRES{11}='\612535_07022020';
allSES{12}='\602101\07022020\602101'; allRES{12}='\602101_07022020';%exclude LFP
allSES{13}='\615883\02052021\145Hz\615883'; allRES{13}='\615883_02052021';
allSES{14}='\615883\03122021\145Hz\615883';allRES{14}='\615883_03122021';
end
if stim_freq==1
ses_sel=[ 2:14];
else
  ses_sel=[ 1:11 13 14];  
end

for ses= ses_sel
    ses
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


stim_vecLFP=zeros(1,size(LFP_data.LFP,1));
for i=1:length(LFP_data.LFP_stim_onset)
timsel=(LFP_data.LFP_stim_onset(i)-0:LFP_data.LFP_stim_offset(i)) ;
stim_vecLFP(timsel)=1;
end

LFP=LFP_data.LFP;
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
aligned_LFP=LFP(delay_frame_LFP:1:end);
stim_vecLFP=stim_vecLFP(delay_frame_LFP:1:end);
%aligned_LFP=aligned_LFP(1:length(moving_period));

v_Intp=interp1(time_vect1,v(1:length(time_vect1)),time_vect2);%  interpolated signal
moving_period=interp1(time_vect1,moving_period(1:length(time_vect1)),time_vect2);%  interpolated signal
aligned_LFP=aligned_LFP(1:length(v_Intp));
stim_vec=interp1(time_vect1,stim_vec(1:length(time_vect1)),time_vect2);%  interpolated signal

v_Intp=v_Intp(1:5:end);
moving_period=moving_period(1:5:end);
aligned_LFP=aligned_LFP(1:5:end);stim_vec=stim_vec(1:5:end);
stim_vecLFP=stim_vecLFP(1:5:end);
%%%%%%%%%%%%%%%%%%%%%%%
if 1

        FS=200;
        if freqband ==1
 Fn = FS/2;FB=[40 100 ];
        else
 Fn = FS/2;FB=[15 30 ];
        end
[B, A] = butter(3, [min(FB)/Fn max(FB)/Fn]);
thetaS= abs(hilbert(filtfilt(B,A,   aligned_LFP)));
 Fn = FS/2;FB=[2 4 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
thetaA= angle(hilbert(filtfilt(B,A,   thetaS)));

 %aligned_LFP= aligned_LFP-thetaS;
   aligned_LFP=zscore(aligned_LFP);
    devianz=( abs(   aligned_LFP)>4);
    lfp_all=[];
    lfp_all.trial{1}(1,:)= zscore(aligned_LFP);%zscore(v-fastsmooth(v,400,1,1));
        lfp_all.trial{1}(2,:)= zscore(v_Intp);%-fastsmooth(v_Intp,10000,1,1));
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP', 'motion'};
      
%%%%%%%%%%%%%%%
deltP1=lfp_all.trial{1}(1,:);
Fn = FS/2;FB=[ 2  4 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
deltP= abs(hilbert(filtfilt(B,A,    deltP1)));
deltPh= angle(hilbert(filtfilt(B,A,    deltP1)));
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

% %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
V1b=deltPh((  (stim_vec==1) &  (moving_period==1)));
V2b=thetaA(( (stim_vec==1) &  (moving_period==1)));
AB=thetaS(( (stim_vec==1) &  (moving_period==1)));
truePLV2=circ_r(circ_dist(V1b,V2b')');
shuf_PLV2=circ_r(circ_dist(V1b',V2b(randperm(length(V2b)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V1=deltPh((  (stim_vec==0) &  (moving_period==1)));
V2=thetaA(( (stim_vec==0) &  (moving_period==1)));
A=thetaS(( (stim_vec==0) &  (moving_period==1)));

llenght=length(V1b);
rdsel=randperm(length(V1));
%V1=V1(rdsel(1:llenght));V2=V2(rdsel(1:llenght));

truePLV=circ_r(circ_dist(V1,V2')');
shuf_PLV=circ_r(circ_dist(V1',V2(randperm(length(V2)))));

%%
 % [ dataFB,p_n]=delta_phase_waveS(AB, V1b );
 % [ dataF,p_n]=delta_phase_waveS(A, V1 );

  mm=mm+1;
  PLVtrue=[PLVtrue,truePLV];  PLVshuf=[PLVshuf,shuf_PLV]; % nostim
    PLVtrue2=[PLVtrue2,truePLV2];  PLVshuf2=[PLVshuf2,shuf_PLV2];%stim
  % aMAP(:,1,mm)=  dataFB-nanmean(dataFB);%stim
  %  aMAP(:,2,mm)=  dataF-nanmean(dataF);%nostim
 aNT=[aNT, stim_freq];
end
end
   end
savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\';
pheight=160

 %  figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
   
V1=(PLVtrue-PLVshuf)';% nostim
V2=(PLVtrue2-PLVshuf2)';%stim
  figure('COlor','w','Position', [ 300 400 130 pheight],'Renderer', 'painters') 
b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
ylim([-0.00 0.12])
 % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'mov_LFP_CFC_Beta_' num2str(stim_freq) '.pdf'])
%
cc=2;
V1=(PLVtrue(aNT==cc)-PLVshuf(aNT==cc))';% nostim
V2=(PLVtrue2(aNT==cc)-PLVshuf2(aNT==cc))';%stim
figure('COlor','w','Position', [ 300 400 190 pheight],'Renderer', 'painters')
line([ 0.8 2.3], [ -0.01 -0.01],'COlor', [ 1 1 1 ],'Linewidth',0.5)
violinplot2(V1,[1.2 ],'ViolinColor', [ 0.3 0.5 0.9; 0 0 0.9])
violinplot2(V2,[1.9 ],'ViolinColor', [ 0.9 0.3 0.7; 0 0 0.9]);axis tight
line([ 0.9 1.5], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
line([ 1.6 2.2], [ 0  0],'COlor', [ 0 0 0 ],'Linewidth',1)
ylim([-0.3 .3]);set(gca,'Xticklabel',[],'Xtick',[])
ylim([-0.01 0.25])
 % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'VIOL_mov_LFP_CFC_' num2str(freqband) '_' num2str(cc) '.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%  STATS %%%%%%%%%%%%%%%
 
cc=2; %effect of baseline vs shuffle, cc=1 (10Hz), cc=2 (145Hz)
%[h,p,ci,stats] = ttest(PLVtrue2(aNT==cc),PLVshuf2(aNT==cc))
[p, h, stats] = signrank(PLVtrue(aNT==cc),PLVshuf(aNT==cc))

cc=2; %effect of stim vs shuffle, cc=1 (10Hz), cc=2 (145Hz)
%[h,p,ci,stats] = ttest(PLVtrue2(aNT==cc),PLVshuf2(aNT==cc))
[p, h, stats] = signrank(PLVtrue2(aNT==cc),PLVshuf2(aNT==cc))

% beta, df=12, p=   2.4414e-04(10Hz), p= 2.4414e-04 (145Hz)
% gamma, df=12, p= 2.4414e-04 (10Hz), p= 2.4414e-04

cc=2; %effect of stim , cc=1 (10Hz), cc=2 (145Hz)
%[h,p,ci,stats] = ttest(PLVtrue(aNT==cc)-PLVshuf(aNT==cc),PLVtrue2(aNT==cc)-PLVshuf2(aNT==cc))
[p, h, stats] = signrank(PLVtrue(aNT==cc)-PLVshuf(aNT==cc),PLVtrue2(aNT==cc)-PLVshuf2(aNT==cc))

%beta, df=12,  p=  0.0327 (10Hz) /  0.2439  (145Hz)
%gamma, df=12,    p=  0.4548(10Hz) / 0.7354(145Hz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc=2;cc1=1; % stim type 10Hz vs. 145Hz
%[h,p,ci,stats] = ttest(PLVtrue2(aNT==cc)-PLVshuf2(aNT==cc),PLVtrue2(aNT==cc1)-PLVshuf2(aNT==cc1))
[p, h, stats] = signrank(PLVtrue2(aNT==cc)-PLVshuf2(aNT==cc),PLVtrue2(aNT==cc1)-PLVshuf2(aNT==cc1))
%beta, df=12,   p = 0.1465
% gamma, df=12, p=0.7354


% cc=2;cc1=1; % effect of stim type
% [h,p,ci,stats] = ttest(PLVtrue(aNT==cc)-PLVshuf(aNT==cc),PLVtrue(aNT==cc1)-PLVshuf(aNT==cc1))
% %Basleine
% %beta, df=12,   0.8372,  p=   0.4188
% %gamma, df=12,   0.3761,  p=    0.7134
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
