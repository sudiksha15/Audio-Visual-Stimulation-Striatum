clear all
close all
%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
addpath(genpath('U:\eng_research_handata\EricLowet'))
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';

stim_freq=2;
aCOHs=[];aSPs=[];aCOHs2=[];aPHs=[];aNT=[];MOT_respALL=[];Stim_respALL=[];aNeuron=[];aPLV=[];   aPph=[];
if stim_freq==1
% 10Hz
allSES{1}='\608448\01062020\608448'; allRES{1}='/608448_01062020'; % no LFP
allSES{2}='\608451\01062020\608451';allRES{2}='/608451_01062020';
allSES{3}='\608452\01062020\608452';allRES{3}='/608451_01062020';
allSES{4}='\608450\03092020\608450';allRES{4}='/608450_03092020';
allSES{5}='\602101\03092020\602101';allRES{5}='/602101_03092020';
allSES{6}='\611311\03092020\611311';allRES{6}='/611311_03092020';
allSES{7}='\608448\01172020\608448';allRES{7}='/608448_01172020';
allSES{8}='\608451\07082020\608451';allRES{8}='/608451_07082020';
allSES{9}='\608452\07082020\608452';allRES{9}='/608452_07082020';
allSES{10}='\611111\07082020\611111';allRES{10}='/611111_07082020';
allSES{11}='\602101\07082020\602101';allRES{11}='/602101_07082020';
allSES{12}='\612535\07082020\612535';allRES{12}='/612535_07082020';
allSES{13}='\615883\02052021\10Hz\615883';allRES{13}='/615883_02052021';
allSES{14}='\615883\03122021\10Hz\615883';allRES{14}='/615883_03122021';
else
% 140Hz
allSES{1}='\608448\01102020\608448';allRES{1}='/608448_01102020'; 
allSES{2}='\608451\01102020\608451';allRES{2}='/608451_01102020';
allSES{3}='\608452\01102020\608452';allRES{3}='/608452_01102020';
allSES{4}='\602101\03112020\602101';allRES{4}='/602101_03112020';
allSES{5}='\611311\03112020\611311';allRES{5}='/611311_03112020';
allSES{6}='\611111\03112020\611111';allRES{6}='/611111_03112020';
allSES{7}='\608450\01172020\608450';allRES{7}='/608450_01172020';
allSES{8}='\608451\07012020\608451';allRES{8}='/608451_07012020';
allSES{9}='\608452\07012020\608452';allRES{9}='/608452_07012020';
allSES{10}='\611111\07012020\611111';allRES{10}='/611111_07012020';
allSES{11}='\612535\07022020\612535';allRES{11}='/612535_07022020';
allSES{12}='\602101\07022020\602101'; allRES{12}='/602101_07022020';%exclude LFP
allSES{13}='\615883\02052021\145Hz\615883'; allRES{13}='/615883_02052021';
allSES{14}='\615883\03122021\145Hz\615883';allRES{14}='/615883_03122021';
end

if stim_freq==1
ses_sel=[13];
else
  ses_sel=[ 13 ];  
end

for ses=ses_sel
if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual\motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual\motion_corrected'])   
end

%cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Sudi_Sridhar\608452\07012020\608452_145Hz_AudioVisual\motion_corrected')



% h5d
motion=h5read('processed_motion.h5','/raw_speed_trace');
tracesF=h5read('processed_trace.h5','/trace');
traces=h5read('processed_trace.h5','/onset_binary_trace');
kk=h5read('processed_neuron_positions.h5','/rois/ROImasks');
load('LFP_ts.mat')
stim_onsets=LFP_data.Stim_onset;
stim_offsets=LFP_data.Stim_offset;
v=motion;
idx=find(isnan(v)==0);
%Remove NANs
v = v(~isnan(v));  % motion signal
% Align traces and motiontraces=traces(idx,:);
traces=traces(idx,:);tracesF=tracesF(idx,:);
if stim_freq==1
path_resp='\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Sudi_Sridhar\Resp_Neurons\'
resp_cells=h5read([ path_resp 'Motion_resp_cells_sustained_10Hz_final.h5'],allRES{ses});
resp_cells=squeeze(resp_cells);
 resp_cells=resp_cells+1;
resp_cellsST=h5read([ path_resp 'Stim_resp_cells_sustained_10Hz_final.h5'],[allRES{ses} '_dec_stn']);
resp_cellsST=squeeze(resp_cellsST);
 resp_cellsST=resp_cellsST+1;
else
  path_resp='\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Sudi_Sridhar\Resp_Neurons\'
resp_cells=h5read([ path_resp 'Motion_resp_cells_sustained_145Hz_final.h5'],allRES{ses});
resp_cells=squeeze(resp_cells);
 resp_cells=resp_cells+1;
 % inc_mov, dec_mov, inc_stn, dec_stn
resp_cellsST=h5read([ path_resp 'Stim_resp_cells_sustained_145Hz_final.h5'],[allRES{ses} '_dec_stn']);
resp_cellsST=squeeze(resp_cellsST);
 resp_cellsST=resp_cellsST+1;  
    
end

 
 
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


LFP=LFP_data.LFP;
%figure(1)
%plot(LFP)
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
aligned_LFP=LFP(delay_frame_LFP:50:end);
aligned_LFP=zscore(aligned_LFP(1:length(moving_period)));
aligned_LFP(abs(zscore(aligned_LFP))>6)=0;
%%%%%%%%%%%%%%%%%%%%%%%
if 1
 
    lfp_all=[];FS=20
    lfp_all.trial{1}(1,:)= (aligned_LFP);%zscore(v-fastsmooth(v,400,1,1));
        lfp_all.trial{1}(2,:)= zscore(v-fastsmooth(v,400,1,1));
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP', 'motion'};
      
%%%%%%%%%%%%%%%
deltLFP=lfp_all.trial{1}(1,:);
deltP1=lfp_all.trial{1}(2,:);
Fn = FS/2;FB=[ 3 4];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
deltP= abs(hilbert(filtfilt(B,A,    deltP1)));
deltPhase= angle(hilbert(filtfilt(B,A,    deltP1)));
deltPhaseLFP= angle(hilbert(filtfilt(B,A,    deltLFP)));

deltsignal= ((filtfilt(B,A,    deltP1)));
deltsignalLFP= ((filtfilt(B,A,    deltLFP)));
thresD= prctile(deltP(moving_period==1),0);   %deltP(moving_period==1)
%                                         
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  cfg = []; %block_type == cfg.blk
% cfg.method ='wavelet';%mtmfft'; %'mvar';
% cfg.output ='fourier';
% cfg.taper='hanning';
% cfg.keeptapers ='yes';
% cfg.keeptrials ='yes';
% cfg.trials='all';cfg.tapsmofrq =6;%
% cfg.channel= ['all']; %chans=cfg.channel;
% cfg.foi= [1:0.5:9.5];    cfg.t_ftimwin =[ones(1,length(cfg.foi))*1];
% cfg.toi= lfp_all.time{1};   cfg.width =6;
% freq2 = ft_freqanalysis(cfg, lfp_all);

[cwt_out,frs]=runcwt(lfp_all, [1 8],FS);
cha=2;
freq2.freq=frs;
wavA = abs(squeeze(cwt_out(1,cha,:,:)));
wavD = angle(squeeze(cwt_out(1,cha,:,:)));

%wavD=squeeze(angle(freq2.fourierspctrm));
shuf=randperm(size(wavD,2));
%wavD=wavD(:,end:-1:1);
wavD2=wavD(:,shuf);
wavD(:,~(deltP>thresD  & (moving_period==1)'))=NaN;
wavD2(:,~(deltP>thresD & (moving_period==1)'))=NaN;

% POW=nanmean(squeeze(abs(freq2.fourierspctrm(:,:,:,deltP>thresD))),2);
% POWB=nanmean(squeeze(abs(freq2.fourierspctrm(:,:,:,moving_period==1))),2);
% figure,plot(freq2.freq, POW)
% hold on,,plot(freq2.freq, POWB)
MOT_resp=[];Stim_resp=[];
mm=0;clear allCOHS allC allCOHS2 allSP allPH1 allNT
for neuron=1:size(tracesF,2)
 [ s_oasis, onsets2 ,onsets]= deconvSudi(tracesF(:,neuron));
minimum_spikes=40;
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
   
    
  if sum(neuron==resp_cells)>0
MOT_resp=[MOT_resp, 1];else; MOT_resp=[MOT_resp, 0];  end
if sum(neuron==resp_cellsST)
Stim_resp=[Stim_resp, 1];else; Stim_resp=[Stim_resp, 0];  end
 
   
    aNeuron=[aNeuron, neuron];aPLV=[aPLV,  circ_r(deltPhase(onsets2))]; 
    aPph=[aPph,  circ_mean(deltPhase(onsets2)')]; end
end
end

if mm>0
    Dx=max(diff(allSP));

calc_transient_perc=0;
aCOHs= [aCOHs, allCOHS(:,Dx>=prctile(Dx,calc_transient_perc))];
aCOHs2= [aCOHs2, allCOHS2(:,Dx>=prctile(Dx,calc_transient_perc))];
aPHs= [aPHs, allPH1(:,Dx>=prctile(Dx,calc_transient_perc))];
aSPs= [aSPs, allSP(:,Dx>=prctile(Dx,calc_transient_perc))];
aNT=[aNT, allNT(Dx>=prctile(Dx,calc_transient_perc))];
MOT_respALL=[MOT_respALL, MOT_resp(Dx>=prctile(Dx,calc_transient_perc))];
Stim_respALL=[Stim_respALL, Stim_resp(Dx>=prctile(Dx,calc_transient_perc))];

end
end
 
savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\'

traceFAC=5;%5
neuron=aNeuron(38), % 2,38
 [ s_oasis, onsets2 ,onsets]= deconvSudi(tracesF(:,neuron));
 sp=zeros(1,length(moving_period));
 sp(onsets2)=1;
figure('COlor','w','Position', [ 300 400 500 180],'Renderer', 'painters')
plot((1:length(tracesF))./1200,tracesF(:,neuron).*traceFAC,'k','Linewidth',1)
%hold on,plot(onsets2./1200,ones(1,length(onsets2)),'.r','Markersize',10)
%hold on,plot(deltPhase)
hold on,plot((1:length(tracesF))./1200,moving_period.*4 -22,'k','LineWidth',1)
hold on,plot((1:length(tracesF))./1200,zscore(deltsignal).*2-12,'Color', [.5 .5 .5])
hold on,plot((1:length(tracesF))./1200,zscore(deltsignalLFP).*1-27,'Color', [.5 .5 .9])
for idx=1:length(onsets2)
line([onsets2(idx)./1200 onsets2(idx)./1200], [-30 10],'COlor', [1 0 0])
end
axis tight
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Example_CA_MOTION_whole_1' '.pdf'])


for idx2=[26 36 ]%20:40
figure('COlor','w','Position', [ 300 400 200 180],'Renderer', 'painters')
plot((1:length(tracesF))./1200,tracesF(:,neuron).*traceFAC,'k','Linewidth',1)
hold on,plot((1:length(tracesF))./1200,zscore(deltsignal).*2-12,'Color', [0.6 0.5 0.5])
hold on,plot((1:length(tracesF))./1200,zscore(deltsignalLFP).*1-23,'Color', [.5 .5 .9])
for idx=1:length(onsets2)
line([onsets2(idx)./1200 onsets2(idx)./1200], [-30 10],'COlor', [1 0 0])
end
axis tight
xlim([onsets2(idx2)./1200-0.02 onsets2(idx2)./1200+0.06])
title(num2str(idx2))
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Example_CA_MOTION_zoom_' num2str(idx2) '.pdf'])

end

figure('COlor','w','Position', [ 300 400 150 90],'Renderer', 'painters')
%subplot(1,2,1)
polarhistogram(deltPhase(sp & moving_period'),14,'Normalization','Probability','FaceColor', [0.6 0.5 0.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Example_CA_MOTION_HIST''.pdf'])
figure('COlor','w','Position', [ 300 400 150 90],'Renderer', 'painters')
polarhistogram(deltPhaseLFP(sp & moving_period'),14,'Normalization','Probability','FaceColor', [0.5 0.5 0.8])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Example_CA_LFP_HIST''.pdf'])


pT=deltPhase>0;
pt=find(diff(pT)>0);
% figure,hold on,plot(zscore(deltsignal));hold on,plot(pt,ones(1,length(pt)),'.r','Markersize',10)
%%%%%%%%%%%%%%%
rast=zeros(1,length(tracesF));rast(onsets2)=1;
 clear dataT2 dataT dataT3
mm2=0;wind=7;
for id=1:2:length(pt)
    if pt(id)>wind  & pt(id)+wind <size(traces,1)  & moving_period(pt(id))>0 
      if min(abs(pt(id)-onsets2))<10
        mm2=mm2+1;
   dataT2(:,mm2)= (tracesF(pt(id)-wind:  pt(id)+wind,neuron  ));
       dataT(:,mm2)= (rast(pt(id)-wind:  pt(id)+wind  ));
         dataT3(:,mm2)= (deltsignal(pt(id)-wind:  pt(id)+wind  ));
      end
    end
    
end

% [b1 b2]=max((dataT))
% [c1 c2]=sort(b2)
% ax=-wind:wind;
% figure,imagesc(-wind:wind,[],(dataT(:,:))')
% for id =1:size(dataT,2)
%    plot(ax(find(dataT(:,id))), ones(1,length(find(dataT(:,id)))) +id,'.k','Markersize',8);hold on,
% end
% hold on,plot(-wind:wind,fastsmooth(nanmean(dataT,2),3,1,1).*300,'k')
% hold on,plot(-wind:wind,zscore(nanmean(dataT3,2)).*10+8,'m')
% colormap(jet);axis xy
% axis tight


figure, imagesc(freq2.freq,[],aCOHs');colormap(jet)


figure('COlor','w','Position', [ 300 400 200 180],'Renderer', 'painters')
polarhistogram(circ_dist(deltPhase',deltPhaseLFP'),14,'Normalization','Probability','FaceColor', [0.5 0.5 0.4])
axis tight
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Example_MOV_LFP_HIST''.pdf'])



t = Tiff('full_projection_m_f_max_min_batch.tif','r');
imageData = read(t);
if 1
figure
imagesc(double(imageData(1:600,1:600))./100)
colormap(gray)
set(gca,'Clim', [-10 80]); hold on,
contour(double(kk(:,:,neuron)),'r')
end
% [z1 z2 z3]=find(kk>0.5)
% 
% figure
% [z1 z2]=contour(double(kk(:,:,neuron)))
% 
nt=20;
   figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
 fill_error_area2(freq2.freq,nanmean(aCOHs(:,aNT>nt),2), nanstd(aCOHs(:,aNT>nt),[],2)./sqrt(size(aCOHs(:,aNT>nt),2)),[ 0.5 0.5 0.5])
  fill_error_area2(freq2.freq,nanmean(aCOHs2(:,aNT>nt),2), nanstd(aCOHs2(:,aNT>nt),[],2)./sqrt(size(aCOHs2(:,aNT>nt),2)),[ 0.5 0.5 0.5])
     plot(freq2.freq,nanmean(aCOHs(:,aNT>nt),2),'r')
plot(freq2.freq,nanmean(aCOHs2(:,aNT>nt),2),'b','Color', [ 0.2 0.2 .2])
axis tight
  %print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LFP_PLV_lock' num2str(stim_freq) '.pdf'])
