clear all
close all
%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
addpath(genpath('U:\eng_research_handata\EricLowet'))
mainroot='U:\eng_research_handata\EricLowet\SUDI';
datapath='U:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\';
stim_freq=1;% 
cha=1; % 1=LFP. 2=ball tracking
aCOHs=[];aSPs=[];aCOHs2=[];aPHs=[];aNT=[];;MOT_respALL=[];Stim_respALL=[];
if stim_freq==1
% 10Hz
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


%%%%%%%%%%%%%%%%%%%%%%%%%
if stim_freq==1
path_resp='C:\Users\sudiksha\Documents\Codes\'
resp_cells=h5read([ path_resp 'Motion_resp_cells_sustained_10Hz_final.h5'],allRES{ses});
resp_cells=squeeze(resp_cells);
 resp_cells=resp_cells+1;
resp_cellsST=h5read([ path_resp 'Stim_resp_cells_sustained_10Hz_final.h5'],[allRES{ses} '_dec_stn']);
resp_cellsST=squeeze(resp_cellsST);
 resp_cellsST=resp_cellsST+1;
else
  path_resp='C:\Users\sudiksha\Documents\Codes\'
resp_cells=h5read([ path_resp 'Motion_resp_cells_sustained_145Hz_final.h5'],allRES{ses});
resp_cells=squeeze(resp_cells);
 resp_cells=resp_cells+1;
 % inc_mov, dec_mov, inc_stn, dec_stn
resp_cellsST=h5read([ path_resp 'Stim_resp_cells_sustained_145Hz_final.h5'],[allRES{ses} '_dec_stn']);
resp_cellsST=squeeze(resp_cellsST);
 resp_cellsST=resp_cellsST+1;  end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Sampling_freq=20;

% Mean across all neuron
mean_traces=zscore(nanmean(traces,2));  %% Here is use median instead of mean as a better population average
% High-speed periods
moving_period=h5read('processed_motion.h5','/moving_period');
moving_period=moving_period(idx);  % moving period in 0 and 1

%2099.95 should be used for some sessions 
time_vect1  = 0:1/Sampling_freq:2099.9;% Original time vector (20Hz)
signal1=v ; % signal
time_vect2 = 0:1/1000:2099.9; % time vector for interpolation (1000Hz)
%signal1_Intp=interp1(time_vect1,signal1,time_vect2);%  interpolated signal


LFP=LFP_data.LFP;
asel=fastsmooth(zscore(fastsmooth(abs(hilbert(LFP)),5,1,1))>4,300,1,1);
LFP(asel>0)=median(LFP)+randn(1, length(find(asel>0))).*std(LFP);
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
aligned_LFP=aligned_LFP(1:length(moving_period));


%%%%%%%%%%%%%%%%%%%%%%%
if 1
 
    lfp_all=[];FS=20
    lfp_all.trial{1}(1,:)= zscore(aligned_LFP);%
        lfp_all.trial{1}(2,:)= zscore(v-fastsmooth(v,400,1,1));
    lfp_all.time{1}= (1:size(lfp_all.trial{1},2))./FS;
    lfp_all.label= {'LFP', 'motion'};
      
%%%%%%%%%%%%%%%
deltP=lfp_all.trial{1}(2,:);
Fn = FS/2;FB=[ 3 4 ];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
deltP= abs(hilbert(filtfilt(B,A,    deltP)));
thresD= prctile(deltP(moving_period==1),30);   %Select only period with motion and at least >30% delta power
%                                         
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  cfg = []; %block_type == cfg.blk
% cfg.method ='wavelet';%mtmfft'; %'mvar';
% cfg.output ='fourier';
% cfg.taper='hanning';
% cfg.keeptapers ='yes';
% cfg.keeptrials ='yes';
% cfg.trials='all';cfg.tapsmofrq =5;%
% cfg.channel= ['all']; %chans=cfg.channel;
% cfg.foi= [1.5:0.5:8];    cfg.t_ftimwin =[ones(1,length(cfg.foi))*1];
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
MOT_resp=[];Stim_resp=[];
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


  if sum(neuron==resp_cells)>0
MOT_resp=[MOT_resp, 1];else; MOT_resp=[MOT_resp, 0];  end
if sum(neuron==resp_cellsST)
Stim_resp=[Stim_resp, 1];else; Stim_resp=[Stim_resp, 0];  end

    end
end
Dx=max(diff(allSP));
end

calc_transient_perc=0; % criterion I am no
aCOHs= [aCOHs, allCOHS(:,Dx>prctile(Dx,calc_transient_perc))];
aCOHs2= [aCOHs2, allCOHS2(:,Dx>prctile(Dx,calc_transient_perc))];
aPHs= [aPHs, allPH1(:,Dx>prctile(Dx,calc_transient_perc))];
aSPs= [aSPs, allSP(:,Dx>prctile(Dx,calc_transient_perc))];
aNT=[aNT, allNT(Dx>prctile(Dx,calc_transient_perc))];
MOT_respALL=[MOT_respALL, MOT_resp(Dx>prctile(Dx,calc_transient_perc))];
Stim_respALL=[Stim_respALL, Stim_resp(Dx>prctile(Dx,calc_transient_perc))];

end
 
savepath='C:\Users\sudiksha\Documents\Codes\Spectral_analysis\Figures\'

   figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
 fill_error_area2(freq2.freq,nanmean(aCOHs(:,MOT_respALL==1),2), nanstd(aCOHs(:,MOT_respALL==1),[],2)./sqrt(size(aCOHs(:,MOT_respALL==1),2)),[ 0.5 0.5 0.5])
 % fill_error_area2(freq2.freq,nanmean(aCOHs2(:,MOT_respALL==1),2), nanstd(aCOHs2(:,MOT_respALL==1),[],2)./sqrt(size(aCOHs2(:,MOT_respALL==1),2)),[ 0.5 0.5 0.5])
    fill_error_area2(freq2.freq,nanmean(aCOHs(:,Stim_respALL==1),2), nanstd(aCOHs(:,Stim_respALL==1),[],2)./sqrt(size(aCOHs(:,Stim_respALL==1),2)),[ 0.5 0.5 0.5])
 %fill_error_area2(freq2.freq,nanmean(aCOHs(:,:),2), nanstd(aCOHs(:,:),[],2)./sqrt(size(aCOHs(:,:),2)),[ 0.5 0.5 0.5])
 %  fill_error_area2(freq2.freq,nanmean(aCOHs2(:,Stim_respALL==1),2), nanstd(aCOHs2(:,Stim_respALL==1),[],2)./sqrt(size(aCOHs2(:,Stim_respALL==1),2)),[ 0.5 0.5 0.5])
% fill_error_area2(freq2.freq,nanmean(aCOHs2(:,:),2), nanstd(aCOHs2(:,:),[],2)./sqrt(size(aCOHs2(:,:),2)),[ 0.5 0.5 0.5])
     plot(freq2.freq,nanmean(aCOHs(:,MOT_respALL==1),2),'r','Linewidth',1.5)
    %     plot(freq2.freq,nanmean(aCOHs(:,:),2),'b-','Linewidth',1.5)
             plot(freq2.freq,nanmean(aCOHs(:,Stim_respALL==1),2),'b','Linewidth',1.5)
%plot(freq2.freq,nanmean(aCOHs2(:,MOT_respALL==1),2),'k','Color', [ 0.2 0.2 .2])
%plot(freq2.freq,nanmean(aCOHs2(:,Stim_respALL==1),2),'k','Color', [ 0.2 0.2 .2])
%plot(freq2.freq,nanmean(aCOHs2(:,:),2),'k','Color', [ 0.2 0.2 .2])

axis tight
title('Mov resp (red)  vs Stim resp')
 % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_lock' num2str(stim_freq) '.pdf'])


   figure('COlor','w','Position', [ 300 400 250 180],'Renderer', 'painters') 
 fill_error_area2(freq2.freq,nanmean(aCOHs(:,MOT_respALL==1),2), nanstd(aCOHs(:,MOT_respALL==1),[],2)./sqrt(size(aCOHs(:,MOT_respALL==1),2)),[ 0.5 0.5 0.5])
  %fill_error_area2(freq2.freq,nanmean(aCOHs2(:,MOT_respALL==1),2), nanstd(aCOHs2(:,MOT_respALL==1),[],2)./sqrt(size(aCOHs2(:,MOT_respALL==1),2)),[ 0.5 0.5 0.5])
 fill_error_area2(freq2.freq,nanmean(aCOHs(:,MOT_respALL==0),2), nanstd(aCOHs(:,MOT_respALL==0),[],2)./sqrt(size(aCOHs(:,MOT_respALL==0),2)),[ 0.5 0.5 0.5])
 % fill_error_area2(freq2.freq,nanmean(aCOHs2(:,MOT_respALL==0),2), nanstd(aCOHs2(:,MOT_respALL==0),[],2)./sqrt(size(aCOHs2(:,MOT_respALL==0),2)),[ 0.5 0.5 0.5])
  plot(freq2.freq,nanmean(aCOHs(:,MOT_respALL==1),2),'r','Linewidth',1.5)
         plot(freq2.freq,nanmean(aCOHs(:,MOT_respALL==0),2),'b-','Linewidth',1.5)
%plot(freq2.freq,nanmean(aCOHs2(:,MOT_respALL==1),2),'k','Color', [ 0.2 0.2 .2])
%plot(freq2.freq,nanmean(aCOHs2(:,MOT_respALL==0),2),'k','Color', [ 0.2 0.2 .2])
axis tight
title('Mov resp (red) vs Mov non-resp')
frsel=find(freq2.freq>3 & freq2.freq<4) 
delt_PLV_mov_resp= nanmean(aCOHs(frsel,MOT_respALL==1));
delt_PLV_non_mov_resp= nanmean(aCOHs(frsel,MOT_respALL==0));
ranksum(delt_PLV_mov_resp,delt_PLV_non_mov_resp)
% %figure,rose(delt_angs)
% 
% [h,p,ci,stats] =ttest(delt_PLV,delt_PLV_SH)
% 
%    figure('COlor','w','Position', [ 300 400 120 90],'Renderer', 'painters') 
% polarhistogram(delt_angs(delt_PLV>0.0),[-pi:0.33:pi],'FaceColor','red','Normalization','probability','Facealpha',0.3);hold on,
%  % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_PHASE' num2str(stim_freq) '.pdf'])
% %   figure,plot(-wind:wind,nanmean(aSPs,2))
% %  fill_error_area2(-wind:wind,nanmean(aSPs,2), nanstd(aSPs,[],2)./sqrt(size(aCOHs,2)),[ 0.5 0.5 0.5])
% 
% % nt=90;
% % figure, imagesc(aCOHs(:,aNT>nt)')
% % colormap(jet)
% % set(gca,'Clim', [-0.05 0.12])
% 
% % figure,imagesc(aSPs')
% % 
% % 
% 
