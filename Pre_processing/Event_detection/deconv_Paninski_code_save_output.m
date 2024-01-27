close all
clear all
addpath(genpath('/home/hanlabadmins/eng_handata/EricLowet/'))
%%  Load motion and mean fluorescence trace
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/612535/07082020/612535_10Hz_AudioVisual/motion_corrected')
mainroot='/home/hanlabadmins/eng_handata/EricLowet/SUDI/';
datapath='/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/';
stim_freq=2;
aSPs=[];aDx=[];
if stim_freq==1
% 10Hz
allSES{1}='/608448/01062020/608448';  % no LFP
allSES{2}='/608451/01062020/608451';
allSES{3}='/608452/01062020/608452';
allSES{4}='/608450/03092020/608450';
allSES{5}='/602101/03092020/602101';
allSES{6}='/611311/03092020/611311';
allSES{7}='/608448/01172020/608448';
allSES{8}='/608451/07082020/608451';
allSES{9}='/608452/07082020/608452';
allSES{10}='/611111/07082020/611111';
allSES{11}='/602101/07082020/602101';
allSES{12}='/612535/07082020/612535';
allSES{13}='/615883/02052021/10Hz/615883';
allSES{14}='/615883/03122021/10Hz/615883';
else
% 140Hz
allSES{1}='/608448/01102020/608448';
allSES{2}='/608451/01102020/608451';
allSES{3}='/608452/01102020/608452';
allSES{4}='/602101/03112020/602101';
allSES{5}='/611311/03112020/611311';
allSES{6}='/611111/03112020/611111';
allSES{7}='/608450/01172020/608450';
allSES{8}='/608451/07012020/608451';
allSES{9}='/608452/07012020/608452';
allSES{10}='/611111/07012020/611111';
allSES{11}='/612535/07022020/612535';
allSES{12}='/602101/07022020/602101'; %exclude LFP
allSES{13}='/615883/02052021/145Hz/615883';
allSES{14}='/615883/03122021/145Hz/615883';
end

for ses=1:length(allSES)
if  stim_freq==1
cd([ datapath  allSES{ses} '_10Hz_AudioVisual/motion_corrected']);
else
 cd([ datapath  allSES{ses} '_145Hz_AudioVisual/motion_corrected'])   
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
%traces=traces(idx,:);tracesF=tracesF(idx,:);
traces=traces;tracesF=tracesF;
Sampling_freq=20;
% Mean across all neuron
mean_traces=zscore(nanmean(traces,2));  %% Here is use median instead of mean as a better population average
% High-speed periods
moving_period=h5read('processed_motion.h5','/moving_period');
moving_period=moving_period(idx);  % moving period in 0 and 1

g
%%%%%%%%%%%%%%%%%%%%%%%
if 1
 


mm=0;clear  allSP 
clear neuron 
event_out=struct();
rising_phase=0;
for neuron=1:size(tracesF,2)
    
    
  %%%%%%%%%%%%%%%%%%%%%% 
    %%%% DECONVULTION CODE %%%%%%%%%%%%%%
 [ s_oasis, onsets2 ,onsets]= deconvSudi(tracesF(:,neuron));
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

s= size(tracesF,1);
event_out(neuron).binary_trace =  zeros(1,s);
event_out(neuron).onset_binary_trace =  zeros(1,s);
event_out(neuron).trace =  zeros(1,s);
wind=20;
peaks=zeros(length(onsets2),1);
for id=1:length(onsets2)
    if onsets2(id)+wind <size(traces,1) 
        [max_amp,peak_idx]=max(tracesF(onsets2(id): onsets2(id)+wind,neuron));  
    end
    rising_phase=[rising_phase,peak_idx-1];
    peaks(id)=onsets2(id)+peak_idx-1;
    %% Store binarized traces 
    if peaks(id)<43200
        event_out(neuron).binary_trace([onsets2(id):peaks(id)])=1;
    end
    event_out(neuron).onset_binary_trace(onsets2(id))=1;
end      
%% Store onsets and peaks
%event_out(neuron).onsets=onsets2;
%event_out(neuron).peaks=peaks;
event_out(neuron).trace=tracesF(:,neuron)';
minimum_spikes=20;
    if length(onsets2) >minimum_spikes
        mm=mm+1

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
Dx=max(diff(allSP));
end

calc_transient_perc=0;
aDx=[aDx, Dx];
aSPs= [aSPs, allSP(:,Dx>prctile(Dx,calc_transient_perc))];
%% Save output 
%save('Binarized_event_detection_final.mat','event_out','-v7.3')
end
rising_phase= rising_phase(2:end);
save('rising_phase_145Hz.mat','rising_phase','-v7.3') 

tim_ax=(-wind:wind)./20;
figure,imagesc(tim_ax,[],aSPs')
colormap(jet)
title('Onset trig averaged  Calc signals')
ylabel('Neurons'); xlabel('Time sec')

figure,imagesc(tim_ax,[],dataT2')
colormap(jet)
title('Onset trig Calc signals SINGLE NEURON')
ylabel('Windows'); xlabel('Time sec')


figure
plot(tracesF(:,neuron))
hold on,%plot(onsets2,ones(1,length(onsets2)),'r.')
plot(onsets2,tracesF(onsets2,neuron),'r.')
plot(peaks,tracesF(peaks,neuron),'b.')

figure
%neuron =5;
plot(tracesF(:,neuron))
hold on,%plot(onsets2,ones(1,length(onsets2)),'r.')
plot(onsets2,tracesF(onsets2,neuron),'r.')
plot(peaks,tracesF(peaks,neuron),'b.')
% 
figure
trace=tracesF(:,neuron);
x_axis = [1:numel(trace)]/20;
plot(x_axis,tracesF(:,neuron))
hold on,%plot(onsets2,ones(1,length(onsets2)),'r.')
%onsets2=event_out(neuron).onsets;
%peaks=event_out(neuron).peaks;
for n=1: length(onsets2)
plot(x_axis(onsets2(n):peaks(n)),trace(onsets2(n):peaks(n)),'r')
end





