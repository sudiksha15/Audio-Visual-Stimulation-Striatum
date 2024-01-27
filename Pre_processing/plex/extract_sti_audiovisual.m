clear all
close all
%function extract_sti
%fileDir='Z:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\615883\02052021\145Hz';
%cd(fileDir);
file = uigetfile('*.plx');%,'MultiSelect','on');

%loading and aligning time frame for imaging, motion and stimulation
[~,frame,~]=plx_event_ts(file, 2);
[~,sti,~]=plx_event_ts(file, 7);
[~,motion,~]=plx_event_ts(file,4);
motion=motion-frame(1);
sti=sti-frame(1);
frame=frame-frame(1);


%loading lfp
[~,names] = plx_adchan_names(file);
[~,adchans] = plx_ad_chanmap(file);
lfp_channel = strmatch('FP01', names);
current_channel_number = adchans(lfp_channel);
lfp_channel_name = strcat(names(lfp_channel,:));
[~, ~, ~, ~, lfp] = plx_ad_v(file, current_channel_number);


%detect beginning and end of stimulation
a=diff(sti);
Sti =sti*20;
[~,a2]=sort(a,'descend');
b=a2(1:4);
a=sort(b);
a=[0;a;numel(Sti)];
plexon_stim_ons=zeros(1,5);
plexon_stim_offs=zeros(1,5);
for i=1:5
    plexon_stim_ons(i)=round(Sti(a(i)+1));
end
for i=1:5
    plexon_stim_offs(i)=round(Sti(a(i+1)));
end

plx =struct('LFP',lfp,'Timestamp_Imaging',frame,'Timestamp_stim',sti,'Timestamp_Motion',motion,'Stim_onset',plexon_stim_ons,'Stim_offset',plexon_stim_offs);
    
%save('Z:\Chengqian Zhou\motion data\608448_01102020_145Hz_plex.mat','plx')