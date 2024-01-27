clear all
close all
%function extract_sti
fileDir='Z:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\608452\11182019';
cd(fileDir);
file = uigetfile('*.plx');%,'MultiSelect','on');

%loading and aligning time frame for imaging, motion and stimulation
[~,frame,~]=plx_event_ts(file, 1);
[~,sti,~]=plx_event_ts(file, 4);
[~,motion,~]=plx_event_ts(file,8);
motion=motion-frame(1);
sti=sti-frame(1);
frame=frame-frame(1);


%lfp not accountable
[~,names] = plx_adchan_names(file);
[~,adchans] = plx_ad_chanmap(file);
lfp_channel = strmatch('FP01', names);
current_channel_number = adchans(lfp_channel);
lfp_channel_name = strcat(names(lfp_channel,:));
[~, ~, ~, ~, lfp] = plx_ad_v(file, current_channel_number);

%set up beginning and end of each stimulation
sn=numel(sti);
plexon_stim_ons=zeros(1,sn);
plexon_stim_offs=zeros(1,sn);
for i=1:sn
     plexon_stim_ons(i)=round(sti(i)*20);
     plexon_stim_offs(i)=round((sti(i)+10)*20);
end

plx =struct('LFP',lfp,'Timestamp_Imaging',frame,'Timestamp_stim',sti,'Timestamp_Motion',motion,'Stim_onset',plexon_stim_ons,'Stim_offset',plexon_stim_offs);
cd('C:\Users\hanlabadmins\Desktop\Motion analysis\movement\mc stimulation')
save('608452_11182019_10Hz_25uA_plex.mat','plx')
