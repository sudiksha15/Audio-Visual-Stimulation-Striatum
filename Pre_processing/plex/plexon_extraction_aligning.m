%%% Wriiten by Sudi
%%% Function to extract all Plexon triggers including timestamps and LFP
%%% Input includes filename along with path and imaging_true=1 if imaging
%%% present in the session else 0 
%%% Output is a struct that has LFP ,timestamps, stim onset and stim offset as fields
%%% It can be saved into mat or hdf5 files  

function plexon_triggers= plexon_extraction_aligning(plx_filename,imaging_true)
%% LFP extraction 
[~,names] = plx_adchan_names(plx_filename);
[~,adchans] = plx_ad_chanmap(plx_filename);
lfp_channel = strmatch('FP01', names);
current_channel_number = adchans(lfp_channel);
lfp_channel_name = strcat(names(lfp_channel,:));

fprintf(['Loading ',lfp_channel_name,'\n'])
[adfreq, ~, ~, ~, lfp] = plx_ad_v(plx_filename, current_channel_number);

%% Timestamps for imaging,stimulation and motion 

fprintf(['Loading Timestamps','\n'])
% Channel numbers - can be input arguments but since fixed for all my plx
% files I am specifying here
channel_imaging=1;
channel_stim=10;
channel_motion=8;

% Timestamp for imaging 
if imaging_true==1
    [~, ts, ~]=  plx_event_ts(plx_filename,channel_imaging);
end
% Timestamp for stimulation 
[~, ts_stim, ~]=  plx_event_ts(plx_filename,channel_stim);
% Timestamp for motion
[~, ts_motion, ~]=  plx_event_ts(plx_filename,channel_motion);
%% Aligning the timestamps

% Align with respect to imaging if imaging 
if imaging_true==1
    % Align with respect to imaging 
    ts_shifted = ts-ts(1);
    ts_motion_shifted = ts_motion-ts(1);
    ts_stim_shifted = ts_stim-ts(1);
else
    % Align with respect to motion if no imaging
    ts_motion_shifted = ts_motion-ts_motion(1);
    ts_stim_shifted = ts_stim-ts_motion(1);
end

%% Find stimulation onset and offset frames 
% Stimulation onset 

stim_ons=7200:7200:36000;
frame_stim_shifted = ts_stim_shifted*20;
plexon_stim_ons = zeros(1,numel(stim_ons));
for i = 1:numel(stim_ons)
    [val,idx]=min(abs(frame_stim_shifted-stim_ons(i)));
    plexon_stim_ons(i)=round(frame_stim_shifted(idx));
end 

% Stimulation offset
plexon_stim_offs = zeros(1,numel(stim_ons));
for i = 1:4
    [val,idx]=min(abs(frame_stim_shifted-stim_ons(i+1)));
    plexon_stim_offs(i)=round(frame_stim_shifted(idx-1));
end 
plexon_stim_offs(5)= round(frame_stim_shifted(end));


%{
a=diff(ts_stim_shifted);
frame_stim_shifted = ts_stim_shifted*20;
[~,a2]=sort(a,'descend');
b=a2(1:4);
a=sort(b);
a=[0;a;numel(ts_stim_shifted)];
plexon_stim_ons=zeros(1,5);
 plexon_stim_offs=zeros(1,5);
for i=1:5
    plexon_stim_ons(i)=round(frame_stim_shifted(a(i)+1));
end
for i=1:5
    plexon_stim_offs(i)=round(frame_stim_shifted(a(i+1)));
end
    
%}



%% Make a struct to save - neat way  
% See what all you want in output - only including shifted timestamps for now - can include original as well 
if imaging_true==1
    plexon_triggers =struct('LFP',lfp,'Timestamp_Imaging',ts_shifted,'Timestamp_stim',ts_stim_shifted,'Timestamp_Motion',ts_motion_shifted,'Stim_onset',plexon_stim_ons,'Stim_offset',plexon_stim_offs);
else
    plexon_triggers =struct('LFP',lfp,'Timestamp_stim',ts_stim_shifted,'Timestamp_Motion',ts_motion_shifted,'Stim_onset',plexon_stim_ons,'Stim_offset',plexon_stim_offs);
end
end