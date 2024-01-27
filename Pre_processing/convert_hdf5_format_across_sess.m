clear all;
close all;

datapath='/home/hanlabadmins/eng_handata/eng_research_handata2/Sudi_Sridhar/';
stim_freq=2;
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

% load binarized. mat
load('Binarized_event_detection_final.mat');
% load plex.mat
% [A,b]=uigetfile('*plex.mat','multiselect','off');
% load(A);
matFiles = dir('*.mat') ; 
matFiles_name = {matFiles.name} ; 
idx = find(~cellfun(@isempty,strfind(matFiles_name,'plex'))) ;
load(matFiles_name{idx})

% Make hdf5 file with all the necessary field from above 
% trace,onset_binary_trace,offset_binary_trace,binary_trace,stim onset,stim
delete processed_trace_final.h5
filename ='processed_trace_final.h5';
%% Add traces and binary traces

inputStruct = event_out;
groups = fieldnames(inputStruct); %Base groups as names of fields in structure
for field = 1:numel(groups) 
    fieldSize = numel([inputStruct.(groups{field})])/numel(inputStruct); %Determine size of field vectors
    dataset = nan(fieldSize, numel(inputStruct)); %Initialize Dataset
    for idx = 1:numel(inputStruct) %Populate Dataset
        dataset(:,idx) = inputStruct(idx).(groups{field});
    end
    h5create(filename, ['/' groups{field}], size(dataset));
    h5write(filename, ['/' groups{field}], dataset);
    descString = 'Post event detection';
            
    h5writeatt(filename, ['/' groups{field}], 'Description', descString);
end



%% Add onset, offset of stim and timestamp of imaging 
dataset = plx.Stim_onset; 
h5create(filename, ['/Stim_onset' ], size(dataset));
h5write(filename, ['/Stim_onset'], dataset);

dataset = plx.Stim_offset; 
h5create(filename, ['/Stim_offset' ], size(dataset));
h5write(filename, ['/Stim_offset'], dataset);
end