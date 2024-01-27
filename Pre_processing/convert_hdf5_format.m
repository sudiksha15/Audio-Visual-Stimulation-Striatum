clear all;
close all;

% Make a function for below input would be mat files, maybe filename? and output is the h5
% file 

% load binarized. mat
load('Binarized_event_detection3.mat');
% load plex.mat
% [A,b]=uigetfile('*plex.mat','multiselect','off');
% load(A);

matFiles = dir('*.mat') ; 
matFiles_name = {matFiles.name} ; 
idx = find(~cellfun(@isempty,strfind(matFiles_name,'plex'))) ;
load(matFiles_name{idx})

% Make hdf5 file with all the necessary field from above 
% trace,onset_binary_trace,offset_binary_trace,binary_trace,stim onset,stim
%delete processed_trace.h5
filename ='processed_trace3.h5';
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

% dataset = plx.Timestamp_Imaging; 
% h5create(filename, ['/Timestamp_imaging' ], size(dataset));
% h5write(filename, ['/Timestamp_imaging'], dataset);

% Add BG traces?


% Also want LFP
% timestamp motion combine with motion stuff or here?

% load in jupyter to check 