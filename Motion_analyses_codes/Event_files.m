close all
clear all
[aA,aa]= uigetfile('*','multiselect','on');  %608452 01062020
for i=1:numel(aA)/3
    i
    b=aA(3*i-2:3*i);
        kk=strfind(b{1},'_');
        name=b{1}(1:kk(3));
    plex=contains(b,'plex');
    motion=contains(b,'motion');
    trace=contains(b,'trace');  %load all file, assign them in three different types
    cd(aa)
    load(b{trace})
    load(b{plex})
    [raw_motion,v]=Clear_Velocity(b{motion},20); %get motion trace
    n=numel(v);  %length of motion trace
    dir=r_out(1).file(1).filename; %get handata2/sudi directory
    N=numel(r_out(1).trace); %length of calcium trace
    NROI=numel(r_out);
    dir=strrep(dir,'Z:','/home/hanlabadmins/eng_handata'); %sometime I couldn't access Z: at calcium computer
    dir=strrep(dir,'\','/');
    [lt,ut]=var_low_high_speed(v,20,0.05,5,2,2);  %get the low and high speed threshold
    MT=plx.Timestamp_Motion(1:n);  %get motion timestamp
    old_IT=plx.Timestamp_Imaging(1:N);
    %sti=plx.Stim_onset;
    %bas=[1,plx.Stim_offset+1200];
    %sti(2,:)=plx.Stim_offset;
    %bas(2,:)=[plx.Stim_onset,n];
    IT=0:0.05:2160-0.05; %ideal imaging timestamp-
    M=motion_onset(v,lt,ut,1);
    LM=motion_onset(v,lt,ut,0);
    M=Timeframe_Conversion(MT,M,IT);
    LM=Timeframe_Conversion(MT,LM,IT);

    %
    [stationary,moving,AP,DP]=motion_event(v,lt,ut);  %get whatever onset/offset indices frame output by the function
   % M=Timeframe_Conversion(MT,M,IT); %convert the indices of MT to that of IT
    %A=Timeframe_Conversion(MT,acceleration,IT);
    S=Timeframe_Conversion(MT,stationary,IT);
    AP=Timeframe_Conversion(MT,AP,IT);
    DP=Timeframe_Conversion(MT,DP,IT);
    moving=Timeframe_Conversion(MT,moving,IT);
    %m=Timeframe_Conversion(MT,moving,IT);
    %}
    V=interp1(MT,v,IT); %do linear interpolation that map motion trace into imaging timestamp
    raw_speed_trace=interp1(MT,raw_motion,IT); 
    %motion_event_data=struct('trace',v','event_idx',moving,'event_time',MT(moving));
    %generate .mat file in motion timestamp
    ti=isnan(V);
    motion_onsets=zeros(1,N);
    Motion_event=zeros(1,N);
    motion_onsets(ti)=nan;
    Motion_event(ti)=nan;
    for j=1:numel(moving)/2
        Motion_event(moving(j,1):moving(j,2))=1;
    end
    for j=1:numel(M)
                motion_onsets(M(j,1))=1;
    end
    broad_motion_onset=zeros(1,N);
    broad_motion_onset(ti)=nan;
    for j=1:numel(LM)
        
        broad_motion_onset(LM(j,1))=1;
    end
    
 
    
    
    stationary_period=zeros(1,N);
    stationary_period(ti)=nan;
    for j=1:numel(S)/2
        stationary_period(S(j,1):S(j,2))=1;
    end
    
   
    
    
    acceleration_period=zeros(1,N);
    acceleration_period(ti)=nan;
    for j=1:numel(AP)/2
        acceleration_period(AP(j,1):AP(j,2))=1;
    end
    deceleration_period=zeros(1,N);
    deceleration_period(ti)=nan;
    for j=1:numel(DP)/2
        deceleration_period(DP(j,1):DP(j,2))=1;
    end
    
    %{
    stimulation_period=zeros(1,N);
    for j=1:numel(sti)/2
        stimulation_period(sti(1,j):sti(2,j))=1;
    end
    baseline_period=zeros(1,N);
    for j=1:numel(bas)/2
        baseline_period(bas(1,j):bas(2,j))=1;
    end
    
    T=zeros(NROI,N);
    for j=1:NROI
   b1=r_out(j).BGtrace;
   a=r_out(j).trace;
   % subtract and detrend 
   a1=a-b1;
   a1=interp1(old_IT,a1,IT);
   
   a1=(a1-mean(a1))/(max(a1)-mean(a1))*100;
   a1=detrend(a1,'linear',1:1000:N);
   %b1=(b1-mean(b1))/(max(b1)-mean(b1))*100;  
   %b1=detrend(b1,'linear',1:1000:N);
   T(j,:)=a1;
   %B(i,:)=b1(con);
    end
    
    cd(dir)
    T_onset=h5read('processed_trace.h5','/onset_binary_trace')';
    T_offset=h5read('processed_trace.h5','/offset_binary_trace')';
    T_trace=h5read('processed_trace.h5','/binary_trace')';
    integration=struct('timestamp',IT,'speed_trace',V,'motion_onset',motion_onset,'motion_event',Motion_event,'m',M,'broad_motion_onset',broad_motion_onset,'broad_motion_event',broad_Motion_event,'lm',LM,'stationary_period',stationary_period,'s',S,'acceleration_period',acceleration_period,'deceleration_period',deceleration_period,'stimulation_period',stimulation_period,'sti',sti,'baseline_period',baseline_period,'calcium_traces',T,'binarized_calcium_traces',T_trace,'calcium_onset',T_onset,'calcium_offset',T_offset);
    cd('Z:\eng_research_handata\Undergrad\Chengqian Zhou\integration')
    save([name,'integration.mat'],'integration')
    %}
    %{
    acceleration_period=zeros(1,N);
    acceleration_period(ti)=nan;
    for j=1:numel(A)/2
        acceleration_period(A(j,1):A(j,2))=true;
    end
    stationary_period=zeros(1,N);
    stationary_period(ti)=nan;
    for j=1:numel(S)/2
        stationary_period(S(j,1):S(j,2))=true;
    end
    moving_period=zeros(1,N);
    moving_period(ti)=nan;
    for j=1:numel(m)/2
        moving_period(m(j,1):m(j,2))=true;
    end
    %}
    %motion_event_out=struct('trace',V,'motion_onset',motion_onset,'acceleration_period',acceleration_period,'stationary_period',stationary_period,'moving_period',moving_period,'low_speed_threshold',lt,'high_speed_threshold',ut);
    %create .h5 file from matlab struct
    new_onset=struct('timestamp',IT,'speed_trace',V,'raw_motion',raw_motion,'raw_speed_trace',raw_speed_trace,'motion_onset_with',motion_onsets,'motion_onset_without',broad_motion_onset,'moving_period',Motion_event,'aceleration_period',acceleration_period,'deceleration_period',deceleration_period,'stationary_period',stationary_period,'low_speed_threshold',lt,'high_speed_threshold',ut,'stationary_starts',S(:,1),'stationary_ends',S(:,2),'moving_starts',moving(:,1),'moving_ends',moving(:,2));
    
    
    cd(dir)
    
    %save('motion_event_detection_output.mat','motion_event_data')
    %save('Binarized_motion_event_detection.mat','new_onset')
    
     
    filename ='processed_motion.h5';
    delete processed_motion.h5
    inputStruct = new_onset;
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



    
  
end
    
        
 
    
    
    
    