close all
clear all
 fileDir='Z:\eng_research_handata2\Sudi_Sridhar\602101\08252020';
 
%basic template of how you extract plx file
        cd(fileDir);
        

       file= uigetfile('*.plx'); 
        
[~,frame,~]=plx_event_ts(file, 1);
[~,sti,~]=plx_event_ts(file, 3);
[~,motion,~]=plx_event_ts(file,4);








