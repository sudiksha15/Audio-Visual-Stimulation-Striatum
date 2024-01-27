close all
clear all
fileDir='Z:\eng_research_handata\eng_research_handata2\Sudi_Sridhar\Motion_on_ball\145 Hz_AV_Stimulation';
cd(fileDir)
A=uigetfile('*.plx','multiselect','on');
%name={'602101_06192020_10Hz','608448_06182020_10Hz','608449_06192020_10Hz','608450_06182020_10Hz','608451_06192020_10Hz','608452_06182020_10Hz','611111_06182020_10Hz','611311_06182020_10Hz','612533_06192020_10Hz','612535_06182020_10Hz'};
name={'608448_07022020_145Hz','611311_07022020_145Hz','612533_07022020_145Hz'};

for i=1:numel(A)
    [~,sti,~]=plx_event_ts(A{i}, 10);
    [~,motion,~]=plx_event_ts(A{i}, 8);
    plx =struct('Timestamp_stim',sti,'Timestamp_Motion',motion);
    cd('C:\Users\hanlabadmins\Desktop\Motion analysis\movement\motion on ball')
    save([name{i},'_plex.mat'],'plx')
    cd(fileDir)
end
A=uigetfile('*.txt','multiselect','on');
for i=1:numel(A)
    copyfile(A{i},['C:\Users\hanlabadmins\Desktop\Motion analysis\movement\motion on ball\',name{i},'_motion.txt'])
end
    