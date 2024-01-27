%session select
clear all;  
close all;
[A,b]= uigetfile('*','MultiSelect','on'); %load all 2 pairs data
a1=contains(A,'145Hz');
a2=contains(A,'10Hz');
%a3=contains(A,'40Hz');
A1=A(a1);
A2=A(a2);
%A3=A(a3);
A={A1,A2};
cd(b);
%{
M{1}=[];
M{2}=[];

Mpair{1}=[];
Mpair{2}=[];
Spair{1}=[];
Spair{2}=[];
Vpair{1}=[];
Vpair{2}=[];
Name{1}={};
Name{2}={};
%}

Ppair{1}=[];
Ppair{2}=[];
%Ppair{3}=[];
x=-1200:4800;
T{1}=zeros(numel(A1)/2,numel(x));
T{2}=zeros(numel(A2)/2,numel(x));

    
for k=1:numel(Ppair)
    for j=1:numel(A{k})/2
        a=A{k}(2*j-1:2*j);
        kk=strfind(a{1},'_');
        name=a{1}(1:kk(3));
        name=strrep(name,'_',' ');
        a1=contains(a,'txt');
        CV =Clear_Velocity(a{a1},20);
        [lt,ut]=var_low_high_speed(CV,20,0.1,2.5,2,2);
        %[~,~,m,~]=motion_event(CV,lt,ut);
        m=motion_onset(CV,lt,ut,0);
        %m=New_motion_event(CV,lt,ut,1);
        % Add lines for onset
        % 
      
        n=numel(CV);
        bm=zeros(1,n);
        if ~isempty(m)
        for i=1:numel(m)/numel(m(1,:))
            
            %bm(m(i,1):m(i,2))=1;
            bm(m(i,1))=1;
        end
        end

        a2=contains(a,'plex');
        load(a{a2});
        MT=plx.Timestamp_Motion;  
        ST=plx.Timestamp_stim;
        sti =sti_extraction(ST);
        sn=numel(sti)/2;
        idx=find_idx(MT,sti);
        bsti=false(n,1);
        bbas=true(n,1);
        for i=1:numel(idx)/numel(idx(:,1))
            bsti(idx(1,i):idx(2,i))=true;
            bbas(idx(1,i):idx(2,i))=false;
            bbas(idx(2,i)+1:idx(2,i)+1200)=false;
        end
        bstim=bm(bsti);

        bbasm=bm(bbas);

 
      Ppair{k}=[Ppair{k};sum(bstim)/sum(bsti)*1200,sum(bbasm)/sum(bbas)*1200];
      
      for i=1:sn
          T{k}(j,:)=T{k}(j,:)+bm(idx(1,i)-1200:idx(1,i)+4800)/sn;
      end
      
      
  
        end
end




%

P145=Ppair{1};
P10=Ppair{2};
%P40=Ppair{3};
%figure(2)
%bi_box_graph (V145,V10,'normalized speed')

figure(1)
bi_box_graph (P145,P10,'% of motion onsets')
%{
figure(2)
subplot(2,2,1)
imagesc(x,1:numel(A1)/2,T{1})
caxis([0,1])
hold on
plot([0,0],[1,numel(A1)/2],'r');
plot([1200,1200],[1,numel(A1)/2],'r');
hold off
title('145Hz binarized speed traces over trials')
ylabel('session index');
subplot(2,2,3)
plot(x,mean(T{1}))
hold on
plot([0,0],[min(mean(T{1})),max(mean(T{1}))],'r');
plot([1200,1200],[min(mean(T{1})),max(mean(T{1}))],'r');
hold off
title('average over sessions')
xlabel('frames after stimulation onsets')
ylabel('p\in motion event')
subplot(2,2,2)
imagesc(x,1:numel(A2)/2,T{2})
caxis([0,1])
hold on
plot([0,0],[1,numel(A2)/2],'r');
plot([1200,1200],[1,numel(A2)/2],'r');
hold off
title('10Hz')
subplot(2,2,4)
plot(x,mean(T{2}))
hold on
plot([0,0],[min(mean(T{2})),max(mean(T{2}))],'r');
plot([1200,1200],[min(mean(T{2})),max(mean(T{2}))],'r');
hold off
%}

%figure(2)
%box_graph (P40,'% aceleration')
%}
%{
figure(2)
bi_box_graph (V145,V10,'cm/s')
figure(2)
bi_box_graph (M145,M10,'cm/s in motion event')
    %}
%figure(3)
%bi_box_graph (-1*M145,-1*M10,'cm/s^2')

        