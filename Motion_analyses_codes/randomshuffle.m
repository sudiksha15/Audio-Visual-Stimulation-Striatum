clear all;  
close all;
[A,b]= uigetfile('*','MultiSelect','on'); %load all 2 pairs data
cd(b)
for j=1:numel(A)/2
    a=A(2*j-1:2*j);
    a1=contains(a,'txt');
        CV =Clear_Velocity(a{a1},20);
        [lt,ut]=var_low_high_speed(CV,20,0.1,2.5,2,2);
        [~,m,~,~]=motion_event(CV,lt,ut);
        %m=motion_onset(CV,lt,ut,0);
        % Add lines for onset
        % 
      
        n=numel(CV);
        bm=zeros(1,n);
        if ~isempty(m)
        for i=1:numel(m)/numel(m(1,:))
            bm(m(i,1):m(i,2))=1;
           %bm(m(i))=1;
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
        p1=mean(bstim)*100;

        bbasm=bm(bbas);
        T=[];
        for i=1:1000
        ransam= randsample(bbasm,numel(bstim));
        T=[T;sum(ransam)/numel(bstim)*100];
        p2=prctile(T,95);
        end
        h=histogram(T,'EdgeColor','none');
        hold on
        h1=plot([p1,p1],[0,max(h.BinCounts)],'r');
        h2=plot([p2,p2],[0,max(h.BinCounts)],'k');
        legend([h1,h2],{'% in sti','95 precentile'})
        hold off
        title(name)
        xlabel('# of shuffle')
        ylabel('% of motion')
        pause
end
        
        
        