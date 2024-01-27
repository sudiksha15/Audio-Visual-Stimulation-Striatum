clear all;
close all;
[A,b]= uigetfile('*','MultiSelect','on');  %load 2 pairs data from motion
cd(b)
for j=1:numel(A)/2
    a=A(2*j-1:2*j);
    kk=strfind(a{1},'_');
    name=a{1}(1:kk(3));
    name=strrep(name,'_',' ');
    a1=find(contains(a,'txt'));
    a2=find(contains(a,'plex'));
    CV =Clear_Velocity(a{a1},20);
    n=numel(CV);
    load(a{a2})



    [lt,ut]=var_low_high_speed(CV,20,0.1,2.5,2,2);

    [~,m,~,~]=motion_event(CV,lt,ut);
    

    %
    bCV=zeros(n,1);
    if ~isempty(m)
    for i=1:numel(m)/numel(m(1,:))
        bCV(m(i,1):m(i,2))=1;
    end
    end
%}
    
     %{
        for i=1:numel(m)/numel(m(1,:))
            bCV(m(i,1):m(i,2))=1;
            plot_motion(1:n,CV,m);
            xlim([m(i,1)-50,m(i,1)+50])
            plot([m(i,1),m(i,1)],[0,max(CV)],'r')
            plot([m(i,1)-50,m(i,1)+50],[lt,lt],'r');
            plot([m(i,1)-50,m(i,1)+50],[ut,ut],'r');
            hold off
            pause
        end
        %}

    
    MT=plx.Timestamp_Motion;
    if numel(MT)<n
        step=mean(diff(MT));
        st=MT(end)+step:step:MT(end)+(n-numel(MT))*step;
        MT=[MT;st'];
    else
        MT=MT(1:n);
    end

    ST=plx.Timestamp_stim;
    sti =sti_extraction(ST);
    sn=numel(sti)/2;
    %sn=numel(ST);
    %sti =ST';
    %sti(2,:)=ST'+10;
    %sti=plx.Stim_onset;
    %sti(2,:)=plx.Stim_offset;
    idx=find_idx(MT,sti);
    oc=[];
    for i=1:sn
        oc=[oc,CV(idx(1,i)-1200:idx(1,i)+4800)];
    end
    oc=mean(oc,2);
    %{
        bsti=false(n,1);
    bbas=true(n,1);
    for i=1:numel(idx)/numel(idx(:,1))
        bsti(idx(1,i):idx(2,i))=true;
        bbas(idx(1,i):idx(2,i))=false;
        bbas(idx(2,i)+1:idx(2,i)+200)=false;
    end
    bstim=bsti&bm;
    bbasm=bbas&bm;
    %}
    %{
    v=CV;
    mlt=ones(numel(v),1)*lt;%mlt is initialized with all points equal to lt
ssv=v>lt;%ssv is the initial binary trace with the initial "islands of ones"


A1=bwlabel(ssv);
M1=[];
for i=1:max(A1)
    B1=find(A1==i);
    M1=[M1;B1(1),B1(end)];
end
nlt=lt;%nlt will be locally adjusted
loop=2;%this process is repeated two times
for k=1:loop
    nlt=nlt+(ut-lt)/10;%increase nlt by 1/10 of the distance between ut and lt    
    A=bwlabel(ssv);%assign indices to each islands
    for i=1:max(A)
        B=find(A==i);%find the location of every island
        if numel(B)>200%only perform adjustment on islands longer than 200 frames
            ssv(B)=v(B)>nlt;%binary trace is locally updated
            mlt(B)=nlt;%local low speed threshold is updated
        end
    end
end
A=bwlabel(ssv);
M2=[];
for i=1:max(A)
    B=find(A==i);
    M2=[M2;B(1),B(end)];
end
    figure(1)
    H1=subplot(2,1,1);
    [h1,~]=plot_motion(MT,CV,M1);
    plot([MT(1),MT(end)],[lt,lt],'r')
    ylabel('speed (cm/s)')
    xlim([min(MT),max(MT)])
    title([name,' fix low speed threshold'])
    legend(h1,'motion events candidate')
    H2=subplot(2,1,2);
    plot_motion(MT,CV,M2);
    h2=plot(MT,mlt,'r');
    xlabel('seconds')
    xlim([min(MT),max(MT)])
    title('locally adjusted low speed threshold')
    legend(h2,'low speed threshold')
    linkaxes([H1,H2],'xy')
    pause
    %}
    %
    figure(1)
    [h1,~]=plot_motion(MT,CV,m);
    plot([MT(1),MT(end)],[lt,lt],'r')
    plot([MT(1),MT(end)],[ut,ut],'r');
    %
    for i=1:sn
        plot([sti(1,i),sti(1,i)],[0,max(CV)],'r')
        plot([sti(2,i),sti(2,i)],[0,max(CV)],'r')
    end
        
    hold off
    ylabel('speed (cm/s)')
    xlabel('seconds')
    xlim([MT(1),MT(end)])
    title(['speed trace of ',name]) 
    legend(h1,'motion event')


    
    %
    figure(2)
    plot(-1200:4800,oc)
    hold on
    plot([0,0],[0,max(oc)],'r')
    plot([1200,1200],[0,max(oc)],'r')
    hold off
    title('speed traces aligned to stimulations')
    ylabel('speed cm/s')
    xlabel('frames after onsets')
    %}
%{
    %H1=subplot(2,4,[1,2,3]);
    plot_motion(MT,CV,m);
    plot([min(MT),max(MT)],[lt,lt],'r')
    plot([min(MT),max(MT)],[ut,ut],'r')
    for i=1:numel(idx)/numel(idx(:,1))
        plot([sti(1,i),sti(1,i)],[0,max(CV)],'r')
        plot([sti(2,i),sti(2,i)],[0,max(CV)],'r')
    end
    hold off
    xlim([min(MT),max(MT)])
    title(name)
    ylabel('speed (cm/s)')
    xlabel('seconds')
    %}
    %{
    figure(2)
    subplot(2,1,1)
    plot(CV(bstim))
    hold on
    plot([0,sum(bstim)],[mean(CV(bstim)),mean(CV(bstim))],'r')
    hold off
    ylim([0,max(CV)])
    xlim([0,sum(bsti)])
    title(['all motion event in sti speed trace ',num2str(mean(CV(bstim)))])
    subplot(2,1,2)
    plot(CV(bsti))
    hold on
    plot([0,sum(bsti)],[mean(CV(bsti)),mean(CV(bsti))],'r')
    hold off
    ylim([0,max(CV)])
    xlim([0,sum(bsti)])
    title(['sti speed trace ',num2str(mean(CV(bsti))),' cm/s ',num2str(sum(bstim)/sum(bsti)*100),'% of time in motion'])
    figure(3)
    subplot(2,1,1)
    plot(CV(bbasm))
    hold on
    plot([0,sum(bbasm)],[mean(CV(bbasm)),mean(CV(bbasm))],'r')
    hold off
    ylim([0,max(CV)])
    xlim([0,sum(bbas)])
    title(['all motion event in bas speed trace ',num2str(mean(CV(bbasm))),' cm/s'])
    subplot(2,1,2)
    plot(CV(bbas))
    hold on
    plot([0,sum(bbas)],[mean(CV(bbas)),mean(CV(bbas))],'r')
    hold off
    ylim([0,max(CV)])
    xlim([0,sum(bbas)])
    title(['bas speed trace ',num2str(mean(CV(bbas))),' cm/s ',num2str(sum(bbasm)/sum(bbas)*100),'% of time in motion'])
    %}
    %{
    H3=subplot(2,4,4);
    for i=1:numel(m)/numel(m(1,:))
        if m(i,1)>40
        plot(-40:40,CV(m(i,1)-40:m(i,1)+40))
        else
        plot(1-m(i,1):40,CV(1:m(i,1)+40))
        end
        hold on
    end
    title('overlay all motion events')
    xlabel('frames')
    ylabel('speed cm/s')
    xlim([-40,40])
    hold off
m=M;
        H2=subplot(2,4,[5,6,7]);
    plot_motion(MT,CV,m);
    plot([min(MT),max(MT)],[lt,lt],'r')
    plot([min(MT),max(MT)],[ut,ut],'r')
    hold off
    xlim([min(MT),max(MT)])
    title('old method')
    xlabel('seconds')
    ylabel('speed (cm/s)')
    H4=subplot(2,4,8);
    for i=1:numel(m)/numel(m(1,:))
        if m(i,1)>40 && m(i,1)+40<=n
        plot(-40:40,CV(m(i,1)-40:m(i,1)+40))
        elseif m(i,1)<=40
        plot(1-m(i,1):40,CV(1:m(i,1)+40))
        elseif m(i,1)+40>n
        plot(-40:n-m(i,1),CV(m(i,1)-40:n))
        end
        hold on
    end
    title('overlay all motion events')
    xlabel('frames')
    ylabel('speed cm/s')
    xlim([-40,40])
    hold off
    linkaxes([H1,H2],'xy')
    linkaxes([H3,H4],'xy')
    
    %}
    pause
  

end

        