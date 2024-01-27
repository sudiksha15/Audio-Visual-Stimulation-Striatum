%function trace_clean_up
%deepdive2

clear all;
close all;

TT=2;% threshold for duration of negative event
Nt=[3,1];% determine the arrangement of subplot (3 by 1)
nt=Nt(1)*Nt(2);
a= uigetfile('*.mat','MultiSelect','on'); 
load(a{1})
load(a{2})

nROI=numel(r_out);
n=numel(r_out(1).trace);

sti=plx.Stim_onset;
sn=numel(sti);
frame=plx.Timestamp_Imaging(1:n);
step=frame(2)-frame(1);
sti(2,:)=plx.Stim_offset;

%{
a= r_out(1).file(1).filename;
b=find(a=='\');

fileDir=a(1:b(end)-1);

cd(fileDir)

 file=uigetfile('*.mat');
        load(file)
        
if exist('roi_list')
    ROI_list=roi_list;
end

imout=get_max_min;
%}


%B=r_out(1).overall_BG_trace;
T1=zeros(nROI,n);
%T2=zeros(nROI,n);
%T3=zeros(nROI,n);
NP=cell(nROI,1); %mean-2std
ML=zeros(nROI,1);
nNP=zeros(nROI,1);
for i=1:nROI
   b=r_out(i).BGtrace;
   a=r_out(i).trace;
   % subtract and detrend 
   c=detrend(a-b,'linear',1000:1000:n);
   T1(i,:)=c;
   %T2(i,:)=detrend(a);
   %T3(i,:)=detrend(b);
   Th=mean(c)-2*std(c);
   bw=c<Th;
   BW=bwlabel(bw);
   a1=max(BW);
   if a1~=0
        d=1:numel(BW);
        e=zeros(a1,2);
        for j=1:a1
            a2=BW==j;
            a3=d(a2);
            e(j,1)=a3(1);
            e(j,2)=a3(end);
        end
        f=e(:,2)-e(:,1)+1;
        g=max(f);
        h=f>TT;
        NP{i}=e(h,:);
        ML(i)=g;
   else
       NP{i}=[];
   end
   [nNP(i),~]=size(NP{i});
 
end


[~,idx]=sort(nNP,'descend');

N=floor(nROI/nt);
remainder=nROI-nt*N;

index=[];
for i=1:N+1
    figure(1)
    if i~=N+1
        for j=1:nt
            subplot(Nt(1),Nt(2),j)
            plot(frame(1:n),T1(idx(j+(i-1)*nt),:))
            xlim([frame(1),frame(n)])
            lt=min(T1(idx(j+(i-1)*nt),:));
            ut=max(T1(idx(j+(i-1)*nt),:));
            ylim([lt-0.2*(ut-lt),ut+0.2*(ut-lt)]);
            title(num2str(nNP(idx(j+(i-1)*nt))))
            
            hold on
            for k=1:sn
                son=frame(sti(1,k));
                soff=frame(sti(2,k));
                plot([son,son],[lt,ut],'r')
                hold on
                plot([soff,soff],[lt,ut],'r')
                hold on
            end
            hold off
            
            set(gca,'tag',num2str(idx(j+(i-1)*nt)))
        end
    else
        for j=1:remainder
            subplot(ceil(remainder/Nt(2)),Nt(2),j)
            plot(frame(1:n),T1(idx(j+(i-1)*nt),:))
            xlim([frame(1),frame(n)])
            lt=min(T1(idx(j+(i-1)*nt),:));
            ut=max(T1(idx(j+(i-1)*nt),:));
            ylim([lt-0.2*(ut-lt),ut+0.2*(ut-lt)]);
            title(num2str(nNP(idx(j+(i-1)*nt))))
            hold on
            for k=1:sn
                son=frame(sti(1,k));
                soff=frame(sti(2,k));
                plot([son,son],[lt,ut],'r')
                hold on
                plot([soff,soff],[lt,ut],'r')
                hold on
            end
            hold off
            set(gca,'tag',num2str(idx(j+(i-1)*nt)))
        end
    end
    
    while true
        w = waitforbuttonpress;
        switch w 
            case 1
                break
            case 0 
                mousept = get(gca,'currentPoint');
                index=[index;str2num(get(gca,'tag'))];
        end
    end    
end

[b,m1,n1] = unique(index,'first');
[c1,d1] =sort(m1);
index = b(d1);


new_r_out=r_out;
new_r_out(index)=[];
            






%{
%at least 5 negative peak are chosen randomly for each ROI
for i=1:nROI
    [N,~]=size(NP{i});
    if N>5
        nn=randsample(N,5,false);
    elseif N~=0
        nn=1:N;
    else
        nn=[];
    end
    if ~isempty(nn)
        nT1=zeros(numel(nn),10*ML(i));
        nT2=nT1;
        nT3=nT2;
        for j=1:numel(nn)
            a=NP{i}(nn(j),:);
            b=10*ML(i)-(a(2)-a(1)+1);
            b1=ceil(b/2);
            b2=floor(b/2);
            c=[a(1)-b1,a(2)+b2];
            if c(1)>0 && c(2)<n
            nT1(j,:)=T1(i,c(1):c(2));
            nT2(j,:)=T2(i,c(1):c(2));
            nT3(j,:)=T3(i,c(1):c(2));
            end
        end
            test=max(nT1,[],2);
            t1=test==0;
            nT1(t1,:)=[];
            nT2(t1,:)=[];
            nT3(t1,:)=[];
            nn(t1)=[];
            
           if ~isempty (nn)
        
        figure(1)
        plot_T({nT1,nT2,nT3},{'trace-BGtrace',['trace  ROI:',num2str(i)],'BGtrace'},[],nn)
        pause
           end
    end
end
            
  %}          
        
        






%{
figure(1)

x=repmat([frame;frame(end)+step]',nROI+1,1);
y=repmat((1:nROI+1)',1,n+1);

cb=sort(reshape(nT1(1:nROI,1:n),[1,nROI*n]),'descend');
    cbt=mean(cb(1:round(numel(cb)/100)));
h=pcolor(x,y,nT1);
title('donut')
    caxis([-1,cbt])
       set(h, 'EdgeColor', 'none');
       colorbar
%}
       
       %{
       figure(3)
       plot(a1)
       title('negative event in donut')
       
       
       
figure(4)
cb=sort(reshape(nT2(1:nROI,1:n),[1,nROI*n]),'descend');
    cbt=mean(cb(1:round(numel(cb)/100)));
h=pcolor(x,y,nT2);
title('overall BG')
    caxis([-1,cbt])
       set(h, 'EdgeColor', 'none');
       colorbar
       r1=sum(nT2,2);
       [a2,I2]=sort(r1,'descend');
       
       
        figure(5)
       plot_traceline(r_out,plx,I2(1:5),[-al,al])
       
       figure(6)
       plot(a2)
       title('negative event in overall BG')
       

%}








%{
top=5;
[~,IDX1]=sort(ns1);
[~,IDX2]=sort(ns2);
figure(1)
plot_traceline(r_out,plx,IDX1(1:top))
figure(2)
if exist('outline')
    overlap_ROI(imout,constract,outline,IDX1(1:top),1)
else
    overlap_ROI(imout,constract,ROI_list,IDX1(1:top),0)
end



figure(3)
plot_traceline(r_out,plx,IDX2(1:top))
figure(4)
constract=[0,5000];
if exist('outline')
    overlap_ROI(imout,constract,outline,IDX2(1:top),1)
else
    overlap_ROI(imout,constract,ROI_list,IDX2(1:top),0)
end

%}