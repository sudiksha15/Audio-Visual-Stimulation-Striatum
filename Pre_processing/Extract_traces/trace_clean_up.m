%function trace_clean_up
%deepdive2
clear all;
close all;

TT=2;% threshold for duration of negative event
nt=20;% display 20 ROIs at once, could adjust it base on size of screen multiple of 2
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
   c=detrend(a-b);
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
            subplot(nt/2,2,j)
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
            subplot(ceil(remainder/2),2,j)
            plot(frame(1:n),T1(idx(j+(i-1)*nt),:))
            xlim([frame(1),frame(n)])
            lt=min(T1(idx(j+(i-1)*nt),:));
            ut=max(T1(idx(j+(i-1)*nt),:));
            ylim([lt-0.5*(ut-lt),ut+0.5*(ut-lt)]);
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
save('608448_01102020_145Hz_clean_trace.mat','new_r_out','-v7.3')
            





