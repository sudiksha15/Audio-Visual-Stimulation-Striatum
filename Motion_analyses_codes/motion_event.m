function [stationary,moving,up_transition,down_transition]=motion_event(v,lt,ut)  %608450_01172020
%v is the speed trace
% is up transition here acceleration period? 
%moving is motion event
n=numel(v);
filter_size=20;
movmean_size=40;
sigmoid_slope=1;
threshold=0.45;
frame_threshold=20;
% Is this used anywhere?
gap_threshold=20;
filter_trace=zeros(n-filter_size+1,1);
for i=1:numel(filter_trace)
    f=diff(v(i:i+filter_size-1));
    % Why sum(f)?
    filter_trace(i)=sum(f);
end
ft=sigmf(filter_trace,[sigmoid_slope,(ut-lt)/3]);
nft=sigmf(-1*filter_trace,[sigmoid_slope,(ut-lt)/3]);
ft=movmean(ft,movmean_size);
nft=movmean(nft,movmean_size);
% Why 0.45
fft=ft>threshold;
nfft=nft>threshold;
A=bwlabel(fft);
B=bwlabel(nfft);


%
up_transition=[];
down_transition=[];
for i=1:max(A)
    a=find(A==i);
    t=v(a(1):a(end)+filter_size-1);
    [~,t1]=min(t);
    [~,t2]=max(t);
    %if t2-t1>frame_threshold
        %if isempty(up_transition) || a(1)-up_transition(end,2)>gap_threshold
            % Why a(1)+t1-1,a(1)+t2-1 ?
            up_transition=[up_transition;a(1)+t1-1,a(1)+t2-1];
        %end
    %end
   
end


for i=1:max(B)
    a=find(B==i);
    t=v(a(1):a(end)+filter_size-1);
    [~,t1]=max(t);
    [~,t2]=min(t);
    %if t2-t1>frame_threshold
        %if isempty(up_transition) || a(1)-up_transition(end,2)>gap_threshold
            % Why a(1)+t1-1,a(1)+t2-1 ?
            down_transition=[down_transition;a(1)+t1-1,a(1)+t2-1];
        %end
    %end
   
end
%}




    
    

%% determine stationary periods

stationary=[];
filter_size=40;
stationary_sigmoid_threshold=0.45;
frame_threshold=40;
sigmoid_slope=0.8;


l=1-sigmf(v,[sigmoid_slope,lt]);
L=movmean(l,filter_size);





L=L>stationary_sigmoid_threshold;
idx=bwlabel(L);
for i=1:max(idx)
       a=find(idx==i);
       N=numel(a);
       if N>frame_threshold 
          stationary=[stationary;a(1),a(end)];
       end
end
%}
%% new definition of moving period  
moving=[];
filter_size=1;
sigmoid_threshold=0.55;
frame_threshold=40;
sigmoid_slope=0.8;


l=sigmf(v,[sigmoid_slope,lt]); %want it to be above ut, just replace lt here with ut and change sigmoid threshold to 0.45
L=movmean(l,filter_size);



L=L>sigmoid_threshold;
idx=bwlabel(L);
for i=1:max(idx)
       a=find(idx==i);
       N=numel(a);
       if N>frame_threshold 
          moving=[moving;a(1),a(end)];
       end
end


% % old definition of moving period - above high speed period 
% moving=[];
% filter_size=1;
% sigmoid_threshold=0.45;
% frame_threshold=80;
% sigmoid_slope=0.8;
% 
% 
% l=sigmf(v,[sigmoid_slope,ut]); %want it to be above ut, just replace lt here with ut and change sigmoid threshold to 0.45
% L=movmean(l,filter_size);
% 
% 
% 
% L=L>sigmoid_threshold;
% idx=bwlabel(L);
% for i=1:max(idx)
%        a=find(idx==i);
%        N=numel(a);
%        if N>frame_threshold 
%           moving=[moving;a(1),a(end)];
%        end
% end



