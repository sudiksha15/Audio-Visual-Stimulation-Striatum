function [T,IT,MT,X1,R1,sti,X,Rsti,Rbas,Msti,Mbas,nROI,cbt]=heatmap(plx,r_out,CV,m,onset,bo1,ao1,bo,ao)    %1 s before 2 s after motion onset occur within stimulation or baseline


nROI=numel(r_out);
n=numel(r_out(1).trace);
sti=plx.Stim_onset;
bas=[1,plx.Stim_offset];
sn=numel(m)/2;
IT=plx.Timestamp_Imaging(1:n);
sti(2,:)=plx.Stim_offset;
bas(2,:)=[plx.Stim_onset,n];
X1=-1*bo1/20:0.05:ao1/20;
R1=zeros(numel(sti)/2,bo1+ao1+1);
for i=1:numel(sti)/2
    R1(i,:)=sti(1,i)-bo1:sti(1,i)+ao1;
end



sti=IT(sti);
bas=IT(bas);
N=numel(CV);
MT=plx.Timestamp_Motion(1:N); 

T=zeros(nROI,n);
X=-1*bo/20:0.05:ao/20;
%B=T;
for i=1:nROI
   b1=r_out(i).BGtrace;
   a=r_out(i).trace;
   % subtract and detrend 
   a1=a-b1;
   
   a1=(a1-mean(a1))/(max(a1)-mean(a1))*100;
   a1=detrend(a1,'linear',1:1000:N);
   %b1=(b1-mean(b1))/(max(b1)-mean(b1))*100;  
   %b1=detrend(b1,'linear',1:1000:N);
   T(i,:)=a1;
   %B(i,:)=b1(con);
end

%really important parameter to adjust
%d=20*10;
%N=1200;






%for sti period
msti=[];
mbas=[];
nsn=numel(onset);
for i=1:nsn
    if sum(sti(1,:)<MT(onset(i)) & MT(onset(i))<sti(2,:))==1
        msti=[msti;onset(i)];
    elseif sum(bas(1,:)<MT(onset(i)) & MT(onset(i))<bas(2,:))==1
        mbas=[mbas;onset(i)];
    end
end
Mbas=mbas;
Msti=msti;
msti=Timeframe_Conversion(MT,msti,IT);
mbas=Timeframe_Conversion(MT,mbas,IT);
n1=numel(msti);
n2=numel(mbas);
Rsti=zeros(n1,bo+ao+1);
Rbas=zeros(n2,bo+ao+1);
for i=1:n1
    Rsti(i,:)=msti(i)-bo:msti(i)+ao;
end
for i=1:n2
    Rbas(i,:)=mbas(i)-bo:mbas(i)+ao;
end
    




%}

   



cbt=zeros(2,1);

 cb=sort(reshape(T,[1,nROI*n]),'descend');
 cbt(1)=mean(cb(1:round(numel(cb)/200)));
 cbt(2)=mean(cb(1:round(numel(cb)/10)));
 %cb2=sort(reshape(B,[1,nROI*N]),'descend');
 %cbt(2)=mean(cb2(1:round(numel(cb2)/200)));

  



