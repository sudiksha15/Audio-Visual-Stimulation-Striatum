function idx2=Timeframe_Conversion(Tf1,idx1,Tf2)
% Tf1-timestamp for idx1, Tf2-timestamp you want to convert idx1 
time=Tf1(idx1);
idx2=idx1;
for i=1:numel(idx1)
    [~,c]=min(abs(Tf2-time(i)));
    idx2(i)=round(c);
end