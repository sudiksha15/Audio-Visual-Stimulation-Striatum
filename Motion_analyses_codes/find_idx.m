function idx=find_idx(Tf,Tb)
idx=Tb;
for i=1:numel(Tb)
    [~,c]=min(abs(Tf-Tb(i)));
    idx(i)=round(c);
end