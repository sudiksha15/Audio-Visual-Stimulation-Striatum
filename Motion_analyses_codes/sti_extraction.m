function sti =sti_extraction(sti_timestamp) %this only work for audiovisual stimulation with 5 trials
a=diff(sti_timestamp);
[~,a2]=sort(a,'descend');
b=a2(1:4);
a=sort(b);
a=[0;a;numel(sti_timestamp)];
sti=zeros(2,5);
for i=1:5
    sti(1,i)=sti_timestamp(a(i)+1);
end
for i=1:5
    sti(2,i)=sti_timestamp(a(i+1));
end