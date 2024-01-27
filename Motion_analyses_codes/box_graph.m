function box_graph (AM,label)
[n1,~]=size(AM);
%[h,p]=ttest(AM(:,1),AM(:,2));    %pair t test  h=1 reject null hypothesis the smaller the p the more significant the difference
[p,h]=signrank(AM(:,1),AM(:,2)); 
boxplot(AM);
hold on
for i=1:n1
    plot([1,2],[AM(i,1),AM(i,2)],'-*','color',[0.5,0.5,0.5])
end
hold off
title(['40Hz h=',num2str(h),' p=',num2str(p)])
xlim([0,3])
xticks([1,2])
xticklabels({'during sti','bas'})
ylabel(label)
