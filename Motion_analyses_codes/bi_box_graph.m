function bi_box_graph (AM,BM,label)
[n1,~]=size(AM);
[n2,~]=size(BM);
subplot(1,2,1)
%[h,p]=ttest(AM(:,1),AM(:,2));    %pair t test  h=1 reject null hypothesis the smaller the p the more significant the difference
[p,h]=signrank(AM(:,1),AM(:,2)); 
%n1=numel(AM)/2; %get the number of pairs
boxplot(AM);
colors =['b';'r'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.2);
end
% hold on
% for i=1:n1
%     plot([1,2],[AM(i,1),AM(i,2)],'-*','color',[0.5,0.5,0.5])
%     hold on
% end
hold off
%title(['145Hz h=',num2str(h),' p=',num2str(p)])
title('145Hz')
xlim([0,3])
xticks([1,2])
xticklabels({'Stim','Baseline'})
ylabel(label)
subplot(1,2,2)
%[h,p]=ttest(BM(:,1),BM(:,2));
[p,h]=signrank(BM(:,1),BM(:,2)); 
%n2=numel(BM)/2;
boxplot(BM);
colors = ['b';'r'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.2);
end
hold on
% for i=1:n2
%     plot([1,2],[BM(i,1),BM(i,2)],'-*','color',[0.5,0.5,0.5])
%     hold on
% end
hold off
%title(['10Hz h=',num2str(h),' p=',num2str(p)])
title('10Hz')
xlim([0,3])
xticks([1,2])
xticklabels({'Stim','Baseline'})