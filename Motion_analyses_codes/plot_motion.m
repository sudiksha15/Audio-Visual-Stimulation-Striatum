function [h1,h2]= plot_motion(MX,CV,m)
v=CV;
if numel(m)<1
   h2=plot(MX,v,'k');
   hold on
   h1=nan;
else
    [~,order]=sort(m(:,1));
    not_moving=m(order,1:2);
    n=numel(CV);
    index=1:n;
    not_moving=[not_moving;index(end),index(end)+1];

   

  
  if not_moving(1,1)>1
      I=1:numel(not_moving)/2-1;
     plot(MX(1:not_moving(1,1)),v(1:not_moving(1,1)),'k')
      hold on
  else
      plot(MX(not_moving(1,1):not_moving(1,2)),v(not_moving(1,1):not_moving(1,2)),'g')
      
      hold on
      plot(MX(not_moving(1,2):not_moving(2,1)),v(not_moving(1,2):not_moving(2,1)),'k')
      hold on
      I=2:numel(not_moving)/2-1;
  end
  
  for i=I
       h1=plot(MX(not_moving(i,1):not_moving(i,2)),v(not_moving(i,1):not_moving(i,2)),'g');
       
      hold on
     h2=plot(MX(not_moving(i,2):not_moving(i+1,1)),v(not_moving(i,2):not_moving(i+1,1)),'k');
      hold on
  end
  %h3=plot([1,MX(end)],[lt,lt],'r');
  %hold on
  %plot([1,MX(end)],[ut,ut],'r')
  

  %xlim([1,MX(end)])
end