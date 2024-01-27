

function fill_error_area2(x,y,errory, shade_color)

x_low = x;
y_low = y -  errory;
x_high = x;
y_high =  y + errory;

 x_high = x_high(end:-1:1);
 y_high = y_high(end:-1:1);

if size(x,1) ~= 1 
x_vec = [x_low; x_high];
y_vec = [y_low; y_high];
else
x_vec = [x_low x_high];
y_vec = [y_low y_high];
end;

x_vec = x_vec(find(~isnan(x_vec)));
y_vec = y_vec(find(~isnan(y_vec)));

figure(gcf); 
hold_stat = ishold;
hold on;

h = fill(x_vec, y_vec, shade_color,'EdgeAlpha',1,'FaceAlpha',0.1);
  
% 
%  hold on; plot(x,y, '-', 'MarkerFaceColor',shade_color,'MarkerSize',5,'Marker','o',...
%     'Color',shade_color, 'Visible','on')
% %hold on; plot(x,y, '-', 'MarkerFaceColor',shade_color,...
 %   'Color',shade_color, 'Visible','on')
% set(h, 'Linewidth' , 2, 'Color', [ 0.1 0.1 0.1])
