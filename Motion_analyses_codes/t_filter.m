function [raw_motion,output]=t_filter(input,range,n)
index=1:numel(input);
reasonable_vel= (input>=range(1)) & (input<=range(2));
input=interp1(index(reasonable_vel),input(reasonable_vel),index);
raw_motion=input';
output = movmean(input,n);
output=output';
