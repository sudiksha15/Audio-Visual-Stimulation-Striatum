function [raw_motion,v] =Clear_Velocity(file_name,n_moving_average) % this function read in the txt file of 
%mouse on a ball motion data with a specific format and output a smoothened
%clean speed trace


 f=dir(file_name);
    Motion_data= importdata(f.name);
    timestamp = Motion_data.data(:,1);
    n=numel(timestamp);
    left_dx = Motion_data.data(:,2);
    left_dy = Motion_data.data(:,3);
    left_dt = Motion_data.data(:,4);
    right_dx = Motion_data.data(:,5);
    water_pin = Motion_data.data(:,8);
    right_dy = Motion_data.data(:,6);
    right_dt = Motion_data.data(:,7);

    ldy=left_dy./left_dt;
    rdy=right_dy./right_dt;
    ldx=left_dx./left_dt;
    rdx=right_dx./right_dt;
    
   
    sensor_angle = 75;
    sr = (sensor_angle/180)*pi;
    
    yl=ldy;
    yr=rdy;
    yl=(yl-yr*cos(sr))/sin(sr);
    v=sqrt(yl.^2+yr.^2)*100;
    [raw_motion,v]=t_filter(v,[0,100],n_moving_average);   %[min,max] represent the desire range of the speed trace
    
    
   
    
    
  