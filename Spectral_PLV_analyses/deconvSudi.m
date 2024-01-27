function   [ s_oasis, onsets2, onsets]= deconvSudi(traces)

sampling_rate=20;
 threshold_min=0;tau_d =1400;  % decay time in ms
tau_r = 350;  % rise time in ms

[c_oasis, s_oasis] = deconvolveCa(traces','thresholded','smin',threshold_min,'optimize_b','window',round(sampling_rate*8));%, 'exp2',[tau_r tau_d  ],'thresholded','smin',threshold_min,'sampling_rate',sampling_rate,'optimize_b','window',round(sampling_rate*4));  %#ok<*ASGLU>
%[c_oasis, s_oasis] = deconvolveCa(traces', 'exp2',[tau_d tau_r  ],'thresholded','smin',threshold_min,'optimize_b','window',round(sampling_rate*14));%, 'exp2',[tau_r tau_d  ],'thresholded','smin',threshold_min,'sampling_rate',sampling_rate,'optimize_b','window',round(sampling_rate*4));  %#ok<*ASGLU>


thres=std(traces-fastsmooth(traces,10,1,1)).*2.5;% THRESHOLD


onsets=find(s_oasis>thres);
onsets2=[];
for id=2:length(onsets)
   if  onsets(id)-onsets(id-1) >3 % onsets(id)-onsets(id-1) >1  &  onsets(id+1)-onsets(id) >1
  onsets2=[onsets2,   onsets(id)];    
   end
end