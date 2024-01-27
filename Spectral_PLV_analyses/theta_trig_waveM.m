

function  [ data]=  delta_phase_wave(POW, trigs,winl,winr)

%size(FT.powspctrm);
%POW=FT.powspctrm;


nt=0;
for trial= 1:size(trigs,1)
    tig = find(trigs(trial,:));tig(tig==0)=[];
 
  
  for tl= 1:length(tig)
    
      if tig(tl)-winl >0  & tig(tl)+ winr <size(POW,2)
          if ~isnan(squeeze(POW(1,tig(tl))))
         nt=nt+1;
         % data(:,:,nt)= Z(:, tig(tl)-winl: tig(tl)+winr);
          data(:,:,nt)= squeeze(POW(:,tig(tl)-winl: tig(tl)+winr)); 
          end
      end
     
  end
    
    
    
end





