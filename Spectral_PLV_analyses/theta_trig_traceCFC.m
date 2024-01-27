

function  [ data]= theta_trig_traceCFC(dat, trigs,  lwin,rwin)

%size(FT.powspctrm);
%POW=FT.powspctrm;
winl=lwin
winr=rwin;
nt=0;
for trial= 1%:size(trigs,2)
    tig = find(trigs);%tig(tig==0)=NaN;
 
  
  for tl= 1:length(tig)
    
      if tig(tl)-winl >0  & tig(tl)+ winr <size(dat,2)
         nt=nt+1;
         % data(:,:,nt)= Z(:, tig(tl)-winl: tig(tl)+winr);
          data(:,nt)= squeeze(dat(tig(tl)-winl: tig(tl)+winr)); 
            
      end
     
  end
    
    
    
end





