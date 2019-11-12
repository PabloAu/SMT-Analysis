%Identify bound tracks

%%The Tracks need to fulfil the condition Rmax for EVERY frame

function Bound_trayectories=identify_boundtracks_V3(trayectories,Bound_threshold);

previous = 0;  % previous = 0 when the previous jump didn't meet 
               % the filter criterium (jump shorter than 3 * locAcc),
               % 1 if it did.
particleIx = ones(1,length(trayectories));

Bound_trayectories={};

for i=1:length(trayectories)
   
    tracktemp = trayectories{1,i};
    tracktemp_length(i) = size(tracktemp,1);
    jd_temp = sqrt((tracktemp(2:end,2)- tracktemp(1:end-1,2)).^2 + (tracktemp(2:end,3)- tracktemp(1:end-1,3)).^2);
    
    for j = 1: length(jd_temp);
     
         if jd_temp(j) < Bound_threshold(1) && ~previous
             
             previous = 1;
             
         elseif jd_temp(j) < Bound_threshold(1)  && previous
           
            particleIx(1,i) = particleIx(i) +1;
            
          
         else
              previous = 0; 
              particleIx(1,i)=1;
         end
     end
     previous = 0;
end



    
 
index=find(particleIx>=(tracktemp_length-1));



for i=1:length(index);
Bound_trayectories{i,1}=trayectories{1,index(i)};

end




end