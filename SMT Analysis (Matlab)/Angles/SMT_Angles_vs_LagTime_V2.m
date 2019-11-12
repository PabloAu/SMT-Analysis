function [Angles, Lag_times] = SMT_Angles_vs_LagTime_V2(tracktemp,Frame_interval,loc_error);

      
%Connecting just the GAPS of 1 frame.
local_length = size(tracktemp,1);
%Connect GAPS of 1 frame
gaps1 = find(abs(diff(tracktemp(:,1)) - 2*Frame_interval) <= eps(0.5));
            if isempty(gaps1);
            else
                for nn=1:size(gaps1);
            local_length = local_length + 1;     
            tracktemp(gaps1(nn)+2:local_length,:) = tracktemp(gaps1(nn)+1:end,:);   

            tracktemp(gaps1(nn)+1,1) = tracktemp(gaps1(nn)+1,1)-Frame_interval; 
            tracktemp(gaps1(nn)+1,2) = ((tracktemp(gaps1(nn)+1,2) + tracktemp(gaps1(nn),2))/2) + (2.*rand(1,1) -1)*loc_error; %Coordinates are in um
            tracktemp(gaps1(nn)+1,3) = ((tracktemp(gaps1(nn)+1,3) + tracktemp(gaps1(nn),3))/2) + (2.*rand(1,1) -1)*loc_error; %Coordinates are in um
                gaps1 = gaps1 + 1;
                end
            end
           
            
 
        dvector = [];
        Lag_times = [];        
        
        %Iterate on all the possible delays
       for i=1:round((length(tracktemp(:,1))/2)-1); 
           
           last_point = floor((length(tracktemp(:,1))-1)/i)*i + 1;     
           dvector = [((tracktemp(i+1:i:last_point,2) - tracktemp(1:i:last_point-i,2))) ((tracktemp(i+1:i:last_point,3) - tracktemp(1:i:last_point-i,3)))]; %In um
           dvector(:,3) = 0; %Make the vector 3D for the cross product. 

                SignTan = [];
                TanTheta_1 = [];
                TanTheta_2 = [];    
                for hh=1:1:size(dvector,1)-1;        
                      TanTheta_1 = [TanTheta_1; norm(cross(dvector(hh,:),dvector(hh+1,:)))];
                      TanTheta_2 = [TanTheta_2; dot(dvector(hh,:),dvector(hh+1,:))];
                      signo = sign(cross(dvector(hh,:),dvector(hh+1,:)));
                      SignTan = [SignTan; signo(3)];  

                end  


            Angles{i} = SignTan.*atan2(TanTheta_1,TanTheta_2); %Angles (sign included)
            Lag_times(end+1) = Frame_interval*i;

       end
       
       
       

end

