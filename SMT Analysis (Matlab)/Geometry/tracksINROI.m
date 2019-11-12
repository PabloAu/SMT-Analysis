%% Extract only the tracks in the subROIs


%Tracks_out_1 are the Tracks inside the subROIs
%Tracks_out_2 are the Tracks outside the subROIs

function [tracks_out_1, tracks_out_2,list_in] = tracksINROI(tracks_in,subROIs,pixel_size,image_names,subROI_names);

tracks_out_1 = tracks_in;
tracks_out_2 = tracks_in;


for i=1:size(tracks_in,2); %Iterate on each cell
    
    idx = find(subROI_names == image_names(i));
    local_subROIs = subROIs{idx};

    for ii=1:size(tracks_in{i},1); %Iterate on each Track

       pos_local_track = tracks_in{i}{ii,1}(:,2:3); 

        for iii=1:size(local_subROIs,2); %Iterate on each subROI
            pos_local_subROIs = (local_subROIs{iii}.mnCoordinates)*pixel_size;
            
                  for iiii = 1:size(pos_local_track,1); %Iterate on each point of each Track

                    in(iiii) = inpoly(pos_local_track(iiii,1:2),pos_local_subROIs);
                  end
                  
                   list_in{i,ii}{iii} = in;
                    if length(find(in == 1)) ==  size(pos_local_track,1);
                    subROI_positive(iii) = 1;

                    else
                    subROI_positive(iii) = 0;
                    end 
                    
                    clear in;
                    clear pos_local_subROIs;
        end
        
                    if isempty(find(subROI_positive == 1));
                    tracks_out_1{i}{ii,1}=[]; %Empty the Tracks 

                    else
                   
                    tracks_out_2{i}{ii,1}=[]; %Empty the Tracks 
                    
                    end
                    
        clear subROI_positive;
        clear pos_local_track;
    end
            clear local_subROIs;

end





end
