%Identify "Butterfly" tracks

%%The Tracks needs to have one jump bigger than mean(jump) +
%%jump_threshold*std(jump)  and  
%%total distance travelled > averagejump*threshold.



function [Butterfly_trayectories,track_segment,Reference_Track_segment,motion_type_segment,TOP_index,BigJumps,BigJump_idx] = identify_Butterfly_tracks_V4(trayectories,jump_threshold,minim_dist,TLOGLOG_fitting_points,R2LIMIT,min_alpha_directed,level,num_points,D0,R0,Conf2JumpThr,Out_percentage,N,P_min_linear,angleTH,Max_Jump_Confined_Butterfly);

fprintf('\n');
fprintf('Detecting Butterfly motion');
fprintf('\n');
fprintf('--------------------------');
fprintf('\n');

%% Add trayectories to @msdanalizer
ma_temp = msdanalyzer(2, 'µm', 's'); 
ma_temp= ma_temp.addAll(trayectories);


%% Second step: identify candidate "butterfly" Tracks.               
BigJumps = zeros(1,length(trayectories));
BigJump_idx = {}; 
Butterfly_trayectories={};


for i=1:length(trayectories);  
    tracktemp = trayectories{1,i};
    [max_dist(i), min_dist(i), avg_dist(i)] = distance_scatter(tracktemp(:,2:3));
    jd_temp{i} = sqrt((tracktemp(2:end,2) - tracktemp(1:end-1,2)).^2 + (tracktemp(2:end,3) - tracktemp(1:end-1,3)).^2);

    
    % 1. Calculate the linearity of segments of N points (N is specified as an input).
    dint = sqrt(((tracktemp(N:1:end,2)-tracktemp(1:1:end-(N-1),2)).^2 + (tracktemp(N:1:end,3)-tracktemp(1:1:end-(N-1),3)).^2));
    % Sum of displacements for first tint points (um)
    sumint=[];
    for k=1:length(jd_temp{i})-(N-2);
        sumint = [sumint;sum(jd_temp{i}(k:k+(N-2)))];
    end
    ratio= dint./sumint;
    %Assign an average ratio to each data point
    P=[];
    for hh=1:length(ratio)+(N-1);
    if hh<N;
        Pc=mean(ratio(1:1:hh));
    elseif hh>length(ratio)
        Pc=mean(ratio(hh-(N-1):1:end));
    else
        Pc=mean(ratio(hh-(N-1):1:hh));
    end
        P =[P; Pc];
    end

    %Identify points that are in a linear segment     
    for hh=1:length(P)
        isdirect_point(hh)= P(hh) >= P_min_linear; % Each data point has a probability assigned between 0 and 1.
    end

    for hh=1:length(jd_temp{i});
        if isdirect_point(hh) == 1 || isdirect_point(hh+1) == 1;
            isdirect_segment{i}(hh) = 1;
        else
            isdirect_segment{i}(hh) = 0;
        end   
    end
    
    
  % 2. Calculate angles between segments:
Sign = []; 
NV = []; 
CosTheta = [];
dvector = [(tracktemp(2:1:end,2)-tracktemp(1:1:end-1,2)) (tracktemp(2:1:end,3)-tracktemp(1:1:end-1,3))];
for hh=1:1:length(dvector);
    NV=[NV;norm(dvector(hh,:))];
end
for hh=1:1:length(dvector)-1;
    CosTheta = [CosTheta; (dot(dvector(hh,:),dvector(hh+1,:)))./(NV(hh)*NV(hh+1))];
    Sign = [Sign; sign(dvector(hh,1)*dvector(hh,2))];
end
ThetaInDegrees = acos(CosTheta)*180/pi; %Angle in degrees (sign not included)
ThetaInDegreesSign = Sign.*acos(CosTheta)*180/pi; %Angle in degrees (sign included)

% Identify directed segments based on angle condition
isdirect_angle{i}(1)=isdirect_segment{i}(1);
for hh=1:length(ThetaInDegrees);
    isdirect_angle{i}(hh+1) = ThetaInDegrees(hh) > angleTH;
end


    
    %Scan through every segment and identify BigJUmps
    for j = 1: length(jd_temp{i});
    jd_temp_mean = mean(jd_temp{i}([1:j-1 j+1:end]));
    jd_temp_std = std(jd_temp{i}([1:j-1 j+1:end]));

   % 3. Check if the Jump is OUTSIDE or INSIDE the polygon defined by previous a posterior segments
    track_temp_pre = tracktemp(1:j,:);
    track_temp_JUMP = tracktemp(j:j+1,:);
    xvals = linspace(track_temp_JUMP(1,2), track_temp_JUMP(2,2), 100);
    yvals = linspace(track_temp_JUMP(1,3), track_temp_JUMP(2,3), 100);
    track_temp_post = tracktemp(j+1:end,:);
    k_pre = boundary(track_temp_pre(:,2),track_temp_pre(:,3));
    k_post = boundary(track_temp_post(:,2),track_temp_post(:,3));
  
        if size(k_pre,1)>2;     
            inside_points_pre{i}{j} = inpolygon(xvals,yvals,track_temp_pre(k_pre,2),track_temp_pre(k_pre,3));
        else
            inside_points_pre{i}{j} = zeros(1,100);
        end

        
        if size(k_post,1)>2;
            inside_points_post{i}{j} = inpolygon(xvals,yvals,track_temp_post(k_post,2),track_temp_post(k_post,3));
        else
            inside_points_post{i}{j} = zeros(1,100);
        end
        
       
           
 % 4. Based on previous calculations, identify the candidates and that
 % perform BigJumps
        if jd_temp{i}(j) >= jd_temp_mean+jump_threshold*jd_temp_std && max_dist(i) >= minim_dist*jd_temp_mean;
        BigJumps(1,i) = BigJumps(1,i) + 1;
        BigJump_idx{i}(j) = 1;
        Jump_displacement{i}(j) = jd_temp{i}(j);
        end  
         
    end      
    
end


%Retrieve the indices of the BigJumps and the general index of the
%"Butterfly trayectories".
index=find(BigJumps>=1);
TOP_index = index;

BigJumps = nonzeros(BigJumps);
BigJump_idx = BigJump_idx(~cellfun('isempty',BigJump_idx));
Jump_displacement = Jump_displacement(~cellfun('isempty',Jump_displacement));



for i=1:length(index);
Butterfly_trayectories{i,1}=trayectories{1,index(i)};
BigJump_idx{i}(size(BigJump_idx{i},2)+1:size(Butterfly_trayectories{i,1},1)-1) = 0;
Jump_displacement{i}(size(Jump_displacement{i},2)+1:size(Butterfly_trayectories{i,1},1)-1) = 0;
end




%% Third step: Segmentate the "butterfly" Tracks
for i=1:size(Butterfly_trayectories,1);    
   
idxxx = find(BigJump_idx{i} == 1);  

    for hh=1:length(isdirect_segment{index(i)});
        if (isdirect_segment{index(i)}(hh) == 1 && isdirect_angle{index(i)}(hh) == 1) && length(find(inside_points_pre{index(i)}{hh} == 0)) >= Out_percentage && length(find(inside_points_post{index(i)}{hh} == 0)) >= Out_percentage;
            isdirect{i}(hh)=1;
        else
            isdirect{i}(hh)=0;
        end
    end

    
     for hh=1:length(isdirect{i});
        if BigJump_idx{i}(hh) == 1;
            directed_segment{i}(hh)=1;
        else if (isdirect{i}(hh) == 1) && ~isempty(find(BigJump_idx{i}(max(1,hh-1):min(length(isdirect{i}),hh+1))) == 1);
            directed_segment{i}(hh)=1;      
            else
            directed_segment{i}(hh)=0;
            end
        end
     end   
    
    %reconnect segments that are "linear" and next to a "directed" segment
    for mm=1:10;
        for hh=1:length(isdirect{i});
            if (isdirect{i}(hh) == 1) && ~isempty(find(directed_segment{i}(max(1,hh-1):min(length(isdirect{i}),hh+1))) == 1);
            directed_segment{i}(hh)=1;      
            end
        end
    end
    
        %Finally, classify as "directed" segments every jump bigger than a
        %threshold.
        for hh=1:length(isdirect{i});
            if jd_temp{index(i)}(hh) > Max_Jump_Confined_Butterfly;
            directed_segment{i}(hh)=1;      
            end
        end
    
  
Segment_cuts{i} = find(diff(directed_segment{i})~=0)+1;
Segment_cuts{i} = [1 Segment_cuts{i} size(Butterfly_trayectories{i,1},1)];



for n=1:length(Segment_cuts{i})-1;   
track_segment{i}{n} = Butterfly_trayectories{i,1}(Segment_cuts{i}(n):Segment_cuts{i}(n+1),:);

if directed_segment{i}(Segment_cuts{i}(n)) == 0;
motion_type_segment{i}(n) = 1; % 1 will correspond to Confined segments
else
motion_type_segment{i}(n) = 0; % 0 will correspond to NOT Confined segments 
end

Reference_Track_segment{i}(n) = i; %Retrieve the Track ID for each segment.
end


end

Butterfly_trajectories_segments = horzcat(track_segment{:}); %Reshape cell array to colapse all the cells in the top level.   
Motion_type_segment = horzcat(motion_type_segment{:}); %Reshape cell array to colapse all the cells in the top level.  
Reference_Track_segment = horzcat(Reference_Track_segment{:}); %Reshape cell array to colapse all the cells in the top level.  

end
