%Code to analyze  tracks from TrackMate (Fiji).
% Pablo Aurelio Gomez Garcia, Matlab 2015a. 2018
clear all
close all
clc
warning off 

%% Some inputs:
%%----------------------------------------------------------
%Define starting directory and input units of the tracks----------------------------------------
%--------------------------------------------------------------
% directory_name 
directory_name = uigetdir('/Volumes');

TrackMate = 1; %If the track file is an .xml file from TrackMate.

Input_space_units = 'micron'; %Specify the input spatial units of the Tracks (With Trackmate I usually use microns) ('pixels' or 'micron')
Input_time_units = 'second'; %Specify the input time units of the Tracks ('frames' or 'second') (With TrackMate I usually use Frames)

%----------------------------------------------------------------
%Binary inputs------------------------------------------------------
%-----------------------------------------------------------------------
dimension = 2; %Dimensionality of the movement. 2 for 2D and 3 for 3D.

Save_Results= 1; %Binary: 1 for save the results.

Save_Figures = 1;
Results_Figure_Folder = '/Result_Figures/All cells';

Filtering = 1; % Binary: 1 for filtering the Tracks based on the specified parameters.

subROIs = 0; %1 for Analyzing subregions (This requires to input the corresponding subROIs from ImageJ). If 0, this analysis won't be performed.
Inside = 0; %Binary: 1 for Analyzing the Tracks INSIDE the subROIS. IF 0, the OUTSIDE Tracks will be analyse.
In_percentage = 0.9; %Normalized percentage of points of a trayectory inside a certain subROI in order to be considered an "INSIDE" track.

%----------------------------------------------------------------------------
%Some Numerical Inputs-------------------------------------------------------------
%-------------------------------------------------------------------------
minimum_track_length = 20;  %Minimum track length in frames to take the track into account

D_free = 0.5; %Estimated value for the Difusion Coeficient of a free diffusing molecule in the nucleus (µm^2/s). This is used to calculate the Bound error probability and to split the Tracks in two populations.
Jump_confined_threshold = 0.1; %Value used to separate confined trajectories into two groups [Average Jump (µm)]. 

Frame_interval = 0.008; %Exposure time in seconds
pixel_size = 0.145;  %Effective Pixel size in µm

TMSD_fitting_points = 3;   %Number of points used of the T-MSD for fitting the Diffussion Coefficient. (LINEAR FITTING) The minimum 3 points for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 3;  %Number of points used of the TE-MSD for fitting the Diffussion Coefficient. (LINEAR FITTING) The minimum 3 points for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; %Number of points used of the TE-MSD for fitting the Diffussion Coefficient. (POWER-LAW FITTING).

max_alpha_conf = 0.8; %Maximum alpha value for the power-fitting of the TE-MSD in order to consider Confined Motion
min_alpha_directed = 1.2; %Minimum alpha value for the power-fitting of the TE-MSD in order to consider Directed Motion
%Everything in the middle will be consider as Brownian Motion

R2LIMIT = 0.6; %Lower limit of the r-squared value for the fitting of the T-MSD and TE-MSD. The Tracks with r-squared lower than that won't be classified.


%------------------------------------------------------------------
%Circle Confined Difussion---------------------------------------------------
%------------------------------------------------------------------
level = 0.01; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in µm]. 
num_points = 8; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (µm^2/s) on the Confined Circle Diffusion Model
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (µm) on the Confined Circle Diffusion Model




%% Initialize Variables----------------------------------------
%%----------------------------------------------------------
%--------------------------------------------------------------
SPACE_UNITS = 'µm'; %This is just for visualization purposes
TIME_UNITS = 's'; %This is just for visualization purposes
ma = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);  %Initialize MSD analizer
ma_AllTracks = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS); 
ma_confined = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined_High_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined_Low_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_directed = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_brownian = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments_confined = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments_directed = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_High_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_Low_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);

MSD=[];
abs_displacement=[];
velocity=[];
total_displacement=[];
Residence=[];


%% Load the Data:
%-------------------------------------------------------------------------
%Tracks-------------------------------------------------------------
%--------------------------------------------------------------------------
d = dir(directory_name);
isub = [d(:).isdir]; 
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

if TrackMate == 1;
%Extract the trayectories and other related information from the results of the tracking  
[list_names,directory_name] = uigetfile(fullfile(directory_name,'*.xml'),'select the files with SMT Tracks from TrackMate','MultiSelect','on');
% [list_names2,directory_name] = uigetfile(fullfile(directory_name,'*.csv'),'select the files with Tracks Statistics from TrackMate','MultiSelect','on');
%%Start importing the TrackMate data.
    if iscell(list_names);
        strB = '';
    for g=1:size(list_names,2);
    [trayectory{g}, metadata{g}] = importTrackMateTracks(fullfile(directory_name,list_names{g}),'clipZ',1);
%     M{g} = readtable(fullfile(strcat(directory_name,'\'),list_names2{g}));
    strOut = sprintf('Loading and extracting the Data: % 4.1f',100*g/(size(list_names,2)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);
    end
   
    
    
    else
    [trayectory{1}, metadata] = importTrackMateTracks(fullfile(directory_name,list_names),'clipZ',1); 

    end
end




%% Convert spatial units to microns if input units were in pixels.    
if Input_space_units == 'pixels';
    for i = 1:size(trayectory,2);
        for j=1:size(trayectory{i},1)
            trayectory{i}{j}(:,2:3) = trayectory{i}{j}(:,2:3)*pixel_size;    
        end
    end
end
   
%From now on, all spatial units are in um.


%% Keep the cell index of each trajectory for cell by cell analysis   
    for i = 1:size(trayectory,2);
        for j=1:size(trayectory{i},1)
            Trayectories_CellReference.AllTracks{i}{j,1} = i;    
        end
    end


 %% Calculate the center, the limits & the area of each nucleus
 %This will performed by calculating the center of mass of all the
 %trajectories on the image (This works if there is only 1 cell per image) 
 for i = 1:size(trayectory,2);
     
     all_points_temp = []; 
     tri_area_cell{i} = [];
     
        for j=1:size(trayectory{i},1);
             all_points_temp = [all_points_temp ; trayectory{i}{j}(:,2:3)];    
        end
        
        Cell_centerOfMass_MaxDist_Area{i}(1:2) = mean(all_points_temp);  %units are in um.   
           
        
  [max_dist_cell{i}, min_dist_cell{i}, avg_dist_cell{i}] = distance_scatter(all_points_temp);
            
         tri = delaunay(all_points_temp(:,1), all_points_temp(:,2));

        for h= 1:size(tri,1);
        %Calculate Area of the triangles
        tri_area(h) = polyarea(all_points_temp(tri(h,:),1),all_points_temp(tri(h,:),2));
        end
        tri_area_cell{i} = sum(tri_area);
            
        Cell_centerOfMass_MaxDist_Area{i}(3) = max_dist_cell{i};
        Cell_centerOfMass_MaxDist_Area{i}(4) = tri_area_cell{i};
        
 end

 
    


%% Filter The Tracks.
%Track Length
if TrackMate == 1;
for ff=1:size(trayectory,2);
    if Filtering == 1;
     
    for u=1:size(trayectory{1,ff},1);   
    Track_length{ff}(u,1) = size(trayectory{1,ff}{u,1},1);   
    end
   
    indices{ff} = find(Track_length{ff}(:,1) >= (minimum_track_length));
    trayectory_filtered{1,ff} = trayectory{1,ff}([indices{ff}]);
    
   Trayectories_CellReference.AllTracksFiltered{1,ff} = Trayectories_CellReference.AllTracks{1,ff}([indices{ff}]);
    
    else
    indices{ff} = 1:1:size(trayectory{1,ff},1);
    trayectory_filtered = trayectory;
    
    Trayectories_CellReference.AllTracksFiltered = Trayectories_CellReference.AllTracks;

    end
end
end

%Filter tracks from SlimFast based on their tracklength
if SlimFast == 1;
    for ff=1:size(trayectory,2);
        if Filtering == 1;
    indices_length{ff} = find(tracklength{ff} >= (minimum_track_length));
    trayectory_filtered{1,ff} = trayectory{1,ff}([indices_length{ff}]);
    
    Trayectories_CellReference.AllTracksFiltered{1,ff} = Trayectories_CellReference.AllTracks{1,ff}([indices_length{ff}]);

        else
    trayectory_filtered = trayectory;  
    Trayectories_CellReference.AllTracksFiltered = Trayectories_CellReference.AllTracks;

    end   
    end
end


%Filter tracks from SimulatedTracks based on their tracklength
if SimulatedTracks == 1;
    for ff=1:size(trayectory,2);
        if Filtering == 1;
    indices_length{ff} = find(tracklength{ff} >= (minimum_track_length));
    trayectory_filtered{1,ff} = trayectory{1,ff}([indices_length{ff}]);
    
    Trayectories_CellReference.AllTracksFiltered{1,ff} = Trayectories_CellReference.AllTracks{1,ff}([indices_length{ff}]);

        else
    trayectory_filtered = trayectory;  
    Trayectories_CellReference.AllTracksFiltered = Trayectories_CellReference.AllTracks;

    end   
    end
end



%% -Load the subROIs from ImageJ corresponding to the Heterochromatin if needed (Just for subROI analysis)
if subROIs == 1;
    
    %Check how many Cells per image there are
    for j = 1:size(list_names,2);
    image_names(j) = str2num(list_names{j}(9:11));
    end
    
    [list_names3,directory_name2] = uigetfile(fullfile(directory_name,'*.zip;*.roi'),'select the files with the subROIs from ImageJ','MultiSelect','on');
    
    if iscell(list_names3);
        
    for j = 1:size(list_names3,2);
    subROI_names(j) = str2num(list_names3{j}(8:10));
    end
    
    idx_member=find(ismember(image_names',subROI_names','rows'));

    if size(image_names,2) == length(idx_member);
    else
        error('One subROI per image is required');
    end
    
    for g=1:size(list_names3,2);
    [sROI{g}] = ReadImageJROI(fullfile(strcat(directory_name2,'\'),list_names3{g}));    
    end   

    else
      sROI = ReadImageJROI(fullfile(strcat(directory_name2,'\'),list_names3));   
      subROI_names = str2num(list_names3(8:10));
       
    
    end
    
%Sort the trayectories inside and outside the heterochromatin regions%
[trayectories_INSIDE, trayectories_OUTSIDE,Trayectory_CellReference_INSIDE,Trayectory_CellReference_OUTSIDE,list_in] = tracksINROI(trayectory_filtered,sROI,pixel_size,image_names,subROI_names,In_percentage,Trayectory_CellReference_filtered);

if Inside == 1;
trayectories = cat(1,trayectories_INSIDE{:});  %CHOOSE ALL HERE
Trayectories_CellReference.AllTracksINSIDE = cat(1,Trayectory_CellReference_INSIDE{:});  %CHOOSE ALL HERE
Trayectories_CellReference.AllTracksINSIDE(find(cellfun(@isempty,Trayectories_CellReference.AllTracksINSIDE))) = []; %Always perform the same operation on the reference of the tracks
else
trayectories = cat(1,trayectories_OUTSIDE{:}); %CHOOSE ALL HERE
Trayectories_CellReference.AllTracksOUTSIDE = cat(1,Trayectory_CellReference_OUTSIDE{:});  %CHOOSE ALL HERE
Trayectories_CellReference.AllTracksOUTSIDE(find(cellfun(@isempty,Trayectories_CellReference.AllTracksOUTSIDE))) = []; %Always perform the same operation on the reference of the tracks
end
  
%Delete the empty rows from the Cell trayectories (just for
%heterochromatin...)
trayectories(find(cellfun(@isempty,trayectories))) = [];

%If the options subROIs is not selected:
else        
trayectories = cat(1,trayectory_filtered{:});  %CHOOSE ALL HERE
Trayectories_CellReference.AllTracksFiltered = cat(1,Trayectories_CellReference.AllTracksFiltered{:});  %CHOOSE ALL HERE
end



%% Multiply the frame by the exposure time (in TrackMate, I track them without specifying the frame interval)
if Input_time_units == 'frames'
for ff=1:size(trayectories,1);
    trayectories{ff,1}(:,1)=trayectories{ff,1}(:,1)*Frame_interval;
end

trayectories = trayectories';
else
trayectories = trayectories';
end
%From now on, trayectories are in um and seconds.




%% Data Analysis
fprintf('\n');     
fprintf('Calculating and plotting......\n');     


%% ------ Create a separate variable for All the Trajectories -----------
ma_AllTracks = ma_AllTracks.addAll(trayectories);
ma_AllTracks = ma_AllTracks.computeMSD;
ma_AllTracks = ma_AllTracks.LogTMSD(TLOGLOG_fitting_points);
ma_AllTracks = ma_AllTracks.TMSD(TMSD_fitting_points);

Trayectories_CellReference.AllTracks = Trayectories_CellReference.AllTracksFiltered;


%% ------------------------------------------------------------------------------
%Add the Trayectories to the msdAnalyzer-----------------------------------
%------------------------------------------------------------------------------   
 if isempty(trayectories)==1;
        error('There are no tracks that fulfill the requirements')
 else
     
ma = ma.addAll(trayectories);
Num_Tracks = size(trayectories,2);

for j=1:Num_Tracks;
Track_Length(j) = size(trayectories{j},1);
Residence_Times(j)=Track_Length(j)*Frame_interval;
end

end
    
%% -------------------------------------------------------------
%Calculate the distance travelled by the Tracks-------------------------
%(This value depends on the Track Length of each Track, but can be informative)--------
%------------------------------------------------------------------
for ggg = 1:size(trayectories,2);
points_coord = trayectories{ggg}(:,2:3);
[max_dist{ggg}, min_dist{ggg}, avg_dist{ggg}] = distance_scatter(points_coord);
end


%% Compute TE-MSD
ma = ma.computeMSD;


%% Preliminary Plotting
figure()
ma_AllTracks.plotTracks;
ma_AllTracks.labelPlotTracks;
set(gca,'Ydir','reverse');
title('All Trajectories')

%figure()
%ma.plotMeanMSD(gca, true)
%mmsd = ma.getMeanMSD;
%temps = mmsd(:,1);
% xs = mmsd(:,2);
% dx_plot = mmsd(:,3) ./ sqrt(mmsd(:,4));
% dxs = mmsd(:,3);

% figure()
% ha = errorbar(temps(2:1:end), xs(2:1:end), dx_plot(2:1:end),'o','MarkerSize',6,'MarkerFaceColor','b','Linewidth',1.4);
% xlim([0 xmax]);
% ylim([0 ymax]);
% 
% 
%   xlabel('Time (seconds)');
%   ylabel('MSD (µm^2)');
%   title('TE-MSD Curves');
%  
%     
% hold on
% drawnow
% errorbar_tick(ha,40);

% [fo, gof,D_lsqFit,lsq_conf_int] = ma.TEMSD(TEMSD_fitting_points);
% ma.labelPlotMSD;
% legend off
% 
% figure()
% hist(Track_Length,100);
% title('Histogram of the Length of the Tracks');
% xlim([0 100]);


%% Compute T-MSD
%This other approach fits every single MSD curve and then the histogram
%of D is very informative
ma = ma.TMSD(TMSD_fitting_points);
good_enough_fit_Ds = find(ma.lfit.r2fit >= R2LIMIT);
Dmean = mean( ma.lfit.a(good_enough_fit_Ds) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit_Ds) ) / 2 / ma.n_dim;
fprintf('**Estimation of the diffusion coefficient from linear fit of the MSD curves (Fitting every MSD curve)**:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, length(good_enough_fit_Ds));

%--------------------------------------------------------------
%CLassify the Tracks based on their Diffussion Coefficients--------
%---------------------------------------------------------------
Ds = ma.lfit.a(good_enough_fit_Ds)/ 2 / ma.n_dim;
idx4 = find(Ds >= D_free);
idx5 = find(Ds < D_free);
%D higher than Db------------------------------
temporal_tray = trayectories(good_enough_fit_Ds);
trayectories_High_D = temporal_tray(idx4);
reference_of_Tracks_High_D = good_enough_fit_Ds(idx4);
ma_High_D = ma_High_D.addAll(trayectories_High_D);
Trayectories_CellReference.High_D = Trayectories_CellReference.AllTracksFiltered(reference_of_Tracks_High_D);

trayectories_Low_D = temporal_tray(idx5);
reference_of_Tracks_Low_D = good_enough_fit_Ds(idx5);
ma_Low_D = ma_Low_D.addAll(trayectories_Low_D);
Trayectories_CellReference.Low_D = Trayectories_CellReference.AllTracksFiltered(reference_of_Tracks_Low_D);



if isempty(ma_High_D.tracks);   
else
% figure()
% ma_High_D.plotMSD;
% title('T-MSD of Tracks with D higher than Db');
% xlim([0 0.1]);
% ylim([0 ymax]);

% figure()
% ma_High_D.plotMeanMSD(gca, true);
% title('TE-MSD of Tracks with D higher than Db');
% xlim([0 0.1]);
% ylim([0 ymax]);

figure()
ma_High_D.plotTracks;
ma_High_D.labelPlotTracks;
title('Trajectories with D higher than D_free','Interpreter', 'none');
end


%% Plotting
figure()
%Take out negative values from the Diffusion Coefficients list
idx2=find(Ds > 0);
histogram(log10(Ds(idx2)));
title('Histogram of diffusion coefficients');
xlabel('Log10(Diffusion Coefficient)')

figure()
%Take out negative values from the Diffusion Coefficients list
idx2=find(Ds > 0);
histogram(Ds(idx2));
title('Histogram of diffusion coefficients');
xlabel('Diffusion Coefficients (um2/s)')
xlim([0,2])

%% Motion type analysis from T-MSD curves (fitting each Track MSD)
ma = ma.LogTMSD(TLOGLOG_fitting_points);

r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;

% Remove bad fits
bad_fits = r2fits < R2LIMIT;
good_enough_fit_alpha = find(r2fits >= R2LIMIT);
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas_filtered = alphas(good_enough_fit_alpha);
%Remove NaN from alphas
%alphas_filtered(find(isnan(alphas_filtered))) = [];

% T-test
[htest, pval] = ttest(alphas_filtered, 1, 0.05, 'left');

if ~htest
    [htest, pval] = ttest(alphas_filtered, 1, 0.05);
end

% Prepare string
str = { [ '\alpha = ' sprintf('%.2f ± %.2f (mean ± std, N = %d)', mean(alphas_filtered), std(alphas_filtered), numel(alphas_filtered)) ] };

if htest
    str{2} = sprintf('Significantly below 1, with p = %.2g', pval);
else
    str{2} = sprintf('Not significantly differend from 1, with p = %.2g', pval);
end


%--------------------------------------------------------------------------------------
%%Separate and plot Tracks based on their alpha value from LogLog fit of the T-MSD---------------
%-------------------------------------------------------------------------------------
temporal_tray = trayectories(good_enough_fit_alpha);
conf_filtered = find(alphas_filtered <= max_alpha_conf);
brownian_filtered = find(alphas_filtered>max_alpha_conf & alphas_filtered<min_alpha_directed);
directed_filtered = find(alphas_filtered >= min_alpha_directed);

reference_of_Tracks_conf = good_enough_fit_alpha(conf_filtered);
reference_of_Tracks_brownian = good_enough_fit_alpha(brownian_filtered);
reference_of_Tracks_directed = good_enough_fit_alpha(directed_filtered);
%Calculate the percentage of coincidence between High Diffusive Tracks and
%the types of motion:
conf_High_D_coincidence = sum(ismember(reference_of_Tracks_High_D,reference_of_Tracks_conf)==1);
brownian_High_D_coincidence = sum(ismember(reference_of_Tracks_High_D,reference_of_Tracks_brownian)==1);
directed_High_D_coincidence = sum(ismember(reference_of_Tracks_High_D,reference_of_Tracks_directed)==1);

High_D_confined_percentage = conf_High_D_coincidence/(conf_High_D_coincidence + brownian_High_D_coincidence + directed_High_D_coincidence);
High_D_brownian_percentage = brownian_High_D_coincidence/(conf_High_D_coincidence + brownian_High_D_coincidence + directed_High_D_coincidence);
High_D_directed_percentage = directed_High_D_coincidence/(conf_High_D_coincidence + brownian_High_D_coincidence + directed_High_D_coincidence);

%---------------------------------------------------------------------------------
%------------------------------------------------------------------------------
Conf_trayectories_filtered = temporal_tray(conf_filtered);
Directed_trayectories_filtered = temporal_tray(directed_filtered);
Brownian_trayectories_filtered = temporal_tray(brownian_filtered);

ma_confined = ma_confined.addAll(Conf_trayectories_filtered);
ma_directed = ma_directed.addAll(Directed_trayectories_filtered);
ma_brownian = ma_brownian.addAll(Brownian_trayectories_filtered);


%Calculate the number of Tracks per type of motion
Percentage_confined = size(ma_confined.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_brownian = size(ma_brownian.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_directed = size(ma_directed.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_butterfly = size(ma_butterfly.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_High_D = size(ma_High_D.tracks,1)/size(ma_AllTracks.tracks,1);


%Compute MSD for all the motion types in separate variables
if isempty(ma_confined.tracks);   
else
ma_confined = ma_confined.computeMSD;
ma_confined = ma_confined.LogTMSD(TLOGLOG_fitting_points);
ma_confined = ma_confined.TMSD(TMSD_fitting_points);

Trayectories_CellReference.Confined = Trayectories_CellReference.AllTracksFiltered(reference_of_Tracks_conf);

end

if isempty(ma_directed.tracks);   
else
ma_directed = ma_directed.computeMSD;
ma_directed = ma_directed.LogTMSD(TLOGLOG_fitting_points);
ma_directed = ma_directed.TMSD(TMSD_fitting_points);

Trayectories_CellReference.Directed = Trayectories_CellReference.AllTracksFiltered(reference_of_Tracks_directed);
end

if isempty(ma_brownian.tracks);   
else
ma_brownian= ma_brownian.computeMSD;
ma_brownian = ma_brownian.LogTMSD(TLOGLOG_fitting_points);
ma_brownian = ma_brownian.TMSD(TMSD_fitting_points);

Trayectories_CellReference.Brownian = Trayectories_CellReference.AllTracksFiltered(reference_of_Tracks_brownian);
end

if isempty(ma_High_D.tracks);  
else
ma_High_D = ma_High_D.computeMSD;
ma_High_D = ma_High_D.LogTMSD(TLOGLOG_fitting_points);
ma_High_D = ma_High_D.TMSD(TMSD_fitting_points);
end

if isempty(ma_Low_D.tracks);  
else
ma_Low_D = ma_Low_D.computeMSD;
ma_Low_D = ma_Low_D.LogTMSD(TLOGLOG_fitting_points);
ma_Low_D = ma_Low_D.TMSD(TMSD_fitting_points);
end

%Separate the Confined Tracks based on their Mean Jump.
for i=1:length(Conf_trayectories_filtered);  
    tracktemp = Conf_trayectories_filtered{i};
    jd_temp_confined_mean(i) = mean(sqrt((tracktemp(2:end,2) - tracktemp(1:end-1,2)).^2 + (tracktemp(2:end,3) - tracktemp(1:end-1,3)).^2));
end
idx6 = find(jd_temp_confined_mean >= Jump_confined_threshold);
idx7 = find(jd_temp_confined_mean < Jump_confined_threshold);

trayectories_confined_High_D = Conf_trayectories_filtered(idx6);
ma_confined_High_D = ma_confined_High_D.addAll(trayectories_confined_High_D);
Trayectories_CellReference.Confined_High_D = Trayectories_CellReference.Confined(idx6);

trayectories_confined_Low_D = Conf_trayectories_filtered(idx7);
ma_confined_Low_D = ma_confined_Low_D.addAll(trayectories_confined_Low_D);
Trayectories_CellReference.Confined_Low_D = Trayectories_CellReference.Confined(idx7);


if isempty(ma_confined_High_D.tracks);  
else
ma_confined_High_D = ma_confined_High_D.computeMSD;
ma_confined_High_D = ma_confined_High_D.LogTMSD(TLOGLOG_fitting_points);
ma_confined_High_D = ma_confined_High_D.TMSD(TMSD_fitting_points);


figure()
ma_confined_High_D.plotTracks;
ma_confined_High_D.labelPlotTracks;
title('Confined Trajectories with Average Jump higher than Jump_threshold','Interpreter', 'none');
set(gca,'Ydir','reverse');

end

if isempty(ma_confined_Low_D.tracks);  
else
ma_confined_Low_D = ma_confined_Low_D.computeMSD;
ma_confined_Low_D = ma_confined_Low_D.LogTMSD(TLOGLOG_fitting_points);
ma_confined_Low_D = ma_confined_Low_D.TMSD(TMSD_fitting_points);


figure()
ma_confined_Low_D.plotTracks;
ma_confined_Low_D.labelPlotTracks;
title('Confined Trajectories with Average Jump lower than Jump_threshold','Interpreter', 'none');
set(gca,'Ydir','reverse');


end

if isempty(ma_butterfly_segments.tracks);  
else
ma_butterfly_segments = ma_butterfly_segments.computeMSD;
ma_butterfly_segments = ma_butterfly_segments.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments = ma_butterfly_segments.TMSD(TMSD_fitting_points);
end




%% Plotting ---------------------------------------------------
%--------------------------------------------------------------------------------------
%--------------------------------------------
%-------------------------------------------------------------------------------------
% figure()
% ma_confined.plotMSD;
% title('T-MSD Confined Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

% figure()
% ma_confined.plotMeanMSD(gca, true);
% title('TE-MSD Confined Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

figure()
ma_confined.plotTracks;
ma_confined.labelPlotTracks;
title('Confined Tracks');
set(gca,'Ydir','reverse');

% figure()
% ma_directed.plotMSD;
% title('T-MSD Directed Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

% figure()
% ma_directed.plotMeanMSD(gca, true);
% title('TE-MSD Directed Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

figure()
ma_directed.plotTracks;
ma_directed.labelPlotTracks;
title('Directed Tracks');
set(gca,'Ydir','reverse');

% figure()
% ma_brownian.plotMSD;
% title('T-MSD Brownian Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

% figure()
% ma_brownian.plotMeanMSD(gca, true);
% title('TE-MSD Brownian Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);
% 
% figure()
% ma_brownian.plotTracks;
% ma_brownian.labelPlotTracks;
% title('Brownian Tracks');
% set(gca,'Ydir','reverse');


%% Save the results from the MSD Analysis
if Save_Results == 1;
fprintf('\n');     
fprintf('Saving...\n'); 
clear ma
ma = ma_AllTracks;
mkdir(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)));
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_AllTracks','.mat')),'ma');
clear ma;
ma = ma_brownian;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Brownian','.mat')),'ma');
clear ma;
ma = ma_confined;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Confined','.mat')),'ma');
clear ma;
ma = ma_directed;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Directed','.mat')),'ma');
clear ma
ma = ma_butterfly;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Butterfly','.mat')),'ma');
clear ma
ma = ma_High_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_High_D_',num2str(D_free),'.mat')),'ma');
clear ma
ma = ma_Low_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Low_D_',num2str(D_free),'.mat')),'ma');
clear ma
ma = ma_confined_High_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_confined_High_Jump_',num2str(Jump_confined_threshold),'.mat')),'ma');
clear ma
ma = ma_confined_Low_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_confined_Low_Jump_',num2str(Jump_confined_threshold),'.mat')),'ma');
clear ma
ma = ma_butterfly_segments_confined;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Butterfly_segments_confined','.mat')),'ma');
clear ma
ma = ma_butterfly_segments_directed;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('MSD_Results_',num2str(Frame_interval),'ms','_Butterfly_segments_directed','.mat')),'ma');

 
% CellReferenceSave = Trayectories_CellReference.AllTracks;
% save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_AllTracks','.mat')),'CellReferenceSave');
CellReference = Trayectories_CellReference.AllTracks;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_AllTracks','.mat')),'CellReference');
CellReference = Trayectories_CellReference.High_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_High_D','.mat')),'CellReference');
CellReference = Trayectories_CellReference.Low_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_Low_D','.mat')),'CellReference');
CellReference = Trayectories_CellReference.Confined;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_Confined','.mat')),'CellReference');
CellReference = Trayectories_CellReference.Directed;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_Directed','.mat')),'CellReference');
CellReference = Trayectories_CellReference.Confined_High_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_Confined_High_D','.mat')),'CellReference');
CellReference = Trayectories_CellReference.Confined_Low_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_Directed_Low_D','.mat')),'CellReference');

if isempty(ma_brownian.tracks);   
else
CellReference = Trayectories_CellReference.Brownian;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('CellReference_Brownian','.mat')),'CellReference');
end



save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('Cells_centerOfMass','.mat')),'Cell_centerOfMass_MaxDist_Area');

%save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(Bound_threshold(1)),'um','_',num2str(minimum_track_length),'min_frames'),strcat('radius_equiv_maxDist_',num2str(Frame_interval),'ms','.mat')),'radius_equiv_maxDist');

Motion_type_percentages{1,2} = 'Confined';
Motion_type_percentages{1,3} = 'Brownian';
Motion_type_percentages{1,4} = 'Directed';
Motion_type_percentages{1,5} = 'All High Diff Coeff';

Motion_type_percentages{2,2} = Percentage_confined;
Motion_type_percentages{2,3} = Percentage_brownian;
Motion_type_percentages{2,4} = Percentage_directed;
Motion_type_percentages{2,5} = Percentage_High_D;

save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames','_subROIs',num2str(subROIs),'_Inside',num2str(Inside)),strcat('Motion_type_percentages','.mat')),'Motion_type_percentages');

end

if Save_Figures == 1;
mkdir(strcat(directory_name,Results_Figure_Folder));
FolderName = strcat(directory_name,Results_Figure_Folder);
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList);
  FigHandle = FigList(iFig);
  FigName   = strcat('Figure_',num2str(iFig));
  savefig(FigHandle, fullfile(FolderName, [FigName, '.fig']));
  saveas(FigHandle, fullfile(FolderName, [FigName, '.tif']));

end

end















