%Script just to plot Tracks from Trackmate with @msdanalyzer

% P. Aurelio, Matlab 2016a. 2018

%Code to analyze tracks from TrackMate (Fiji).
clear all
close all
clc

%% Some inputs:
%%----------------------------------------------------------
%Initialize Variables----------------------------------------
%--------------------------------------------------------------
dimension = 2; %Dimensionality of the movement
SPACE_UNITS = 'µm'; %This is just for visualization purposes
TIME_UNITS = 's'; %This is just for visualization purposes

directory_name = uigetdir('C:\Users\pgomez\Desktop\Temporal\FoxA1\JON SMT data\FOXA1 WT\TrackMate');

[list_names,directory_name] = uigetfile(fullfile(directory_name,'*.xml'),'select the files with SMT Tracks from TrackMate','MultiSelect','on');

 if iscell(list_names);
        strB = '';
    for g=1:size(list_names,2);
    [trayectory{g}, metadata{g}] = importTrackMateTracks(fullfile(strcat(directory_name,'\'),list_names{g}),'clipZ',1);
    strOut = sprintf('Loading and extracting the Data: % 4.1f',100*g/(size(list_names,2)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);
    end
      
    else
    [trayectory{1}, metadata] = importTrackMateTracks(fullfile(strcat(directory_name,'\'),list_names),'clipZ',1); 
  
 end
 
 
if iscell(list_names);
    
    for g=1:size(list_names,2);
        
ma = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);  %Initialize MSD analizer        
ma = ma.addAll(trayectory{g}');
 
figure()
ma.plotTracks;
ma.labelPlotTracks;
title(strcat('Tracks corresponding to_',g));
 
clear ma
    end
   
else
  
ma = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);  %Initialize MSD analizer        
ma = ma.addAll(trayectory{1});
 
figure()
ma.plotTracks;
ma.labelPlotTracks;
title(strcat('Tracks'));
 
end 
     
     
     
