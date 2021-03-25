#Written by Pablo Aurelio Gomez Garcia (PhD Student at ICFO) (January, 2020)
#Scrip for TrackMate batch tracking over multiple images and multiple subROIs

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.oldlap import SimpleLAPTrackerFactory
from fiji.plugin.trackmate.tracking.sparselap import SimpleSparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.action import ExportStatsToIJAction
from fiji.plugin.trackmate.action import ExportAllSpotsStatsAction
from fiji.plugin.trackmate.action import TrackBranchAnalysis
from fiji.plugin.trackmate.graph import GraphUtils
from ij.plugin import HyperStackConverter
from ij.measure import ResultsTable
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
from ij import IJ
from ij import ImagePlus
from ij.plugin.frame import RoiManager
from ij.gui import Roi
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import sys
import os
from java.io import File
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer

import fiji.plugin.trackmate.io.TmXmlReader as TmXmlReader
import fiji.plugin.trackmate.action.ExportTracksToXML as ExportTracksToXML
import fiji.plugin.trackmate.io.TmXmlWriter as TmXmlWriter
import fiji.plugin.trackmate.features.ModelFeatureUpdater as ModelFeatureUpdater
import fiji.plugin.trackmate.features.SpotFeatureCalculator as SpotFeatureCalculator
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzerFactory as SpotContrastAndSNRAnalyzerFactory
import fiji.plugin.trackmate.features.FeatureAnalyzer as FeatureAnalyzer
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzer as SpotContrastAndSNRAnalyzer
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory as SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzer as SpotIntensityAnalyzer
import fiji.plugin.trackmate.features.track.TrackSpeedStatisticsAnalyzer as TrackSpeedStatisticsAnalyzer
import fiji.plugin.trackmate.features.TrackFeatureCalculator as TrackFeatureCalculator
import fiji.plugin.trackmate.features.SpotFeatureCalculator as SpotFeatureCalculator
import fiji.plugin.trackmate.action.ExportAllSpotsStatsAction as ExportAllSpotsStatsAction
import fiji.plugin.trackmate.features.track.TrackSpeedStatisticsAnalyzer as TrackSpeedStatisticsAnalyzer
import fiji.plugin.trackmate.features.edges.EdgeTargetAnalyzer as EdgeTargetAnalyzer
import fiji.plugin.trackmate.features.edges.EdgeTimeLocationAnalyzer as EdgeTimeLocationAnalyzer
import fiji.plugin.trackmate.util.TMUtils as TMUtils


#................................................................................................................#
#......................................................INPUTS....................................................#
#................................................................................................................#

#Directory
Folder = '//freenas01.linux.crg.es/pcosma/pgomez/180201_PerfectDay_H2B_500ms/Live/ESC Clone2/560nm/';
OutFolfer = '//freenas01.linux.crg.es/pcosma/pgomez/180201_PerfectDay_H2B_500ms/Live/ESC Clone2/560nm/Int.150 Jump.400 Close.200 Gap.2/';


extension = '.tif'
Folder_ROIs = Folder + 'ROIs/';

#Experimental values
#Pixel_size = 0.160; #Space units are in microns.
#Frame_interval = 0.015; #Time units are in seconds.
Pixel_size = 0.160; #Space units are in microns.
Frame_interval = 0.5; #Time units are in seconds.

Initial_frames = 0; #Exclude this number of initital frames from the stacks.
Final_frames = 0; #Exclude this number of final frames from the stacks.

Close_images = 1; #Binary decision. 1 for closing the images after tracking. This saves RAM memory.
Analyze_tracks = 1; #Binary decision. 1 for inclusing the feature analysis of the tracks to the module.
Export_spot_statistics = 1; #Binary decision. 1 for exporting the spot statistics table as a .csv file. **If you want to export the spot statistics, Analyze_tracks variable must be equal to 1**

#Localization
Subpixel_localization = True;
Median_filtering = True;
Spot_radius = 0.4; #Units are in microns
Intensity_threshold = 150.0;

#Tracking
Track_splitting = False;
Track_merging = False;
Gap_closing = True;
Max_jump = 0.4; #Max jump from frame to frame. Units are in microns.
Gap_closing_distance = 0.2; #Maximum spatial distance when closing a GAP. Units are in microns.
Max_frame_gap = 2; #Maximum number of frame to close the GAP.

#................................................................................................................#
#................................................................................................................#


#................................................................................................................#
#......................................................PROCESSING................................................#
#................................................................................................................#
ROI_idx = 0;
for filename in os.listdir(Folder):
	if filename.endswith(extension):
	
		print(filename)
	
		imp = IJ.openImage(Folder + filename)
		imp.show()
	
		Calib = imp.getCalibration();
	 	Calib.frameInterval = Frame_interval;
		Calib.pixelDepth = Pixel_size;
		Calib.pixelWidth = Pixel_size;
		Calib.pixelHeight = Pixel_size;
		Calib.unit = "microns";	
		imp.setCalibration(Calib);
	
		num_frames = imp.getNFrames();
		num_z_stacks = imp.getNSlices();
		
		if num_frames < num_z_stacks:
			imp.setDimensions(1,num_frames,num_z_stacks);
				
		    
		    		       
		ROI_list = os.listdir(Folder_ROIs);

		Roi_File = Folder_ROIs + ROI_list[ROI_idx];
		print(ROI_list[ROI_idx]);
		rm = RoiManager();
		rm = RoiManager.getInstance();
		rm.runCommand("Open", Roi_File);
		total_rois = rm.getCount();

		for idx in range(total_rois):
		
			#----------------------------
			# Create the model object now
			#----------------------------
			    
			# Some of the parameters we configure below need to have
			# a reference to the model at creation. So we create an
			# empty model now.	    
			model = Model();		    
			# Send all messages to ImageJ log window.
			model.setLogger(Logger.IJ_LOGGER);
			model.setPhysicalUnits(Calib.getUnit(), Calib.getTimeUnit());

			#------------------------
			# Prepare settings object
			#------------------------	       
			settings = Settings();
					
			rm.select(imp,idx);
			actual_roi = rm.getRoi(idx);
			imp.setRoi(actual_roi);
			settings.setFrom(imp);
			settings.tstart = Initial_frames;
			settings.tend = imp.getNFrames() - Final_frames -1;
			#settings.roi = rm.getSelectedRoisAsArray();
			       
			# Configure detector - We use the Strings for the keys
			settings.detectorFactory = LogDetectorFactory();
			settings.detectorSettings = { 
			    'DO_SUBPIXEL_LOCALIZATION' : Subpixel_localization,
			    'RADIUS' : Spot_radius,
			    'TARGET_CHANNEL' : 1,
			    'THRESHOLD' : Intensity_threshold,
			    'DO_MEDIAN_FILTERING' : Median_filtering,
			};  
			    
			# Configure spot filters - Classical filter on quality
			#filter1 = FeatureFilter('QUALITY', 0, True);
			#settings.addSpotFilter(filter1);
			     
			# Configure tracker - We want to allow merges and fusions
			#settings.trackerFactory = SimpleLAPTrackerFactory();
			settings.trackerFactory = SimpleSparseLAPTrackerFactory();
			settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap(); # almost good enough
			settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = Track_splitting;
			settings.trackerSettings['ALLOW_TRACK_MERGING'] = Track_merging;
			settings.trackerSettings['LINKING_MAX_DISTANCE'] = Max_jump;
			settings.trackerSettings['ALLOW_GAP_CLOSING'] = Gap_closing;
			settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = Gap_closing_distance;
			settings.trackerSettings['MAX_FRAME_GAP'] = Max_frame_gap;
			    
		
			#----------------
			# Configure track features analyzers
			#----------------
			if Analyze_tracks:
				settings.addEdgeAnalyzer(EdgeTargetAnalyzer());
				settings.addEdgeAnalyzer(EdgeTimeLocationAnalyzer());
				settings.addSpotAnalyzerFactory(SpotIntensityAnalyzerFactory());
				settings.addSpotAnalyzerFactory(SpotContrastAndSNRAnalyzerFactory());
				settings.addTrackAnalyzer(TrackSpeedStatisticsAnalyzer());
				settings.addTrackAnalyzer(TrackDurationAnalyzer());
			    
			#----------------
			# Configure track filters
			#----------------		
			# the bottom right of the image. Track displacement must be above 10 pixels.		    
			#filter2 = FeatureFilter('TRACK_DISPLACEMENT', 0, True);
			#settings.addTrackFilter(filter2);
			    
			    
			#-------------------
			# Instantiate plugin
			#-------------------
			    
			trackmate = TrackMate(model, settings);
			       
			#--------
			# Process
			#--------
			    
			trackmate.checkInput()
			print("Input checked")
			#if not ok:
			    #sys.exit(str(trackmate.getErrorMessage()))
			    
			trackmate.process()
			print("Process finished")
			#if not ok:
			    #sys.exit(str(trackmate.getErrorMessage()))
			    
			       
			#----------------
			# Display results
			#----------------
			     
			selectionModel = SelectionModel(model)
			displayer =  HyperStackDisplayer(model, selectionModel, imp)
			displayer.render()
			displayer.refresh()
			
			# Echo results with the logger we set at start:
			model.getLogger().log(str(model))


			#----------------
			# export needed files for further analysis
			#----------------
			filename_only = os.path.splitext(filename)[0];
			out_filename = filename_only + '_Roi' + str(idx+1) + '_Tracks.xml';		
			outFile = File(OutFolfer, out_filename);
				
			#out_filename_stats = filename_only + '_Roi' + str(idx+1) + '_Track_Statistics.csv';
			out_filename_spots = filename_only + '_Roi' + str(idx+1) + '_Spots_Statistics.csv';
			#outFile_stats = File(OutFolfer, out_filename_stats);
			#outFile_spots = File(OutFolfer, out_filename_spots);
			outFile_spots = OutFolfer + out_filename_spots;	
			
			ExportTracksToXML.export(model,settings,outFile);

			if Export_spot_statistics:
				ExportStatsToIJAction(selectionModel).execute(trackmate);
				IJ.selectWindow("Track statistics");
				IJ.run("Close");
				IJ.selectWindow("Links in tracks statistics");
				IJ.run("Close");
				IJ.selectWindow("Spots in tracks statistics");		
				IJ.saveAs("Results",outFile_spots);
				IJ.run("Close");
			
		rm.runCommand("Select All");
		rm.runCommand("delete");
	
	 	if Close_images:
	 		imp.close();
	    
	   	ROI_idx = ROI_idx + 1;
	    
		continue
		
	else:
		continue



