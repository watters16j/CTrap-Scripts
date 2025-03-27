# -*- coding: utf-8 -*-
"""
This script is designed to take in .tdms and .h5 kymograph file inputs and perform batch processing of foci
for downstream analysis process. Lines are tracked using two of Lumicks' line tracking algorithm: information 
about which can be found at - https://lumicks-pylake.readthedocs.io/en/stable/tutorial/kymotracking.html 

This script automatically exports the lines tracked for each of the colors requested (time data in second,
coordinate data in nm) as well as the kymotracking settings for each call of the kymotracker in a .xlsx file.
Metadata from the files is stored in a separate .csv file.

There are built in options to extract the photon counts around the foci (as determined by the line_width
variable) and the distance between a chosen foci color and the other color's foci.
All additional data is stored in the .xlsx file.

See keyword arguments and tracking algorithm parameters for more information.

---------------------------------------------------
Optional keyword arguments to skip the input steps:
To input keywords when you call the script input them as a space separated list of variable=value
Do not include apostrophes as shown in these definitions
Ex: tracking_method=greedy line_width=4, colors_to_track=RG

tracking_method : greedy or lines (as a string)
    This variable determines which of the lumicks tracking methods the script uses. The "track_greedy" method
    or the "track_lines" method
    
file_name : string
    This variable defines the base file to save the line data/kymotracker setting for file_name_summary.xlsx 
    and the metadata file as file_name + _metadata.csv
    
colors_to_track : string of any combination of R and/or G and/or B
    This variable defines which channels are going to be put through the kymotracker algorithm
    
opt_to_extract_intensities : string of "yes" or "no"
    This variable allows the user to control whether or not the script also outputs the sum of the photon counts
    around each foci. The range over which the sum is measured is based off of the line_width variable
        num_frames = math.ceil(line_width/2)

opt_to_extract_distance_between_foci: string of "yes" or "no"
    The variable goes through each of the "base" lines (as defined by colors_to_track_distance) and finds the closest point
    in each of the other color traces (ex: the closest red line to a green base line) and reports the distance.
    Negative distance means the tracked line is below the baseline, and positive distance means the tracked line is above the baseline

opt_for_area_selection : integer, 1 or 2
    This variable determine the method used to select the area fed into the kymotracker algorithm.
    Method 1 just involves choosing a point above and below the region of interest - and the algorithm uses a rectangle with those as the vertical points
    Method 2 involves choosing the top left then top right point and then every point of the bottom area where you want to define.
        The bottom line will be a horizontal line until the first clicked point, then will track linearly to the next point and so on until the end of the kymograph
        From the last point, the last defined slope will will continue until the end of the kymograph

line_width : integer
    This variable sets the default line_width to whatever integer value you set

color_to_track_distance : string of "R" "G" or "B"
    Defines the base line for distance extract.

--------------------------------------------------
Parameters (description taken from Lumicks Script):
track_greedy:
    line_width : float
        Expected line width in pixels.
    pixel_threshold : float
        Intensity threshold for the pixels. Local maxima above this intensity level will be designated as a line
        origin.
    window : int
        Number of kymograph lines in which the particle is allowed to disappear (and still be part of the same line).
    sigma : float or None
        Uncertainty in the particle position. This parameter will determine whether a peak in the next frame will be
        linked to this one. Increasing this value will make the algorithm tend to allow more positional variation in
        the lines. If none, the algorithm will use half the line width.
    vel : float
        Expected velocity of the traces in the image. This can be used for non-static particles that are expected to
        move at an expected rate (default: 0.0).
    diffusion : float
        Expected diffusion constant (default: 0.0). This parameter will influence whether a peak in the next frame
        will be connected to this one. Increasing this value will make the algorithm allow more positional variation
        in.
    sigma_cutoff : float
        Sets at how many standard deviations from the expected trajectory a particle no longer belongs to this trace.
        Lower values result in traces being more stringent in terms of continuing (default: 2.0).
    filter_line_length : int (this one is normally tied into a different function in Lumicks' scripts)
        Sets a threshold for a filter where any lines shorter than this threshold will be removed from the data set.

track_lines:
    line_width:float
        Expected line width in pixels.
    max_lines:int
        Maximum number of lines to trace.
    start_threshold:float
        Threshold for the value of the derivative.
    continuation_threshold:float
        Derivative threshold for the continuation of a line.
    angle_weight: float
        Factor which determines how the angle between normals needs to be weighted relative to distance. High values push for straighter lines. Weighting occurs according to distance + angle_weight * angle difference
    filter_line_length : int (this one is normally tied into a different function in Lumicks' scripts)
        Sets a threshold for a filter where any lines shorter than this threshold will be removed from the data set.
"""

import lumicks.pylake as lk
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from nptdms import TdmsFile as td
import math
import glob
import pandas as pd
from tkinter import filedialog
from tkinter import *
import sys
import os

def extract_lines_data(filepath,dict_kymotracking_method_storage, color_list, h5_kymo_object=""):        
    string_size_inside_loop = 74
    def extract_TDMS_channel_data(channelData):
        count = 0
        photonCountArray = np.empty([pixelsPerLine,numLines])
        for lineNum in range(0,numLines):
            for pixelNum in range(0,pixelsPerLine):
                photonCountArray[pixelNum,lineNum] = channelData[count]
                count += 1
        return photonCountArray
    
    def get_correct_input_data(photon_count_array,rescaled_photon_count_array):                
        fig, ax = plt.subplots(nrows=1,ncols=1)

        ax.imshow(rescaled_photon_count_array,aspect="auto")    
        
        ax.set_title(f'{filepath}',weight='bold',size=16)
        plt.tight_layout()
        
        analysis_area_mapping=np.zeros(photon_count_array.shape)
        
        if opt_for_area_selection == 1:
            print('\n' + "#"*string_size_inside_loop)
            print('Input area selections of\n[1] Top Point of Region of Interest\n[2] Bottom Point of Region of Interest')
            print("#"*string_size_inside_loop + '\n')
            points = plt.ginput(2,timeout=45)
            plt.close()
            
            top_bead_data = math.floor(points[0][1])
            bottom_bead_data = math.ceil(points[1][1])
            
            analysis_area_mapping=np.zeros(photon_count_array.shape)
        
            analysis_area_mapping[top_bead_data,:]=1
            analysis_area_mapping[bottom_bead_data,:]=1
        # option to select a region above
        elif opt_for_area_selection == 2:
            print('\n' + "#"*string_size_inside_loop)
            print('Assumed Location of the stationary bead is at the top\nInput area selections of\n[1] Below top bead near the start of the kymograph\n[2]Below top bead near the end of the kymograph')
            print('[3] Bottom left point above the bead\n[4+] Any point where the bottom bead trajectory is changed\n[Last] Bottom right point above the bottom bead')
            print("#"*string_size_inside_loop + '\n')
            points = plt.ginput(-1,timeout=20)
            plt.close()
            top_bead_top_first = math.ceil(points[0][1])
            top_bead_top_end = math.ceil(points[1][1])
            
            top_bead_data = min(top_bead_top_end,top_bead_top_first)
            
            length_of_list = len(points)
            bottom_bead_list = []
            for i in range(2,length_of_list):
                bottom_bead_list.append(points[i])
            
            analysis_area_mapping[max(top_bead_top_first,top_bead_top_end),:] = 1 # take maximum input of the 
            
            analysis_area_mapping[min(math.floor(bottom_bead_list[0][1]),math.floor(bottom_bead_list[1][1])),0:math.ceil(bottom_bead_list[1][0])] = 1
            
            previous_time_val = math.ceil(bottom_bead_list[1][0])
            previous_position_val = math.floor(bottom_bead_list[1][1])
            
            array_for_bottom_bead_filter = np.zeros(analysis_area_mapping.shape[1])
            
            array_for_bottom_bead_filter[:math.ceil(bottom_bead_list[1][1])] = min(math.floor(bottom_bead_list[0][1]),math.floor(bottom_bead_list[1][1]))
            
            for numSteps in range(2,len(bottom_bead_list)):
                current_time_val = math.ceil(bottom_bead_list[numSteps][0])
                current_position_val = math.floor(bottom_bead_list[numSteps][1])
                
                slope = (current_position_val - previous_position_val) / (current_time_val - previous_time_val)
                for i in range(previous_time_val,current_time_val):
                    analysis_area_mapping[math.floor(previous_position_val+(i-previous_time_val)*slope),i] = 1
                    array_for_bottom_bead_filter[i] = math.floor(previous_position_val+(i-previous_time_val)*slope)
                
                previous_time_val = current_time_val
                previous_position_val = current_position_val
            
            #mapping last slope line
            for i in range(previous_time_val,numLines):
                analysis_area_mapping[math.floor(previous_position_val+(i-previous_time_val)*slope),i] = 1
                array_for_bottom_bead_filter[i] = math.floor(previous_position_val+(i-previous_time_val)*slope)

        fig, ax = plt.subplots(nrows=1,ncols=1)
        ax.imshow(rescaled_photon_count_array,alpha=0.9,aspect="auto")
        ax.imshow(analysis_area_mapping,cmap="binary",alpha = 0.2,aspect="auto")
        ax.set_title(f"Mapping Area of Analysis\n{filepath}",weight='bold',size=16)
        ax.axis('off')
        plt.tight_layout()
        plt.show()
        
        if opt_for_area_selection == 1:
            return top_bead_data, bottom_bead_data
        elif opt_for_area_selection == 2:
            return top_bead_data, array_for_bottom_bead_filter
    
    def filter_image_data(photon_count_array,top_bead_data,bottom_bead_filter):
        if opt_for_area_selection == 1:
            area_of_analysis = photon_count_array[top_bead_data:bottom_bead_filter,:]
        else:
            max_distance_to_analyze = np.amax(bottom_bead_filter)
            bottom_bead_filter = bottom_bead_filter - top_bead_data
            
            area_of_analysis = photon_count_array[top_bead_data:int(max_distance_to_analyze),:]
            
            
            for i in range(0,photon_count_array.shape[1]):
                area_of_analysis[int(bottom_bead_filter[i]):,i] = 0

        return area_of_analysis
    
    def call_kymotracker(photon_count_area_of_analysis,kymotracker_dict_values,color_lines_tracked):
        # call the kymotracker prototype to find the lines
        happy_with_kymotracking = "no"
        
        max_of_photon_count_area = np.amax(photon_count_area_of_analysis)
        rescaled_area_of_analysis = photon_count_area_of_analysis * (255 / max_of_photon_count_area)
        rescaled_area_of_analysis = rescaled_area_of_analysis.astype(int)
        
        while happy_with_kymotracking != "yes":
            if tracking_method == 1:
                lines_tracked = lk.track_greedy(photon_count_area_of_analysis,
                                                line_width=kymotracker_dict_values["line_width"],
                                                pixel_threshold=kymotracker_dict_values["pixel_threshold"],
                                                window = kymotracker_dict_values["window"],
                                                sigma = kymotracker_dict_values["sigma"],
                                                vel = kymotracker_dict_values["vel"],
                                                diffusion = kymotracker_dict_values["diffusion"],
                                                sigma_cutoff = kymotracker_dict_values["sigma_cutoff"])
            else:
                lines_tracked = lk.track_lines(photon_count_area_of_analysis,    
                                               line_width = kymotracker_dict_values["line_width"],
                                               max_lines = kymotracker_dict_values["max_lines"],
                                               start_threshold = kymotracker_dict_values["start_threshold"],
                                               continuation_threshold = kymotracker_dict_values["continuation_threshold"],
                                               angle_weight = kymotracker_dict_values["angle_weight"])
                
            filtered_tracked_lines = lk.filter_lines(lines_tracked,kymotracker_dict_values["filter_line_length"])
            
            # add refine lines centroid method once it is added to lumicks
            
            # take the extracted lines and add it to a makeshift plot for quality control purposes
            lines_of_photon_counts = np.zeros(photon_count_area_of_analysis.shape)

            print(f"{len(filtered_tracked_lines)} lines tracked")
        
            # map the tracked lines in a pixelated manner
            for line in filtered_tracked_lines:
                time_vals = line.time_idx
                coordinate_vals = line.coordinate_idx
                len_time_vals = len(time_vals)
                for index in range(0,len_time_vals):
                    if tracking_method ==1:
                        lines_of_photon_counts[int(coordinate_vals[index]),time_vals[index]] = 1 #not subpixel mapping but it will be close
                    elif tracking_method == 2:
                        lines_of_photon_counts[int(coordinate_vals[index]),int(time_vals[index])] = 1 #not subpixel mapping but it will be close
            # format the plots to assess quality of kymotracking
            matplotlib.pyplot.close('all')
            
            fig, ax = plt.subplots(nrows=2,ncols=1)
            ax[0].imshow(photon_count_area_of_analysis,cmap=color_lines_tracked,aspect="auto")
            ax[1].imshow(lines_of_photon_counts,cmap="binary",aspect="auto")
            ax[0].get_xaxis().set_ticks([])
            ax[0].get_yaxis().set_ticks([])
            ax[1].get_xaxis().set_ticks([])
            ax[1].get_yaxis().set_ticks([])
            ax[0].set_title(f"{color_lines_tracked[:-1]} Area of Analysis",weight='bold',size=16)
            ax[1].set_title("Mapped Lines",weight='bold',size=16)
            ### need to refine this math to make it look better - calling plt.tight_layout() ten times
            plt.tight_layout()
            
            date_header = filepath.split(" ")[0] + " "
            title_of_figure = filepath.strip(date_header)
            title_of_figure = title_of_figure.split(".tdms")[0] + "_"
            fig.savefig(title_of_figure+color_lines_tracked[:-1]+"_Channel.tiff",bbox_inches="tight", pad_inches = 0)
            plt.show()
            
            happy_with_kymotracking = input("Are you happy with the kymotracking?\n").lower()
            
            if happy_with_kymotracking != "yes":
                print('\n'+"#"*string_size_inside_loop)
                print('Would you like to change the parameters of the line fitting?')
                print('Write a comma separated list of the new values you would like to re-define')
                print("Ex: 'bead_line_width=8,line_width=8")
                print("Previous values are displayed below")
                print("#"*string_size_inside_loop)
                for key in kymotracker_dict_values:
                    print(f"{key} --- {kymotracker_dict_values[key]}")
                print("#"*string_size_inside_loop+'\n')
                string_of_variables_to_change = input("Input comma separated list of variables to change:\n")
                print("#"*string_size_inside_loop)
                
                string_of_variables_to_change = string_of_variables_to_change.split(",")
                try:
                    for value in string_of_variables_to_change:
                        value = value.split("=")
                        if "-" in value[0] or "+" in value[0] or " " in value[0]: # handling errors that deal with hitting - or + instead of the equals sign
                            print("Non-text value (+,-, ) was found in the list of variables to change - was skipped")
                            continue
                        
                        if value[0] in kymotracker_dict_values.keys():
                            try:
                                kymotracker_dict_values[value[0].strip(" ")] = int(value[1])
                            except:
                                try:
                                    kymotracker_dict_values[value[0]] = float(value[1])
                                except:
                                    kymotracker_dict_values[value[0]] = value[1]
                except:
                    print("No values to change were found")
                    pass
                
        return filtered_tracked_lines, lines_of_photon_counts, kymotracker_dict_values, kymotracker_dict_values["line_width"]
    
    # extract necessary data for analysis
    if ".tdms" in filepath:
        #open the data file            
        tdms_file = td.open(filepath)
        data = tdms_file['Data']
        time_data_ms = data['Time (ms)'][:]
        
        totalNumPixels = 0 #placeholder to reduce redundant calculations if multiple channels are being examined

        #extract necessary infromation from the metadata
        metadata = td.read_metadata(filepath)
        metadataDict = metadata.properties
        pixelsPerLine = metadataDict["Pixels per line"]
        pixel_size_nm = int(float(metadataDict["Scan Command.Scan Command.scanning_axes.0.pix_size_nm"]))    
    
        #extract data from each channel that is desired           
        pixel_ch1 = data['Pixel ch 1'][:]
        pixel_ch2 = data['Pixel ch 2'][:]
        pixel_ch3 = data['Pixel ch 3'][:]
        totalNumPixels = pixel_ch2.shape[0]
    
        numLines = int(totalNumPixels / pixelsPerLine)

        delta_line_time = max(time_data_ms) / (1000 * numLines)

        pixel_size_um = pixel_size_nm / 1000 # define pixel size for future distance measurements
    
        red_photon_counts = extract_TDMS_channel_data(pixel_ch1)
        green_photon_counts = extract_TDMS_channel_data(pixel_ch2)
        blue_photon_counts = extract_TDMS_channel_data(pixel_ch3)
    elif ".h5" in filepath:
        h5_file_object = lk.File(filepath) # build in an option to do this for each kymograph file? - do this before entering this function or after? before is probably easiest
        kymo_object = h5_file_object.kymos[h5_kymo_object]
        
        #extract photon count data
        red_photon_counts = kymo_object.red_image
        green_photon_counts = kymo_object.blue_image
        blue_photon_counts = kymo_object.green_image
        
        pixel_size_um = kymo_object.pixelsize_um
        pixel_size_nm = pixel_size_um * 1000
        
        delta_line_time = (kymo_object.timestamps[0,1]-kymo_object.timestamps[0,0]) / 1000000000 
        
        numLines = red_photon_counts.shape[1]        
        metadataDict = h5_file_object.description
    else:
        print('The filepath being analyzed was neither a .h5 or .tdms file.')
        exit()

    #define the specific kymotracker method from lumicks
    kymotracker_dict_values = {}    
    
    # set defaults by changing these variable declarations
    if tracking_method == 1: #track greedy
        kymotracker_dict_values["line_width"] = 5
        kymotracker_dict_values["pixel_threshold"] = 1
        kymotracker_dict_values["window"] = 8
        kymotracker_dict_values["sigma"] = None
        kymotracker_dict_values["vel"] = 0.0
        kymotracker_dict_values["diffusion"] = 0
        kymotracker_dict_values["sigma_cutoff"] = 1.0
        kymotracker_dict_values["filter_line_length"] = 100
    else: #track lines
        kymotracker_dict_values["line_width"] = 5
        kymotracker_dict_values["max_lines"] = 10
        kymotracker_dict_values["start_threshold"] = 0.005 
        kymotracker_dict_values["continuation_threshold"] = 0.005
        kymotracker_dict_values["angle_weight"] = 10.0
        kymotracker_dict_values["filter_line_length"] = 50
    
    if def_line_width != 0:
        kymotracker_dict_values["line_width"] = def_line_width
        
    #loop until user is happy with the area selection - add logical selection for complicated systems       
    rgb_image = np.dstack((red_photon_counts,green_photon_counts,blue_photon_counts))
    
    max_photon_count_value = np.amax(rgb_image)  
    rgb_image_modified = rgb_image * (255 / max_photon_count_value)
    rgb_image = rgb_image.astype(int)
    rgb_image_modified = rgb_image_modified.astype(int)
    
    #partition the specific areas
    logical_for_inputs = "no"
    while logical_for_inputs.lower() != "yes":
        top_bead_data, bottom_bead_filter = get_correct_input_data(rgb_image,rgb_image_modified)
        logical_for_inputs=input("Are you happy with the area selections?\n")
    
    
    """
    Go through the iterative kymotracker calling until the user is happy with the line tracing method for all of the colors desired 
    """
    if "R" in color_list:
        red_area_of_analysis = filter_image_data(red_photon_counts,top_bead_data,bottom_bead_filter)
        red_lines, red_line_array, kymotracker_dict_values, red_line_width = call_kymotracker(red_area_of_analysis, kymotracker_dict_values, color_lines_tracked = "Reds")
        
        for key in kymotracker_dict_values:
            dict_kymotracking_method_storage[key].append(kymotracker_dict_values[key])
        dict_kymotracking_method_storage["color_tracked_list"].append("Red")
        
    if "G" in color_list:
        green_area_of_analysis = filter_image_data(green_photon_counts,top_bead_data,bottom_bead_filter)
        green_lines, green_line_array, kymotracker_dict_values, green_line_width = call_kymotracker(green_area_of_analysis, kymotracker_dict_values, color_lines_tracked = "Greens") 
        
        for key in kymotracker_dict_values:
            dict_kymotracking_method_storage[key].append(kymotracker_dict_values[key])
        dict_kymotracking_method_storage["color_tracked_list"].append("Green")
    
    if "B" in color_list:
        blue_area_of_analysis = filter_image_data(blue_photon_counts,top_bead_data,bottom_bead_filter)
        blue_lines, blue_line_array, kymotracker_dict_values, blue_line_width = call_kymotracker(blue_area_of_analysis,kymotracker_dict_values, color_lines_tracked = "Blues")
            
        for key in kymotracker_dict_values:
            dict_kymotracking_method_storage[key].append(kymotracker_dict_values[key])
        dict_kymotracking_method_storage["color_tracked_list"].append("Blue")
    
    """
    Output basic line data for analysis
    Option to output the sum of photon counts across a region as defined by the line_width
    """
    dict_for_traces = {}  
    if "R" in color_list:
        red_trace_count=1  
        if opt_to_extract_intensities == "yes":
            red_num_pixels = math.ceil(red_line_width / 2)
        
        for red_line in red_lines:
            red_line_time_vals = red_line.time_idx
            red_line_coordinate_vals = red_line.coordinate_idx
            dict_for_traces['Red Trace ' + str(red_trace_count) + " Time (s)"] = np.asarray(red_line_time_vals) * delta_line_time
            dict_for_traces['Red Trace ' + str(red_trace_count) + " Position (nm)"] = (np.asarray(red_line_coordinate_vals)+ top_bead_data) * pixel_size_nm
            if opt_to_extract_intensities == "yes":
               summed_intensity_values = red_line.sample_from_image(num_pixels = red_num_pixels)
               dict_for_traces['Red Trace ' + str(red_trace_count) + " Summed Photon Counts"] = np.asarray(summed_intensity_values)
               
            red_trace_count += 1
    
    if "G" in color_list:
        green_trace_count=1  
        if opt_to_extract_intensities == "yes":
            green_num_pixels = math.ceil(green_line_width / 2)
        
        for green_line in green_lines:
            green_line_time_vals = green_line.time_idx
            green_line_coordinate_vals = green_line.coordinate_idx
            dict_for_traces['Green Trace ' + str(green_trace_count) + " Time (s)"] = np.asarray(green_line_time_vals) * delta_line_time
            dict_for_traces['Green Trace ' + str(green_trace_count) + " Position (nm)"] = (np.asarray(green_line_coordinate_vals) + top_bead_data) * pixel_size_nm
            green_trace_count += 1
            if opt_to_extract_intensities == "yes":
               summed_intensity_values = green_line.sample_from_image(num_pixels = green_num_pixels)
               dict_for_traces['Green Trace ' + str(green_trace_count) + " Summed Photon Counts"] = np.asarray(summed_intensity_values)

    if "B" in color_list:
        blue_trace_count=1  
        
        if opt_to_extract_intensities == "yes":
            blue_num_pixels = math.ceil(blue_line_width/ 2)
        
        for blue_line in blue_lines:
            blue_line_time_vals = blue_line.time_idx
            blue_line_coordinate_vals = blue_line.coordinate_idx
            dict_for_traces['Blue Trace ' + str(blue_trace_count) + " Time (s)"] = np.asarray(blue_line_time_vals) * delta_line_time
            dict_for_traces['Blue Trace ' + str(blue_trace_count) + " Position (nm)"] = (np.asarray(blue_line_coordinate_vals) + top_bead_data) * pixel_size_nm
            if opt_to_extract_intensities == "yes":
               summed_intensity_values = blue_line.sample_from_image(num_pixels = blue_num_pixels)
               dict_for_traces['Blue Trace ' + str(blue_trace_count) + " Summed Photon Counts"] = np.asarray(summed_intensity_values)
            blue_trace_count += 1
    
    # export closest line data based on options
    if opt_to_extract_distance_between_foci == "yes" and len(color_list)>1:
        def filter_closest_line_to_base_line(single_base_line_time, single_base_line_coord, length_of_single_base_line, temp_line_object):
            array_for_time_values = []
            array_for_distance_values = []
        
            for line_time_index in range(0,length_of_single_base_line):
                temp_base_time = single_base_line_time[line_time_index]
                temp_base_coordinate = single_base_line_coord[line_time_index]
                best_distance = base_area_of_analysis.shape[0]
                best_distance_no_count = base_area_of_analysis.shape[0]
                
                for temp_line in temp_line_object:
                    temp_line_time_vals = temp_line.time_idx
                    temp_line_coordinate_vals = temp_line.coordinate_idx
                    try:
                        index_of_temp_point = temp_line_time_vals.index(int(temp_base_time))
                        temp_coordinate = temp_line_coordinate_vals[index_of_temp_point]
                        temp_distance = (temp_base_coordinate - temp_coordinate)
                        if abs(temp_distance) < abs(best_distance):
                            best_distance=temp_distance
                    except:
                        pass
                    
                if best_distance != best_distance_no_count:
                    array_for_time_values.append(temp_base_time * delta_line_time)
                    array_for_distance_values.append(best_distance * pixel_size_nm)
        
            #convert the values to one amenable to pandas dataframe
            array_for_time_values = np.asarray(array_for_time_values)
            array_for_distance_values = np.asarray(array_for_distance_values)
            return array_for_time_values, array_for_distance_values
        
        # define base line
        removed_base_line_color = color_list   
        if "R" in color_to_track_distance:
            base_lines = red_lines
            base_area_of_analysis = red_area_of_analysis
            removed_base_line_color = removed_base_line_color.replace("R","")
        elif "G" in color_to_track_distance:
            base_lines = green_lines
            base_area_of_analysis = green_area_of_analysis
            removed_base_line_color = removed_base_line_color.replace("G","")
        else:
            base_lines = blue_lines
            base_area_of_analysis = blue_area_of_analysis
            removed_base_line_color = removed_base_line_color.replace("B","")
        
        base_line_trace_count = 0
        for base_line in base_lines:
            base_line_trace_count += 1
            base_time_vals = base_line.time_idx
            base_coordinate_vals = base_line.coordinate_idx
            length_base_time_vals = len(base_time_vals)
            
            if "R" in removed_base_line_color:
                time_array_values, distance_array_values = filter_closest_line_to_base_line(base_time_vals,base_coordinate_vals,length_base_time_vals,red_lines)
                dict_for_traces[f"Closest Distance Base[{color_to_track_distance}] Trace "+str(base_line_trace_count)+" Red Time (s)"] = time_array_values
                dict_for_traces[f"Closest Distance Base[{color_to_track_distance}] Trace "+str(base_line_trace_count)+" Red Distance (nm)"] = distance_array_values
            
            if "G" in removed_base_line_color:
                time_array_values, distance_array_values = filter_closest_line_to_base_line(base_time_vals,base_coordinate_vals,length_base_time_vals,green_lines)
                dict_for_traces["Closest Distance Base[{color_to_track_distance}] Trace "+str(base_line_trace_count)+" Green Time (s)"] = time_array_values
                dict_for_traces["Closest Distance Base[{color_to_track_distance}] Trace "+str(base_line_trace_count)+" Green Distance (nm)"] = distance_array_values
        
            if "B" in removed_base_line_color:
                time_array_values, distance_array_values = filter_closest_line_to_base_line(base_time_vals,base_coordinate_vals,length_base_time_vals,blue_lines)
                dict_for_traces["Closest Distance Base[{color_to_track_distance}] Trace "+str(base_line_trace_count)+" Blue Time (s)"] = time_array_values
                dict_for_traces["Closest Distance Base[{color_to_track_distance}] Trace "+str(base_line_trace_count)+" Blue Distance (nm)"] = distance_array_values
    elif opt_to_extract_distance_between_foci == "yes":
        print('Option to extract closest distance between foci was stated but only one option was in the color list')
        
    # send data to the spreadsheet
    pd_data_frame_for_filepath = pd.DataFrame.from_dict(dict_for_traces,orient='index')
    pd_data_frame_for_filepath = pd_data_frame_for_filepath.transpose()
    
    date_header = filepath.split(" ")[0] + " "
    name_for_sheet = filepath.strip(date_header)
    name_for_sheet = name_for_sheet.split(".tdms")[0]
    
    if len(name_for_sheet) >=31:
        name_for_sheet = name_for_sheet[-30:]
        print(f"File name was too long to title an excel sheet name - shortened to {name_for_sheet}")
    
    pd_data_frame_for_filepath.to_excel(writer,sheet_name=name_for_sheet,index=False,header=True)
    
    return dict_kymotracking_method_storage, metadataDict

# call tkinter dialog box to let the user navigate to the desired folder
root = Tk()
root.withdraw()
folder_selected = filedialog.askdirectory()
os.chdir(folder_selected)

# collect candidate files
filelist_1 = glob.glob("*.tdms")
filelist_2 = glob.glob("*.h5")
filelist = filelist_1 + filelist_2

#pre-define variables
tracking_method = 0
file_name = 0
colors_to_track = 0
opt_to_extract_intensities = 0
opt_to_extract_distance_between_foci = 0
opt_for_area_selection = 0
def_line_width =0
color_to_track_distance = 0

#add an option to manually define the answers
if len(sys.argv) > 1:
    for string_input in sys.argv:
        if "tracking_method" in string_input:
            if "greedy" in string_input:
                tracking_method=1
            if "lines" in string_input:
                tracking_method=2
        if "line_width" in string_input:
            def_line_width = int((string_input.split("="))[-1])
        if "file_name" in string_input:
            file_name = (string_input.split("="))[-1]
        if "colors_to_track" in string_input:
            colors_to_track = (string_input.split("="))[-1].upper()
        if "opt_to_extract_intensities" in string_input:
            opt_to_extract_intensities = (string_input.split("="))[-1].lower()
            if opt_to_extract_intensities != "yes" and opt_to_extract_intensities != "no":
                opt_to_extract_intensities=0
        if "opt_to_extract_distance_between_foci" in string_input:
            opt_to_extract_distance_between_foci = (string_input.split("="))[-1].lower()
            if opt_to_extract_distance_between_foci != "yes" and opt_to_extract_distance_between_foci != "no":
                opt_to_extract_distance_between_foci=0
        if "opt_for_area_selection" in string_input:
            opt_for_area_selection = int((string_input.split("="))[-1])
            if opt_for_area_selection != 1 and opt_for_area_selection != 2:
                opt_to_extract_distance_between_foci=0
        if "color_to_track_distance" in string_input:
            color_to_track_distance = (string_input.split("="))[-1].upper()
            if color_to_track_distance != "R" and color_to_track_distance != "G" and color_to_track_distance != "B":
                color_to_track_distance = 0
    
# add option to name summary files - printing out the file names to make sure user is happy with what is in the folder
max_filepath_length = len(max(filelist,key=len))
max_separator_string = 78
if max_filepath_length > max_separator_string:
    max_separator_string = max_filepath_length


print('-'*max_separator_string)
print('List of files to be analyzed')
print('#'*max_separator_string)
for file in filelist:
    print(file)
print('#'*max_separator_string)
if file_name == 0:
    file_name = input("Please sort the files so that only similar experiments are in the same folder:\nInput file name to save (do not add any file extension)?\n")
print('-'*max_separator_string)

if tracking_method == 0:
    tracking_method = int(input("\n" + '-'*max_separator_string + "\nPlease indicate desired tracking method from the lumicks.pylake options\n[1] track_greedy\n[2] track_lines\nInput integer number of the correct method:\n"))

#allow for mistakes in the user input
if tracking_method != 1 and tracking_method != 2:
    while tracking_method != 1 and tracking_method !=2:
        tracking_method = int(input("\n" + '#'*max_separator_string + "\nPlease re-indicate desired tracking method\n[1] track_greedy\n[2] track_lines\n\nInput integer number of the correct method\nDo not include brackets in the input\n" + '#'*max_separator_string + "\n"))

#call correct dictionary
dict_kymotracking_method_storage = {}
if tracking_method == 1:
    dict_kymotracking_method_storage["line_width"] = []
    dict_kymotracking_method_storage["pixel_threshold"] = []  
    dict_kymotracking_method_storage["window"] = []  
    dict_kymotracking_method_storage["sigma"] = []  
    dict_kymotracking_method_storage["vel"] = []  
    dict_kymotracking_method_storage["diffusion"] = []  
    dict_kymotracking_method_storage["sigma_cutoff"] = []
    dict_kymotracking_method_storage["filter_line_length"] = []  
    dict_kymotracking_method_storage["color_tracked_list"] = []
elif tracking_method == 2:
    dict_kymotracking_method_storage["line_width"] = []
    dict_kymotracking_method_storage["max_lines"] = []  
    dict_kymotracking_method_storage["start_threshold"] = []  
    dict_kymotracking_method_storage["continuation_threshold"] = []  
    dict_kymotracking_method_storage["angle_weight"] = []  
    dict_kymotracking_method_storage["filter_line_length"] = []  
    dict_kymotracking_method_storage["color_tracked_list"] = []
else: #exit script in a controlled manner if user does not input a correct number
    print("Correct input was not detected in the tracking method input\nPlease input the number without brackets next time\nEnding Program")
    exit()

# get user inputs
if colors_to_track == 0:
    print('-'*max_separator_string)
    colors_to_track = input("\n" + '-'*max_separator_string + "\nChoose colors to track for all files\nInput RGB values, Ex: RG or RGB:\n").upper()
    print('-'*max_separator_string+"\n")
if opt_to_extract_intensities == 0:
    print('-'*max_separator_string)
    opt_to_extract_intensities = input("Would you like to extract the photon counts sum of the lines?\nInput yes/no:\n").lower()
    print('-'*max_separator_string+"\n")

if len(colors_to_track) > 1:
    if opt_to_extract_distance_between_foci == 0:
        print('-'*max_separator_string)
        opt_to_extract_distance_between_foci = input("Would you like to extract the distance between tracked lines?\nInput yes/no:\n").lower()
        if opt_to_extract_distance_between_foci == "yes":
            color_to_track_distance = input("What color would you like to use as the base for extracting this distance?\nInput R/G/B\n")
        print('-'*max_separator_string+"\n")
    elif opt_to_extract_distance_between_foci == "yes" and color_to_track_distance == 0:
        print('-'*max_separator_string)
        color_to_track_distance = input("What color would you like to use as the base for extracting this distance?\nInput R/G/B\n")
        print('-'*max_separator_string+"\n")
    
if opt_for_area_selection == 0:
    print('-'*max_separator_string)
    opt_for_area_selection = int(input("Choose option of how to manually input the area of analysis:\n[1] Manually define top and bottom positions of a rectangle to analyze\n[2] Manually define a more complex region (containing pulls and relaxes)\nInput integer number of the correct method\n"))
    print('-'*max_separator_string+"\n")

writer= pd.ExcelWriter(file_name+"_summary.xlsx", engine = "xlsxwriter")
output_file = open(file_name + "_metadata_doc.csv","w")

#write the metadata dictionary in a separate file
output_file.write("Notes:\n")
if tracking_method == 1:
    output_file.write("Lumicks' track_greedy alogrithim is used to track lines in the trace and extract different data types\n")
elif tracking_method == 2:
    output_file.write("Lumicks' track_lines alogrithim is used to track lines in the trace and extract different data types\n")
output_file.write(f"Metadata for {file_name}_summary.xlsx - all traces:\n\n")

tdms_count = 1
h5_count = 1
for filepath in filelist:
    if ".tdms" in filepath:        
        #call the line extraction function
        dict_kymotracking_method_storage, metadataDict = extract_lines_data(filepath,dict_kymotracking_method_storage, color_list=colors_to_track)
        
        #write out metadata information to see difference in file type
        if tdms_count == 1: # write headers
            tdms_count += 1
            output_file.write(",")
            # write column headings
            for key in metadataDict:
                output_file.write(f"{str(key)},")
            output_file.write("\n")
            for key in metadataDict:
                temp_string_to_write = str(metadataDict[key]).replace('\n',' ')
                output_file.write(temp_string_to_write)
                output_file.write(",")
            output_file.write("\n")
        else:
            output_file.write(f"{filepath},")
            for key in metadataDict:
                temp_string_to_write = str(metadataDict[key]).replace('\n',' ')
                output_file.write(temp_string_to_write)
                output_file.write(",")
            output_file.write("\n")
    
    elif ".h5" in filepath:
        h5_file_object = lk.File(filepath)
        kymo_obj_list = list(h5_file_object.kymos)
        
        for kymo_obj in kymo_obj_list:
            dict_kymotracking_method_storage, metadataDict = extract_lines_data(filepath,dict_kymotracking_method_storage, color_list=colors_to_track,h5_kymo_object=kymo_obj)
            
            metadataString = metadataDict.replace("\n",",")            
            if h5_count == 1: # write headers
                h5_count +=1    
                output_file.write("Filepath,Kymograph Pointer,Description by Channel,")
                output_file.write("\n")
                output_file.write(f"{filepath},")
                output_file.write(f"{kymo_obj},")
                output_file.write(f"{metadataString},")
                output_file.write("\n")
            else: # write metadata
                output_file.write(f"{filepath},")
                output_file.write(f"{kymo_obj},")
                output_file.write(f"{metadataString},")
                output_file.write("\n")

# convert the dictionary of kymotracker settings used for this folder as the last sheet in the summary document
pd_data_frame = pd.DataFrame.from_dict(dict_kymotracking_method_storage)
pd_data_frame.to_excel(writer,sheet_name="Kymotracker Settings",index=False,header=True)

# properly save and close both the .xlsx summary document and the metadata .csv file
writer.save()
output_file.close()