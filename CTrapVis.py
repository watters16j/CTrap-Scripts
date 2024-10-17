"""
#################################################################################################
BSD 2-Clause License
Copyright (c) 2021, John Watters
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################################


This script is the second version of a GUI designed to dynamically look through
.h5 files from a Lumicks C-Trap instrument. Users are able to choose what plots
they want (Trap Position vs. Time, Force vs. Distance, Force vs. Time, RGB image, 
or a combination of RGB and one of the force/trap position options) with options to customize 
the plot to add a title, scale RGB image values, shift axis of all available axis. 

If both plots are dependent on time the x-axes are linked by the manual entries,
but the individual plots can be zoomed in/out of using the interactive toolbar.
The save name for both the image and the metadata files can be customized by the 
"File Name to Save" entry and the "Image Format" option. The "Draw Plot" button
pulls in all of the GUI information and draws/plots the correct plot back on 
the interface. The "Quit" button destroyed the tkinter GUI and quits the python execution. 
The "Export Force" button pulls a separate Tkinter window to allow the user
to export desired force data to a .csv file for further analysis. The "Export
Image No Axis" button exports the scan image without axis/extra space being 
plotted in the GUI window as .png files. The "Export ImageJ Montage" button 
exports scan image in such a way that the axis are shown as well as a timestamp
of the image for easy transformation into an ImageJ montage. The "Save GUI Image"
re-draws the plot and takes a screenshot of the GUI image and saves it as a 
file type and name inputted in the "File Name" selection. The "Open KymoTracker"
button opens a new window used to track lines and extract data using the algorithms
in the pylake library.

Tested in Python 3.8.6
Modules Used: numpy, lumicks.pylake, matplotlib, tkinter, pandas, glob, os, tifffile
If any of the modules are not avaliable run "pip install moduleName" in terminal
or "conda install moduleName"

Author: John Watters
* Force Extract Button Code was taken from unpublished code by Dr. Michael Wasserman
"""

# import necessary modules
import numpy as np
import lumicks.pylake as lk
import matplotlib
matplotlib.use("Agg")
import os
import glob
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import pandas as pd
import math
import tifffile as tiff

# test version info of lumicks.pylake to make sure the image reconstruction bug is avoided
if int((lk.__version__).split(".")[1]) < 8:
    raise ValueError("Please update your lumicks.pylake package (with the command 'pip install --upgrade lumicks.pylake') to update to at least version 0.8.1 to fix an image reconstruction bug/have KymoTracking functionalities")


class KymoTrackerGUI():
    def __init__(self,kt_master,filepath,typePointer,RGB_Data,mod_RGB_Data):
        self.kt_master = kt_master
        print("\nTo select a region to analyze, either drag a rectangle on the RGB image or hit the 'Define Region of Interest'\nand follow those instructions.\n")
        kymoPointer = lk.File(filepath).kymos[typePointer]
        red_channel_data = RGB_Data[:,:,0]
        green_channel_data = RGB_Data[:,:,1]
        blue_channel_data = RGB_Data[:,:,2]
        
        # set variables for max and min positions of a rectangle drag
        num_timestamps = RGB_Data[:,:,0].shape[1]
        pixel_line_length = RGB_Data[:,:,0].shape[0]

        dx = kymoPointer.pixelsize_um[0] * 1000 #pixel size in nm
        dt = kymoPointer.line_time_seconds #scan time in s <-- this has been different in previous programs
        print(dt)
        
        fig, ax = plt.subplots(1,1,constrained_layout=True)
        fig.set_dpi(130)
        #,extent=[0, maxTrueTime, 0, maxTrueDist]
        ax.imshow(mod_RGB_Data, aspect="auto")
        ax.axis('off')
        
        frameForKTCanvas = tk.ttk.Frame(kt_master,relief=tk.FLAT)
        frameForKTCanvas.grid(row=0,rowspan=10,column=0,columnspan=1,sticky="nw",padx=0,pady=0)
            
        canvas=FigureCanvasTkAgg(fig,master=frameForKTCanvas)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvas.draw()
        
        toolbar = NavigationToolbar2Tk(canvas, frameForKTCanvas)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        """
        The next two functions are used in drawing the rectangle as you draw it on the plot
        """
        def get_rect_dimensions(eclick,erelease):
            global basic_area
            x, y = eclick.inaxes.transData.inverted().transform((eclick.x, eclick.y))
            x2, y2 = erelease.inaxes.transData.inverted().transform((erelease.x, erelease.y))
            
            min_x = int(min( x , x2 ))
            min_y = int(min( y , y2))
            
            delta_x = int(abs(x2 - x))
            delta_y = int(abs(y2 - y))
            
            basic_area = [ [min_x , min_x + delta_x] , [ min_y , min_y + delta_y] ]
            return
        
        def draw_temp_rectangle(event):
            return
        
        """
        This function does the work of actually tracking the line objects.
        The offset terms are used to define the region of interest/plot the correct
        position and time values.
        """
        def call_track_lines(event):
            global offset_x
            global offset_y
            global filtered_red_lines
            global filtered_green_lines
            global filtered_blue_lines
            
            filtered_red_lines = ""
            filtered_green_lines = ""
            filtered_blue_lines = ""
            
            """
            Call lumicks.pylake tracking algorithm
            """
            def track_lines_one_color(one_channel_data,tracking_style):
                if tracking_style == "Greedy":
                    lines_tracked = lk.track_greedy(one_channel_data,
                                                line_width=int(entryLineWidthGreedy.get()),
                                                pixel_threshold=int(entryPixelThresholdGreedy.get()),
                                                window = int(entryWindow.get()),
                                                sigma = float(entrySigma.get()),
                                                vel = float(entryVel.get()),
                                                diffusion = float(entryDiffusion.get()),
                                                sigma_cutoff = float(entrySigmaCutoff.get()))
                else:
                    lines_tracked = lk.track_lines(one_channel_data,    
                                               line_width = int(entryLineWidthLines.get()),
                                               max_lines = int(entryMaxLines.get()),
                                               start_threshold = float(entryStartThreshold.get()),
                                               continuation_threshold = float(entryContinuationThreshold.get()),
                                               angle_weight = float(entryAngleWeight.get()))
                
                filtered_tracked_lines = lk.filter_lines(lines_tracked,int(entryLineLenGreedy.get()))
                    
                return filtered_tracked_lines
            
            """
            If the user wants to use the custom area selection, the image data is filtered
            through this function.
            """
            def filter_custom_area(colorArray,listOfCoords):
                global custom_x_max
                offset_x = listOfCoords[0][1]
                offset_y = 0 #for the custom area method the window is defined through the whole kymograph
                
                list_of_x_vals = [item[1] for item in listOfCoords[1:]]
                custom_x_max = max(list_of_x_vals)
                
                previous_time_val = listOfCoords[1][0]
                previous_position_val = listOfCoords[1][1]
            
                filtered_color_array = colorArray[offset_x:custom_x_max,:]
                filtered_color_array[previous_position_val:,0:previous_time_val] = 0
                
                for numSteps in range(2,len(listOfCoords)):
                    current_time_val = math.ceil(listOfCoords[numSteps][0])
                    current_position_val = math.floor(listOfCoords[numSteps][1])

                    slope = (current_position_val - previous_position_val) / (current_time_val - previous_time_val)

                    for i in range(previous_time_val,current_time_val):
                        filtered_color_array[math.floor(previous_position_val+(i-previous_time_val)*slope):,i] = 0
                        
                    previous_time_val = current_time_val
                    previous_position_val = current_position_val

                return filtered_color_array, offset_x, offset_y
            
            """
            If the user wants to use the basic area selection, the image data is filtered
            through this function.
            """
            def filter_basic_area(colorArray,basicAreaCoords):
                offset_x = basicAreaCoords[1][0]
                offset_y = basicAreaCoords[0][0]
                
                filtered_color_array = colorArray[offset_x:basicAreaCoords[1][1],offset_y:basicAreaCoords[0][1]]
                return filtered_color_array, offset_x, offset_y
            
            """
            This function plots the tracked lines either on the image plot or on 
            the additional plot.
            """
            def plot_tracked_lines(line_obj,string_for_color,offset_x,offset_y):
                for line in line_obj:
                    time_vals = line.time_idx
                    time_vals = [y+offset_y for y in time_vals]
                    coordinate_vals = line.coordinate_idx
                    coordinate_vals = [x+offset_x for x in coordinate_vals]
                    
                    axForTraces.plot(time_vals,coordinate_vals,linewidth=1,color=string_for_color)
                    
                return
            
            """
            Filtering the area of analysis using either the custom defintion (requires you to enter
            the TopLevel feature and click on the correct places and have the complexAreaOption selected)
            or by using the last drawn rectangle in this window.
            """
            if complexAreaOption.state() == ('selected',):
                try:
                    filtered_red_channel_data , offset_x, offset_y = filter_custom_area(red_channel_data,custom_area_pointers)
                    filtered_green_channel_data , offset_x, offset_y = filter_custom_area(green_channel_data,custom_area_pointers)
                    filtered_blue_channel_data , offset_x, offset_y = filter_custom_area(blue_channel_data,custom_area_pointers)
                except:
                    try:
                        print("\nCustom area has not been defined yet, defaulting to tracking the basic_area parameter.")
                        filtered_red_channel_data , offset_x, offset_y = filter_basic_area(red_channel_data,basic_area)
                        filtered_green_channel_data , offset_x, offset_y = filter_basic_area(green_channel_data,basic_area)
                        filtered_blue_channel_data , offset_x, offset_y = filter_basic_area(blue_channel_data,basic_area)
                    except:
                        print("Both potential area options (custom using click method or dragging a rectangle) have not been defined. Defaulting to tracking whole image.")
                        filtered_red_channel_data = red_channel_data
                        filtered_green_channel_data = green_channel_data
                        filtered_blue_channel_data = blue_channel_data
                        offset_x=0
                        offset_y=0
            else:
                try:
                    filtered_red_channel_data , offset_x, offset_y= filter_basic_area(red_channel_data,basic_area)
                    filtered_green_channel_data , offset_x, offset_y= filter_basic_area(green_channel_data,basic_area)
                    filtered_blue_channel_data , offset_x, offset_y= filter_basic_area(blue_channel_data,basic_area)
                except:
                        print("No previous rectangle dragged to define the region of interest. Defaulting to tracking whole image.")
                        filtered_red_channel_data = red_channel_data
                        filtered_green_channel_data = green_channel_data
                        filtered_blue_channel_data = blue_channel_data
                        offset_x=0
                        offset_y=0
            
            #now generate the plot
            if separatePlotOpt.state() == ('selected',):
                kt_fig, (axRGB,axForTraces) = plt.subplots(nrows=2,ncols=1,constrained_layout=True,sharex=True,sharey=True)
                axRGB.imshow(mod_RGB_Data, aspect="auto")
                axRGB.axis('off')
                axForTraces.axis('off')
            else:
                kt_fig, axForTraces = plt.subplots(nrows=1,ncols=1)
                axForTraces.imshow(mod_RGB_Data, aspect="auto")
                axForTraces.axis('off')
            
            kt_fig.set_dpi(130)
            
            tracking_method = comboboxMethod.get()
            
            if redLinesVar.state() == ('selected',):
                filtered_red_lines = track_lines_one_color(filtered_red_channel_data, tracking_method)
                if refineLinesOptCB.state() == ('selected',):
                    try:
                        filtered_red_lines = lk.refine_lines_centroid(filtered_red_lines,line_width=int(entryLineWidthLines.get()))
                    except:
                        print("No red lines were tracked")
                        
                plot_tracked_lines(filtered_red_lines,"red", offset_x, offset_y)
                
            if greenLinesVar.state() == ('selected',):
                filtered_green_lines = track_lines_one_color(filtered_green_channel_data, tracking_method)
                if refineLinesOptCB.state() == ('selected',):
                    try:
                        filtered_green_lines = lk.refine_lines_centroid(filtered_green_lines,line_width=int(entryLineWidthLines.get()))
                    except:
                        print("No green lines were tracked")
                plot_tracked_lines(filtered_green_lines,"green", offset_x, offset_y)
                
            if blueLinesVar.state() == ('selected',):
                filtered_blue_lines = track_lines_one_color(filtered_blue_channel_data, tracking_method)
                if refineLinesOptCB.state() == ('selected',):
                    try:
                        filtered_blue_lines = lk.refine_lines_centroid(filtered_blue_lines,line_width=int(entryLineWidthLines.get()))
                    except:
                        print("No blue lines were tracked")
                plot_tracked_lines(filtered_blue_lines,"blue", offset_x, offset_y)
            
            #plot the areas
            if showRegionOpt.state() == ('selected',):
                if complexAreaOption.state() != ('selected',):
                    delta_x = basic_area[0][1]-basic_area[0][0]
                    delta_y = basic_area[1][1]-basic_area[1][0]
            
                    rect_basic_area = matplotlib.patches.Rectangle((offset_y,offset_x),delta_x,delta_y,linewidth=1,edgecolor='gray',facecolor='none')
                    axForTraces.add_patch(rect_basic_area)
                else:
                    if len(custom_area_pointers) == 2:
                        axForTraces.axhline(offset_x,color="gray",linewidth=1)
                        axForTraces.axhline(custom_x_max,color="gray",linewidth=1)
                    else:
                        axForTraces.plot([0 , red_channel_data.shape[1]],[ offset_x,offset_x ],color="gray",linewidth=1)
                
                        previous_time_val = custom_area_pointers[1][0]
                        previous_position_val = custom_area_pointers[1][1]
                         
                        axForTraces.plot([0,previous_time_val],[previous_position_val,previous_position_val],color="gray",linewidth=1)
                        
                        for numSteps in range(2,len(custom_area_pointers)):
                            current_time_val = math.ceil(custom_area_pointers[numSteps][0])
                            current_position_val = math.floor(custom_area_pointers[numSteps][1])
                            
                            axForTraces.plot([previous_time_val,current_time_val],[previous_position_val,current_position_val],color="gray",linewidth=1)
                            
                            slope_for_last = (current_position_val - previous_position_val) / (current_time_val - previous_time_val)
                            
                            previous_time_val = current_time_val
                            previous_position_val = current_position_val
                        
                        step_til_end = red_channel_data.shape[1]-previous_time_val
                        axForTraces.plot([previous_time_val,red_channel_data.shape[1]] , [previous_position_val,previous_position_val+step_til_end*slope_for_last] ,color="gray",linewidth=1)
                        
            frameForKTCanvas = tk.ttk.Frame(kt_master,relief=tk.FLAT)
            frameForKTCanvas.grid(row=0,rowspan=10,column=0,columnspan=1,sticky="nw",padx=0,pady=0)
            
            canvas=FigureCanvasTkAgg(kt_fig,master=frameForKTCanvas)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvas.draw()
        
            toolbar = NavigationToolbar2Tk(canvas, frameForKTCanvas)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            
            #add the rectangle selector functionality again
            if separatePlotOpt.state() == ('selected',):
                draw_temp_rectangle.RS = matplotlib.widgets.RectangleSelector(axRGB, get_rect_dimensions, drawtype='box',rectprops=rectproperties)
                kt_fig.canvas.mpl_connect('button_press_event', draw_temp_rectangle)
            else:
                draw_temp_rectangle.RS = matplotlib.widgets.RectangleSelector(axForTraces, get_rect_dimensions, drawtype='box',rectprops=rectproperties)
                kt_fig.canvas.mpl_connect('button_press_event', draw_temp_rectangle)
            
            return
        
        """
        This function allows for the custom defintion of a region of interest,
        which is useful if the area you are looking at contains pulling/relaxing
        regions of the kymograph and you need to exclude the bead area.
        """
        def define_area_of_analysis(event):
            global custom_area_pointers
            custom_area_pointers = []            
            
            if complexAreaOption.state() != ("selected",):
                complexAreaOption.invoke()
            
            def escape_top_level(event):
                customAreaMaster.destroy()
                return
                
            def return_coordinates(event):
                if len(custom_area_pointers) < 2:
                    print('Only one point selected')
                else:
                    print('Multiple points selected')
                customAreaMaster.destroy()
                return
            
            def add_point(event):
                x, y = event.inaxes.transData.inverted().transform((event.x , event.y))
                x = int(x)
                y = int(y)
                custom_area_pointers.append([x,y])
                return
            
            customAreaMaster = tk.Toplevel()
            customAreaMaster.title(f"Choose custom area for: {typePointer}")
            if complexAreaOption.state() != ('selected',):
                print('For the custom area you define here to take effect make sure to click the "Custom Area Selection?" option!')
        
            fig_CAM, ax_CAM = plt.subplots(1,1,constrained_layout=True)
            ax_CAM.imshow(mod_RGB_Data,aspect="auto")
            ax_CAM.axis('off')

            frameForCAM = tk.ttk.Frame(customAreaMaster,relief=tk.FLAT)
            frameForCAM.grid(row=0,rowspan=10,column=0,columnspan=1,sticky="nw",padx=0,pady=0)
                
            canvasCAM=FigureCanvasTkAgg(fig_CAM,master=frameForCAM)
            canvasCAM.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvasCAM.draw()
        
            toolbarCAM = NavigationToolbar2Tk(canvasCAM, frameForCAM)
            toolbarCAM.update()
            canvasCAM.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            
            stringForDescription = '********Assumption: Top bead is stationary********\nThis will analyze the kymograph for all timepoints\n\nSupply clicks in a specific order to define a custom region\n\n'
            stringForDescription1 = 'Click Number\n[ 1 ] Defines the top limit of the custom area\n[ 2 ]Define starting point for the bottom line of the area \n * Implies a horizontal line to this point\n[ 3+ ] Define any point where you want to change\n * Implies a straight line from previous point to this point'
            tk.Label(customAreaMaster,text=stringForDescription,font=('Helvetica', 10, 'bold'), justify="left",anchor="s").grid(row=4,column=1,sticky="sw",pady=0)
            tk.Label(customAreaMaster,text=stringForDescription1,font=('Helvetica', 10), justify="left",anchor="s").grid(row=5,column=1,sticky="nw",pady=0)
            fig_CAM.canvas.mpl_connect('button_release_event', add_point)
            
            customAreaMaster.bind("<Escape>",escape_top_level)
            customAreaMaster.bind("<Return>",return_coordinates)
            customAreaMaster.mainloop()
            return
        
        """
        Accessory function for the copy_data and extract_data functions that 
        inputs a line object and outputs it to the desired line properties packaged
        in a dictionary.
        """
        def append_data_to_dict(dict_obj,lines_obj,string_descriptor):
            for count, line in enumerate(lines_obj):
                #extract necessary data
                time_vals = line.time_idx
                coordinate_vals = line.coordinate_idx
                
                #add the correct offset value
                time_vals = [y+offset_y for y in time_vals]
                coordinate_vals = [x+offset_x for x in coordinate_vals]
                
                #scale data using the metadata for dx (in nm) dt (in s)
                time_vals = [t*dt for t in time_vals]
                coordinate_vals = [x*dx for x in coordinate_vals]
                
                dict_obj[string_descriptor+" #" +str(count+1) +" Time(s)"] = time_vals
                dict_obj[string_descriptor+" #" +str(count+1) +" Position(nm)"] = coordinate_vals
                if extractIntensitiesOpt.state() == ("selected",):
                    summed_intensity_values = line.sample_from_image(num_pixels = math.ceil(float(entryLineWidthGreedy.get())))
                    dict_obj[string_descriptor+ " #" + str(count+1) + " Summed Photon Counts"] = np.asarray(summed_intensity_values)
            return dict_obj
        
        """
        This function lets the user copy the tracked lines data to the clipboard for use
        in other applications. The major difference between this and the extract button is
        that the copy function does not separate the color of the lines into different
        sheets in the excel library.
        """
        def copy_kt_data(event):
            dict_for_copy = {}
            if filtered_red_lines != "":
                dict_for_copy = append_data_to_dict(dict_for_copy, filtered_red_lines, "Red Line")
            if filtered_green_lines != "":
                dict_for_copy = append_data_to_dict(dict_for_copy, filtered_green_lines, "Green Line")
            if filtered_blue_lines != "":
                dict_for_copy = append_data_to_dict(dict_for_copy, filtered_blue_lines, "Blue Line")
            
            pd_data_frame_dict = pd.DataFrame.from_dict(dict_for_copy,orient='index')
            pd_data_frame_dict = pd_data_frame_dict.transpose()
            pd_data_frame_dict.to_clipboard(sep='\t',index=False)
            return
        
        """
        This function lets the user extract the tracked lines data to a .xlsx file.
        The major difference between this and the extract button is that the copy 
        function does not separate the color of the lines into different sheets in
        the .xlsx document.
        """
        def extract_data(event):
            writer = pd.ExcelWriter(filepath[:-3].replace(" ","_")+"_tracked_lines.xlsx")
            
            if filtered_red_lines != "":
                dict_for_red_extract = {}
                dict_for_red_extract = append_data_to_dict(dict_for_red_extract, filtered_red_lines, "Red Line")
                
                pd_data_frame_dict = pd.DataFrame.from_dict(dict_for_red_extract,orient='index')
                pd_data_frame_dict = pd_data_frame_dict.transpose()
                pd_data_frame_dict.to_excel(writer,sheet_name="Red Lines",index=False,header=True)
            if filtered_green_lines != "":
                dict_for_green_extract = {}
                dict_for_green_extract = append_data_to_dict(dict_for_green_extract, filtered_green_lines, "Green Line")
            
                pd_data_frame_dict = pd.DataFrame.from_dict(dict_for_green_extract,orient='index')
                pd_data_frame_dict = pd_data_frame_dict.transpose()
                pd_data_frame_dict.to_excel(writer,sheet_name="Green Lines",index=False,header=True)
            if filtered_blue_lines != "":
                dict_for_blue_extract = {}
                dict_for_blue_extract = append_data_to_dict(dict_for_blue_extract, filtered_blue_lines, "Blue Line")
                
                pd_data_frame_dict = pd.DataFrame.from_dict(dict_for_blue_extract,orient='index')
                pd_data_frame_dict = pd_data_frame_dict.transpose()
                pd_data_frame_dict.to_excel(writer,sheet_name="Blue Lines",index=False,header=True)
                
            writer.save()
            print('Extraction completed!')
            return
        
        """
        This function swaps the parameters options when the user switches tracking methods
        """
        def swap_parameters(event):
            try:
                frameForGreedyMethod.grid_forget()
            except:
                pass
                
            try:
                frameForLinesMethod.grid_forget()
            except:
                pass
            
            if comboboxMethod.get() == "Greedy":
                frameForGreedyMethod.grid(row=1,column=2,sticky="nw")
            else:
                frameForLinesMethod.grid(row=1,column=2,sticky="nw")
            return
        
        
        def quitKymotracker(event):
            kt_master.destroy()
            return
        
        # define method option frames
        frameForMethodOption = tk.ttk.Frame(kt_master)
        frameForMethodOption.grid(row=0,column=1,columnspan=3)
        
        pad_vertical=2
        tk.ttk.Label(frameForMethodOption,text="Choose Tracking Method: ",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2)
        comboboxMethod= tk.ttk.Combobox(frameForMethodOption,values=['Greedy','Lines'],width=10)
        comboboxMethod.set(['Greedy'])
        comboboxMethod.grid(row=1,column=0,columnspan=2,padx=2)
        comboboxMethod.bind("<<ComboboxSelected>>",swap_parameters)
                
        tk.ttk.Label(frameForMethodOption,text="Track Red Lines?").grid(row=2,column=0,sticky="nw")
        tk.ttk.Label(frameForMethodOption,text="Track Green Lines?").grid(row=3,column=0,sticky="nw")
        tk.ttk.Label(frameForMethodOption,text="Track Blue Lines?").grid(row=4,column=0,sticky="nw")
        
        maxPhotonValues = np.amax(RGB_Data)
        
        redLinesVar = tk.ttk.Checkbutton(frameForMethodOption)
        redLinesVar.grid(row=2,column=1,sticky="nw")
        redLinesVar.invoke()
        
        if np.amax(RGB_Data[:,:,0]) != maxPhotonValues:
            redLinesVar.invoke()
        
        greenLinesVar = tk.ttk.Checkbutton(frameForMethodOption)
        greenLinesVar.grid(row=3,column=1,sticky="nw")
        greenLinesVar.invoke()
        
        if np.amax(RGB_Data[:,:,1]) != maxPhotonValues:
            greenLinesVar.invoke()

        blueLinesVar = tk.ttk.Checkbutton(frameForMethodOption)
        blueLinesVar.grid(row=4,column=1,sticky="nw")
        blueLinesVar.invoke()
        
        if np.amax(RGB_Data[:,:,2]) != maxPhotonValues:
            blueLinesVar.invoke()
        
        # write the separate grid options for tracking
        frameForGreedyMethod = tk.ttk.Frame(kt_master)
        frameForGreedyMethod.grid(row=1,column=2,sticky="nw")
    
        tk.ttk.Label(frameForGreedyMethod,text="Pixel Threshold: ").grid(row=0,column=0,sticky="w",pady=pad_vertical)
        entryPixelThresholdGreedy = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entryPixelThresholdGreedy.grid(row=0,column=1,sticky="w",pady=pad_vertical)
        entryPixelThresholdGreedy.insert(0,"1")
        
        tk.ttk.Label(frameForGreedyMethod,text="Line Width: ").grid(row=1,column=0,sticky="w",pady=pad_vertical)
        entryLineWidthGreedy = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entryLineWidthGreedy.grid(row=1,column=1,sticky="w",pady=pad_vertical)
        entryLineWidthGreedy.insert(0,"6")
        
        tk.ttk.Label(frameForGreedyMethod,text="Window: ").grid(row=2,column=0,sticky="w",pady=pad_vertical)
        entryWindow = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entryWindow.grid(row=2,column=1,sticky="w",pady=pad_vertical)
        entryWindow.insert(0,"8")
        
        tk.ttk.Label(frameForGreedyMethod,text="Sigma: ").grid(row=3,column=0,sticky="w",pady=pad_vertical)
        entrySigma = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entrySigma.grid(row=3,column=1,sticky="w",pady=pad_vertical)
        entrySigma.insert(0,"2")
        
        tk.ttk.Label(frameForGreedyMethod,text="Vel: ").grid(row=4,column=0,sticky="w",pady=pad_vertical)
        entryVel = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entryVel.grid(row=4,column=1,sticky="w",pady=pad_vertical)
        entryVel.insert(0,"0")
        
        tk.ttk.Label(frameForGreedyMethod,text="Diffusion: ").grid(row=5,column=0,sticky="w",pady=pad_vertical)
        entryDiffusion = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entryDiffusion.grid(row=5,column=1,sticky="w",pady=pad_vertical)
        entryDiffusion.insert(0,"0")
        
        tk.ttk.Label(frameForGreedyMethod,text="Sigma Cutoff: ").grid(row=6,column=0,sticky="w",pady=pad_vertical)
        entrySigmaCutoff = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entrySigmaCutoff.grid(row=6,column=1,sticky="w",pady=pad_vertical)
        entrySigmaCutoff.insert(0,"1")
        
        tk.ttk.Label(frameForGreedyMethod,text="Min. Line Length: ").grid(row=7,column=0,sticky="w",pady=pad_vertical)
        entryLineLenGreedy = tk.ttk.Entry(frameForGreedyMethod,width=5)
        entryLineLenGreedy.grid(row=7,column=1,sticky="w",pady=pad_vertical)
        entryLineLenGreedy.insert(0,"10")

        #define trackLinesMethod parameters
        frameForLinesMethod = tk.ttk.Frame(kt_master)
        
        tk.ttk.Label(frameForLinesMethod,text="Line Width: ").grid(row=0,column=0,sticky="w",pady=pad_vertical)
        entryLineWidthLines = tk.ttk.Entry(frameForLinesMethod,width=8)
        entryLineWidthLines.grid(row=0,column=1,sticky="w",pady=pad_vertical)
        entryLineWidthLines.insert(0,"6")
        
        tk.ttk.Label(frameForLinesMethod,text="Max Number Lines: ").grid(row=1,column=0,sticky="w",pady=pad_vertical)
        entryMaxLines = tk.ttk.Entry(frameForLinesMethod,width=8)
        entryMaxLines.grid(row=1,column=1,sticky="w",pady=pad_vertical)
        entryMaxLines.insert(0,"10")
        
        tk.ttk.Label(frameForLinesMethod,text="Start Threshold: ").grid(row=2,column=0,sticky="w",pady=pad_vertical)
        entryStartThreshold = tk.ttk.Entry(frameForLinesMethod,width=8)
        entryStartThreshold.grid(row=2,column=1,sticky="w",pady=pad_vertical)
        entryStartThreshold.insert(0,"0.005")

        tk.ttk.Label(frameForLinesMethod,text="Continuation Threshold: ").grid(row=3,column=0,sticky="w",pady=pad_vertical)
        entryContinuationThreshold = tk.ttk.Entry(frameForLinesMethod,width=8)
        entryContinuationThreshold.grid(row=3,column=1,sticky="w",pady=pad_vertical)
        entryContinuationThreshold.insert(0,"0.005")

        tk.ttk.Label(frameForLinesMethod,text="Angle Weight: ").grid(row=4,column=0,sticky="w",pady=pad_vertical)
        entryAngleWeight = tk.ttk.Entry(frameForLinesMethod,width=8)
        entryAngleWeight.grid(row=4,column=1,sticky="w",pady=pad_vertical)
        entryAngleWeight.insert(0,"10")

        tk.ttk.Label(frameForLinesMethod,text="Min. Line Length: ").grid(row=5,column=0,sticky="w",pady=pad_vertical)
        entryLineLenLines = tk.ttk.Entry(frameForLinesMethod,width=8)
        entryLineLenLines.grid(row=5,column=1,sticky="w",pady=pad_vertical)
        entryLineLenLines.insert(0,"50")
    
        # write additional options frame
        # if you want to change the default settings, delete the second .invoke() if you want it selected or add a second .invoke() to deselect
        frameForAdditionalOpt = tk.ttk.Frame(kt_master)
        frameForAdditionalOpt.grid(row=2,column=1,columnspan=3)
        tk.ttk.Label(frameForAdditionalOpt,text="Custom Area Selection?").grid(row=0,column=0)    
        
        complexAreaOption = tk.ttk.Checkbutton(frameForAdditionalOpt)
        complexAreaOption.grid(row=0,column=1)
        complexAreaOption.invoke()
        complexAreaOption.invoke()
        
        tk.ttk.Label(frameForAdditionalOpt,text="Add Refine Lines Method?").grid(row=1,column=0)    
        
        refineLinesOptCB = tk.ttk.Checkbutton(frameForAdditionalOpt)
        refineLinesOptCB.grid(row=1,column=1)
        refineLinesOptCB.invoke()
        
        tk.ttk.Label(frameForAdditionalOpt,text="Extract Line Intensities?").grid(row=2,column=0)    
        
        extractIntensitiesOpt = tk.ttk.Checkbutton(frameForAdditionalOpt)
        extractIntensitiesOpt.grid(row=2,column=1)
        extractIntensitiesOpt.invoke()
        extractIntensitiesOpt.invoke()
        
        tk.ttk.Label(frameForAdditionalOpt,text="Show on separate plot?").grid(row=3,column=0)    
        
        separatePlotOpt = tk.ttk.Checkbutton(frameForAdditionalOpt)
        separatePlotOpt.grid(row=3,column=1)
        separatePlotOpt.invoke()
        separatePlotOpt.invoke()
        
        tk.ttk.Label(frameForAdditionalOpt,text="Show ROI?").grid(row=4,column=0)    
        
        showRegionOpt = tk.ttk.Checkbutton(frameForAdditionalOpt)
        showRegionOpt.grid(row=4,column=1)
        showRegionOpt.invoke()
        showRegionOpt.invoke()
        
        # write button frame
        ktButtonFrame = tk.ttk.Frame(kt_master)
        ktButtonFrame.grid(row=3,column=1,columnspan=3,sticky="n")
        
        redefineComplexAreaButton = tk.ttk.Button(ktButtonFrame,text="Define Region of Interest",width=25)
        redefineComplexAreaButton.pack(side="top",padx=4,pady=3)
        runKymotrackerButton = tk.ttk.Button(ktButtonFrame,text="Run KymoTracker",width=25)
        runKymotrackerButton.pack(side="top",padx=4,pady=3)
        extractDataAndQuitButton = tk.ttk.Button(ktButtonFrame,text="Extract Data to .xlsx",width=25)
        extractDataAndQuitButton.pack(side="top",padx=5,pady=3)
        tk.ttk.Label(ktButtonFrame,text="Keyboard Shortcuts",justify="left",font=('Helvetica', 10,'bold')).pack(side="top",anchor="nw",pady=2)
        tk.ttk.Label(ktButtonFrame,text="Enter - Run KymoTracker\nCtrl+C - Copy Data to Clipboard\nCtrl+D - Define Custom Area of Analysis\nCtrl+E - Extract Data to .xlsx\nEsc - Quit KymoTracker GUI",justify="left",font=('Helvetica', 8)).pack(side="top",anchor="nw",pady=4)
        
        # bind buttons and functions
        rectproperties = dict(facecolor='cyan', edgecolor = 'blue',alpha=0.2, fill=True)
        draw_temp_rectangle.RS = matplotlib.widgets.RectangleSelector(ax, get_rect_dimensions, drawtype='box',rectprops=rectproperties)
        fig.canvas.mpl_connect('button_press_event', draw_temp_rectangle)
        runKymotrackerButton.bind("<ButtonRelease-1>",call_track_lines)
        redefineComplexAreaButton.bind("<ButtonRelease-1>",define_area_of_analysis)
        extractDataAndQuitButton.bind("<ButtonRelease-1>",extract_data)
        
        #bind keyboard shortcuts
        kt_master.bind("<Escape>",quitKymotracker)
        kt_master.bind("<Return>",call_track_lines)
        kt_master.bind("<Control-c>",copy_kt_data)
        kt_master.bind("<Control-C>",copy_kt_data)
        kt_master.bind("<Control-e>",extract_data)
        kt_master.bind("<Control-E>",extract_data)
        kt_master.bind("<Control-d>",define_area_of_analysis)
        kt_master.bind("<Control-D>",define_area_of_analysis)

# Using a class to describe the GUI interface and provide the different functions
class CTrapGUI(): 
    #this global definition is necessary for the line extract functionality
    global line_list_for_line_extract
    line_list_for_line_extract = []
    
    """
    Upon initialization, build the core GUI system of labels and buttons with default values
    the "change folder" button can then be used to change the file directory to the file of your choice.
    """
    def __init__(self, master):         
        self.__version__ = "1.0.1"
        self.__versionDate__ = "3/1/2021"
        self.__cite__ = "Watters, J.W. (2020) C-Trap .h5 Visualization GUI. Retrieved from https://harbor.lumicks.com/"

        def extract_trap_position_data(h5file,timestampsForIndexing=('',''),timestampsForScanIndexing=('',''),multiScanShading=0):
            amtToDS = float(entryDownSample.get())
            forceString = forceChannelPulldown.get()
            downsampleOpt = checkValueDownsampleOpt.get()
            
            #get sample rate to determine the amount to downsample
            sample_rate = h5file['Force HF']['Force 1x'].sample_rate
            
            trapOpt = whichTrapPosValue.get()
            
            descriptorString = "Trap Position " + trapOpt + " (nm)"
            if downsampleOpt == 1:
                descriptorString += " Downsampled to " + str(amtToDS) + " Hz"
            
            if multiScanShading == 0: #option being zero means that the user does not want the individual scan image highlighted
                if downsampleOpt == 1:    
                    yData = h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                    xData = (h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                    yMin = np.amin(yData)
                    yData = (yData - yMin) * 1000
                    
                    maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                    xMin = np.min(xData)
                    xMax = np.max(xData)
                    xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                else:
                    yData = h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                    yMin = np.amin(yData)
                    yData = (yData - yMin) * 1000
                    
                    xData = (h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps)
                    maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                    xMin = np.min(xData)
                    xMax = np.max(xData)
                    xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                return yData, xData, descriptorString
            else:
                if downsampleOpt == 1:
                    yDataFull = h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                    xDataFull = (h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                    yMin = np.amin(yData)
                    yData = (yData - yMin) * 1000
                    
                    xMin = np.min(xDataFull)
                    xMax = np.max(xDataFull)
                    xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))                            
                    
                    #index into the Force file for the partial data to highlight
                    scanY = h5file["Trap position"][trapOpt][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                    scanY = (scanY - yMin) * 1000
                    scanX = (h5file["Trap position"][trapOpt][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                    scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                else:
                    yDataFull = h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                    yMin = np.amin(yData)
                    yData = (yData - yMin) * 1000
                    
                    xDataFull = (h5file["Trap position"][trapOpt][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                    xMin = np.min(xDataFull)
                    xMax = np.max(xDataFull)
                    xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))                            
                    
                    #index into the Force file for the partial data to highlight
                    scanY = h5file["Trap position"][trapOpt][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                    scanY = (scanY - yMin) * 1000
                    scanX = (h5file["Trap position"][trapOpt][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                    scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                return yDataFull, xDataFull, scanY, scanX, descriptorString

        """
        Function that when called takes all of the force option/plotting options from the interface
        and extracts the correct force and time/distance values from the .h5 file
        Used a different function for FD curves to limit confusion
        """
        def extract_force_data(h5file,timestampsForIndexing=('',''),timestampsForScanIndexing=('',''),multiScanShading=0):
            amtToDS = float(entryDownSample.get())
            forceString = forceChannelPulldown.get()
            distanceVarString = whichDistanceValue.get()
            downsampleOpt = checkValueDownsampleOpt.get()
            stringForceChannel = 'Force ' + forceString
            
            #get sample rate to determine the amount to downsample
            sample_rate = h5file['Force HF']['Force 1x'].sample_rate
            
            if multiScanShading == 0:
                """
                The combobox logical gate is to determine Force vs Time or Force vs. Distance.
                The len > 3 logical gate is to see if the user wants the weighted average of the X and Y forces on the bead or a singular channel.
                If downsampleOpt = 1 then the downsampled data will be extracted instead of the HF data.
                For the FD option --> only low frequency force data is acquired.
                """
                if comboboxForNonRGB.get() == "Force-Time":
                    #force time data collection
                    if len(forceString) < 3:
                        #one channel option
                        if downsampleOpt == 1:
                            yData = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            xData = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                        else:
                            yData = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            xData = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps)
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                    else:
                        #average of channels for one of the beads
                        if downsampleOpt == 1:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yData = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            xData = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                        else:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yData = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            xData = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                            #assert 316 == 1
                else: #fd option
                    if len(forceString) < 3:
                        #get force data for a specific channel
                        yData = h5file["Force LF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        #get distance data that corresponds to the distance
                        xData = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                    else:
                        #get distance data
                        xComponent = h5file["Force LF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yComponent = h5file["Force LF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yData = np.sqrt(np.square(xComponent) + np.square(yComponent))
                        #get force data that corresponds to the distance
                        xData = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data                
                
                if forceString == '2x':
                    yData = yData * -1
                return yData, xData
            
            #logical gate for multiple scan images - to highlight what force regime you are in
            if multiScanShading != 0:
                maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                #same FT vs FD gate
                if comboboxForNonRGB.get() == "Force-Time":
                    #force time data collection
                    if len(forceString) < 3:
                        #one channel option
                        if downsampleOpt == 1:
                            #generate full data
                            yDataFull = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            xDataFull = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))                            
                            
                            #index into the Force file for the partial data to highlight
                            scanY = h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            scanX = (h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                            xMin = np.min(scanX)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                        else:
                            yDataFull = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            xDataFull = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps)
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))                            
                            
                            #index into the Force file for the partial data to highlight
                            scanY = h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                            scanX = (h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].timestamps)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                    else:
                        #average of channels for one of the beads
                        if downsampleOpt == 1:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yDataFull = np.sqrt(np.square(xComponent) + np.square(yComponent))

                            xDataFull = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))

                            #partial index
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yComponent = h5file["Force HF"]["Force " + forceString[0]+"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            scanY = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            scanX = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS))
                            scanX = scanX.timestamps
                            xMin = np.min(scanX)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                        else:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yDataFull = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            xDataFull = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))
                            
                            #partial index
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                            scanY = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            scanX = h5file["Force HF"]["Force 1x"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].timestamps
                            xMin = np.min(scanX)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                            
                else: #fd option
                    if len(forceString) < 3:
                        #get force data for a specific channel
                        yDataFull = h5file["Force LF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        #get distance data that corresponds to the distance
                        xDataFull = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        
                        #get partial data
                        scanY = h5file["Force LF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        scanX = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                    else:
                        #get force data
                        xComponent = h5file["Force LF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yComponent = h5file["Force LF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yDataFull = np.sqrt(np.square(xComponent) + np.square(yComponent))
                        #get distance data that corresponds to the force
                        xDataFull = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        
                        #get partial data
                        xComponent = h5file["Force LF"]["Force " + forceString[0] +"x"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                        yComponent = h5file["Force LF"]["Force " + forceString[0] +"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                        scanY = np.sqrt(np.square(xComponent) + np.square(yComponent))
                        scanX = h5file["Distance"][distanceVarString][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                
                if forceString == '2x':
                    yDataFull = yDataFull * -1
                    try:
                        scanY = scanY * -1
                    except:
                        pass
                    
                return yDataFull, xDataFull, scanY, scanX
        
        
        #code to modify RGB values
        def modify_rgb_image(RGB_code):
            bightnessAddition = int(entryBrightness.get())
            
            grayscaleOption = grayscaleOpt.get()
            
            if grayscaleOption == "No":
                mod_RGB = RGB_code[:,:,:] * [float(entryRed.get()), float(entryGreen.get()), float(entryBlue.get())]
                if bightnessAddition != 0:
                    mod_RGB = (mod_RGB + bightnessAddition)
            else:
                if grayscaleOption == "R":
                     mod_RGB = RGB_code[:,:,0] * float(entryRed.get())
                elif grayscaleOption == "G":
                    mod_RGB = RGB_code[:,:,1] * float(entryGreen.get())
                else:
                    mod_RGB = RGB_code[:,:,2] * float(entryBlue.get())
                
            mod_RGB = mod_RGB.astype(int)
            return mod_RGB
        
        """
        Code to reset the boundary entry values and to remove previous labels describing the maximums
        - max values are set to '-' to help other code recognize when a new plot is being made and new maxes need to be defined
        - min values are set to 0 - if you are analyzing a Force-Distance plot then this zero will be used as a logical gate to set
        the distance minimum value to the distance minimum boundary
        """
        def reset_axis_items(diffFileOpt):
            if diffFileOpt != 0:
                    entryTimeMin.delete(0, "end")
                    entryYForceMin.delete(0, "end")
                    entryDistMin.delete(0, "end")
                    entryYRGBMin.delete(0, "end")
                    entryTimeMax.delete(0, "end")
                    entryDistMax.delete(0,"end")
                    entryYForceMax.delete(0, "end")
                    entryYRGBMax.delete(0, "end")
                    entryScanWidthMin.delete(0,"end")
                    entryScanWidthMax.delete(0,"end")
                    entryTrapPosMax.delete(0,"end")
                    entryTimeMin.insert(0,'0')
                    entryYForceMin.insert(0,'0')
                    entryDistMin.insert(0,'0')
                    entryYRGBMin.insert(0,'0')
                    entryScanWidthMin.insert(0,"0")
                    entryTimeMax.insert(0,'-')
                    entryYForceMax.insert(0,'-')
                    entryDistMax.insert(0,'-')
                    entryYRGBMax.insert(0,'-')
                    entryScanWidthMax.insert(0,"-")
                    entryTrapPosMax.insert(0,"-")
                    
                    try:
                        labelTimeMax.grid_forget()
                    except:
                        pass
                    try:
                        labelDistMax.grid_forget()
                    except:
                        pass
                    try:    
                        labelYRGBMax.grid_forget()
                    except:   
                        pass
                    try:
                        labelYForceMax.grid_forget()
                    except:
                        pass
                    try:
                        labelScanWidthMax.grid_forget()
                    except:
                        pass
                    try:
                        labelTrapPosMax.grid_forget()
                    except:
                        pass
            return
        
        # three methods to extract and plot the different objects
        #Function to extract data type that is specific to the Kymo data type
        def extractAndPlotKymo(h5file,resetBoundsOpt,extract_photons_only="",extract_other_data_only=""):                
            global labelTimeMax
            global labelDistMax
            global labelYRGBMax
            global labelYForceMax
            global labelScanWidthMax
            global labelTrapPosMax
                
            kymoString = typePulldown.get()
            splitKymoString = kymoString.split('-')
            kymoNumber = "-".join(splitKymoString[1:])
            kymoPointer = h5file.kymos[kymoNumber]
            
            reset_axis_items(resetBoundsOpt)
            
            minTimeIndex = kymoPointer.timestamps[0,0].astype(np.int64)
            maxTimeIndex = np.max(kymoPointer.timestamps).astype(np.int64)
            
            RGB_unaltered = saved_color_data
            
            if extract_photons_only != "":
                return RGB_unaltered
            
            if extract_other_data_only != "":
                if comboboxForNonRGB.get() == "Force-Time" or comboboxForNonRGB.get() == "Force-Distance":
                    yData, xData = extract_force_data(h5file,timestampsForIndexing=(minTimeIndex, maxTimeIndex))
                    descriptor = "Force " + forceChannelPulldown.get() + " (pN)"
                    if checkValueDownsampleOpt.get() == 1:
                        descriptor += " Downsampled to " + str(float(entryDownSample.get())) + " Hz"
                else:
                    yData, xData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                return xData, yData, descriptor
            
            if plottingOpt.get() == "Both":
                fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True)
                
                dx = kymoPointer.pixelsize_um[0] * 1000 #pixel size in nm
                dt = (kymoPointer.timestamps[0,1]-kymoPointer.timestamps[0,0]) / 1000000000 #scan time in s
                RGB_altered = modify_rgb_image(RGB_unaltered)
                
                maxTime = len(kymoPointer.timestamps[0,:])
                numberPixels = len(RGB_unaltered)
                maxTrueTime = maxTime*dt
                maxTrueDist = dx*numberPixels/1000 #convert to uM
                
                if entryTimeMax.get() == '-':
                    entryTimeMax.delete(0, "end")
                    entryTimeMax.insert(0,str(round(maxTrueTime,2)))
                    labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueTime,2)))
                    labelTimeMax.grid(row=1,column=3)
                if entryYRGBMax.get() == '-':
                    entryYRGBMax.delete(0, "end")
                    entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                    labelYRGBMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueDist,2)))
                    labelYRGBMax.grid(row=3,column=3)
                    
                if aspectOptionVar.get():
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect="auto")
                else:
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (maxTrueTime / maxTrueDist))
                
                
                if grayscaleOpt.get() == "No":
                    ax1.imshow(RGB_altered, **{**default_kwargs})
                else:
                    ax1.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())

                ax1.set_ylabel(u'Position(\u03bcm)')
                ax1.set_xlabel('Time(s)')
                ax1.set_xlim([float(entryTimeMin.get()),float(entryTimeMax.get())])
                ax1.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                    
                if comboboxForNonRGB.get() == "Force-Distance":
                    forceData, distData = extract_force_data(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                    ax2.plot(distData,forceData)
                        
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                        
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax2.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax2.set_ylabel('Force(pN)')
                    ax2.set_xlabel(u'Distance(\u03bcm)')
                        
                    forceString = forceChannelPulldown.get()
                
                    if len(forceChannelPulldown.get()) > 2:
                        bead_ID = forceString[0]
                        combinedString = 'Bead ' + bead_ID + ' vs. Distance'
                    else:
                        combinedString = 'Channel ' + forceString + ' vs. Distance'
                        ax2.set_title(combinedString)
                elif comboboxForNonRGB.get() == "Trap Pos.-Time":
                    trapPosData, timeData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                    ax2.plot(timeData,trapPosData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0, "end")
                        entryTimeMin.insert(0,"0")
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryTrapPosMax.get() == '-':
                        entryTrapPosMin.delete(0, "end")
                        entryTrapPosMin.insert(0,"0")
                        entryTrapPosMax.delete(0, "end")
                        entryTrapPosMax.insert(0,str(round(max(trapPosData),1)))
                        labelTrapPosMax = tk.ttk.Label(frameForAxis,text=str(round(max(trapPosData),1)))
                        labelTrapPosMax.grid(row=6,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryTrapPosMin.get()),float(entryTrapPosMax.get()))
                    ax2.set_ylabel('Trap Position(nm)')
                    ax2.set_xlabel('Time(s)')
                    ax2.set_title(descriptor)
                else:
                    forceData, timeData = extract_force_data(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                    ax2.plot(timeData,forceData)

                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax2.set_ylabel('Force(pN)')
                    ax2.set_xlabel('Time(s)')
                        
                    forceString = forceChannelPulldown.get()
                        
                    if len(forceChannelPulldown.get()) > 2:
                        bead_ID = forceString[0]
                        if checkValueDownsampleOpt.get() == 1:
                            combinedString = 'Bead ' + bead_ID + ' Downsampled to ' + entryDownSample.get() + ' Hz'
                        else:
                            combinedString = 'Bead ' + bead_ID
                    else:
                        if checkValueDownsampleOpt.get() == 1:
                            combinedString = forceString + ' Channel Downsampled to ' + entryDownSample.get() + ' Hz'
                        else:
                            combinedString = forceString + ' Channel'
                    ax2.set_title(combinedString)

            elif plottingOpt.get() == "Non-RGB Only":
                fig, ax = plt.subplots(constrained_layout=True)
                
                if comboboxForNonRGB.get() == "Force-Time":
                    
                    minTimeIndex = kymoPointer.timestamps[0,0].astype(np.int64)
                    maxTimeIndex = np.max(kymoPointer.timestamps).astype(np.int64)
                    forceData, timeData = extract_force_data(h5file,timestampsForIndexing= (minTimeIndex,maxTimeIndex))
                    ax.plot(timeData,forceData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                        
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel('Time(s)')
                    
                elif comboboxForNonRGB.get() == "Trap Pos.-Time":
                    trapPosData, timeData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                    ax.plot(timeData,trapPosData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0, "end")
                        entryTimeMin.insert(0,"0")
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryTrapPosMax.get() == '-':
                        entryTrapPosMin.delete(0, "end")
                        entryTrapPosMin.insert(0,"0")
                        entryTrapPosMax.delete(0, "end")
                        entryTrapPosMax.insert(0,str(round(max(trapPosData),1)))
                        labelTrapPosMax = tk.ttk.Label(frameForAxis,text=str(round(max(trapPosData),1)))
                        labelTrapPosMax.grid(row=6,column=3)
                        
                    ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax.set_ylim(float(entryTrapPosMin.get()),float(entryTrapPosMax.get()))
                    ax.set_ylabel('Trap Position(nm)')
                    ax.set_xlabel('Time(s)')
                else:
                    forceData, distData = extract_force_data(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                    ax.plot(distData,forceData)
                        
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel(u'Distance(\u03bcm)')
            
            else:
                dx = kymoPointer.pixelsize_um[0] * 1000 #pixel size in nm
                dt = (kymoPointer.timestamps[0,1]-kymoPointer.timestamps[0,0]) / 1000000000 #scan time in s
                RGB_altered = modify_rgb_image(RGB_unaltered)
                maxTime = len(kymoPointer.timestamps[0,:])
                numberPixels = len(RGB_unaltered)
                maxTrueTime = maxTime*dt
                maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                if entryTimeMax.get() == '-':
                    entryTimeMax.delete(0, "end")
                    entryTimeMax.insert(0,str(round(maxTrueTime,2)))
                    labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueTime,2)))
                    labelTimeMax.grid(row=1,column=3)
                if entryYRGBMax.get() == '-':
                    entryYRGBMax.delete(0, "end")
                    entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                    labelYRGBMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueDist,2)))
                    labelYRGBMax.grid(row=3,column=3)
                    
                fig, ax = plt.subplots()
                
                if aspectOptionVar.get():
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect="auto")
                else:
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (maxTrueTime / maxTrueDist))
                if grayscaleOpt.get() == "No":
                    ax.imshow(RGB_altered, **{**default_kwargs})
                else:
                    ax.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                ax.set_ylabel(u'Position(\u03bcm)')
                ax.set_xlabel('Time(s)')
                ax.set_xlim([float(entryTimeMin.get()),float(entryTimeMax.get())])
                ax.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
            return fig
            
        """
        This function is an offshoot of the extractAndPlotScans if the RGB data is 3 dimensional (a stack)
        The multiScanShading option tells the program whether or not to highlight the force data from the different scan image (vs. the stack)
        Currently multiScanShading option has not been tested for Force-Distance curves
        """
        def extractAndPlotMultipleScans(h5file,scanPointer,resetBoundsOpt,extract_photons_only="",extract_other_data_only=""):
            global labelTimeMax
            global labelDistMax
            global labelYRGBMax
            global labelYForceMax
            global labelScanWidthMax
            global labelTrapPosMax
            
            scanString = typePulldown.get()
            splitScanString = scanString.split('-')
            scanNumber = "-".join(splitScanString[1:])
            scanPointer = h5file.scans[scanNumber]
                
            #plot options are already reset in the extractAndPlotStack header
            #see if the stack being observed is a new stack object
            if resetBoundsOpt == 1:
                scanNumber = scanPointer.num_frames
                scaleForStack.config(from_=1, to=scanNumber)
                scaleForStack.set(1)
                    
            #scanPointerFrame = scanPointer(frame=2)
            timestampArray = scanPointer.timestamps
            scanTimeStamp = (timestampArray[int(scaleForStack.get())-1,:,:]) 
            stackRGB = saved_color_data[int(scaleForStack.get())-1,:,:,:]
            if extract_photons_only != "":
                return stackRGB
            
            minTimeScanObj = np.min(timestampArray[0,0,0]).astype(np.int64)
            maxTimeScanObj = np.amax(timestampArray[int(scanPointer.num_frames)-1,:,:]).astype(np.int64)
            minTimeScan = scanTimeStamp[0,0].astype(np.int64)
            maxTimeScan = np.max(scanTimeStamp[:,:]).astype(np.int64)
            
            if extract_other_data_only != "":
                if comboboxForNonRGB.get() == "Force-Time" or comboboxForNonRGB.get() == "Force-Distance":
                    yData, xData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                    descriptor = "Force " + forceChannelPulldown.get() + " (pN)"
                    if checkValueDownsampleOpt.get() == 1:
                        descriptor += " Downsampled to " + str(float(entryDownSample.get())) + " Hz"
                else:
                    yData, xData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                return xData, yData, descriptor
            
            if plottingOpt.get() == "Both":
                fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True)
                RGB_unaltered = stackRGB
                dx = scanPointer.pixelsize_um[0] * 1000 #pixel size in nm
                totalScanWidth = scanPointer.scan_width_um[0]
                dt = (scanPointer.timestamps[0,0,1]-scanPointer.timestamps[0,0,0]) / 1000000000 #scan time in s
                RGB_altered = modify_rgb_image(RGB_unaltered)
                    
                #maxTime = len(timestampArray[0,:])
                numberPixels = len(RGB_unaltered)
                #maxTrueTime = maxTime*dt
                maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                if entryScanWidthMax.get() == '-':
                    entryScanWidthMax.delete(0,"end")
                    entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                    labelScanWidthMax = tk.ttk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                    labelScanWidthMax.grid(row=5,column=3)
                
                if entryYRGBMax.get() == '-':
                    entryYRGBMax.delete(0, "end")
                    entryYRGBMax.insert(0,str(round(maxTrueDist,3)))
                    labelYRGBMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                    labelYRGBMax.grid(row=3,column=3)
                
                if aspectOptionVar.get():
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect="auto")
                else:
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                    
                if grayscaleOpt.get() == "No":
                    ax1.imshow(RGB_altered, **{**default_kwargs})
                else:
                    ax1.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                ax1.set_ylabel(u'Position(\u03bcm)')
                ax1.set_xlabel(u'Width(\u03bcm)')
                ax1.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                ax1.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                    
                if comboboxForNonRGB.get() == "Force-Time":
                    if multiScanPlotOpt.get() == 1:
                        forceData, timeData, scanForceData, scanTimeData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj), timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        
                        ax2.plot(timeData,forceData)
                        ax2.plot(scanTimeData, scanForceData)
                    else:
                        forceData, timeData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                        ax2.plot(timeData,forceData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0,"end")
                        entryTimeMin.insert(0,'0')
                        entryTimeMax.delete(0,"end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                            
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0,"end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax2.set_ylabel('Force(pN)')
                    ax2.set_xlabel('Time(s)')
                        
                    forceString = forceChannelPulldown.get()
                        
                    if len(forceChannelPulldown.get()) > 2:
                        bead_ID = forceString[0]
                        if checkValueDownsampleOpt.get() == 1:
                            combinedString = 'Bead ' + bead_ID + ' Downsampled to ' + entryDownSample.get() + ' Hz'
                        else:
                            combinedString = 'Bead ' + bead_ID
                    else:
                        if checkValueDownsampleOpt.get() == 1:
                            combinedString = forceString + ' Channel Downsampled to ' + entryDownSample.get() + ' Hz'
                        else:
                            combinedString = forceString + ' Channel'
                    ax2.set_title(combinedString)
                    
                elif comboboxForNonRGB.get() == "Trap Pos.-Time":
                    if multiScanPlotOpt.get() == 1:
                        trapPosData, timeData, partialTrapPos, partialTimeData,descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj),timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        ax2.plot(timeData,trapPosData)
                        ax2.plot(partialTimeData, partialTrapPos)
                    else:
                        trapPosData, timeData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                        ax2.plot(timeData,trapPosData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0, "end")
                        entryTimeMin.insert(0,"0")
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryTrapPosMax.get() == '-':
                        entryTrapPosMin.delete(0, "end")
                        entryTrapPosMin.insert(0,"0")
                        entryTrapPosMax.delete(0, "end")
                        entryTrapPosMax.insert(0,str(round(max(trapPosData),1)))
                        labelTrapPosMax = tk.ttk.Label(frameForAxis,text=str(round(max(trapPosData),1)))
                        labelTrapPosMax.grid(row=6,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryTrapPosMin.get()),float(entryTrapPosMax.get()))
                    ax2.set_ylabel('Trap Position(nm)')
                    ax2.set_xlabel('Time(s)') 
                    ax2.set_title(descriptor)
                else:
                    if multiScanPlotOpt.get() == 1:
                        forceData, distData, scanForceData, scanDistData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj),timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        ax2.plot(distData,forceData)
                        ax2.plot(scanDistData, scanForceData)
                    else:
                        forceData, distData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                        ax2.plot(distData,forceData)
                        
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax2.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax2.set_ylabel('Force(pN)')
                    ax2.set_xlabel(u'Distance(\u03bcm)')
                        
                    forceString = forceChannelPulldown.get()
                        
                    if len(forceChannelPulldown.get()) > 2:
                        bead_ID = forceString[0]
                        combinedString = 'Bead ' + bead_ID + ' vs. Distance'
                    else:
                        combinedString = 'Channel ' + forceString + ' vs. Distance'
                    ax2.set_title(combinedString)

                        
            elif plottingOpt.get() == "Non-RGB Only":
                fig, ax = plt.subplots(constrained_layout=True)
                if comboboxForNonRGB.get() == "Force-Time":
                    if multiScanPlotOpt.get() == 1:
                        forceData, timeData, scanForceData, scanTimeData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj), timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                    
                        ax.plot(timeData,forceData)
                        ax.plot(scanTimeData, scanForceData)
                    else:
                        forceData, timeData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                        ax.plot(timeData,forceData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                        
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=4)
                        
                    ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel('Time(s)')
                
                elif comboboxForNonRGB.get() == "Trap Pos.-Time":
                    if multiScanPlotOpt.get() == 1:
                        trapPosData, timeData, partialTrapPos, partialTimeData,descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj),timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        ax2.plot(timeData,trapPosData)
                        ax2.plot(partialTimeData, partialTrapPos)
                    else:
                        trapPosData, timeData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                        ax2.plot(timeData,trapPosData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0, "end")
                        entryTimeMin.insert(0,"0")
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryTrapPosMax.get() == '-':
                        entryTrapPosMin.delete(0, "end")
                        entryTrapPosMin.insert(0,"0")
                        entryTrapPosMax.delete(0, "end")
                        entryTrapPosMax.insert(0,str(round(max(trapPosData),1)))
                        labelTrapPosMax = tk.ttk.Label(frameForAxis,text=str(round(max(trapPosData),1)))
                        labelTrapPosMax.grid(row=6,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryTrapPosMin.get()),float(entryTrapPosMax.get()))
                    ax2.set_ylabel('Trap Position(nm)')
                    ax2.set_xlabel('Time(s)')
                else:
                    if multiScanPlotOpt.get() == 1:
                        forceData, distData, scanForceData, scanDistData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj),timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        ax.plot(distData,forceData)
                        ax.plot(scanDistData, scanForceData)
                    else:
                        forceData, distData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                        ax.plot(distData,forceData)
                        
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel(u'Distance(\u03bcm)')

            else:
                RGB_unaltered = stackRGB
                dx = scanPointer.pixelsize_um[0] * 1000 #pixel size in nm
                dt = (scanPointer.timestamps[0,1]-scanPointer.timestamps[0,0]) / 1000000000 #scan time in s
                totalScanWidth = scanPointer.scan_width_um[0]

                RGB_altered = modify_rgb_image(RGB_unaltered)
                #maxTime = len(scanPointer.timestamps[0,:])
                numberPixels = len(RGB_unaltered)
                #maxTrueTime = maxTime*dt
                maxTrueDist = dx*numberPixels/1000 #convert to uM
                
                if entryScanWidthMax.get() == '-':
                    entryScanWidthMax.delete(0, "end")
                    entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                    labelScanWidthMax = tk.ttk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                    labelScanWidthMax.grid(row=5,column=3)
                if entryYRGBMax.get() == '-':
                    entryYRGBMax.delete(0, "end")
                    entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                    labelYRGBMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                    labelYRGBMax.grid(row=3,column=3)
                    
                fig, ax = plt.subplots()
                
                if aspectOptionVar.get():
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect="auto")
                else:
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                
                if grayscaleOpt.get() == "No":
                    ax.imshow(RGB_altered, **{**default_kwargs})
                else:
                    ax.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                ax.set_ylabel(u'Position(\u03bcm)')
                ax.set_xlabel(u'Width(\u03bcm)')
                ax.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                ax.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])      
            return fig
            
            
        #Function to extract data type that is specific to the Scan data type
        def extractAndPlotScan(h5file,resetBoundsOpt,extract_photons_only="",extract_other_data_only=""):              
            scanString = typePulldown.get()
            splitScanString = scanString.split('-')
            scanNumber = "-".join(splitScanString[1:])
            scanPointer = h5file.scans[scanNumber]
                
            reset_axis_items(resetBoundsOpt)
            
            #logical operator to split off if the scan is not one image (aka it is a stack)
            shapeLenRGB = len(saved_color_data.shape)
            if shapeLenRGB > 3:
                if extract_other_data_only  != "":
                    xData, yData, descriptor = extractAndPlotMultipleScans(h5file,scanPointer,resetBoundsOpt,extract_other_data_only=extract_other_data_only)
                    return xData, yData, descriptor
                elif extract_photons_only != "":
                    photonData = extractAndPlotMultipleScans(h5file,scanPointer,resetBoundsOpt,extract_photons_only=extract_photons_only)
                    return photonData
                else:
                    fig = extractAndPlotMultipleScans(h5file,scanPointer,resetBoundsOpt,extract_photons_only=extract_photons_only,extract_other_data_only=extract_other_data_only)
                    return fig
            global labelTimeMax
            global labelDistMax
            global labelYRGBMax
            global labelYForceMax
            global labelScanWidthMax
            global labelTrapPosMax
            
            minTimeScanObj = np.min(scanPointer.timestamps[0,:]).astype(np.int64)
            maxTimeScanObj = np.max(scanPointer.timestamps[0,:]).astype(np.int64)
            
            RGB_unaltered = saved_color_data
            if extract_photons_only != "":
                return RGB_unaltered
            
            if extract_other_data_only != "":
                if comboboxForNonRGB.get() == "Force-Time" or comboboxForNonRGB.get() == "Force-Distance":
                    yData, xData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                    descriptor = "Force " + forceChannelPulldown.get() + " (pN)"
                    if checkValueDownsampleOpt.get() == 1:
                        descriptor += " Downsampled to " + str(float(entryDownSample.get())) + " Hz"
                else:
                    yData, xData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                return xData, yData, descriptor
            
            if plottingOpt.get() == "Both":
                fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True)
                
                dx = scanPointer.pixelsize_um[0] * 1000 #pixel size in nm
                #dt = (scanPointer.timestamps[0,1]-scanPointer.timestamps[0,0]) / 1000000000 #scan time in s
                totalScanWidth = RGB_unaltered.shape[1] * dx / 1000
                RGB_altered = modify_rgb_image(RGB_unaltered)
                #maxTime = len(scanPointer.timestamps[0,:])
                numberPixels = len(RGB_unaltered)
                #maxTrueTime = maxTime*dt
                maxTrueDist = dx*numberPixels/1000 #convert to uM
                
                if entryScanWidthMax.get() == '-':
                    entryScanWidthMax.delete(0, "end")
                    entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                    labelScanWidthMax = tk.ttk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                    labelScanWidthMax.grid(row=5,column=3)
                if entryYRGBMax.get() == '-':
                    entryYRGBMax.delete(0, "end")
                    entryYRGBMax.insert(0,str(round(maxTrueDist,3)))
                    labelYRGBMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                    labelYRGBMax.grid(row=3,column=3)
                
                if aspectOptionVar.get():
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect="auto")
                else:
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                
                if grayscaleOpt.get() == "No":
                    ax1.imshow(RGB_altered, **{**default_kwargs})
                else:
                    ax1.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                ax1.set_ylabel(u'Position(\u03bcm)')
                ax1.set_xlabel(u'Width(\u03bcm)')
                ax1.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                ax1.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                    
                if comboboxForNonRGB.get() == "Force-Time":
                    forceData, timeData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                    ax2.plot(timeData,forceData)
                    
                    if entryTimeMax.get() == '-':
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax2.set_ylabel('Force(pN)')
                    ax2.set_xlabel('Time(s)')
                        
                    forceString = forceChannelPulldown.get()
                    
                    if len(forceChannelPulldown.get()) > 2:
                        bead_ID = forceString[0]
                        if checkValueDownsampleOpt.get() == 1:
                            combinedString = 'Bead ' + bead_ID + ' Downsampled to ' + entryDownSample.get() + ' Hz'
                        else:
                            combinedString = 'Bead ' + bead_ID
                    else:
                        if checkValueDownsampleOpt.get() == 1:
                            combinedString = forceString + ' Channel Downsampled to ' + entryDownSample.get() + ' Hz'
                        else:
                            combinedString = forceString + ' Channel'
                    
                    ax2.set_title(combinedString)
                    
                elif comboboxForNonRGB.get() == "Trap Pos.-Time":
                    trapPosData, timeData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                    ax2.plot(timeData,trapPosData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0, "end")
                        entryTimeMin.insert(0,"0")
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryTrapPosMax.get() == '-':
                        entryTrapPosMin.delete(0, "end")
                        entryTrapPosMin.insert(0,"0")
                        entryTrapPosMax.delete(0, "end")
                        entryTrapPosMax.insert(0,str(round(max(trapPosData),1)))
                        labelTrapPosMax = tk.ttk.Label(frameForAxis,text=str(round(max(trapPosData),1)))
                        labelTrapPosMax.grid(row=6,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryTrapPosMin.get()),float(entryTrapPosMax.get()))
                    ax2.set_ylabel('Trap Position(nm)')
                    ax2.set_xlabel('Time(s)')
                    ax2.set_title(descriptor)
                else:
                    forceData, distData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                    ax2.plot(distData,forceData)
                            
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax2.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax2.set_ylabel('Force(pN)')
                    ax2.set_xlabel(u'Distance(\u03bcm)')
                        
                    forceString = forceChannelPulldown.get()
                    
                    if len(forceChannelPulldown.get()) > 2:
                        bead_ID = forceString[0]
                        combinedString = 'Bead ' + bead_ID + ' vs. Distance'
                    else:
                        combinedString = 'Channel ' + forceString + ' vs. Distance'
                    ax2.set_title(combinedString)

            elif plottingOpt.get() == "Non-RGB Only":
                fig, ax = plt.subplots(constrained_layout=True)
                if comboboxForNonRGB.get() == "Force-Time":
                    forceData, timeData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                    ax.plot(timeData,forceData)
                    
                    if entryTimeMax.get() == '-':
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                        
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=4)
                        
                    ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel('Time(s)')
                        
                    forceString = forceChannelPulldown.get()
                    
                elif comboboxForNonRGB.get() == "Trap Pos.-Time":
                    trapPosData, timeData, descriptor = extract_trap_position_data(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                    ax2.plot(timeData,trapPosData)
                        
                    if entryTimeMax.get() == '-':
                        entryTimeMin.delete(0, "end")
                        entryTimeMin.insert(0,"0")
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(max(timeData),2)))
                        labelTimeMax = tk.ttk.Label(frameForAxis,text=str(round(max(timeData),2)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryTrapPosMax.get() == '-':
                        entryTrapPosMin.delete(0, "end")
                        entryTrapPosMin.insert(0,"0")
                        entryTrapPosMax.delete(0, "end")
                        entryTrapPosMax.insert(0,str(round(max(trapPosData),1)))
                        labelTrapPosMax = tk.ttk.Label(frameForAxis,text=str(round(max(trapPosData),1)))
                        labelTrapPosMax.grid(row=6,column=3)
                        
                    ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                    ax2.set_ylim(float(entryTrapPosMin.get()),float(entryTrapPosMax.get()))
                    ax2.set_ylabel('Trap Position(nm)')
                    ax2.set_xlabel('Time(s)')
                    
                else:
                    forceData, distData = extract_force_data(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                    ax.plot(distData,forceData)
                    
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0,"end")
                        entryYForceMin.insert(0,str(round(max(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel(u'Distance(\u03bcm)')
                        
                    forceString = forceChannelPulldown.get()
            else:                
                dx = scanPointer.pixelsize_um[0] * 1000 #pixel size in nm
                #dt = (scanPointer.timestamps[0,1]-scanPointer.timestamps[0,0]) / 1000000000 #scan time in s
                totalScanWidth = scanPointer.scan_width_um[0]
                RGB_altered = modify_rgb_image(RGB_unaltered)
                #maxTime = len(scanPointer.timestamps[0,:])
                numberPixels = len(RGB_unaltered)
                #maxTrueTime = maxTime*dt
                maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    
                if entryScanWidthMax.get() == '-':
                    entryScanWidthMax.delete(0, "end")
                    entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                    labelScanWidthMax = tk.ttk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                    labelScanWidthMax.grid(row=5,column=3)
                    
                if entryYRGBMax.get() == '-':
                    entryYRGBMax.delete(0, "end")
                    entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                    labelYRGBMax = tk.ttk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                    labelYRGBMax.grid(row=3,column=3)
                    
                    
                fig, ax = plt.subplots()
                
                if aspectOptionVar.get():
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect="auto")
                else:
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                
                if grayscaleOpt.get() == "No":
                    ax.imshow(RGB_altered, **{**default_kwargs})
                else:
                    ax.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())

                ax.set_ylabel(u'Position(\u03bcm)')
                ax.set_xlabel(u'Width(\u03bcm)')
                ax.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                ax.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
            return fig

            
        def extractAndPlotFD(h5file,resetBoundsOpt,extract_other_data_only=""):             
            global labelTimeMax
            global labelDistMax
            global labelYRGBMax
            global labelYForceMax
            global labelScanWidthMax
            global labelTrapPosMax
            
            fdString = typePulldown.get()
            splitFDString = fdString.split('-')
            fdNumber = splitFDString[-1]
            fdPointer = h5file.fdcurves[fdNumber]
                
            reset_axis_items(resetBoundsOpt)
                
            if plottingOpt.get() == "Non-RGB Only" or plottingOpt.get() == "Both":
                if plottingOpt.get() == "Both" or plottingOpt.get() == "RGB Only":
                    print('FD Curves are only currently supported for Force-Distance plotting only\nIf RGB Data is in the file select the kymo/scan object associated with it to get Force and/or RGB Images')                         
                    
                forceString = forceChannelPulldown.get()
                distanceVarString = whichDistanceValue.get()
                if len(forceString) < 3:
                    modifiedFDPointer = fdPointer.with_channels(force=forceString,distance=distanceVarString[-1])
                    forceData = modifiedFDPointer.f.data
                    if forceString == "2x":
                        forceData = -1*forceData
                    distData = modifiedFDPointer.d.data
                else:
                    modifiedFDPointer = fdPointer.with_channels(force=forceString[0],distance=distanceVarString[-1])
                    forceData = modifiedFDPointer.f.data
                    distData = modifiedFDPointer.d.data
                
                if extract_other_data_only != "":
                    return distData,forceData, forceString
                
                fig, ax = plt.subplots(constrained_layout=True)
                ax.plot(distData,forceData)
                    
                if entryDistMax.get() == '-':
                    entryDistMin.delete(0, "end")
                    entryDistMin.insert(0,str(round(min(distData),2)))
                    entryDistMax.delete(0, "end")
                    entryDistMax.insert(0,str(round(max(distData),2)))
                    labelDistMax = tk.ttk.Label(frameForAxis,text=str(round(max(distData),2)))
                    labelDistMax.grid(row=2,column=3)
                        
                if entryYForceMax.get() == '-':
                    entryYForceMin.delete(0, "end")
                    entryYForceMin.insert(0,str(round(min(forceData),2)))
                    entryYForceMax.delete(0, "end")
                    entryYForceMax.insert(0,str(round(max(forceData),2)))
                    labelYForceMax = tk.ttk.Label(frameForAxis,text=str(round(max(forceData),2)))
                    labelYForceMax.grid(row=4,column=3)
                        
                ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                ax.set_ylabel('Force(pN)')
                ax.set_xlabel(u'Distance(\u03bcm)')
            else:
                print('FD Curves are only supported for Force-Distance plotting only\nIf RGB Data is in the file select the kymo/scan object associated with it to get Force-Time or RGB Images')
                fig, ax = plt.subplots(constrained_layout=True)
            return fig
        
        """
        Function that will be linked to both the Save Image and Draw Plot feature.
        This function resets the channel data and writes the channel data (writeMetadat(file)) 
            if the new file/component type is different from the last one plotted.
        The typeComponent GUI input is then scanned and then output is sent to either
            extractAndPlotKymo, extractAndPlotScan (and its subset extractAndPlotStack),
            and extractAndPlotFD.
        """
        def generateFigure(extract_photon_count_options="",extract_other_data=""):
            #remove previous metadata and add metadata description from the .h5 file
            def writeMetadata(file):
                global metadataLabel
                
                #remove the label of previous metadata
                try:
                    metadataLabel.grid_forget()
                except:
                    pass
                
                metadataText = file.description
                metadataLabel = tk.ttk.Label(frameForMetadata,text=metadataText,justify=tk.LEFT)
                metadataLabel.grid(row=1,column=0,rowspan=1,columnspan=1,sticky='w',padx=2)
                return
            
            #logic to determine if it is the first time this is being plotted
            figureDrawRequest = [directoryPulldown.get(),typePulldown.get(),forceChannelPulldown.get(),whichDistanceValue.get(),comboboxForNonRGB.get(),checkValueDownsampleOpt.get(),whichTrapPosValue.get()]
            
            currentFile = lk.File(figureDrawRequest[0])
            
            #if resetPlotOpt = 0 then the maximums will not be reset
            resetPlotOpt = 0
            
            #if it is the same file - dont update any exterior parameters
            if figureDrawRequest == lastPlottedFig:
                pass
            #if it is in the same H5 file - dont update metadata just update graph maximums
            elif figureDrawRequest[0] == lastPlottedFig[0]:
                lastPlottedFig[1] = figureDrawRequest[1]
                lastPlottedFig[2] = figureDrawRequest[2]
                lastPlottedFig[3] = figureDrawRequest[3]
                lastPlottedFig[4] = figureDrawRequest[4]
                lastPlottedFig[5] = figureDrawRequest[5]
                lastPlottedFig[6] = figureDrawRequest[6]
                resetPlotOpt = 2
            #if it is in a different H5 file - update metadata and update graph maximums
            else:
                lastPlottedFig[0] = figureDrawRequest[0]
                lastPlottedFig[1] = figureDrawRequest[1]
                lastPlottedFig[2] = figureDrawRequest[2]
                lastPlottedFig[3] = figureDrawRequest[3]
                lastPlottedFig[4] = figureDrawRequest[4]
                lastPlottedFig[5] = figureDrawRequest[5]
                lastPlottedFig[6] = figureDrawRequest[6]
                resetPlotOpt = 1
                writeMetadata(currentFile)
            
            fileComponent = typePulldown.get()
            splitFileType= fileComponent.split('-')
            filetype = splitFileType[0]
            
            #logical filter to determine what h5 component type is analyzed
            if filetype == "kymos":
                if extract_photon_count_options == "" and extract_other_data == "":
                    figureReturned = extractAndPlotKymo(currentFile,resetPlotOpt)
                elif extract_photon_count_options != "":
                    not_transformed_RGB = extractAndPlotKymo(currentFile,resetPlotOpt,extract_photons_only=extract_photon_count_options)
                    return not_transformed_RGB
                else:
                    x_data, other_data, descriptor = extractAndPlotKymo(currentFile,resetPlotOpt,extract_other_data_only=extract_other_data)
                    return x_data, other_data, descriptor
            elif filetype == "scans":
                if extract_photon_count_options == "" and extract_other_data == "":
                    figureReturned = extractAndPlotScan(currentFile,resetPlotOpt)
                elif extract_photon_count_options != "":
                    not_transformed_RGB = extractAndPlotScan(currentFile,resetPlotOpt,extract_photons_only=extract_photon_count_options)
                    return not_transformed_RGB
                else:
                    x_data, other_data, descriptor = extractAndPlotScan(currentFile,resetPlotOpt,extract_other_data_only=extract_other_data)
                    return x_data, other_data, descriptor
            else:
                if extract_photon_count_options == "" and extract_other_data == "":
                    figureReturned = extractAndPlotFD(currentFile,resetPlotOpt)
                elif extract_photon_count_options != "":
                    print("No photon count data is collected in fdcurve objects")
                    return "no_RGB_found"
                else:
                    x_data, other_data, descriptor = extractAndPlotFD(currentFile,resetPlotOpt,extract_other_data_only=extract_other_data)
                    return x_data, other_data, descriptor
            
            #Make title of the graph
            figureReturned.suptitle(entryPlotTitle.get(), fontsize=16,va='top')
            #figureReturned.set_size_inches(8,6.4)
            figureReturned.set_dpi(110)
            #make frame, canvas and toolbar objects for the figure
            frameForCanvas = tk.ttk.Frame(master,relief=tk.FLAT)
            frameForCanvas.grid(row=0,rowspan=20,column=0,columnspan=2,sticky="nw",padx=0,pady=0)
            
            canvas=FigureCanvasTkAgg(figureReturned,master=frameForCanvas)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvas.draw()
            
            toolbar = NavigationToolbar2Tk(canvas, frameForCanvas)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            return figureReturned
        
        #Combobox bound function that lists the different .h5 files in that folder, 
        #different file components, and different distance options
        def changeFileComponents(event):
            file = lk.File(directoryPulldown.get())
            listOfFileTypes = []
            
            for key in file.kymos.keys():
                listOfFileTypes.append('kymos-'+key)
            
            for key in file.scans.keys():
                listOfFileTypes.append('scans-'+key)
            
            for key in file.fdcurves.keys():
                listOfFileTypes.append('fdcurves-'+key)

            typePulldown['values']= listOfFileTypes
            typePulldown.current('0')
            
            #now pull in the distance options and write it to the corresponding combobox pulldown
            try:
                listOfDistOptions = []
                for i in file['Distance']:
                    listOfDistOptions.append(i)
                whichDistanceValue['values'] = listOfDistOptions
                whichDistanceValue.current('0')
            except:
                whichDistanceValue['values'] = []
                whichDistanceValue.set('')
                pass
            
            try:
                listOfTrapPos = []
                for i in file['Trap position']:
                    listOfTrapPos.append(i)
                whichTrapPosValue['values'] = listOfTrapPos
                whichTrapPosValue.current('0')
            except:
                whichTrapPosValue['values'] = []
                whichTrapPosValue.set('')
                pass

            # load color data
            global saved_color_data
            splitTypePulldown = (typePulldown.get()).split('-')[0]
            
            if splitTypePulldown == "fdcurves":
                saved_color_data = 0   
            elif splitTypePulldown == "kymos":
                saved_color_data = file.kymos["-".join(typePulldown.get().split('-')[1:])].rgb_image
            else:
                saved_color_data = file.scans["-".join(typePulldown.get().split('-')[1:])].rgb_image
                
            fileNameToSave = typePulldown.get()
            entrySaveFile.delete(0, "end")
            entrySaveFile.insert(0,fileNameToSave)
            preload_RGB_and_changeSaveName(1)
            return
        
        #function to reset file comboboxes upon selecting a new folder        
        def changeH5FileDir(event):
            pointerToDir = filedialog.askdirectory(parent=master, title='Please select a directory to analyze')
        
            os.chdir(pointerToDir)
            h5FileList= glob.glob("*.h5")
            if len(h5FileList) > 0:
                h5FileListTuple = tuple(h5FileList)
                directoryPulldown['values'] = h5FileListTuple
                maxFileNameLength = len(max(h5FileList,key=len))
                
                if maxFileNameLength > 60:
                    directoryPulldown['width'] = maxFileNameLength + 5
                else:
                    directoryPulldown['width'] = 60
                directoryPulldown.current(0)
            else:    
                print('No .h5 files were found in the folder selected')
                directoryPulldown['values'] = []
                directoryPulldown.set('')
                typePulldown['values'] = []
                typePulldown.set('')
                entrySaveFile.delete(0, "end")
                return
                
            
            #same commands as the changeFileComponents() function
            file = lk.File(directoryPulldown.get())
            listOfFileTypes = []
            
            for key in file.kymos.keys():
                listOfFileTypes.append('kymos-'+key)
            
            for key in file.scans.keys():
                listOfFileTypes.append('scans-'+key)
            
            for key in file.fdcurves.keys():
                listOfFileTypes.append('fdcurves-'+key)
            
            typePulldown['width'] = len(max(listOfFileTypes,key=len)) + 5
            typePulldown['values']= listOfFileTypes
            typePulldown.current('0')
            
            entrySaveFile['width'] = len(max(listOfFileTypes,key=len)) + 8
            entryPlotTitle['width'] =  len(max(listOfFileTypes,key=len)) + 8
            
            fileNameToSave = typePulldown.get()
            entrySaveFile.delete(0, "end")
            entrySaveFile.insert(0,fileNameToSave)
            
            #now pull in the distance options
            try:
                listOfDistOptions = []
                for i in file['Distance']:
                    listOfDistOptions.append(i)
                whichDistanceValue['values'] = listOfDistOptions
                whichDistanceValue.current('0')
            except:
                whichDistanceValue['values'] = []
                whichDistanceValue.set('')
                pass
            #now pull in the trap position options
            try:
                listOfTrapPos = []
                for i in file['Trap position']:
                    listOfTrapPos.append(i)
                whichTrapPosValue['values'] = listOfTrapPos
                whichTrapPosValue.current('0')
            except:
                whichTrapPosValue['values'] = []
                whichTrapPosValue.set('')
                pass
            
            # load color data
            preload_RGB_and_changeSaveName(1)
            return
        
        #define RGB Constants - change these values based on your normal images
        def define_Global_Defaults():
            global defaultDict
            defaultDict = {}
            defaultDict['red']= 10
            defaultDict['green']= 10
            defaultDict['blue']= 10
            return
        
        #updatePlot button bound event to generate the figure
        def buildPlot(event):
            global figureToSave
            figureToSave = generateFigure()
            return
        
        #saveImageButton bound event to generate a figure, refresh the image,
        #   and save image/metadata
        def saveFigure(event):
            figureToSave = generateFigure()
            
            imageStringPrefix = ((entrySaveFile.get()).replace(" ","_")).replace("-","_")
            imageSuffix = imageFormatOption.get()
            plt.savefig(imageStringPrefix + '.' +imageSuffix,bbox_inches="tight")
            
            #split metadata from the h5 file
            fileName = lk.File(directoryPulldown.get())
            
            metaData = fileName.description
            
            #open,write, and save metadata File
            metaDataFileString = imageStringPrefix.replace(' ','_')+ '_desc' +'.txt'
            metaDataFile = open(metaDataFileString,'w')
            metaDataFile.write(f'Metadata for image analysis of {directoryPulldown.get()} --- {typePulldown.get()}\n')
            
            metaDataFile.write(metaData)
            
            metaDataFile.write('------------------------------------------\n')
            metaDataFile.write(f'Image Scaling Factors\nR\t{entryRed.get()}\nG\t{entryGreen.get()}\nB\t{entryBlue.get()}\n')
            metaDataFile.write(f'Brightness Shift: {entryBrightness.get()}\n')
            
            arrayForKymoOrScan = (typePulldown.get()).split('-')
            
            if  arrayForKymoOrScan[0] == "scans":
                scanNumber = arrayForKymoOrScan[-1]
                scanPointer = fileName.scans[scanNumber]
                
                if scanPointer.num_frames == 1:
                    timestampArray = scanPointer.timestamps
                    
                    metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0]/1e9} seconds\n")
                    metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[:,:]) - timestampArray[0,0])/1e9} seconds\n")
                else:
                    timestampArray = scanPointer.timestamps
                    metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0,0]/1e12} seconds\n")
                    metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[0,:,:]) - timestampArray[0,0,0])/1e9} seconds\n")

            return
        
        """
        This function has been borrowed from Michael Wasserman's code to extract the image data
        for the current file being analyzed. This code is different from the batch extract code (h5_image_extract.py)
        is that no force vs. RGB summary plot are made. This functionality is still availiable using the save plot
        button to get force-correlated RGB data (or Non-RGB Only data)
        This makes the image extract button much faster because you do not have to read in large force datasets
        """
        def extractImageCTrap(event):
            temp_file_name = directoryPulldown.get()
            
            # extract and save experimental description
            def save_exp_desc(exp_desc,filename_without_extension):
                imageStringPrefix = filename_without_extension
                
                fileName = lk.File(directoryPulldown.get())
                
                metaDataFileString = imageStringPrefix.replace(' ','_')+ '_desc' +'.txt'
                metaDataFile = open(metaDataFileString,'w')
                metaDataFile.write(f'Metadata for image analysis of {directoryPulldown.get()} --- {typePulldown.get()}\n')
                
                metaDataFile.write(exp_desc)
            
                metaDataFile.write('------------------------------------------\n')
                metaDataFile.write(f'Image Scaling Factors\nR\t{entryRed.get()}\nG\t{entryGreen.get()}\nB\t{entryBlue.get()}\n')
                metaDataFile.write(f'Brightness Shift: {entryBrightness.get()}\n')
            
                arrayForKymoOrScan = (typePulldown.get()).split('-')
            
                if  arrayForKymoOrScan[0] == "scans":
                    scanNumber = arrayForKymoOrScan[-1]
                    scanPointer = fileName.scans[scanNumber]
                
                    if scanPointer.num_frames == 1:
                        timestampArray = scanPointer.timestamps
                        metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0]/1e9} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[:,:]) - timestampArray[0,0])/1e9} seconds\n")
                    else:
                        timestampArray = scanPointer.timestamps
                        metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0,0]/1e9} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[0,:,:]) - timestampArray[0,0,0])/1e9} seconds\n")
            
            # save scans without labeled axes. Important for conserving pixel num etc
            def save_image_for_CTrapViewer(image_num, image_type):
                if image_type == 'kymo':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.kymos[s].pixelsize_um[0] * 1000 #pixel size in nm
                    #dt = tempfile.kymos[s].json['scan volume']['scan axes'][0]['scan time (ms)'] #scan time in ms
                    # Note: the above dt extraction is incorrect if you manually add a wait time between line scans
                    dt = (tempfile.kymos[s].timestamps[0,1]-tempfile.kymos[s].timestamps[0,0])/1e6 #scan time in ms
                    # save kymograph without labeled axes. Important for conserving pixel num etc
                    filename_png = filename_without_extension + "_dx_" + str(dx) + "nm_dt_" + str(dt) + "ms" + ".png"
                    modified_RGB = modify_rgb_image(tempfile.kymos[s].rgb_image)
                    maxTime = len(tempfile.kymos[s].timestamps[0,:])
                    numberPixels = len(modified_RGB)
                    maxTrueTime = maxTime*dt /1000
                    maxTrueDist = dx*numberPixels/1000 
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect=(modified_RGB.shape[0] / modified_RGB.shape[1]) * (maxTrueTime / maxTrueDist))
                    fig, ax = plt.subplots()
                    ax.imshow(modified_RGB,**{**default_kwargs})
                    ax.axis('off')
                    plt.savefig(filename_png,bbox_inches="tight",pad_inches=0)
                    plt.close()
                elif image_type == 'scan':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.scans[t].pixelsize_um[0] * 1000 #pixel size in nm
                    if tempfile.scans[t].num_frames == 1: #single 2D scan
                        dt = (tempfile.scans[t].timestamps[0,1]-tempfile.scans[t].timestamps[0,0])/1e6 #line time in ms
                        image_slice = tempfile.scans[t].rgb_image
                        image_slice_mod = modify_rgb_image(image_slice)
                        heightScan = image_slice.shape[1] * dx / 1e3
                        totalScanWidth =  tempfile.scans[t].scan_width_um[0]
                        default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_mod.shape[0] / image_slice_mod.shape[1]) * (totalScanWidth/ heightScan))
                            
                        filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms" + ".png"
                        fig, ax = plt.subplots(constrained_layout=True)
                        ax.imshow(image_slice_mod, **{**default_kwargs})
                        plt.axis('off')
                        plt.savefig(filename_png,bbox_inches="tight", pad_inches = 0)
                        
                    elif tempfile.scans[t].num_frames > 1: #image stack
                        dt = (tempfile.scans[t].timestamps[0,1,0]-tempfile.scans[t].timestamps[0,0,0])/1e6 #line time in ms

                        for i in range(tempfile.scans[t].num_frames):
                            frame_num = i + 1
                            filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms f" + str(frame_num) + ".png"
                            image_slice = tempfile.scans[t].rgb_image[i,:,:,:]
                            heightScan = image_slice.shape[1] * dx / 1e3
                            image_slice_int = modify_rgb_image(image_slice) #convert to uint8 for saving as .png <- JWW changed this, might need to change back
                            totalScanWidth =  tempfile.scans[t].scan_width_um[0]
                            
                            #listOfInterpols = ['none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric','catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
                            default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_int.shape[0] / image_slice_int.shape[1]) * (totalScanWidth/ heightScan),interpolation="nearest")
                            
                            fig, ax = plt.subplots(constrained_layout=True)
                            ax.imshow(image_slice_int, **{**default_kwargs})
                            ax.axis('off')
                            #filename_png = item + filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms f" + str(frame_num) + ".png"
                            plt.savefig(filename_png,bbox_inches = 'tight', pad_inches = 0)
          
            tempfile = lk.File(temp_file_name) # load h5 file
            filename_without_extension = tempfile.h5.filename.replace(".h5", "")  # "file.h5" -> "file"
            exp_desc = tempfile.description        
            save_exp_desc(exp_desc, filename_without_extension) # save .txt with experimental description
                
            # loop to extract each kymo 
            for s in list(tempfile.kymos): 
                image_type = 'kymo'                                 
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_CTrapViewer(s, image_type)
               
            # loop to extract each 2D scan 
            for t in list(tempfile.scans): 
                image_type = 'scan'            
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_CTrapViewer(t, image_type)
            return
        
        """
        This is the same basic functionality as the extractImageCTrap() function.
        The most important changes are that axis are labeled for and multi-image
        scans are labeled with a timestamp. One image scans are formatted the same
        way as multi-image scans so that they can be put into an ImageJ montage easily.
        """
        def extractImageImageJ(event):
            # extract and save experimental description
            def save_exp_desc(exp_desc,filename_without_extension):
                imageStringPrefix = filename_without_extension
                
                fileName = lk.File(directoryPulldown.get())
                
                metaDataFileString = imageStringPrefix.replace(' ','_')+ '_desc' +'.txt'
                metaDataFile = open(metaDataFileString,'w')
                metaDataFile.write(f'Metadata for image analysis of {directoryPulldown.get()} --- {typePulldown.get()}\n')
                
                metaDataFile.write(exp_desc)
            
                metaDataFile.write('------------------------------------------\n')
                metaDataFile.write(f'Image Scaling Factors\nR\t{entryRed.get()}\nG\t{entryGreen.get()}\nB\t{entryBlue.get()}\n')
                metaDataFile.write(f'Brightness Shift: {entryBrightness.get()}\n')
            
                arrayForKymoOrScan = (typePulldown.get()).split('-')
                if  arrayForKymoOrScan[0] == "scans":
                    scanNumber = arrayForKymoOrScan[-1]
                    scanPointer = fileName.scans[scanNumber]
                
                    
                    if scanPointer.num_frames == 1:
                        timestampArray = scanPointer.timestamps
                        firstTimestamp = timestampArray[0,0]/1e9
                        timestepBetweenImages = (np.max(timestampArray[:,:]) - timestampArray[0,0])/1e9
                        metaDataFile.write(f"First Timestamp Value: {firstTimestamp} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {timestepBetweenImages} seconds\n")
                    else:
                        timestampArray = scanPointer.timestamps
                        firstTimestamp = timestampArray[0,0,0]/1e9
                        timestepBetweenImages = (np.max(timestampArray[0,:,:]) - timestampArray[0,0,0])/1e9
                        metaDataFile.write(f"First Timestamp Value: {firstTimestamp} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {timestepBetweenImages} seconds\n")
                
                return
            
            # save scans without labeled axes. Important for conserving pixel num etc
            def save_image_for_ImageJ(image_num, image_type):
                if image_type == 'kymo':
                    print('No special settings are applied to kymographs when using the "Export Image For ImageJ Montage" Button - it just shows the image with the axis labeled')
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.kymos[s].pixelsize_um[0] * 1000 #pixel size in nm
                    #dt = tempfile.kymos[s].json['scan volume']['scan axes'][0]['scan time (ms)'] #scan time in ms
                    # Note: the above dt extraction is incorrect if you manually add a wait time between line scans
                    dt = (tempfile.kymos[s].timestamps[0,1]-tempfile.kymos[s].timestamps[0,0])/1e6 #scan time in ms
                    # save kymograph without labeled axes. Important for conserving pixel num etc
                    filename_png = filename_without_extension + "_dx_" + str(dx) + "nm_dt_" + str(dt) + "ms" + ".png"
                    modified_RGB = modify_rgb_image(tempfile.kymos[s].rgb_image)
                    maxTime = len(tempfile.kymos[s].timestamps[0,:])
                    numberPixels = len(modified_RGB)
                    maxTrueTime = maxTime*dt / 1000
                    maxTrueDist = dx*numberPixels/1000 #different calculation because dt is calculated differently elsewhere
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect=(modified_RGB.shape[0] / modified_RGB.shape[1]) * (maxTrueTime / maxTrueDist))
                    fig, ax = plt.subplots()
                    ax.imshow(modified_RGB,**{**default_kwargs})
                    ax.set_ylabel(u'Position(\u03bcm)')
                    ax.set_xlabel(u'Time(s)')
                    plt.savefig(filename_png,bbox_inches="tight")
                    plt.close()
                elif image_type == 'scan':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.scans[t].pixelsize_um[0] * 1000 #pixel size in nm

                    if tempfile.scans[t].num_frames == 1: #single 2D scan
                        dt = (tempfile.scans[t].timestamps[0,1]-tempfile.scans[t].timestamps[0,0])/1e6 #line time in ms
                        image_slice = tempfile.scans[t].rgb_image
                        image_slice_mod = modify_rgb_image(image_slice)
                        heightScan = image_slice.shape[0] * dx / 1e3
                        #totalScanWidth =  tempfile.scans[t].scan_width_um
                        totalScanWidth = image_slice.shape[1] * dx / 1e3
                        default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_mod.shape[0] / image_slice_mod.shape[1]) * (totalScanWidth/ heightScan))
                            
                        filename_png = filename_without_extension + "_ImageJ_dx_" + str(dx) + "nm_dt_" + str(dt) + "ms" + ".png"
                        fig, ax = plt.subplots(constrained_layout=True)
                        ax.imshow(image_slice_mod, **{**default_kwargs})
                        ax.set_ylabel(u'Position(\u03bcm)')
                        ax.set_xlabel(u'Scan Width(\u03bcm)')
                        fig.suptitle(' ', y = -0.04) #empty title is made so that the dimensions can match a separate stack if one image was taken as a pre-FRAP image
                        plt.savefig(filename_png,bbox_inches="tight", pad_inches = 0)
                        
                    elif tempfile.scans[t].num_frames > 1: #image stack
                        dt = (tempfile.scans[t].timestamps[0,1,0]-tempfile.scans[t].timestamps[0,0,0])/1e6 #line time in ms
                        
                        initialTimeTotal = tempfile.scans[t].timestamps[0,0,0]
                        
                        for i in range(tempfile.scans[t].num_frames):
                            frame_num = i + 1
                            filename_png = filename_without_extension + "ImageJ_dx_" + str(dx) + "nm_dt_" + str(dt) + "ms_f" + str(frame_num) + ".png"
                            image_slice = tempfile.scans[t].rgb_image[i,:,:,:]
                            initialTimeScan = tempfile.scans[t].timestamps[i,0,0]
                            heightScan = image_slice.shape[0] * dx / 1e3
                            image_slice_int = modify_rgb_image(image_slice) #convert to uint8 for saving as .png <- JWW changed this, might need to change back
                            #totalScanWidth =  tempfile.scans[t].scan_width_um
                            totalScanWidth = image_slice.shape[1] * dx / 1e3
                            #listOfInterpols = ['none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric','catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
                            default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_int.shape[0] / image_slice_int.shape[1]) * (totalScanWidth/ heightScan),interpolation="nearest")
                            
                            fig, ax = plt.subplots(constrained_layout=True)
                            ax.imshow(image_slice_int, **{**default_kwargs})
                            ax.set_ylabel(u'Position(\u03bcm)')
                            ax.set_xlabel(u'Scan Width(\u03bcm)')
                            fig.suptitle(f'{round((initialTimeScan - initialTimeTotal) / 1e9,2)} s', y = -0.04)
                            plt.savefig(filename_png,bbox_inches = 'tight', pad_inches = 0)
                            plt.close()
            
            temp_file_name = directoryPulldown.get()
            tempfile = lk.File(temp_file_name) # load h5 file
            filename_without_extension = (tempfile.h5.filename.replace(".h5", "")).replace(" ","_")  # "file.h5" -> "file"
            exp_desc = tempfile.description        
            save_exp_desc(exp_desc, filename_without_extension) # save .txt with experimental description
                
            # loop to extract each kymo 
            for s in list(tempfile.kymos): 
                image_type = 'kymo'                                 
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_ImageJ(s, image_type)
               
            # loop to extract each 2D scan 
            for t in list(tempfile.scans): 
                image_type = 'scan'            
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_ImageJ(t, image_type)
            return
        
        """
        This function has been borrowed from Michael Wasserman's code to extract the desired force data
        for the current file being analyzed.
        This task is computationally expensive if you want the HF data - be wary of this
        """
        def extractForceMethod(event):
            def extractForceCommand(settings_list):
                # Extract force data
                filename = directoryPulldown.get()
                temp_file = lk.File(filename)
                filename_no_extension = filename.replace(".h5",'')
                force1xHF = temp_file['Force HF']['Force 1x']
                force1yHF = temp_file['Force HF']['Force 1y']
                force2xHF = temp_file['Force HF']['Force 2x']
                force2yHF = temp_file['Force HF']['Force 2y']
                sample_rate = force1xHF.sample_rate # Hz
                downsampled_rate = float(entryDownSample.get())
                
                force1x_downsamp = force1xHF.downsampled_by(int(sample_rate/downsampled_rate))
                force1y_downsamp = force1yHF.downsampled_by(int(sample_rate/downsampled_rate))
                force2x_downsamp = force2xHF.downsampled_by(int(sample_rate/downsampled_rate))
                force2y_downsamp = force2yHF.downsampled_by(int(sample_rate/downsampled_rate))
                pooled_force_data = [force1xHF, force1yHF, force2xHF, force2yHF, force1x_downsamp,
                                     force1y_downsamp, force2x_downsamp, force2y_downsamp]
                
                #extract HF and downsampled time stamps
                time = force1xHF.timestamps/1e9 # time traces (seconds)
                time = time - time[0]
                time_downsamp = force1x_downsamp.timestamps/1e9
                time_downsamp = time_downsamp - time_downsamp[0]

                #indices for what force channels to take                
                force_ind = [i for i, v in enumerate(settings_list) if v[:][1] == True]
                if len(force_ind) > 0:
                    all_column_names, force_bool = zip(*settings_list) #split pairs
                    column_names =  []
                    forces_to_save = []
                if any(force_bool[0:3]): # if any HF data is to be saved
                    forces_to_save.append(time)
                    column_names.append('Time (s)')
                    for i in range(4):
                        if force_bool[i] == True:
                            forces_to_save.append(pooled_force_data[i].data)
                            column_names.append(all_column_names[i])
                if any(force_bool[4:7]): # if any downsampled data is to be saved
                    forces_to_save.append(time_downsamp)
                    column_names.append('Time downsamp (s)')               
                    for j in range(4,8):
                        if force_bool[j] == True:
                            forces_to_save.append(pooled_force_data[j].data)
                            column_names.append(all_column_names[j]) 
                else:
                    print('No data channels were selected')
                    return
                
                #save metadata to a .txt file
                exp_desc = temp_file.description
                filename_meta = filename_no_extension + "_desc.txt"
                if len(exp_desc) > 0 and not os.path.exists(filename_meta):
                    meta_file = open(filename_meta, "w")
                    meta_file.write(exp_desc)
                    meta_file.close()
                
                ## preallocate pd DataFrame with #indices = max length of force data
                grouped_data = pd.DataFrame(index=range(len(forces_to_save[0])))
                for i in range(len(forces_to_save)):
                    # if shorter than index length, append NaN: 
                    if len(forces_to_save[i]) < len(forces_to_save[0]):
                        forces_to_save[i] = np.append(forces_to_save[i], np.repeat(np.nan, len(forces_to_save[0])-len(forces_to_save[i])))
                    #transpose force arrays for saving in columns
                    forces_to_save[i]=forces_to_save[i].reshape(-1,1)
                    grouped_data.insert(i, column_names[i], forces_to_save[i])
                
                grouped_data = grouped_data.fillna('') # fill NaN with ''
                filename_force = filename_no_extension + "_force.csv"
                grouped_data.to_csv(filename_force, sep=',', index=False, header=True)
                return
            
            def extractAndDestroyForceOrigin():
                settings_list = [('Force1x HF', force1xHF_var.get()), 
                     ('Force1y HF', force1yHF_var.get()), 
                     ('Force2x HF', force2xHF_var.get()), 
                     ('Force2y HF', force2yHF_var.get()),                      
                     ('Force1x downsamp', force1x_downsamp_var.get()), 
                     ('Force1y downsamp', force1y_downsamp_var.get()), 
                     ('Force2x downsamp', force2x_downsamp_var.get()), 
                     ('Force2y downsamp',force2y_downsamp_var.get())]
                extractForceCommand(settings_list)
                forceRoot.destroy()
                return
            
            #Pop Up Menu
            forceRoot = tk.Toplevel()
            forceRoot.wm_title("Force Extract Option")
            
            # Initialize variables for save option checkboxes
            force1xHF_var = tk.BooleanVar()
            force1yHF_var = tk.BooleanVar()
            force2xHF_var = tk.BooleanVar()
            force2yHF_var = tk.BooleanVar()    
            force1x_downsamp_var = tk.BooleanVar()
            force1x_downsamp_var.set(True) # most common save option, so selected by default
            force1y_downsamp_var = tk.BooleanVar()
            force2x_downsamp_var = tk.BooleanVar()
            force2y_downsamp_var = tk.BooleanVar()
            
            # Establish checkboxes for each save option
            tk.ttk.Checkbutton(forceRoot, text="Force1x_HF",variable=force1xHF_var).grid(row=1, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force1y_HF",variable=force1yHF_var).grid(row=2, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force2x_HF",variable=force2xHF_var).grid(row=3, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force2y_HF",variable=force2yHF_var).grid(row=4, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force1x_downsamp",variable=force1x_downsamp_var).grid(row=5, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force1y_downsamp",variable=force1y_downsamp_var).grid(row=6, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force2x_downsamp",variable=force2x_downsamp_var).grid(row=7, sticky="w")
            tk.ttk.Checkbutton(forceRoot, text="Force2y_downsamp",variable=force2y_downsamp_var).grid(row=8, sticky="w")
            tk.ttk.Button(forceRoot, text='Extract data', command=extractAndDestroyForceOrigin).grid(row=10, sticky="w", pady=8)
            return        
        
        """
        This function extracts the photon counts of all three channels to a .xlsx file.
        """
        def extract_photon_counts(event):
            three_color_photon_data = saved_color_data
            
            if extractPhotonCountsOpt.get() == ".xlsx":
                if len(three_color_photon_data.shape) > 3:
                    numFrame = int(scaleForStack.get())
                    three_color_photon_data = three_color_photon_data[numFrame-1,:,:,:]
                    print("For scan files with multiple frames, only output the current frame number")
                    writer = pd.ExcelWriter((((typePulldown.get()).replace(" ","_")).replace("-","_"))+"_frame_"+str(numFrame)+"_extracted_photon_counts.xlsx")
                else:
                    writer = pd.ExcelWriter((((typePulldown.get()).replace(" ","_")).replace("-","_"))+"_extracted_photon_counts.xlsx")
                
                df_red = pd.DataFrame(three_color_photon_data[:,:,0],dtype=int)
                df_green = pd.DataFrame(three_color_photon_data[:,:,1],dtype=int)
                df_blue = pd.DataFrame(three_color_photon_data[:,:,2],dtype=int)
                
                
                df_red.to_excel(writer, sheet_name='Red Photon Counts',header=False,index=False)
                df_green.to_excel(writer, sheet_name='Green Photon Counts',header=False,index=False)
                df_blue.to_excel(writer, sheet_name='Blue Photon Counts',header=False,index=False)
                writer.save()
            elif extractPhotonCountsOpt.get() == ".tif":
                tiff.imsave((((typePulldown.get()).replace(" ","_")).replace("-","_"))+".tif",three_color_photon_data)
            return
        
        """
        The three functions below are activated upon the "Extract Line Scans" button
        being hit. The user can click points of interest on the RGB plot and when the
        user's cursor exits the window the line scans are extracted to a .xlsx file.
        The option for vertical or horizontal lines is shown in the GUI as well.
        """
        def leave_figure_after_extract_line_scans(event):
            buttonToExtractLineScans['state'] = tk.NORMAL
            comboboxForLineScan['state'] = tk.NORMAL
            figureToSave.canvas.mpl_disconnect(leave_figure_event)
            figureToSave.canvas.mpl_disconnect(add_point_to_line_extract)
            
            writer = pd.ExcelWriter((((typePulldown.get()).replace(" ","_")).replace("-","_"))+"_extracted_line_scans_"+ comboboxForLineScan.get().replace(".","") + ".xlsx")
            dict_for_red_line_scans = {}
            dict_for_green_line_scans = {}
            dict_for_blue_line_scans = {}
            
            h5_filepath = directoryPulldown.get()
            
            stringType = typePulldown.get().split("-")
            typePointer = "-".join(stringType[1:])
            stringType = stringType[0]

            h5_file = lk.File(h5_filepath)
            
            color_data = saved_color_data
            if stringType == "kymos":
                kymoPointer = h5_file.kymos[typePointer]
                max_y = kymoPointer.pixelsize_um[0] * 1000 * saved_color_data.shape[0] / 1000 #pixel size in nm
                max_x = ((kymoPointer.timestamps[0,1]-kymoPointer.timestamps[0,0]) / 1000000000) * saved_color_data.shape[1] #scan time in s
            elif stringType == "scans":
                scanPointer = h5_file.scans[typePointer]
                max_y = scanPointer.pixelsize_um[0] * 1000 * saved_color_data.shape[1] / 1000 #pixel size in nm
                max_x = scanPointer.scan_width_um[0]
                if len(color_data.shape) > 3:
                    color_data = saved_color_data[int(scaleForStack.get())-1,:,:,:]
            else:
                print('fdcurves do not have color data, and therefore no line scans were extracted')
                return
            
            #time/position is used in the most common case - for analyzing kymographs
            if comboboxForLineScan.get() == 'Vert.':
                list_of_x_vals = [item[0] for item in line_list_for_line_extract]
                dict_for_red_line_scans["Y Range"] = np.linspace(0,max_y,num=color_data.shape[0])
                dict_for_green_line_scans["Y Range"] = np.linspace(0,max_y,num=color_data.shape[0])
                dict_for_blue_line_scans["Y Range"] = np.linspace(0,max_y,num=color_data.shape[0])
                for counter, x_pos in enumerate(list_of_x_vals):
                    index = int( (x_pos / max_x) *color_data.shape[1] )
                    red_signal = color_data[:,index,0]
                    green_signal = color_data[:,index,1]
                    blue_signal = color_data[:,index,2]
                    
                    fractional_x_pos = round(x_pos/max_x*100,2)
                    dict_for_red_line_scans[f"Red Scan #{counter+1}\nFractional Position - {fractional_x_pos}% of X-Axis"] = red_signal
                    dict_for_green_line_scans[f"Green Scan #{counter+1}\nFractional Position - {fractional_x_pos}% of X-Axis"] = green_signal
                    dict_for_blue_line_scans[f"Blue Scan #{counter+1}\nFractional Position - {fractional_x_pos}% of X-Axis"] = blue_signal
            else:
                list_of_pos_vals = [item[1] for item in line_list_for_line_extract]
                dict_for_red_line_scans["X Range"] = np.linspace(0,max_x,num=color_data.shape[1])  
                dict_for_green_line_scans["X Range"] = np.linspace(0,max_x,num=color_data.shape[1])
                dict_for_blue_line_scans["X Range"] = np.linspace(0,max_x,num=color_data.shape[1])
                for counter, position in enumerate(list_of_pos_vals):
                    index = int( (position / max_y)*color_data.shape[0] )
                    red_signal = color_data[color_data.shape[0]-index,:,0]
                    green_signal = color_data[color_data.shape[0]-index,:,1]
                    blue_signal = color_data[color_data.shape[0]-index,:,2]
                    
                    fractional_y_pos = round(position/max_y*100,2)
                    dict_for_red_line_scans[f"Red Scan #{counter+1}\nFractional Position - {fractional_y_pos}% of Y-Axis"] = red_signal
                    dict_for_green_line_scans[f"Green Scan #{counter+1}\nFractional Position - {fractional_y_pos}% of Y-Axis"] = green_signal
                    dict_for_blue_line_scans[f"Blue Scan #{counter+1}\nFractional Position - {fractional_y_pos}% of Y-Axis"] = blue_signal
            
            data_frame_red = pd.DataFrame.from_dict(dict_for_red_line_scans)
            data_frame_green = pd.DataFrame.from_dict(dict_for_green_line_scans)
            data_frame_blue = pd.DataFrame.from_dict(dict_for_blue_line_scans)
            
            data_frame_red.to_excel(writer,sheet_name="Red",header=True,index=False)
            data_frame_green.to_excel(writer,sheet_name="Green",header=True,index=False)
            data_frame_blue.to_excel(writer,sheet_name="Blue",header=True,index=False)
            writer.save()
            
            line_list_for_line_extract.clear()
            return
        
        def add_point_to_line_extract(event):
            x, y = event.inaxes.transData.inverted().transform((event.x , event.y))
            line_list_for_line_extract.append([x,y])
            return
        
        def extract_line_scans_by_click(event):
            global leave_figure_event
            global click_to_add
            
            print("Click points on the plot that you want to extract the vertical/horizontal lines of. The file will save the .xlsx file once your cursor leaves the figure window.")
            
            buttonToExtractLineScans['state'] = tk.DISABLED
            comboboxForLineScan['state'] = tk.DISABLED
            leave_figure_event = figureToSave.canvas.mpl_connect('figure_leave_event', leave_figure_after_extract_line_scans)
            click_to_add = figureToSave.canvas.mpl_connect('button_press_event', add_point_to_line_extract)
            return
        
        """
        This functions extracts the color data of the file type whenever the kymo/scan
        pointer in the file is changed and auto-fills the fileNameToSave entry to the 
        key for the specific file object.
        """
        def preload_RGB_and_changeSaveName(event):
            fileNameToSave = typePulldown.get()
            entrySaveFile.delete(0, "end")
            entrySaveFile.insert(0,fileNameToSave)
            
            global saved_color_data
            h5_filepath = directoryPulldown.get()
            splitTypePulldown = (typePulldown.get()).split('-')[0]

            if splitTypePulldown == "fdcurves":
                saved_color_data = 0   
            elif splitTypePulldown == "kymos":
                h5_filepath = directoryPulldown.get()
                h5file = lk.File(h5_filepath)
                saved_color_data = h5file.kymos["-".join(typePulldown.get().split('-')[1:])].rgb_image
            else:
                h5_filepath = directoryPulldown.get()
                h5file = lk.File(h5_filepath)
                saved_color_data = h5file.scans["-".join(typePulldown.get().split('-')[1:])].rgb_image
            return
        
        """
        This function allows for the Ctrl+C shortcut to copy Non-RGB data to the clipboard
        for other analysis.
        
        *Using high frequency,not-downsampled data can cause the program to freeze
        while it extracts/copies all of the data.
        """
        def copy_data(event):
            x_data, y_data, y_data_string = generateFigure(extract_other_data="yes")
            
            y_data_string = (typePulldown.get()).split('-')[1] + "\n" +y_data_string
            
            if comboboxForNonRGB.get() == "Force-Distance":
                x_data_string = "Distance (uM)"
            else:
                x_data_string = "Time (s)"
            
            df_to_copy = pd.DataFrame({x_data_string:x_data,y_data_string:y_data})
            df_to_copy.to_clipboard(sep='\t',index=False)
            return
        
        """
        This function is used to provide visualization of the frame number since the 
        tk.ttk version of the slider does not have a natural way of displaying that
        number.
        """
        def print_scale_output(event):
            entryScaleSlider.configure(state='normal')
            entryScaleSlider.delete(0,"end")
            entryScaleSlider.insert(0,str(int(scaleForStack.get())))
            entryScaleSlider.configure(state='readonly')
            return
        
        """
        This function opens a separate GUI window to do apply lumicks.pylake algorithms
        to track lines and extract desired properties (time,sub-pixel position, and intensity).
        Information about how the algorithm extracts the data is avaliable at:
            https://lumicks-pylake.readthedocs.io/en/stable/tutorial/kymotracking.html
        
        Most of the design elements of this window were taken from the KymoWidget developed
        by the lumicks team as well.
        """
        def callKymotracker(event):
            stringType = typePulldown.get().split("-")
            typePointer = "".join(stringType[1:])
            stringType = stringType[0]
            if stringType == "kymos": 
                modified_RGB_data = modify_rgb_image(saved_color_data)
                
                kymoTrackerRoot = tk.Tk()
                kymoTrackerRoot.config(bg="gray94")
                kymoTrackerRoot.title(f'KymoTracker -- {directoryPulldown.get()} -- {typePulldown.get()}')
                kymoTracker_gui = KymoTrackerGUI(kymoTrackerRoot,directoryPulldown.get(),typePointer,saved_color_data,modified_RGB_data) #call the Application class
                kymoTrackerRoot.mainloop()
            else:
                print(f"{stringType} file type detected. Only kymograph files are applicable to use in the kymotracker functionality")
            return
        """
        quitButton bound event to destory the tkinter window and exit out of python
        """
        def totalQuit(event):
            master.destroy()
            quit()
            return
        
        self.master = master
        #build the simple GUI
        define_Global_Defaults()

        #add directory system - values to be assigned dynamically later
        frameForFileAccess = tk.ttk.Frame(master)
        frameForFileAccess.grid(row=0,column=2,columnspan=2,rowspan=2,sticky='nw',pady=2,padx=2)
        tk.ttk.Label(frameForFileAccess,text="H5 Files in Directory",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=1,sticky='nw')
        
        directoryPulldown = tk.ttk.Combobox(frameForFileAccess,width=60)
        directoryPulldown.grid(row=1,column=0,columnspan=4,sticky="nw")
        tk.ttk.Label(frameForFileAccess,text="File Components",font=('Helvetica', 10, 'bold')).grid(row=2,column=0,columnspan=1,sticky='w')
        
        typePulldown = tk.ttk.Combobox(frameForFileAccess,width=17)
        typePulldown.grid(row=2,column=1,columnspan=2,sticky="nw")
        #Labels and entries for Plot Title and Image/metadata
        tk.ttk.Label(frameForFileAccess,text="Plot Title",font=('Helvetica', 10, 'bold')).grid(row=3,column=0,sticky='nw')
        tk.ttk.Label(frameForFileAccess,text="File Name",font=('Helvetica', 10, 'bold')).grid(row=4,column=0,sticky='nw')
        entryPlotTitle = tk.ttk.Entry(frameForFileAccess)
        entryPlotTitle.grid(row=3,column=1,columnspan=2,sticky='nw')
        entrySaveFile = tk.ttk.Entry(frameForFileAccess)
        entrySaveFile.grid(row=4,column=1,columnspan=2,sticky='nw')
        
        imageFormatOption = tk.ttk.Combobox(frameForFileAccess,values=['PNG','PDF','TIFF','JPEG','SVG'],width=5)
        imageFormatOption.grid(row=4,column=3,padx=5,sticky="w")
        imageFormatOption.current('0')    
        
        ###add color scaling options
        #Labels
        frameForColorOpt = tk.ttk.Frame(master)
        frameForColorOpt.grid(row=3,column=3,rowspan=1,columnspan=1,sticky='nw',pady=2)
        tk.ttk.Label(frameForColorOpt,text="Photon Count Multiplier",font=('Helvetica', 10, 'bold')).grid(row=0 ,column=0,columnspan=2,sticky='nw',pady=2)
        #RGBscalingFrame = tk.ttk.Frame(frameForColorOpt)
        #RGBscalingFrame.grid(row=1,column=1)
        tk.ttk.Label(frameForColorOpt,text="Red Scalar").grid(row=1, column=0,sticky='w')
        tk.ttk.Label(frameForColorOpt,text="Green Scalar").grid(row=2, column=0,sticky='w')
        tk.ttk.Label(frameForColorOpt,text="Blue Scalar").grid(row=3, column=0,sticky='w')
        tk.ttk.Label(frameForColorOpt,text="Brightness").grid(row=4, column=0,sticky='w')
            
        #Entries for color analysis
        entryRed = tk.ttk.Entry(frameForColorOpt,width=6)
        entryRed.insert(0,str(defaultDict['red']))
        entryRed.grid(row=1,column=1)
        entryGreen = tk.ttk.Entry(frameForColorOpt,width=6)
        entryGreen.insert(0,str(defaultDict['green']))
        entryGreen.grid(row=2,column=1)
        entryBlue = tk.ttk.Entry(frameForColorOpt,width=6)
        entryBlue.insert(0,str(defaultDict['blue']))
        entryBlue.grid(row=3,column=1)
        entryBrightness = tk.ttk.Entry(frameForColorOpt,width=6)
        entryBrightness.grid(row=4,column=1)
        entryBrightness.insert(0,'0')
        tk.ttk.Label(frameForColorOpt,text="Greyscale Options").grid(row=5 ,column=0,columnspan=1,sticky='nw')
        grayscaleOpt = tk.ttk.Combobox(frameForColorOpt,values=['No','R','G','B'],width=3)
        grayscaleOpt.current('0')
        grayscaleOpt.grid(row=5, column=1)
        
        buttonToExtractPhotonCounts = tk.ttk.Button(frameForColorOpt,text="Extract Photon Counts")
        buttonToExtractPhotonCounts.grid(row=7,column=0,columnspan=1,pady=6)
        extractPhotonCountsOpt = tk.ttk.Combobox(frameForColorOpt,values=[".tif",".xlsx"],width=5)
        extractPhotonCountsOpt.grid(row=7,column=1,columnspan=1,padx=2,pady=2)
        extractPhotonCountsOpt.current('0')
        buttonToExtractLineScans = tk.ttk.Button(frameForColorOpt,text="Extract Line Scans")
        buttonToExtractLineScans.grid(row=6,column=0,columnspan=1,pady=2)
        comboboxForLineScan = tk.ttk.Combobox(frameForColorOpt,values=['Vert.','Horiz.'],width=5)
        comboboxForLineScan.grid(row=6,column=1,columnspan=1,padx=2,pady=2)
        comboboxForLineScan.current('0')
        
            
        ###add force options
        frameForForceOptions = tk.ttk.Frame(master)
        frameForForceOptions.grid(row=2,column=3,rowspan=1,columnspan=1,sticky='nw',pady=2)
        tk.ttk.Label(frameForForceOptions,text="Force Options",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='w',pady=2)
        tk.ttk.Label(frameForForceOptions,text="Sampling(Hz)",).grid(row=1, column=0,sticky='w')
        tk.ttk.Label(frameForForceOptions,text="Channel",).grid(row=2, column=0,sticky='w')
        
        entryDownSample = tk.ttk.Entry(frameForForceOptions,width=9)
        entryDownSample.grid(row=1,column=1,sticky='e')
        entryDownSample.insert(0,'100')
        forceChannelPulldown = tk.ttk.Combobox(frameForForceOptions,values=['1x','1y','1-Both','2x','2y','2-Both','3x','3y','3-Both','4x','4y','4-Both'],width=6)
        forceChannelPulldown.grid(row=2,column=1,columnspan=2,sticky='e')
        forceChannelPulldown.current('3')
        checkValueDownsampleOpt = tk.IntVar(value=1)
        tk.ttk.Label(frameForForceOptions,text="Downsample?").grid(row=3,column=0,columnspan=1,sticky="nw") #,sticky="nw"
        forceOptions = tk.ttk.Checkbutton(frameForForceOptions,var=checkValueDownsampleOpt)
        forceOptions.grid(row=3,column=1,columnspan=2)
            
        ###add plotting options
        pady_plotting_opts = 3
        padx_plotting_opts = 2
        frameForPlotting = tk.ttk.Frame(master)
        frameForPlotting.grid(row=2,column=2,rowspan=1,columnspan=1,sticky='nw',pady=2)
        tk.ttk.Label(frameForPlotting,text="Plotting Options",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='nw',pady=2)
        tk.ttk.Label(frameForPlotting,text="Which Plot(s)").grid(row=1,column=0,sticky='nw',padx=padx_plotting_opts,pady=pady_plotting_opts)
        tk.ttk.Label(frameForPlotting,text="Non-RGB Plot").grid(row=2,column=0,sticky='nw',padx=padx_plotting_opts,pady=pady_plotting_opts)
        tk.ttk.Label(frameForPlotting,text="Distance Data").grid(row=3,column=0,sticky='nw',padx=padx_plotting_opts,pady=pady_plotting_opts)
        tk.ttk.Label(frameForPlotting,text="Trap Position").grid(row=4,column=0,sticky='nw',padx=padx_plotting_opts,pady=pady_plotting_opts)
        plottingOpt = tk.ttk.Combobox(frameForPlotting,values=['Both','RGB Only','Non-RGB Only'],width=15)
        plottingOpt.grid(row=1,column=1,padx=padx_plotting_opts,pady=pady_plotting_opts)
        plottingOpt.current('1')
        comboboxForNonRGB = tk.ttk.Combobox(frameForPlotting,values=['Force-Time','Force-Distance','Trap Pos.-Time'],width=15)
        comboboxForNonRGB.grid(row=2,column=1,padx=padx_plotting_opts,pady=pady_plotting_opts,sticky="nw")
        comboboxForNonRGB.current('0')
        whichDistanceValue = tk.ttk.Combobox(frameForPlotting,width=15)
        whichDistanceValue.grid(row=3,column=1,padx=padx_plotting_opts,pady=pady_plotting_opts)
        whichTrapPosValue = tk.ttk.Combobox(frameForPlotting,width=15)
        whichTrapPosValue.grid(row=4,column=1,padx=padx_plotting_opts,pady=pady_plotting_opts)
        
        
        ###axis scaling interface
        frameForAxis = tk.ttk.Frame(master)
        frameForAxis.grid(row=3,column=2,rowspan=1,columnspan=1,sticky='nw',pady=2)
        tk.ttk.Label(frameForAxis,text="Axis Options",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=4,sticky='nw',pady=2)
        tk.ttk.Label(frameForAxis,text="Time").grid(row=1,column=0,sticky='w')
        tk.ttk.Label(frameForAxis,text="Distance").grid(row=2,column=0,sticky='nw')
        tk.ttk.Label(frameForAxis,text="Y RGB").grid(row=3,column=0,sticky='nw')
        tk.ttk.Label(frameForAxis,text="Y Force").grid(row=4,column=0,sticky='nw')
        tk.ttk.Label(frameForAxis,text="Scan Width").grid(row=5,column=0,sticky='nw')
        tk.ttk.Label(frameForAxis,text="Trap Position").grid(row=6,column=0,sticky='nw')
        
        axisEntryWidth = 6
        padAxisEntry = 2
        #Minimum axis entries - set by default to zero
        entryTimeMin = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryTimeMin.grid(row=1,column=1,padx=padAxisEntry)
        entryTimeMin.insert(0,'0')
        entryDistMin = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryDistMin.grid(row=2,column=1,padx=padAxisEntry)
        entryDistMin.insert(0,'0')
        entryYRGBMin = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryYRGBMin.grid(row=3,column=1,padx=padAxisEntry)
        entryYRGBMin.insert(0,'0')
        entryYForceMin = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryYForceMin.grid(row=4, column=1,padx=padAxisEntry)
        entryYForceMin.insert(0,'0')
        entryScanWidthMin = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryScanWidthMin.grid(row=5,column=1,padx=padAxisEntry)
        entryScanWidthMin.insert(0,'0')
        entryTrapPosMin = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryTrapPosMin.grid(row=6,column=1,padx=padAxisEntry)
        entryTrapPosMin.insert(0,'0')
        
        #Maximum axis entries - not defined by any default since it depends on the file
        entryTimeMax = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryTimeMax.grid(row=1,column=2,padx=padAxisEntry)
        entryTimeMax.insert(0,'-')
        entryDistMax = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryDistMax.grid(row=2,column=2,padx=padAxisEntry)
        entryDistMax.insert(0,'-')
        entryYRGBMax = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryYRGBMax.grid(row=3,column=2,padx=padAxisEntry)
        entryYRGBMax.insert(0,'-')
        entryYForceMax = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryYForceMax.grid(row=4, column=2,padx=padAxisEntry)
        entryYForceMax.insert(0,'-')
        entryScanWidthMax = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryScanWidthMax.grid(row=5,column=2,padx=padAxisEntry)
        entryScanWidthMax.insert(0,'-')
        entryTrapPosMax = tk.ttk.Entry(frameForAxis,width=axisEntryWidth)
        entryTrapPosMax.grid(row=6,column=2,padx=padAxisEntry)
        entryTrapPosMax.insert(0,'-')
        
        aspectOptionVar = tk.BooleanVar()
        aspectOptionVar.set(True)
        tk.ttk.Label(frameForAxis,text="Re-Scale RGB to Window").grid(row=7,column=0,columnspan=2,sticky='w')
        tk.ttk.Checkbutton(frameForAxis,variable=aspectOptionVar).grid(row=7, column=2)
            
        buttonWidth = 25
        #Buttons to link to commands
        buttonFrame = tk.ttk.Frame(master)
        buttonFrame.grid(row=0,column=4,rowspan=20,sticky="nw",pady=2)
        updateH5File= tk.ttk.Button(frameForFileAccess,text="Change Folder")
        updateH5File.grid(row=0,column=0,columnspan=5,pady=1,sticky="e")
        updatePlot = tk.ttk.Button(buttonFrame,text="Draw Plot",width=buttonWidth)
        updatePlot.pack(side="top",padx=4,pady=2)
        exportForceButton = tk.ttk.Button(buttonFrame,text="Export Force",width=buttonWidth)
        exportForceButton.pack(side="top",padx=4,pady=2)
        exportImageButton = tk.ttk.Button(buttonFrame,text="Export Image No Axis",width=buttonWidth)
        exportImageButton.pack(side="top",padx=4,pady=2)
        exportForImageJ = tk.ttk.Button(buttonFrame,text="Export ImageJ Montage",width=buttonWidth)
        exportForImageJ.pack(side="top",padx=4,pady=2)
        openKymotrackerButton = tk.ttk.Button(buttonFrame,text="Open KymoTracker",width=buttonWidth)
        openKymotrackerButton.pack(side="top",padx=4,pady=2)
        saveImageButton = tk.ttk.Button(buttonFrame,text="Save GUI Image",width=buttonWidth)
        saveImageButton.pack(side="top",padx=4,pady=2)
        quitButton = tk.ttk.Button(buttonFrame,text="Quit?",width=buttonWidth) #button to quit Tkinter GUI
        quitButton.pack(side="top",padx=4,pady=2)
        tk.ttk.Label(buttonFrame,text="Keyboard Shortcuts:",font=('Helvetica', 10, 'bold'),justify="left").pack(side="top",anchor="w",padx=4)
        tk.ttk.Label(buttonFrame,text="Enter - Build Plot\nCtrl+O - Change Directory\nCtrl+C - Copy Data to Clipboard\nCtrl+R - Extract Photon Counts\nCtrl+S - Save GUI Image\nCtrl+K - Open KymoTracker\nEsc - Quit the GUI",justify="left",font=('Helvetica', 8)).pack(side="top",anchor="nw",padx=4)
        
        #Inputs and Labels for metadata
        frameForMetadata = tk.ttk.Frame(master)
        frameForMetadata.grid(row=4,column=2,rowspan=1,columnspan=1,sticky='nw',pady=2)
        tk.ttk.Label(frameForMetadata,text="Metadata",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='nw',pady=2)
        
        lastPlottedFig = ['','','','','','','']
        
        #Inputs and labels to configure plotting through a multiple image scan
        frameForSlider = tk.ttk.Frame(master)
        frameForSlider.grid(row=4,column=3,rowspan=1,columnspan=1,sticky="nw",pady=2)
        tk.ttk.Label(frameForSlider,text="Scan Image Frame",font=('Helvetica', 10, 'bold'),justify="left").grid(row=0,column=0,sticky="nw",pady=2)
        
        sliderIntVar = tk.IntVar()
        sliderIntVar.set(1)
        scaleForStack = tk.ttk.Scale(frameForSlider,from_=1,to=10,orient=tk.HORIZONTAL,command=print_scale_output)
        scaleForStack.grid(row=2,column=0,columnspan=3,pady=2)
        
        tk.ttk.Label(frameForSlider,text="Frame Number: ").grid(row=3,column=0,columnspan=1,sticky="nw")
        entryScaleSlider = tk.ttk.Entry(frameForSlider,width=3,justify=tk.CENTER)
        entryScaleSlider.insert(0,str(scaleForStack.get()))
        entryScaleSlider.grid(row=3,column=1,pady=2,sticky="nw")
        entryScaleSlider.configure(state='readonly')
        
        multiScanPlotOpt = tk.IntVar(value=1)
        tk.ttk.Label(frameForSlider,text="Highlight Scan Range?").grid(row=4,column=0,columnspan=1,sticky="nw")
        plotScanTime = tk.ttk.Checkbutton(frameForSlider,var=multiScanPlotOpt)
        #plotScanTime.select()
        plotScanTime.grid(row=4,column=1,pady=2,sticky="e")
                
        #bind the proper buttons to the desired functions
        quitButton.bind("<ButtonRelease-1>",totalQuit)
        updateH5File.bind("<ButtonRelease-1>",changeH5FileDir)
        exportForceButton.bind("<ButtonRelease-1>",extractForceMethod)
        exportImageButton.bind("<ButtonRelease-1>",extractImageCTrap)
        exportForImageJ.bind("<ButtonRelease-1>",extractImageImageJ)
        directoryPulldown.bind("<<ComboboxSelected>>",changeFileComponents)
        updatePlot.bind("<ButtonRelease-1>",buildPlot)
        typePulldown.bind("<<ComboboxSelected>>",preload_RGB_and_changeSaveName)
        openKymotrackerButton.bind("<ButtonRelease-1>",callKymotracker)
        buttonToExtractPhotonCounts.bind("<ButtonRelease-1>",extract_photon_counts)
        saveImageButton.bind("<ButtonRelease-1>",saveFigure)
        buttonToExtractLineScans.bind("<ButtonPress-1>",extract_line_scans_by_click)
        
        #bind keyboard shortcuts
        master.bind("<Control-c>",copy_data)
        master.bind("<Control-C>",copy_data)
        master.bind("<Control-r>",extract_photon_counts)
        master.bind("<Control-R>",extract_photon_counts)
        master.bind("<Control-o>",changeH5FileDir)
        master.bind("<Control-O>",changeH5FileDir)
        master.bind("<Control-s>",saveFigure)
        master.bind("<Control-S>",saveFigure)
        master.bind("<Control-k>",callKymotracker)
        master.bind("<Control-K>",callKymotracker)
        master.bind("<Return>",buildPlot)
        master.bind("<Escape>",totalQuit)
    
#call commands to open the gui class and loop continuosly
root = tk.Tk()
root.title('C-TrapVis v1.0.1')
root.config(bg="gray94")
my_gui = CTrapGUI(root) #call the Application class
root.mainloop()