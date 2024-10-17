"""
#################################################################################################
BSD 2-Clause License
Copyright (c) 2022, John Watters
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

This script is used to extract the sum of photon counts in a line scan
from chosen areas of a kymograph in a .h5 file. This analysis method is useful
when you want to measure raw intensity values of fluorescent molecules to quantify
how many molecules are present without worrying about potential biases introduced
by using line tracking methods.

The user will be asked to  navigate to the .h5 file of interest and  apply a multiplier
value to photon counts if desired (for display purposes only monitoring the method 
used in CTrapVis.py). The user will then chose how they want to define the area to 
analyze either by [1] Explicitly defined dimensions [2] Drag a window on the kymograph
in a pop-up window. Then the user will be asked to click on other regions to extract
where the click defines the center left point of the box (in the cartoon below,  '-'
is the number of timepoints being analyzed, '|' is the number of pixels being analyzed,
and 'X' marks the spot where the user clicks).


----------
|        |
|        |
X        |
|        |
|        |
----------

Once the points are selected on the plot, then you can exit out of the plot. This will
show you a separate image showing the regions the script will extract. If you used the 
dragging a kymograph window method of defining the area dimensions, then you will see
the original box drawn in blue with the additionally clicked boxes drawn in orange.

If the user is happy with the point selection, they can confirm that through user input
and the script will be extracted the data to a .xlsx file. There are two types of sheets in the
output .xlsx. The first sheet type shows the sum of the columns of each region (at each time
point) and the extracted simple statistics from each region (average and standard 
deviation) for each channel you want to extract. The second sheet records some metadata
that might be useful if you need to remake any plots/redo any analysis. In addition,
an image of the extracted regions is also outputted as a .png for a quick resource 
on which region #'s align with other regions of the plot.

Enjoy!
"""

from matplotlib import pyplot as plt
import numpy as np
import lumicks.pylake as lk
from tkinter import filedialog
from tkinter import *
import os
from pandas import DataFrame as df
from pandas import ExcelWriter
import matplotlib
from math import ceil, floor
import datetime

##############################################################################
# Formatting Parameters
number_padding = 50
rectproperties = dict(facecolor='cyan', edgecolor = 'blue',alpha=0.2, fill=True)

# Output Parameters
opt_to_show_regions_extracted = "yes"

##############################################################################
#Function defintion section
#############################################################################
# function to extract the necessary data from .h5 files
def extract_image_data(filepath):
    h5_file = lk.File(filepath)
    kymo_obj = list(h5_file.kymos)[0]
    
    line_time_seconds = h5_file.kymos[kymo_obj].line_time_seconds
    
    im_g = h5_file.kymos[kymo_obj].green_image
    im_r = h5_file.kymos[kymo_obj].red_image
    im_b = h5_file.kymos[kymo_obj].blue_image
    
    pixel_size = float(h5_file.kymos[kymo_obj].pixelsize_um[0]) * 1000
    
    return im_g, im_r, im_b, line_time_seconds, pixel_size

# get maximum photon counts for rescaling later
def extract_max_counts(list_to_extract_maxes):
    len_list_to_extract = len(list_to_extract_maxes)
    
    if len_list_to_extract == 0:
        print('List is empty!')
        return
    
    list_to_export = []
    
    for i in range(0,len_list_to_extract):
        list_to_export.append(np.amax(list_to_extract_maxes[i]))    
    
    return list_to_export

"""
The next two functions are used in drawing the rectangle as you draw it on the plot
"""
def get_rect_dimensions(eclick,erelease):
    global basic_area
    global rect_one
    global rect_two
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

def plot_on_fig(center_point,edge_color='orange',text_to_label=""):
    patch_to_plot = matplotlib.patches.Rectangle((center_point[1], center_point[0] - ceil(dims_to_extract[0]/2)),dims_to_extract[1],dims_to_extract[0],linewidth=1,edgecolor=edge_color,facecolor='none')
    axRGB.add_patch(patch_to_plot)
    
    if text_to_label != "":
        axRGB.annotate(text_to_label,(center_point[1] + ceil(dims_to_extract[1]/2), center_point[0]), color='white', weight='bold', fontsize=10, ha='center', va='center')
    
    return

def add_point_to_extract(event):
    x, y = event.inaxes.transData.inverted().transform((event.x , event.y))
    list_of_points.append([floor(y),floor(x)])
    #plot_on_fig(list_of_points[-1])
    return

def send_dict_to_excel(dict_obj,average_list,std_list,writer_obj,sheet_name):
    dict_obj[''] = np.array([])
    dict_obj['Region #'] = np.asarray(list_of_region_number)
    dict_obj['# Timepoints'] = np.asarray(list_of_num_time_steps)
    dict_obj['Average'] = np.asarray(average_list)
    dict_obj['Standard Deviation'] = np.asarray(std_list)
    
    #now send dict to excel sheet
    pd_data_frame_dict = df.from_dict(dict_obj,orient='index')
    pd_data_frame_dict = pd_data_frame_dict.transpose()
    pd_data_frame_dict.to_excel(writer,sheet_name=sheet_name,index=False,header=True)
    return

##############################################################################
# Import Data

print("#"*number_padding+'\n')
print('Please navigate to and select the .h5 file of interest')
root = Tk()
root.withdraw()
file_selected = filedialog.askopenfilename()
folder_selected = "/".join(file_selected.split('/')[:-1])
os.chdir(folder_selected)

which_color_option = input("What colors would you like to output?\nExample Inputs ('RG' for just Red and Green or 'RGB' for all):\n")

green_photon_counts, red_photon_counts, blue_photon_counts, line_time, pixel_size_nm = extract_image_data(file_selected)
max_rgb_counts = extract_max_counts([red_photon_counts,green_photon_counts,blue_photon_counts])

time_array = np.linspace(0,green_photon_counts.shape[1] * line_time,green_photon_counts.shape[1])
assert len(time_array) == green_photon_counts.shape[1], 'Time array generation is not working'
assert time_array[0] == 0, 'Making sure the time array starts at zero'

maxTrueTime = line_time * green_photon_counts.shape[1]
maxTrueDist = pixel_size_nm * green_photon_counts.shape[0]

print('File import is complete!')
print("#"*number_padding+'\n')

scaling_opt = input('Would you like to apply a multiplier to the scaling (to mimic CTrapVis)?\nIf no, just hit enter. If yes enter the R-G-B multiplier you want\nEx: 0-25-0 to multiply Green by 25\n')
if scaling_opt != '':
    scaling_opt = list(map(float,scaling_opt.split('-')))
    print(f"\nMultiplying photon counts by:\nRed - {scaling_opt[0]}\nGreen - {scaling_opt[1]}\nBlue - {scaling_opt[2]}\n")


today = datetime.datetime.now()
date_string = today.strftime("_%m_%d_%Y")
##############################################################################
# define initial region to extract
opt_to_define_area = 2
opt_to_define_area = int(input("Would you like to define an area using:\n[1] Explicitly defined dimensions\n[2] Drag a window on the kymograph\nPlease type in '1' or '2' below and hit enter:\n"))

list_of_points = [] #pre-define list of points to start adding to - only will be added to here if they use the drag method
print("#"*number_padding+'\n')
if opt_to_define_area == 1:
    print(f'The dimensions of the kymograph are {int(green_photon_counts.shape[0])} pixels in each line scan and {int(green_photon_counts.shape[1])} time points')
    print('Please enter the dimensions of the area you want to extract\nExample: 11x10 for 11 pixels in the line scan direction at 10 time points')
    print('If the number of pixels in the line scan axis is not odd, 1 will be added to make it odd:')
    dims_to_extract = input('Enter dimension:\n').split('x')
    dims_to_extract  = list(map(int, dims_to_extract))
    
    if dims_to_extract[0] % 2 == 0:
        dims_to_extract[0]= dims_to_extract[0] + 1
    
    assert len(dims_to_extract) == 2,'Dimensions to extract were not entered properly, follow the format "11x10"'
    
if opt_to_define_area == 2:
    fig, axRGB = plt.subplots(nrows=1,ncols=1,constrained_layout=True)
    figManager = plt.get_current_fig_manager()

    if scaling_opt == '':
        axRGB.imshow(np.dstack((red_photon_counts,green_photon_counts,blue_photon_counts)) / int(max(max_rgb_counts)), aspect="auto")
    elif len(scaling_opt) == 3:
        axRGB.imshow(np.dstack((red_photon_counts*scaling_opt[0],green_photon_counts*scaling_opt[1],blue_photon_counts*scaling_opt[2])).astype(int), aspect="auto")
    else:
        print('Import of scaling factor for imaging was done incorrectly. Ending program')
        exit()
    axRGB.axis('off')

    print('Draw a box to define the area to extract pixels from!')
    plt.title('Draw a box to define the area to extract pixels from!')
    draw_temp_rectangle.RS = matplotlib.widgets.RectangleSelector(axRGB, get_rect_dimensions, drawtype='box',rectprops=rectproperties)
    fig.canvas.mpl_connect('button_press_event', draw_temp_rectangle)
    
    plt.show()
    plt.close()
    
    dims_to_extract = [basic_area[1][1]-basic_area[1][0], basic_area[0][1]-basic_area[0][0]]
    
    if dims_to_extract[0] % 2 == 0:
        dims_to_extract[0]= dims_to_extract[0] + 1
    
    list_of_points.append([ceil((basic_area[1][1]+basic_area[1][0])/2),basic_area[0][0]])

print(f'Box dimensions\n' + '-'*number_padding + f'\nPixels {dims_to_extract[0]}\nTimepoints {dims_to_extract[1]}')

##############################################################################
print("#"*number_padding+'\n')
print("Now you will select different regions to analyze!")

user_is_happy = 'no'

while user_is_happy == 'no':
    print("Click on the point you want the box to be centered on")
    fig, axRGB = plt.subplots(nrows=1,ncols=1,constrained_layout=True)
    figManager = plt.get_current_fig_manager()

    click_to_add = fig.canvas.mpl_connect('button_press_event', add_point_to_extract)
    
    if scaling_opt == '':
        axRGB.imshow(np.dstack((red_photon_counts,green_photon_counts,blue_photon_counts)) / int(max(max_rgb_counts)), aspect="auto")
    elif len(scaling_opt) == 3:
        axRGB.imshow(np.dstack((red_photon_counts*scaling_opt[0],green_photon_counts*scaling_opt[1],blue_photon_counts*scaling_opt[2])).astype(int), aspect="auto")
    else:
        print('Import of scaling factor for imaging was done incorrectly. Ending program')
        exit()
    axRGB.axis('off')

    if opt_to_define_area == 2:
        plot_on_fig(list_of_points[0],edge_color='blue')
        
    plt.title('Click on the point(s) you want the box to be centered on')
    plt.show()
    print("Close out of the plot to see the regions selected and confirm those regions work for you!")    
    plt.close()
    
    fig, axRGB = plt.subplots(nrows=1,ncols=1,constrained_layout=True)
    figManager = plt.get_current_fig_manager()
    
    click_to_add = fig.canvas.mpl_connect('button_press_event', add_point_to_extract)
    
    if scaling_opt == '':
        axRGB.imshow(np.dstack((red_photon_counts,green_photon_counts,blue_photon_counts)) / int(max(max_rgb_counts)), aspect="auto")
    elif len(scaling_opt) == 3:
        axRGB.imshow(np.dstack((red_photon_counts*scaling_opt[0],green_photon_counts*scaling_opt[1],blue_photon_counts*scaling_opt[2])).astype(int), aspect="auto")
    else:
        print('Import of scaling factor for imaging was done incorrectly. Ending program')
        exit()
    axRGB.axis('off')
    
    if opt_to_define_area == 2:
        plot_on_fig(list_of_points[0],edge_color='blue')
        
    num_points = len(list_of_points)
    
    if opt_to_define_area == 1:
        for i in range(0,num_points):
            plot_on_fig(list_of_points[i],edge_color='orange')
    elif opt_to_define_area == 2:
        for i in range(1,num_points):
            plot_on_fig(list_of_points[i],edge_color='orange')
    
    plt.title('Showing regions to extract')

    plt.show()
    user_is_happy = input("Do these regions work for you? If not, you can restart the picking process by inputting 'No' here\nInput 'Yes' or 'No:\n").lower()    
    plt.close()
    
    if user_is_happy == 'yes':
        print('Fantastic! Extracting data to .xlsx spreadsheet now')
    elif user_is_happy == 'no':
        print('Better luck next time! Clearing list of points and restarting.')
        if opt_to_define_area == 1:
            list_of_points = []
        elif opt_to_define_area == 2:
            list_of_points = [list_of_points[0]]
    else:
        print('User input could not be determined. Restarting picking process... Please enter "Yes" or "No" next time.')
        user_is_happy = "no"
        if opt_to_define_area == 1:
            list_of_points = []
        elif opt_to_define_area == 2:
            list_of_points = [list_of_points[0]]


#extract data and output to a .xlsx spreadsheet
##############################################################################
len_list_of_points= len(list_of_points)
red_photon_count_dict = {}
green_photon_count_dict = {}
blue_photon_count_dict = {}
array_of_center_points = {}

one_time_array = np.arange(0,dims_to_extract[1]*line_time,line_time)
red_photon_count_dict['Time (s)'] = one_time_array
green_photon_count_dict['Time (s)'] = one_time_array
blue_photon_count_dict['Time (s)'] = one_time_array

array_of_center_points['Filepath'] = file_selected.split('/')[-1].split('.')[0]
array_of_center_points['Number Pixels Each Region'] = dims_to_extract[0]
array_of_center_points['Number Timepoints Each Region'] = dims_to_extract[1]

list_of_region_number = []
list_of_center_pixel = []
list_of_first_timestamp = []
list_of_num_time_steps = []

red_regions_average_list = []
green_regions_average_list = []
blue_regions_average_list = []

red_regions_std_list = []
green_regions_std_list = []
blue_regions_std_list = []

for i in range(0,len_list_of_points):
    #sum_of_trace = np.sum(top_polymerase_photon_counts,axis=0)
    current_point = list_of_points[i]
    
    list_of_center_pixel.append(current_point[0])
    list_of_first_timestamp.append(current_point[1])
    list_of_region_number.append(i+1)
    list_of_num_time_steps.append(dims_to_extract[1])
    
    #(center_point[1], center_point[0] - ceil(dims_to_extract[0]/2)),dims_to_extract[1],dims_to_extract[0]
    green_area_of_analysis = green_photon_counts[current_point[0] - ceil(dims_to_extract[0]/2):current_point[0] + ceil(dims_to_extract[0]/2),current_point[1]:current_point[1]+dims_to_extract[1]]
    red_area_of_analysis = red_photon_counts[current_point[0] - ceil(dims_to_extract[0]/2):current_point[0] + ceil(dims_to_extract[0]/2),current_point[1]:current_point[1]+dims_to_extract[1]]
    blue_area_of_analysis = blue_photon_counts[current_point[0] - ceil(dims_to_extract[0]/2):current_point[0] + ceil(dims_to_extract[0]/2),current_point[1]:current_point[1]+dims_to_extract[1]]
    
    red_summed_lines = np.sum(red_area_of_analysis,axis=0)
    green_summed_lines = np.sum(green_area_of_analysis,axis=0)
    blue_summed_lines = np.sum(blue_area_of_analysis,axis=0)

    assert len(red_summed_lines) == len(green_summed_lines), "Red and green line objects aren't the same length"
    assert len(red_summed_lines) == len(blue_summed_lines), "Red and blue line objects aren't the same length"
    assert len(one_time_array) == len(blue_summed_lines), "Time array and line arrays aren't the same length"
    
    #export summed lines to dictionary
    red_photon_count_dict[f'Region {i+1}'] = red_summed_lines
    green_photon_count_dict[f'Region {i+1}'] = green_summed_lines
    blue_photon_count_dict[f'Region {i+1}'] = blue_summed_lines
    
    red_regions_average_list.append(np.average(red_summed_lines))
    green_regions_average_list.append(np.average(green_summed_lines))
    blue_regions_average_list.append(np.average(blue_summed_lines))
    
    red_regions_std_list.append(np.std(red_summed_lines))
    green_regions_std_list.append(np.std(green_summed_lines))
    blue_regions_std_list.append(np.std(blue_summed_lines))
    
# now append the list of points to the metadata sheet
array_of_center_points['Region #'] = np.asarray(list_of_region_number)
array_of_center_points['Center Pixel'] = np.asarray(list_of_center_pixel)
array_of_center_points['First Frame (in pixels)'] = np.asarray(list_of_first_timestamp)  

# optional argument to output the region mapping
if opt_to_show_regions_extracted.lower() == 'yes':
    fig, axRGB = plt.subplots(nrows=1,ncols=1,constrained_layout=True)
    figManager = plt.get_current_fig_manager()
        
    if scaling_opt == '':
        axRGB.imshow(np.dstack((red_photon_counts,green_photon_counts,blue_photon_counts)) / int(max(max_rgb_counts)), aspect="auto")
    elif len(scaling_opt) == 3:
        axRGB.imshow(np.dstack((red_photon_counts*scaling_opt[0],green_photon_counts*scaling_opt[1],blue_photon_counts*scaling_opt[2])).astype(int), aspect="auto")
    else:
        print('Import of scaling factor for imaging was done incorrectly. Ending program')
        exit()
    axRGB.axis('off')
    
    if opt_to_define_area == 2:
        plot_on_fig(list_of_points[0],edge_color='blue',text_to_label='1')
        
    num_points = len(list_of_points)
    
    if opt_to_define_area == 1:
        for i in range(0,num_points):
            plot_on_fig(list_of_points[i],edge_color='orange',text_to_label=str(i+1))
    elif opt_to_define_area == 2:
        for i in range(1,num_points):
            plot_on_fig(list_of_points[i],edge_color='orange',text_to_label=str(i+1))
    
    plt.savefig((file_selected.split('/')[-1].split('.')[0])+"_extracted_intensity_regions"+ date_string +".png")
    plt.close()

# Now sending the data to a .xlsx spreadsheet
##############################################################################    
# export data to .xlsx sheet
writer = ExcelWriter(str(file_selected.split('/')[-1].split('.')[0])+"_extracted_photon_counts" + date_string + ".xlsx")

if 'R' in which_color_option:
    send_dict_to_excel(red_photon_count_dict,red_regions_average_list,red_regions_std_list,writer,'Red')
if 'G' in which_color_option:
    send_dict_to_excel(green_photon_count_dict,green_regions_average_list,green_regions_std_list,writer,'Green')
if 'B' in which_color_option:
    send_dict_to_excel(blue_photon_count_dict,blue_regions_average_list,blue_regions_std_list,writer,'Blue')

#export region settings
pd_data_frame_dict = df.from_dict(array_of_center_points,orient='index')
pd_data_frame_dict = pd_data_frame_dict.transpose()
pd_data_frame_dict.to_excel(writer,sheet_name='Metadata',index=False,header=True)

#call the final command to save the writer
writer.save()

print('Data succesfully exported!')
