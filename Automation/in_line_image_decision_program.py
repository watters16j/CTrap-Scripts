import bluelake
from bluelake import trap1, trap2, microstage, shutters, fluidics, pause, power, timeline, reset_force, confocal, excitation, force_feedback, mirror1
import time
import datetime
from winsound import Beep, MessageBeep
from skimage.measure import label, regionprops
from skimage.morphology import closing, remove_small_objects
from math import sqrt
from tifffile import imsave
import numpy as np

# user defined inputs for experiment image processing
def general_experiment_protocol(scan_method, kymograph_method):
	time_of_report = time.strftime('%X').replace(":",'_')
	timeline.mark_begin(time_of_report)

	# experiment is set up as a try finally loop so that if an error occurs it will still end the marker file
	#   Denoting marker files is ideal because it keeps the kymograph image paired with the corresponding FD curve
	try:
		# generate a 2D scan image to see if a polymerase is loaded
		time_pre_scan = timeline.current_time
		pause(0.1)
		scan_confocal_scan(scan_method)
		pause(0.1)
		time_post_scan = timeline.current_time

		# reconstruct the scan image inline and apply processing logic
		outputted_decision = reconstruct_scan_and_make_imaging_decision(time_pre_scan,time_post_scan,time_of_report)
		# option for imaging processing --> 'yes' to use the automation, 'no' to manually curate or just image everything
		opt_to_use_auto_imaging_processing = 'yes'


		if opt_to_use_auto_imaging_processing == 'yes':
			if outputted_decision == 'not single pol':
				print('Imaging processing being applied and this scan did not meet the necessary criteria --> restating loop now')
				timeline.mark_end(export=False)
				return #ending the program before the tether is kymograph imaged

		#buffer flow in all channels to refresh protein in tube
		#start_flow_all_channels(0.1)
		#pause(5)
		#stop_flow()

		# experiment imaging
		print(f'Starting Kymograph Imaging - {outputted_decision}')
		logical_if_tether_broke = kymo_confocal_scan(kymograph_method,time_to_image=500) #starting current image collection settings, timing is in seconds
		print('Kymograph done')

	finally:
		timeline.mark_end(export=False) #it can autoexport the files if you want if you change the export variable to True
	
	#wash the channels again briefly so beads do not stick
	start_flow_all_channels(0.1)
	shutters.clear(1,2,delay_ms = 200)
	pause(5)
	stop_flow()
	return

def reconstruct_scan_and_make_imaging_decision(point1,point2,export_string):
	def reconstruct_scan_image(time1,time2):
		def process_timestamps_and_output_width(time_one,time_two):
			#take the input times and find the corect time points to reconstruct
			confocal_diagnostics_pointer = timeline['Confocal diagnostics']['Frame clock'][time_one:time_two]
			output_diagnostics = confocal_diagnostics_pointer.data
			diagnostics_time = confocal_diagnostics_pointer.timestamps

		
			#forward search to make sure we are not starting at a middle point 
			# need to add in a portion to automatically detect the lines
			time_to_start = 0
			num_line_scans = 0
			for i in range(0,len(output_diagnostics)-1):
				if output_diagnostics[i] == 0 and output_diagnostics[i+1] == 1:
					if time_to_start == 0:
						#print(f"Index of forward search is {diagnostics_time[i]} nanoseconds")
						time_to_start = diagnostics_time[i]
						break

			#reverse search to make sure we are not ending in the middle of a frame
			time_to_end = 0
			for i in reversed(range(1,len(output_diagnostics))):
				if output_diagnostics[i] == 0 and output_diagnostics[i-1] == 1:
					#print(f"Index of reverse search is {diagnostics_time[i]} nanoseconds")
					time_to_end = diagnostics_time[i]
					break

			return time_to_start, time_to_end

		def calculate_len_line_scan(map_array):
			count_of_pixels_in_line = 0
		
			for index in range(0,len(map_array)-1):
				if map_array[index] == 2:
					count_of_pixels_in_line += 1
				if map_array[index] == 2 and map_array[index+1] == 0:
					break

			return count_of_pixels_in_line

		# account for edge cases where the number of timelines are weird
		t1,t2 = process_timestamps_and_output_width(time1,time2)

		# extract data to reconstruct image --> using Cy5 fluoresence in this processing
		sensor=timeline['Photon count']['Red'][t1:t2].data
		mapping=timeline['Info wave']['Info wave'][t1:t2].data

		num_pixels_in_line = calculate_len_line_scan(mapping)
		#print(f"{num_pixels_in_line} pixels in single line")

		#calculate how many pixels are in this file
		assert len(sensor)==len(mapping)
		total_num_of_pixels = len(mapping[mapping==2]) #this assumption is fine given a fixed area
		#print(f"{total_num_of_pixels} pixels detected")
	
		#define width and heighth
		width = num_pixels_in_line
		height = (total_num_of_pixels / width)
		assert height - int(height) == 0, "This must be zero or else something is going wrong with the timestamp processing"
		height = int(height)

		#indexing to extract the photon count values
		pixels=np.zeros(total_num_of_pixels)
		p=0
		ind=0
		for i in np.arange(len(mapping)):
			if mapping[i]<2:
				p=p+sensor[i]*mapping[i]
			else:
				pixels[ind]=(p+sensor[i])
				p=0
				ind+=1

		# reconstruct the image with the given pixel dimensions
		#   You might need to change the ordering if you scan across the dimension differently (x vs y dimension)
		final_image = pixels.reshape((height,width)) # F order is necessary to fill each vertical line then go to the next vertical line
		return final_image

	# You will need to make changes to this function - specifically in the logical checks portion to make it work for your protein of interest
	def process_image_output_good_or_not(im):
		min_threshold = 1
		max_threshold = 6
		binary_map = np.logical_and(im > min_threshold, im < max_threshold)

		bw = closing(im > min_threshold)
		bw = remove_small_objects(bw,min_size=6)
		labeled_img = label(bw)

		#print(f"{np.max(labeled_img)} features detected")
		props = regionprops(labeled_img,im)

		if len(props) == 0:
			return 'not single pol' #just in case it is a blank scan image
		elif len(props) == 1:
			index_of_closest_point = 0
		else: # multiple regions were found --> find the region that is closest to the middle
			index_of_closest_point = 0
			center_point = [int(im.shape[0]/2),int(im.shape[1]/2)]
			closest_distance=  sqrt(im.shape[0]**2 + im.shape[1]**2) #max distance corner to corner
			count = 0
			for region in props:
				location_region = [region['Centroid'][0],region['Centroid'][1]]
				distance_from_center = sqrt((region['Centroid'][0]-center_point[0])**2 +(region['Centroid'][1]-center_point[1])**2)
				if distance_from_center < closest_distance:
					index_of_closest_point = count
					closest_distance = distance_from_center
				count += 1

		# logical checks based on intensity, area, and distance from the center
		# change this is you want to make the pol detection more or less stringent 
		distance_threshold_pixels = 20
		area_threshold = 60
		max_intensity_threshold = 15
	
		if closest_distance > distance_threshold_pixels or props[index_of_closest_point]['area'] > area_threshold or props[index_of_closest_point]['max_intensity'] > max_intensity_threshold:
			if closest_distance > distance_threshold_pixels:
				print(f"Not a single pol because the closest distance was {closest_distance} when it needed to be less than {distance_threshold_pixels}")
			if props[index_of_closest_point]['area'] > area_threshold:
				print(f"Not a single pol because the foci area was {props[index_of_closest_point]['area']} when it needed to be less than {area_threshold}")
			if props[index_of_closest_point]['max_intensity'] > max_intensity_threshold:
				print(f"Not a single pol because the max intensity was {props[index_of_closest_point]['max_intensity']} when it needed to be less than {max_intensity_threshold}")
			return 'not single pol' # something was detected but it does not meet the criteria of a good single polymerase
		else:
			return 'correct single pol'
	
	#call the respective functions
	reconstructed_img = reconstruct_scan_image(point1,point2)
	outputted_decision = process_image_output_good_or_not(reconstructed_img)

	# save the file as a .tif for future analysis
	try:
		imsave(export_string+" "+outputted_decision+".tif",reconstructed_img)
	except:
		print('Image output did not work')\

	print(f'Image processing decision --> {outputted_decision}')
	return outputted_decision

# for this script take an image
def scan_confocal_scan(method):
	py_auto_gui_click_and_move(click_image_scan) #ensure it is in kymo mode --> important because bluelake does not remember this in the imaging presets
	
	# check if the method exists and if it doesnt just perform a kymograph with the standard settings
	try:
		confocal.start_scan(method)              
	except:
		confocal.start_scan() # start the active configuration
		print('Specified confocal method is not contained in this user profile -- starting active configuration')

	# wait until the confocal is done with the scan (if your preset is in continuous mode it will get stuck here!)
	while confocal.is_scanning:
		pause(1)   
	
	confocal.abort_scan() 
	return

# for this script to take a kymograph image --> need to customize by the parameters of your experiment
def kymo_confocal_scan(method,time_to_image):
	py_auto_gui_click_and_move(click_kymo_button) #ensure it is in kymo mode --> important because bluelake does not remember this in the imaging presets
	check_confocal_time = time.time()
	confocal.start_scan(method)

	pause(10)         
	# move tether into a new channel
	microstage.move_to(name_junction_channel)
	microstage.move_to(name_ch4)	  
	#loop to keep scanning until we reach the end
	
	
	while confocal.is_scanning and (time.time() - check_confocal_time) < time_to_image: #checking to make sure the user hasn't ended the scan manually
		pause(1)

	confocal.abort_scan()
	return False

# This function executes the imaging commmands and moves the stage according to users instructions. Must change the code in the function to change the imaging method
general_experiment_protocol(scan_method,kymograph_method)
print("Experiment done!")
		