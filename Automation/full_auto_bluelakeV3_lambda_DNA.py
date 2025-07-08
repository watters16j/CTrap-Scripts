import bluelake
from bluelake import trap1, trap2, microstage, shutters, fluidics, pause, power, timeline, reset_force, confocal, excitation_lasers, force_feedback, mirror1
import time
import datetime
from winsound import Beep, MessageBeep

# Maximum number of tethers to try before ending the script
#   If I can connect the computer to the internet safely then I would like to be able to e-mail the user when the maximum number of tethers are reached - but I am not sure if that is possible. 
max_tethers = 50                     # maximum number of tethers to try in total

### variables to change based on the specifics of your DNA tether
# using fractional lengths to target where we should start fishing (fractional_contour_length_fishing)
estimated_contour_length = 16.49 #in microns!
fractional_contour_length_for_fish = 0.65                # 0 --> 1 This variable is the starting point to start fishing (the beads will get sequentially closer together each fishing cycle)
fractional_contour_length_for_fishing_extension = 1.05   # >1 by a little bit - when calling force feedback loop to extend the tether if the bead reaches this distance it will confirm that a tether has not been caught and will start the next fishing cycle
fractional_contour_length_for_extension = 0.74           # 0 --> 1 This variable is called when setting the tether in 
fractional_contour_length_for_FD = 1.5                   # >1 however long you want to pull the tether for the FD curve
distance_for_single_tether_rejection = 10                # distance in microns that if you do not reach when pulling to verify single tether, it will reject it as a multiple tether

# Variables for Fishing
fishing_force_target = 15                                # Fishing force (pN) to target with a force feedback loop to see if a tether has formed
force_to_pull_to_image = 2                                                                              # Force to pull to before imaging
step_size_fishing = 0.1                                                                                 # Amount to shorten the tether each fishing cycle (microns)
fishing_move_speed = 16                                                                                  # Speed to move move trap 1 (right bead in this set-up) back to the minimum distance each cycle
# Automatically calculate the distances using the contour length and fractional distances provdied earlier
fishing_start_point = estimated_contour_length*fractional_contour_length_for_fish                       # starting point for the first fishing cycle (will get closer each fishing cycle)
fishing_distance_to_cycle = estimated_contour_length*fractional_contour_length_for_fishing_extension    # maximum distance which the force feedback should hit before restarting the cycle
extension_start_point_var = estimated_contour_length*fractional_contour_length_for_extension            # starting point to stretch DNA before imaging/making the FD curve post imaging
FD_max_distance = estimated_contour_length*fractional_contour_length_for_FD
# The two lists are the settings the Force Feedback loop during either fishing or extension right before imaging (if you require different settings for the two different modes)
force_feedback_variables_fishing = [3, 0, 0, 2000,False]                                                #[kp, ki, kd, max_step,reversed] parameters list for PID
force_feedback_variables_extension = [2, 0, 0, 1000,False]                                              #[kp, ki, kd, max_step,reversed] parameters list for PID 

# Variables for Bead Catching
#   Why this is important: The bluelake software isnt able to guess where each trap is so you have to manually define a place where you are holding Trap 2 (the bead on the left)
#   For this set-up, I keep Trap 2 (on the LEFT) left of the 50 micron marker
#   If you do not have this set up correctly it will practically never catch beads correctly, but the previous Lumicks solution to this issue (not manually hard-coding in a location)
#        resulted in much slower bead catching
centered_of_x_scan = 23                         # place where the left and right bead should NOT cross! this is useful for the bead picking algorithm
x_pos_to_pull_to_for_bead_catching = 25      # x-distance (as noticed in the waypoints in the blulake software) to where we should move mirror 1 to make sure it doesnt overlay with mirror 2
                                                #   This is important to (i) keep the right bead far away from the center line and (ii) far enough away to not have interference during the force calibration 
bead_threshold_score = 90                       # minimum threshold score for beads (will reject anything < 0.94 --> this is scaled up 100x in the backend which is why this value is 94 instead of 0.94)

# It is best practice to manually catch two beads and then re-align and re-save the confocal areas to ensure that you will get good imaging
# The names  that I have in here are linked to the methods I keep in my experimental profile
scan_method = "Alignment_Scan" #0.1 s pixel time, 2D scan of the beads and the tether in between -> not used by gabby so kept the name the same
kymograph_method = "Kymo" #0.1 s pixel time slice through the middle of the bead (need to calibrate this location before)

#Channel names as in the UI
name_bead_channel = "beads"
name_dna_channel = "DNA"
name_buffer_only = "J1" #making buffer only J2 because that matches what gabby's ctrap set-up is
name_junction_channel = "J1"
name_ch4 = "Ch4"

# Set up variables that will get real time data from the instrument
# match_scoreX refers to the bead matching score based on the template
distance = timeline["Distance"]["Distance 1"]
force = timeline["Force LF"]["Trap 2"]
match_score1 = timeline["Tracking Match Score"]["Bead 1"]
match_score2 = timeline["Tracking Match Score"]["Bead 2"]
bead_1_position_x = timeline["Bead position"]["Bead 1 X"]

# Some additional flow parameters
flow_pressure = 0.1
time_to_wait_for_flow = 1.0 #in seconds
flow_timer_to_shut_off_ch4and_ch5 = 20 #seconds

# user defined inputs for experiment image processing
# for this script take an image
def scan_confocal_scan(method):
    try:
        confocal.start_scan(method)              
    except:
        confocal.start_scan() # start the active configuration
        print('Specified confocal method is not contained in this user profile -- starting active configuration')

    while confocal.is_scanning:
        pause(1)   
    
    confocal.abort_scan() 
    return

def kymo_confocal_scan(method,match_threshold):
    try:
        confocal.start_scan(method)              
    except:
        confocal.start_scan() # start the active configuration
        print('Specified confocal method is not contained in this user profile -- starting active configuration') 
    
    # wait for scan to finish (not required)
    #microstage.move_to(name_ch4)
    pause_counter = 0
    time_to_move = 4                    #seconds
    time_to_incubate_post_move = 400     #seconds
    time_to_image_post_move2 = 10      #seconds

    #summing the different time moves
    time_to_move_from_protein_channel = time_to_move + time_to_incubate_post_move
    time_to_stop_scanning = time_to_move + time_to_incubate_post_move + time_to_image_post_move2
  

    while confocal.is_scanning: #checking to make sure the user hasn't ended the scan manually
        pause(1)
        pause_counter += 1
        if pause_counter == time_to_move:
            microstage.move_to(name_ch4)
            print(f'correct move at {pause_counter} pause counter')
        if pause_counter > time_to_move_from_protein_channel: 
            microstage.move_to(name_buffer_only)
        if pause_counter > time_to_stop_scanning:
            break

    confocal.abort_scan()
    return False #since not pulling too much in this script the tether should always still be there


def experiment_imaging(kymograph_method,match_threshold):
    #Beep(2500,200) insert new tone to play
    print('Starting Kymograph Imaging')
    logical_if_tether_broke = kymo_confocal_scan(kymograph_method,match_threshold) #starting current image collection settings (probably kymograph)
    print('Kymograph done')
    
    return logical_if_tether_broke #this logical is false if we still have a tether, true if there is no tether attached


def throw_if_beads_lost(match_threshold):
    """Raise an exception if we lose the beads."""
    if match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        print("Beads lost, restarting protocol")
        bluelake.force_feedback.enabled= False
        return True
    else:
        return False

def catch_beads(match_threshold, flow_timer, pressure=flow_pressure,initial_opt='no'):
    def determine_bead_pos(current_one_bead_pos,x_threshold):
        if current_one_bead_pos - x_threshold > 0: #bead is on the right
            return "right"
        elif current_one_bead_pos - x_threshold < 0:
            return "left"
        else:
            print("The single bead is right on the center line - erroring out!")
            exit()

    """Starts the flow and attempts to catch beads. Toggles the shutters when match score is too low."""
    start_flow(pressure)
    intial_time_for_flow = time.time()
    if initial_opt != 'yes':
        print(f"Waiting {time_to_wait_for_flow} sec for flow to begin.")
        pause(time_to_wait_for_flow)

    print("Moving to bead channel.")
    mirror1.move_to(x=x_pos_to_pull_to_for_bead_catching,speed=5)
    microstage.move_to(name_bead_channel)

    #clearing any beads in the trap
    trap_delay_ms = 200 #I changed this trap delay to be faster - since I fixed the cycling issue
    trap_delay_ms_long = trap_delay_ms

    if initial_opt == "yes":
        shutters.clear(1,2,delay_ms = trap_delay_ms_long)
    else: #logical loop to check bead quality
        if match_score2.latest_value > 0 and match_score2.latest_value < match_threshold: #this means there must be two beads - one of them is bad because its still in the while loop
            if 0 < match_score2.latest_value < match_threshold: #check if bead on the right (by definition bead 2, on trap 1) is bad
                shutters.clear(1,delay_ms=trap_delay_ms)

            if 0 < match_score1.latest_value < match_threshold: #check if bead on the left is bad (if we removed the bead on the right this should still be labeled as bead 1)
                shutters.clear(2,delay_ms=trap_delay_ms)

        if 0 < match_score1.latest_value < match_threshold and match_score2.latest_value < match_threshold: #this means that there must be one bead, and the bead is bad
            bead_opt = determine_bead_pos(bead_1_position_x.latest_value,centered_of_x_scan) #find where the bad bead is
            
            if bead_opt == "left": #bad bead is on the left (2) so clear trap 2
                shutters.clear(2,delay_ms=trap_delay_ms)
            elif bead_opt == "right": #bad bead is on the right (1) so clear trap 1
                shutters.clear(1,delay_ms=trap_delay_ms)

    print("Trapping beads.")
    start_time = time.time() #initialize clock for long delay trap clear

    max_time_before_quitting = 60*10 #seconds


    while match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        """Drop beads that do not fulfill the template"""
        if time.time() - intial_time_for_flow > max_time_before_quitting: #haven't caught a bead in 10 minutes - kill the program to save time
            print(f'Performing safe exit of the script because have not caught beads in {max_time_before_quitting/60} minutes!')
            stop_flow()
            fluidics.start_venting() #john added short distance before
            exit()

        if time.time() - intial_time_for_flow > flow_timer: #logical gate to stop flow in protein channel
            close_protein_channels() #this function should close channel 4 and 5


        if match_score2.latest_value > 0 and match_score2.latest_value < match_threshold: #this means there must be two beads - one of them is bad because its still in the while loop
            if 0 < match_score2.latest_value < match_threshold: #check if bead on the right (by definition bead 2, on trap 1) is bad
                shutters.clear(1,delay_ms=trap_delay_ms)

            if 0 < match_score1.latest_value < match_threshold: #check if bead on the left is bad (if we removed the bead on the right this should still be labeled as bead 1)
                shutters.clear(2,delay_ms=trap_delay_ms)

        if 0 < match_score1.latest_value < match_threshold and match_score2.latest_value < match_threshold: #this means that there must be one bead, and the bead is bad
            bead_opt = determine_bead_pos(bead_1_position_x.latest_value,centered_of_x_scan) #find where the bad bead is
            
            if bead_opt == "left": #bad bead is on the left (2) so clear trap 2
                shutters.clear(2,delay_ms=trap_delay_ms)
            elif bead_opt == "right": #bad bead is on the right (1) so clear trap 1
                shutters.clear(1,delay_ms=trap_delay_ms)

        """If it's taking too long, maybe something is stuck in the trap. Clear traps that currently don't have beads."""
        if time.time() - start_time > 5:
            if match_score1.latest_value > match_threshold:
                bead_opt = determine_bead_pos(bead_1_position_x.latest_value,centered_of_x_scan)
            
                #now we want to clear the trap that DOESNT have a bead (we have found the trap that does have a bead)
                if bead_opt == "left": #good bead is on the left, clear right trap
                    shutters.clear(1,delay_ms=trap_delay_ms_long)
                elif bead_opt == "right": #bad bead is on the right
                    shutters.clear(2,delay_ms=trap_delay_ms_long)
            else:
                shutters.clear(1,2,delay_ms=trap_delay_ms_long)
            
            #reset start time
            start_time = time.time()

        pause(0.25)

    print("Got beads!")
    time_under_flow = time.time() - intial_time_for_flow
    return time_under_flow

def force_feedback_call(target_pN,max_distance,mode):
    force_feedback_logical = 1
    force_feedback.set_target(target_pN)
    force_feedback.set_device("1")
    force_feedback.set_angle(0)
    force_feedback.set_lock_motion_angle()
    force_feedback.set_detector("Trap 2") # Force detector
    force_feedback.set_frequency(200)

    # access the different force feedback variable between the fishing and extension modes
    #   Ex: force_feedback_variables_extension = [0.4, 0, 0, 100,False]  [kp, ki, kd, max_step,reversed] parameters list for PID 
    if mode == "fishing": 
        force_feedback.set_pid_settings(kp=force_feedback_variables_fishing[0], ki=force_feedback_variables_fishing[1], kd=force_feedback_variables_fishing[2], max_step=force_feedback_variables_fishing[3],reversed=force_feedback_variables_fishing[4])
    else:
        force_feedback.set_pid_settings(kp=force_feedback_variables_extension[0], ki=force_feedback_variables_extension[1], kd=force_feedback_variables_extension[2], max_step=force_feedback_variables_extension[3],reversed=force_feedback_variables_extension[4])

    # reset the force at the start of the extension
    if mode == "extension":
        pause(2)
        reset_force()
        pause(2)

    force_feedback.enabled = True
    
    previous_distance = distance.latest_value
    
    while force_feedback_logical == 1: 
        pause(0.2)
        current_distance = distance.latest_value
        
        if distance.latest_value > max_distance:
            max_distance_opt = 2
            #print('No DNA found, fishing again')
            #print(f"{distance.latest_value-previous_distance} microns")
            bluelake.force_feedback.enabled= False
            return -1
        
        if current_distance - previous_distance < 0: #assumption is that the feedback is fast enough that it wont hit less than zero before 
            force_feedback.enabled = False
            return 1
        
        previous_distance = current_distance

    force_feedback.enabled = False

    return force_feedback_logical


def goto_distance(target, match_threshold, speed=1, tolerance=0.1):
    """Move trap 1 until it reaches the `target` distance from trap 2.

    Note: This throws an error if the beads are lost (since we will not have a reliable distance then either)"""
    dx = target - distance.latest_value
    
    restart_option = throw_if_beads_lost(match_threshold)
    if restart_option:    
        return restart_option

    while abs(dx) > tolerance:  # um
        if dx > 0:
            trap1.move_by(dx=+0.1, speed=speed)
        else:
            trap1.move_by(dx=-0.1, speed=speed)

        dx = target - distance.latest_value
        
        restart_option = throw_if_beads_lost(match_threshold)
        if restart_option:    
            return restart_option
    
    return False



def check_if_should_close_protein_channels(time1,time2,flow_timer):
    if time1 + time2 > flow_timer:
        close_protein_channels()

    return

#function to take trapped beads in the bead channel to the DNA channel and catch a DNA tether
# to change the number of fishing cycles change max_retries
def catch_dna(min_distance, max_distance, match_threshold, fishing_force_threshold, fishing_speed, step_size_each_cycle, min_distance_extension,time_flow_was_on,flow_timer,max_retries=20):
    """Moves to the DNA channel and starts oscillating the trap until a prescribed force threshold is reached."""
    print("Moving to DNA channel")
    microstage.move_to(name_dna_channel)

    initial_time_for_flow_counting = time.time() #time to check with flow
    check_if_should_close_protein_channels(time.time()-initial_time_for_flow_counting,time_flow_was_on,flow_timer) #line to make sure we aren't wasting protein

    #go to min distance to start fishing
    restart_opt = goto_distance(min_distance, match_threshold, speed=6) #if you dont reset the force before this does it error out <- is this true? JWW note 10/19/2022

    #spacer segment to let the force normalize and the check to make sure the beads are still there
    pause(1)
    if restart_opt:
        return restart_opt
    pause(1)
    check_if_should_close_protein_channels(time.time()-initial_time_for_flow_counting,time_flow_was_on,flow_timer) #line to make sure we aren't wasting protein

    #set initial parameters for fishing
    attempts = 0
    restart_opt = False
    fishing_logical = -1 #-1 means no tethers, 1 means tether
    
    reset_force()
    pause(1)
    #fishing while loop that uses force feedback (ONLY COMPATIBLE WITH BLUELAKE V2)
    while fishing_logical != 1:
        check_if_should_close_protein_channels(time.time()-initial_time_for_flow_counting,time_flow_was_on,flow_timer) #line to make sure we aren't wasting protein

        restart_opt = goto_distance(min_distance - (step_size_each_cycle*(attempts)), 0.1, speed=fishing_speed)
        pause(0.25)
        if restart_opt:
            return restart_opt

        restart_opt = goto_distance(min_distance, 0.1, speed=fishing_speed)
        pause(0.5)
        reset_force()
        pause(0.5)

        # method to catch DNA using force feedback
        # is_dna_caught = -1 if the DNA is NOT caught, and 1 if it is potentially caught (or multiple tether)    
        is_dna_caught = force_feedback_call(target_pN=fishing_force_threshold,max_distance=max_distance,mode='fishing')
        if is_dna_caught == -1:
            if attempts > max_retries:
                print(f"Max retries {max_retries} reached.")
                restart_opt = True
                return restart_opt
            attempts += 1

            continue
        elif is_dna_caught == 1:
            fishing_logical = 1
            continue

    
    #stop flow and then prepare for stretching
    fluidics.start_venting() #john added short distance before
    stop_flow()
    restart_opt = goto_distance(min_distance_extension, match_threshold, speed=3) #John added this code so that the tether doesnt rip when it moves
    if restart_opt:
        return restart_opt
    
    #micromicrostage control to get the tether ready for imaging
    print("Moving to buffer channel.")

    pause(1)
    microstage.move_to(name_buffer_only)
    pause(1)
    return False #This program sets up the beads with hopefully tether in the junction to be stretched to force to image with a zeroed force



def prepare_tether_for_imaging(target_pN,max_distance,distance_to_reject_if_multiple,error_threshold,distance_to_relax_if_multiple):
    #extend to target force and measure the length of the tether
    #then do quality control on if that distance makes sense for a single tether
    is_dna_caught = force_feedback_call(target_pN,max_distance,mode='extension') #hard coded in 
    current_distance = distance.latest_value
    restart_loop_opt = False
    
    if is_dna_caught == 1:
        Beep(2500,200)
        pause(0.1)
        Beep(2500,200)
    elif is_dna_caught == -1:
        print(f"No DNA tether detected, restarting protocol")
        restart_loop_opt = True
        return restart_loop_opt, is_dna_caught
    
    error_percentage = (abs(estimated_distance - current_distance) / estimated_distance) * 100
    print(f"Tether stretched to approximately {force_to_pull_to_image} pN at {current_distance} microns")
    print(f"Error % - {error_percentage} for estimated tether length of {estimated_distance}")
    
    # logic to remove tether
    if current_distance > distance_to_reject_if_multiple:
        print('This is a single tether!')
        return restart_loop_opt, 1
    else:
        print("Multiple tethers, restarting loop")
        return True, 2

def make_fd_curve(min_distance, max_distance, match_threshold, fd_speed):
    """Measure an F, d curve."""
    restart_opt = goto_distance(min_distance, match_threshold, speed=3)
    if restart_opt:
        return restart_opt
    
    pause(1)
    reset_force()
    pause(1)
    
    #call the extension method
    goto_distance(max_distance, match_threshold, speed=fd_speed, tolerance=0.2)
    
    pause(1)
    return

def set_pressure(target):
    """Increase pressure until we are above a certain target level"""
    while fluidics.pressure < target:
        fluidics.increase_pressure()
        pause(1.0)

def start_flow(pressure):
    print("Starting Flow.")
    set_pressure(pressure)
    fluidics.open(1, 2, 3, 4, 5, 6)

def stop_flow():
    print("Stopping flow.")
    fluidics.close(1, 2, 3, 4, 5, 6)

def high_pressure_flush(high_pressure,waiting_time_seconds,final_pressure=flow_pressure):
    #ensure the correct channels are open
    fluidics.open(1,2,3,4,5,6)

    # faster pause time to get up to pressure faster
    while fluidics.pressure < high_pressure:
        fluidics.increase_pressure()
        pause(0.01)

    pause(waiting_time_seconds)
    
    while fluidics.pressure > final_pressure:
        fluidics.start_venting()
        pause(0.1)

    fluidics.open(1,2,3,4,5,6) #the channels should re-open in the bead catching program but this is just in case they do not


def close_protein_channels():
    #fluidics.close(4,5)
    #skipping this because it is a C2 chip
    return

def workflow(match_threshold, dna_fishing_speed, starting_distance_fishing, contour_length_max, extension_start_point, max_length_FD_curve,fishing_distance_max,
                fishing_force_target, fishing_step_size, force_to_pull_to_image, max_tethers,flow_timer,distance_to_reject_if_multiple):

    tether_num = 0                      # Tether number counter
    timer_for_force_calibration = 0     #
    restart_loop_opt = False            # This is a boolean logical variable that is referenced frequently to make sure different quality control markers are met (aka still have beads that match the template)
                                        # If this value becomes true at any point in one of the functions it will exit the function almsot immediately and restart the whole loop (i.e will try to catch beads again)
    last_cal_was_bad = False            # variable to catch if the last calibration was bad
    previous_time_flushed = time.time() # not flushing on cycle one -> can incorporate if it is not doing well

    time_in_seconds_to_wait_for_flush = 60*5 #seconds times number of minutes to wait to do the bead flushing protocol (will only check at the start of each bead cycle) 

    # This while loop continues until the maximum number of tethers has been reached
    while tether_num <= max_tethers:
        tether_num += 1
        timer_for_force_calibration += 1
        print(f"Tether no. {tether_num} out of {max_tethers}")
        
        # Check to see if the timer has been triggered to flush the beads 
        if time.time() - previous_time_flushed > time_in_seconds_to_wait_for_flush or tether_num == 1:
            print('Performing high pressure flush to get beads off of the surface')
            print(f'It has been {(time.time() - previous_time_flushed)/60} minutes since the last flush')
            high_pressure_flush(high_pressure=1.0,waiting_time_seconds=10)
            previous_time_flushed = time.time()
            print('High pressure wash completed')

        # Entering the loop to catch beads, make sure the left trap (trap 2!) is to the left of the waypoint defined earlier
        #   Currently, each cycle both beads are 
        time_flow_was_on = catch_beads(match_threshold=match_threshold,flow_timer=flow_timer,initial_opt='yes')
        

        if restart_loop_opt: #exit condition if (i) a bead is lost or (ii) force calibration fails
            restart_loop_opt = False
            continue


        # This program moves the beads to the DNA channel and iteratively fish for a DNA tether
        #   For fishing, I found the best way to do this is to use the force feedback module
        #   Every cycle the force feedback loop is turned on (with a target force of fishing_force_target), the beads will move further apart, and if the beads distance exceeds the max_distance variable it will 
        #       turn off the force feedback and go back to the initial distance (minus a small step size*number of cycles to increase fishing likelihood) and restart the cycle.
        #   To detect if a tether is caught, the change in distance is measured during each force feedback loop:
        #       If no tether, this value will always be positive, because the distance between the beads will increase with each second until it reaches the maximum and the force feedback loop is turned off.
        #       If there is a tether, this value will approach and fluctuate around zero - because the feedback has reached the correct force and is only making minute changes to the mirror position
        #       Once the 'instantaneous' change in distance becomes negative (due to the fluctuations around zero) the program will exit the fishing cycle and relax the tether back to the original fishing point 
        #            before stopping the flow and moving the beads to the juntion channel  
        #   If no DNA is caught in 20 different fishing cycles we reset the loop using the restart_loop_opt boolean
        restart_loop_opt = catch_dna(min_distance=starting_distance_fishing, max_distance=fishing_distance_max,
                                     match_threshold=match_threshold, fishing_force_threshold=fishing_force_target,fishing_speed=dna_fishing_speed,
                                     step_size_each_cycle=fishing_step_size, min_distance_extension=extension_start_point,time_flow_was_on=time_flow_was_on, flow_timer=flow_timer)
        
        if restart_loop_opt: #check for exit condition if the beads don't match or if max number of retries were reached
            restart_loop_opt = False
            continue
        
        #now stretch the DNA to the imaging force in the junction channel, and checking for double tether
        error_range = 20 #%
        restart_loop_opt, is_single_tether = prepare_tether_for_imaging(target_pN=force_to_pull_to_image,max_distance=fishing_distance_max,distance_to_reject_if_multiple=distance_to_reject_if_multiple,error_threshold=error_range,distance_to_relax_if_multiple=starting_distance_fishing)
        if restart_loop_opt: #check for exit condition if we either lost beads or do not have a single tether
            restart_loop_opt = False
            continue
        
        #tether has been pulled to the correct force --> starting experiment
        logical_tether_broke = experiment_imaging(kymograph_method,match_threshold)
        print("Imaging done!")
 
workflow(match_threshold=bead_threshold_score,                  # Minimal template match threshold
            dna_fishing_speed=fishing_move_speed,               # DNA fishing speed [-]
            starting_distance_fishing=fishing_start_point,      # Minimal distance when fishing for DNA [um]
            contour_length_max=estimated_contour_length,        # Maximal distance when fishing for DNA [um] --> contour length
            extension_start_point = extension_start_point_var,  # Distance to relax the tether once caught before moving to buffer channel and preparing for 
            max_length_FD_curve = FD_max_distance,              # Maximum distance to pull to during the force-extension cycle after imaging
            fishing_distance_max = fishing_distance_to_cycle,   # Maximum distance that when hit during a fishing cycle the bead automatically returns to the minimum distance
            fishing_force_target= fishing_force_target,         # Force to target each fishing cycle (2X direction)
            fishing_step_size=step_size_fishing,                # Distance change each fishing cycle (gets negated in the script) so that the beads get closer with each failed fishing cycle
            force_to_pull_to_image = force_to_pull_to_image,    # Force threshold to extend to when 
            max_tethers=max_tethers,                            # Maximal number of attempted tethers [#]
            flow_timer=flow_timer_to_shut_off_ch4and_ch5,
            distance_to_reject_if_multiple=distance_for_single_tether_rejection)       # Time (in s) to wait for flow to turn off in channels 4 and 5       
 

# once you have debugged your system, comment out the above line and make it so that 
"""
try:        
    workflow(match_threshold=bead_threshold_score,              # Minimal template match threshold
            dna_fishing_speed=fishing_move_speed,               # DNA fishing speed [-]
            starting_distance_fishing=fishing_start_point,      # Minimal distance when fishing for DNA [um]
            contour_length_max=estimated_contour_length,        # Maximal distance when fishing for DNA [um] --> contour length
            extension_start_point = extension_start_point_var,  # Distance to relax the tether once caught before moving to buffer channel and preparing for 
            max_length_FD_curve = FD_max_distance,              # Maximum distance to pull to during the force-extension cycle after imaging
            fishing_distance_max = fishing_distance_to_cycle,   # Maximum distance that when hit during a fishing cycle the bead automatically returns to the minimum distance
            fishing_force_target= fishing_force_target,         # Force to target each fishing cycle (2X direction)
            fishing_step_size=step_size_fishing,                # Distance change each fishing cycle (gets negated in the script) so that the beads get closer with each failed fishing cycle
            force_to_pull_to_image = force_to_pull_to_image,    # Force threshold to extend to when 
            max_tethers=max_tethers,                            # Maximal number of attempted tethers [#]
            flow_timer=flow_timer_to_shut_off_ch4and_ch5)       # Time (in s) to wait for flow to turn off in channels 4 and 5
except:
    fluidics.start_venting()
    stop_flow()
    try:
        force_feedback.enabled = False
    except:
        print('Check if you should turn off force feedback')
"""