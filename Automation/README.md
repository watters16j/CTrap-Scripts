Hello!

Attached are my attempts to make life easier by automating operation of the CTrap during data collection. This does not include any passivation or washing!

All of these example scripts were built off of a similar backbone, but are custom built for various tethers and experimental set-ups, so you will have to tune it to fit your instrument but you will save a lot of time in long run!

# General Automation Pipeline
## Set-Up
* Set up your optics as you would normally
* Catch a pair of beads (the assumed geometry is trap 2 on the left and trap 1 on the right for these scripts)
* Perform force calibration and make sure you are in a good z-plane for imaging

## Script Operation
Once you set up the parameters and hit run the script operates in this order
1. Move to bead channel and catch beads
2. Move to DNA channel and fish for DNA
3. Verify single tether formation
4. Start imaging method and perform whatever stage movements or tether movements are needed
5. End image collection
6. Restart

If the process fails at any of these steps, the whole process will restart until you reach the number of max_tethers you want to attempt.

### A few notes on development
These are development notes that I think were critical to making this system work robustly for our lab. I want to state them explicitly in case helps you during your tinkering with this automation code.
1. Bead tracking - huge note that it took me a long time to realize the bead "Tracking Match Score" score is agnostic to trap position (it could either be in trap 1 or 2)
2. DNA fishing - specifically when you want to use low forces during fishing (to keep your a protein loaded complex intact for example) the published ways in the lumicks example did not work well, and in general that method was very slow because it was moving the traps at a linear rate. Incorporation of force feedback makes fishing more robust to transient fluctutations in force that might come from aggregates or additional DNA molecules hitting the bead and can be efficiently tuned to be much faster in terms of fishing cycles per minutes.
3. Image reconstruction - this is only present in the full_auto_bluelakeV3_with_image_decision, but in this I demonstrate how to reconstruct a Scan image from the raw output of the instrument. This is  useful to makeing imaging based decisions on if you should continue on with an experiment or not. In this case it was used to optimze if the loaded protein was present so that we did not waste precious instrument time and protein samples imaging an empty DNA tether. --> I decided to put this as a standalone program that you would have to plug and play manually in the script for simplicities sake.

## Parameters to Change
Most of these features are at the start of the code
### Fishing Parameters (specifics of the tether)
* CRITICAL - change centered_of_scan to be the value (as shown in the x-scan center) that your beads 1 and 2 will NEVER cross. This is crucial to make the bead catching process fast.
* estimated_contour_length, fractional_contour_length_for_fish, fractional_contour_length_for_fishing_extension, fractional_countour_length_for_extension, fractional_contour_length_for_FD -> change re: tether length
* fishing_force_target, force_to_pull_to_image, step_size_fishing, fishing_move_speed, force_feedback_variables_fishing -> to change fishing process
* force_feedback_variables_extension -> to change how slow/fast you want this extension to happen
* x_pos_to_pull_to_for_bead_catching -> change this depending on how separate you want your beads when catching beads
* bead threshold score
* section #Channel names as in the UI -> edit channel names as you use
* scan_method, kymograph_method -> Imaging methods to call during the process
* flow_pressure, time_to_wait_for_flow, flow_timer_to_shut_off_ch4and_ch5 -> options to change how flow is established
* high_pressure_flush -> a useful method that is called ~every 5 minutes in this version of the script. Established high flow rate to re-establish laminar flow

# Additional Notes
Please reach out to me at johnwatters97@gmail.com if you have any technical issues! There will likely be a few pain points getting this to work best for your system but it will save a lot of time in the long run.
