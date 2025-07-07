Hello!

Attached are my attempts to make life easier by automating operation of the CTrap druing data collection. This does not include any passivation or washing!

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
3. Image reconstruction - this is only present in the full_auto_bluelakeV3_with_image_decision, but in this I demonstrate how to reconstruct a Scan image from the raw output of the instrument. This is  useful to makeing imaging based decisions on if you should continue on with an experiment or not. In this case it was used to optimze if the loaded protein was present so that we did not waste precious instrument time and protein samples imaging an empty DNA tether.

## Parameters to Change
Most of these features are at the start of the code
### Fishing Parameters (specifics of the tether)
* estimated_contour_length, fractional_contour_length_for_fish, fractional_contour_length_for_fishing_extension, fractional_countour_length_for_extension, fractional_contour_length_for_FD
* fishing_force_target, fishing_start_point, fishing 
