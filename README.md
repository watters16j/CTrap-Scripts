# CTrap-Scripts
Scripts for data analysis/automation using the Lumicks CTrap instrument. There are two different scripts here. "CTrapVis.py" is best used for general looking through kymograph images one by one and extracting desired features. "Area_photon_count_extractor.py" is best used to batch extract photon counts from static molecules.  

## CTrapVis
This script is the first version of a GUI designed to dynamically look through .h5 files from a Lumicks C-Trap instrument as part of a project from Professor Shixin Liu's Laboratory of Nanoscale Biophysics and Biochemistry at The Rockefeller University. This script allows users to take full advantage of the tools in the pylake library to visualize and extract data of interest without python scripting knowledge or the need for a computer with Bluelake software. Script users are able to choose what plots they want (Force vs. Distance, Force vs. Time, RGB image, or a combination of RGB and one of the force options) with options to customize the plot to add a title, scale RGB image values, shift axis of all available axis. This script can analyze kymograph, scans/multiple scans, and force-distance objects contained within a .h5 file

Comparable axes are linked by the manual entries, but the individual plots can be zoomed in/out of using the interactive toolbar. The save name for both the image and the metadata files can be customized by the "File Name to Save" entry and the "Image Format" option. The "Draw Plot" button pulls in all of the GUI information and draws/plots the correct plot back on the interface (shortcut - enter key). The "Show RGB Histogram" button shows the distribution of red, green, and blue pixel intensities after any image manipulation. The "Quit" button destroyed the tkinter GUI and quits the python execution (shortcut - escape key). The "Export Force For Origin" button pulls a separate Tkinter window to allow the user to export desired force data to a .csv file for further analysis in software packages like Origin. The "Export Image For CTrapViewer" exports .png files compatible with Dr. Ioddo Heller's Lab CTrapViewer software. The "Export Image For ImageJ" button exports scan image in such a way that the axis are shown as well as a timestamp of the image for easy transformation into an ImageJ montage.

Python Modules Used: numpy, lumicks.pylake, matplotlib, tkinter, pandas, glob, os

To Run:
1.	Download the CTrapVis.py script and put it in a folder.
2.	Open up terminal and change the working directory to the folder containing the correct script.
3.	Type the command "python CTrapVis.py" and hit enter
4.	After a few seconds, the GUI interface should show up (Note: this may look different on a Mac becuase tkinter widgets match the OS)
  ![image](https://github.com/user-attachments/assets/bf1c73f7-6e09-43b6-a459-c5a2ab29bd60)
5.	Hit the change folder button and navigate to the folder containing the .h5 files that you want to visualize
6.	Choose the file in the "H5 Files in Directory" pulldown menu, as well as the specifc object in the .h5 file to visualize by changing the option in the "File Components" pulldown menu.
  ![image](https://github.com/user-attachments/assets/3b1a261d-2631-4f94-8386-8fd9a4676430)
7.	Other parameters in the Plotting Options, Force Options, or the Photon Count Multipliers can be changed before hitting the "Draw Plot" button to draw the desired plot on the canvas to the left of the option panel (see example Outputs of GUI For Different Data Objects)
8.	Save desired data/image using the export/save buttons on the right of the GUI

***Example Outputs of GUI For Different Data Objects:***

File name and metadata are redacted from these example outputs at the individual's request.

*Kymograph*

![image](https://github.com/user-attachments/assets/e3983404-068c-443b-adb1-d396fc6daaaa)

*Scan*

![image](https://github.com/user-attachments/assets/ec8d8d28-6f34-4a3b-bf17-e503dac7b898)

*Fdcurve*

![image](https://github.com/user-attachments/assets/e99e50a6-ef93-4d7c-bd2a-ad7223e843e8)

**Notes of Potential Interest:**
* "Both" force channel data outputs the magnitude of the force on the bead
* Force 2x data if automatically inverted - because that is standard practice for the Liu lab, this might have been changed based on your version of Bluelake.
  - To remove this, search for the text "if forceString == '2x':" and delete the contents of that if statement
* The default image showing up is RGB only because it is better to only load in the RGB data to test which "Photon Count Multiplier" values give the best image. After this, one can switch to plotting both
* For large kymographs (> 1 GB) - loading high-frequency force data over such a large time window can cause the program to freeze while it completes the calculation (>1 minute)
* Export Image For ImageJ button is uniquely suited for droplet fusion/FRAP experiments where you want to export similar images with both time and position data 
* Doesn't apply any additional functionality for kymograph objects/just scan objects with multiple frames
* For scans with multiple frames - the scan image frame slider will become active and let you toggle through the images. The highlight scan option will add an additional trace covering the range of the force regime that is represented at the same time as the scan image being displayed
* If the GUI window is too large for your screen you can change this by lowering the .set_dpi() parameter from 110 until it doesn't exceed your screen limits (search for "figureReturned.set_dpi(110)")
* The "Fix Image Reconstruction?" option is a vestigial function that would only apply to a user if they are using a version of lumicks.pylake < v0.6.0
  - More info in the changelog: https://lumicks-pylake.readthedocs.io/en/stable/changelog.html
* Images are being reconstructed by directly relating photon count (for each pixel for each color) to an RGB value (after multiplying by the photon count multiplier value and adding the brightness addition value -> always extract the raw .tiff images and adjust contour in something like ImageJ to make sure you are not artifically changing the relative brightness between different foci. 

Versions Tested In:
* Python - 3.7.6
* NumPy - 1.16.5
* lumicks.pylake - 0.8.1 (This script will not work with later lumicks.pylake versions - especially anything v1 after)
* Matplotlib - 3.2.1
* pandas - 1.0.5

## Area Photon Count Extractor
This script is used to extract the sum of photon counts in a line scan from chosen areas of a kymograph in a .h5 file. This analysis method is useful when you want to measure raw intensity values of fluorescent molecules to quantify how many molecules are present without worrying about potential biases introduced by using line tracking methods (most relevant when you are trying to determine the stochiometry in the 1/2/3 molecule range).

The user will be asked via text inputs to navigate to the .h5 file of interest and apply a multiplier value to photon counts if desired (for display purposes only monitoring the method used in CTrapVis.py). The user will then choose how they want to define the area to analyze either by [1] Explicitly defined dimensions or [2] dragging a window on the kymograph in a pop-up window. Then the user will be asked to click on other regions to extract where the click defines the center-left point of the box (in the cartoon below,  '-' is the number of time points being analyzed, '|' is the number of pixels being analyzed, and 'X' marks the spot where the user clicks).

![image](https://github.com/user-attachments/assets/2e321bd3-ff3b-4e4c-905b-54b4cf0082be)

Once the points are selected on the plot, then you can exit the plot. This will show you a separate image showing the regions the script will extract. If you used the  'dragging a kymograph window' method of defining the area dimensions, then you will see the original box drawn in blue with the additionally clicked boxes drawn in orange.

If the user is happy with the point selection, they can confirm that through user input and the script will be extracted the data to a .xlsx file. There are two types of sheets in the output .xlsx. The first sheet-type shows the sum of the columns of each region (at each time point) and the extracted simple statistics from each region (average and standard deviation) for each channel you want to extract. The second sheet-type records some metadata that might be useful if you need to remake any plots/redo any analysis. In addition, an image of the extracted regions is also outputted as a .png for a quick resource on which region #'s align with other regions of the plot.

*Example .png output*

![image](https://github.com/user-attachments/assets/5bb63271-7049-4384-93ac-de9ef7e71516)

*Example .xlsx output*

![image](https://github.com/user-attachments/assets/69020d61-db7f-49e1-a27f-ff37505dc85e)

![image](https://github.com/user-attachments/assets/03b6c37a-6ac0-4420-bdf3-84b236c5c753)

Version Information:
* Lumicks.pylake - 0.8.1
* pandas - 1.4.1
* matplotlib - 3.5.1
* numpy - 1.22.3

## Feedback/Questions/Concerns
Please direct any feedback/issues/constructive criticism/correspondence to jwatters@rockefeller.edu
