PeakSeeker
==========

What is it?
-----------
PeakSeeker is a project created by Jonathan Lu, under the invaluable supervision of Dr. Michael J. Trnka and Dr. Shenheng Guan in the Burlingame Lab at UCSF. PeakSeeker is a comprehensive and flexible program to process, detect peaks, and assign charge states in complicated native electrospray ionization protein complex mass spectra. It is command line based, and runs in an automatic or manual mode. Raw spectra are read in as tab-delimited text files, as copy-pasted from XCalibur. Deconvoluted envelopes and their mass information can be saved as picture and text files, respectively. Manual mode uses interactive figure plotting via matplotlib and iPython.

How does it work?
-----------------
The details of the algorithm will be available in an upcoming paper.

Installation
------------
* You must have Python 2.7, scipy, numpy, and matplotlib installed. You have one of the following options:
  * the Anaconda Distribution's iPython (http://continuum.io/downloads) 
    * If using Anaconda, run the program through the iPython qtconsole, not Spyder. 
  * the Enthought Canopy Environment (https://store.enthought.com/downloads/). 
  * manual installation of libraries

* Download and save both PeakSeeker.py and OPTIONS.txt

* Specify the directory path of OPTIONS.txt in the program code of PeakSeeker.py:
  * Control+f "opener = open('./OPTIONS.txt')" This is at line 1378.
  * Change './OPTIONS.txt' to your OPTIONS file path, still in apostrophes
  
* Similarly, Change the first line of OPTIONS.txt to the directory path of the mass spectrum .txt file, as copy-pasted from XCalibur.
    * Do NOT remove the header lines of the text file with the mass spectrum. The program will remove these automatically.
    * The first line should look like "directory	(directory path of text file)"

Practice Manual Mode Runthrough
-------------------
1. Save the SampleTextFile.txt
2. Change the first line of OPTIONS.txt to the path of SampleFile.txt. Save.
3. Run PeakSeeker.py
4. A plot will pop up. This lists 5 peaks with dots and centroids. These are the 5 tallest peaks. The plot will lose interactivity after a few seconds.
5. DO NOT CLOSE THE PLOT.  Return to the command line. If you would like to continue interacting with the plot, enter /. This will close the plot and reopen it for a brief interactivity time. Otherwise, enter 1. This is the number of the tallest peak. The program will now iterate charge states to this peak.
6. Another plot will pop up. This lists several peak series with different colored dots. The dot height represents the height of the simulated peak in the simulated charge envelope. Peaks with the same colored dots are in the same envelope. The number next to the colored dots in the key is the corresponding charge state of the central peak.
7. DO NOT CLOSE THE PLOT. Return to the command line. It should list the charges in the plot, the corresponding masses, and the list of scores. A lower score means a better fit. The charges are already ordered from lowest to highest score.
8. Enter 48. This will save the envelope and mark off the peak as simulated. Note that 48 fits better than 50, for example, because the dots are closer to the heights of the real peaks.
9. Another plot will pop up. This lists the 5 tallest peaks that haven't been simulated, similar to step 4.
10. DO NOT CLOSE THE PLOT. Return to the command line. Enter 1. The program will now iterate charge states to this peak.
11. Another plot will pop up, again listing the possible charge states.
12. DO NOT CLOSE THE PLOT. Return to the command line. Enter 50. Note that this has a higher score than 25. However, 25 gives a mass that is half of the true mass (the envelope has every other peak). The scoring algorithm will sometimes prefer charge states with less peaks, because these are easier to fit well.
13. Another plot will pop up with the remaining peaks. There are no more major charge envelopes, so enter n.
14. The final plots will display all of the simulated envelopes, the sum of the simulated envelopes, and the subtracted spectra. You can save these plots.


The program isn't working! What do I do?
----------------------------------------
See Strategies for Using PeakSeeker.txt

To-Do List
----------
* Add a capability for simulating envelopes based on user's input parameters.
* Document and clean up code.
