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
* You must have Python 2.7, scipy, numpy, and matplotlib installed. You can do any one of the following:
  * Download Anaconda (http://continuum.io/downloads) (HIGHLY RECOMMENDED)
    * Chrome will block the package as "malicious." Go to Downloads and allow the download.
  * Download Canopy Express (https://store.enthought.com/downloads/). (HIGHLY RECOMMENDED)
  * Download any of the other options in the Scipy Stack (http://www.scipy.org/install.html).
  * Download plain old Python 2.7 (https://www.python.org/downloads/), and manually install the libraries (http://www.scipy.org/install.html, scroll down to "Individual Binary and Source Packages")

* Download the PeakSeeker pacakge in a zip. Save to a convenient directory.

* Open PeakSeeker.py in a text editor.
  * If using Anaconda, launch Spyder and open PeakSeeker.py.
  * If using Canopy, launch Canopy and open PeakSeeker.py
  * If using Python, launc your typical code editor.
 
* Specify the directory path of OPTIONS.txt in the code ofPeakSeeker.py:
  * Control+f "opener = open('./OPTIONS.txt')" This is at line 1378.
  * Change './OPTIONS.txt' to your OPTIONS file path, still in apostrophes
  * Save.
  
* Open OPTIONS.txt. You may need to use Word for the file to appear in columns.
    * Change the first line of OPTIONS.txt to the file path of SampleTextFile.txt
      * For other spectra, copy the spectrum from Xcalibur and save in a text file. Change the first line to the new file path.
        * Do NOT remove the header lines; the program will remove these automatically.
    * The first line should look like "directory	(file path of text file)"
    * Save.

Run through Manual Mode once to test your installation.
-------------------------------------------------------
1. Run PeakSeeker.py
  * If using Anaconda, launch the iPython-qtconsole. Change the directory to PeakSeeker.py's folder. Then, run PeakSeeker. 

2. A plot (as below) will pop up. This lists 5 peaks with dots and centroids. These are the 5 tallest peaks. You should be able to interact with this plot (try clicking the magnifying glass to zoom in on certain parts). Depending on the backend of your matplotlib, the plot will A) be displayed for 5 seconds before it loses interactivity, or B) stay interactive for as long as you want, and then once you close the plot, the command line will wait 5 seconds before displaying.
   * If the plot is not interactive at all, matplotlib may be using the wrong backend. Kill the program, then go to lines 23 and 25 in PeakSeeker.py to delete the '#' sign. This should switch to the correct backend, if you are using Anaconda or Canopy.  
   * ^ If you have A), DO NOT CLOSE THE PLOT. It will close automatically once you enter in the command line.   
   * ^ If you have B), you need to close the plot in order to return to the command line. You will have to wait 5 seconds. You can adjust the wait time by changing __PeakDispTime__ and __ChargeDispTime__ in OPTIONS.txt

  ![Peak Display 1](http://i.imgur.com/UM7x7WN.png)


3. Return to the command line^. If you would like to continue interacting with the plot, enter /. This will close the plot and reopen it for a brief interactivity time. Otherwise, enter 1. This is the number of the tallest peak. The program will now iterate charge states to this peak.

  ![Command Line Peak](http://i.imgur.com/32hCfTD.png)


4. Another plot will pop up. This lists several peak series with different colored dots. The dot height represents the height of the simulated peak in the simulated charge envelope. Peaks with the same colored dots are in the same envelope. The number next to the colored dots in the key is the corresponding charge state of the central peak. The top number is the best fit, according to the scoring function. Note that 48 fits better than 50, for example, because the dots are closer to the heights of the real peaks.
   ![Charge Display 1](http://i.imgur.com/5cSYEsb.png)

7. Return to the command line^. It should list the charges in the plot, the corresponding masses, and the list of scores. A lower score means a better fit. The charges are already ordered from lowest to highest score.
8. Enter 48. This will save the envelope and mark off the peak as simulated.

   ![Command Line Charge 1](http://i.imgur.com/wLhtBso.png)


9. Another plot will pop up. This lists the 5 tallest peaks that haven't been simulated, similar to step 4.
   ![Peak Display 2](http://i.imgur.com/nyclcn6.png)
10. Return to the command line^. Enter 1. The program will now iterate charge states to this peak.

11. Another plot will pop up, again listing the possible charge states.
   ![Charge Display 2](http://i.imgur.com/61C91nE.png)

12. Return to the command line^. Enter 50. Note that 25 has a better score. However, 25 gives a mass that is half of the true mass (the red dots only have every other peak). The scoring algorithm will sometimes prefer charge states with less peaks, because these are easier to fit well.

   ![Command Line Charge 2](http://i.imgur.com/Ym4IOHz.png)


13. Another plot will pop up with the remaining peaks. There are no more major charge envelopes, so enter n.
14. The final plots will display all of the simulated envelopes, the sum of the simulated envelopes, and the subtracted spectra. You can save these plots.

![Final Display 1](http://i.imgur.com/0H6Gql6.png)
![Final Display 2](http://i.imgur.com/qy9h1wp.png)
![Final Output](http://i.imgur.com/u02UBVH.png)


User Manual
===========

### STARTING PARAMETERS
* __Automatic__
	* The program tries to find charge envelopes, starting from the most intense peak.
	* It can make mistakes, based on the scoring function. See the section "Charge State Assignment" below. 
	
* __Manual__
	* The user can choose which peak to deconvolute, and which charge state is the best answer from a list of charge states.
	* The program works by displaying the spectrum with a list of options. See Step 4 in the manual installation above.
		* __PeakDispTime__ Time in seconds the peak plot will display.
		* __ChargeDispTime__ Time in seconds the charge plot will display.
		* Return to the command line to make your choice^.

* Limit the m/z window, charge state range, mass range, and peak width to those which are desired. This can improve runtimes.


### PROCESSING

- Make sure only the desired processing (background remove, smoothing, or Savitzky-Golay) is turned on.

- The Savitzky-Golay filter preserves peak shape better than the Smooth. 

### PEAK DETECTION

- __PeakFindWay__
	- i simply finds local maxima. This should be used with a filter in case of noisy spectra. This works for most spectra. 
	- c performs a continuous wavelet transform. This is useful for very noisy spectra. See scipy.signal.find_peaks_cwt
	- m uses Massign's fixed and adjusted threshold. See the SI of Morgner, N.; Robinson, C. Anal. Chem. 2012, 84, 2939-2948. 

- __JustPeaks__ Use this to see which peaks are detected before running the whole program workflow.

- __Threshold__ Use to filter small noisy peaks out.

- __FindOverlaps__ If you notice odd peak shapes, there may be overlapping signals.
  - __WeightWidths__ If your overlap fitting is imperfect, the starting guessses for the peak parameter may be off. Turn this on to weight the starting widths by the signal's heights, instead of giving all peaks equal widths.

### CHARGE STATE ASSIGNMENT

- __MassTolerance__ Set this to equal half of the FWHM of a typical peak in the spectrum.
	- Make big enough so that you're finding the peaks close to the simulated peak.
	- Make small enough so that you don't find too many wrong peaks close to the simulated peak.
	- EXPLANATION: when we iterate charge states, we calculate the m/zs at other charge states. Then, we search for peaks that correspond to those other m/zs. If the m/z is within the Mass Tolerance away from a peak's centroid, then that peak is considered. If there are multiple peaks found, the one that has the least deviation from the simulated peak is deemed correct. This deviation is determined by the scoring function.

- __ScoreLimit, Mass ErrWeight, HeightErrWeight, WidthErrWeight__
  - The Scoring Function determines:
	  1. which real peak is closest to the simulated peak in the simulated charge envelope.
	  2. how well the simulated peak approximates the real peak, and therefore how well the simulated charge envelope fits the spectrum.
	  3. the correct charge state in Automatic Mode.

  - THE SCORING FUNCTION MAKES MISTAKES.
	  - It can give preference to charge states that find less numbers of peaks, because it is easier for them to fit into a simulated charge 	envelope.
		 - This can happen with charge states that are 1 off the true charge state, or charge states that are half the true charge state.
		  - Use min_peak_num to filter those out.

- Prefer mass error over height error for spectra with many different peaks.
	- ScoreLimit: 50000
	- MassErrWeight: 50000
	- HeightErrWeight: 100
	- WidthErrWeight: Doesn't really matter

- Prefer height error over mass error for spectra with overlapping peaks.
	- ScoreLimit: 50000
	- MassErrWeight: 1000
	- HeightErrWeight: 500
	- WidthErrWeight: Doesn't really matter

- __AttemptLim__ The number of charge state options the user can choose from.

- __MaxSimulations__ The program will try to find as many envelopes as in MaxSimulations BEFORE SUBTRACTING.
	- The peaks in a given charge envelope can be selected from previous charge envelopes, since there is the possibility of overlapping peaks. 
	- Possible problems:
		- Later envelopes may erroneously find peaks that have already been accounted for.
		- If MaxSimulations is too large and many envelopes share the same peaks, the program may come up with negative or very big 				abundances, since it will try to fit all of the envelopes at the same time.
	- Use Manual Mode to choose the correct peaks.
	- Use RepeatSearch so you don't have to fit all of them at the same time. This subtracts the found envelopes and repeat searching.

### DISPLAYS AND SAVING SUBTRACTED SPECTRUM, MASS INFORMATION
- __DISPLAYS:__
  - Displays can be zoomed in (magnifying glass image) and out.
  - The window can be shifted (the cross with arrows image). 
  - The pictures can be saved in one of several formats (save file image).
  - You can go back to the original display. (home image)
  - See matplotlib for ways to manipulate displays.

- __SaveSubtract__ saves the subtracted spectrum as a text file in the same directory as the original spectrum, with the same header.

- __SaveMasses__ saves the mass information as a text file in the same directory as the original spectrum.

To-Do List
----------
* Add a capability for simulating envelopes based on user's input parameters.
* Document and clean up code.
