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
* You must have Python 2.7, scipy, numpy, and matplotlib installed. I highly recommend using iPython, either in the Anaconda Distribution (http://continuum.io/downloads) or Enthought Canopy Environment (https://store.enthought.com/downloads/). If using Anaconda, run the program through the iPython qtconsole, not Spyder. 

* Download and save both PeakSeeker.py and PEAKSEEKER OPTIONS.txt

* In order to use the program, you must specify the directory path of PEAKSEEKER OPTIONS.txt in the program code of PeakSeeker.py. Control+f "opener = open('C:\Python27\PEAKSEEKER OPTIONS.txt')" This is around line 1353.

* Similarly, PEAKSEEKER OPTIONS.txt must specify, in the first line, the directory path of a text file with the mass spectra (i.e. copy-pasted from XCalibur) 

* Do NOT remove the header lines of the text file. The program will remove these automatically.

The program isn't working! What do I do?
----------------------------------------
See Strategies for Using PeakSeeker.txt

To-Do List
----------
* Add a capability for simulating envelopes based on user's input parameters.
* Document and clean up code.
