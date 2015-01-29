""" Copyright 2013, 2014 Jonathan Lu, Michael J. Trnka, Shenheng Guan.
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>."""

import time
import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import stats
from scipy import signal
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
#plt.switch_backend('Qt4Agg')
from math import factorial
import operator
from itertools import groupby
from matplotlib.ticker import MultipleLocator


def Smooth(array, smooth_window):
    """This function smooths a list of values by replacing each value with the
    average of all of the values in a window of size 2*smooth_window +1 centered
    at that point."""
    a = len(array)
    temp_intensity = array[:] #don't use the changed intensity values
    for i in range(a):
        if smooth_window<=i and i <= a-smooth_window-1:
            for j in range(i-smooth_window,i):
                array[i]+= temp_intensity[j]
            for j in range(i+1, i+smooth_window+1):
                array[i]+= temp_intensity[j]
            array[i] /= smooth_window*2 +1
        elif i < smooth_window:
            for j in range(i):
                array[i] += temp_intensity[j]
            for j in range(i+1, 2*i+1):
                array[i] += temp_intensity[j]
            array[i]/=(2*i+1)
        elif i >= a-smooth_window - 1:
            for j in range(i+1, a):
                array[i] += temp_intensity[j]
            for j in range(2*i-a +1 , i):
                array[i]+= temp_intensity[j]
            array[i] /=  (2*a -2*i-1)
        else: print "ERROR in Smooth: out of range"
    return array

def BackgroundRemove(array, step_window, smooth_window):
    """This function subtracts the background from the spectrum. Background is
    defined as the minimum within a window of size 2*step_window +1 centered at
    each value in the list. This background is then smoothed, and then subtracted from each
    point."""
    global mzarray
    background = np.zeros(len(array))
    for i in range(len(array)):
        if step_window<=i and i<=len(array)-step_window-1:
            background[i] += min(array[i-(step_window):i+(step_window)+1])
        elif i<step_window:
            background[i] += min(array[0:2*i+1])
        elif i >len(array)-step_window-1:
            background[i] += min(array[2*i-len(array)+1: len(array)])
        else: print "ERROR in BackgroundRemove: out of range"
    Smooth(background, smooth_window)
    array = np.array(array)- background
    #plt.plot(mzarray, background)
    return(array)

class Peaks:
    """The Peak object stores the m/z and intensity data points of the peak, as
        well as calculated values such as the peak's height, centroid, and
        whether or not it overlaps other peaks. Deconvolution is performed on
        the peak objects rather than the raw spectrum."""
    def _init_(self):
        self.simulated = False     #marks whether peak is already in a previously identified charge envelope
        self.mz = []             #list of m/z values of peak's data points
        self.intensity = []     #list of intensity values of peak's data points
        self.indices= []         #list of indices of peak's data points 
        self.height = 0         #intensity of highest point in peak
        self.centroid = 0         #weighted average of m/z values of peak
        self.fit_center = 0         #center of Gaussian fitted to peak
        self.fit_height = 0     #height of Gaussian fitted to peak
        self.std = 5
        self.fit_std = 5         #standard deviation of Gaussian fitted to peak
        self.left = None         #Peak index of peak adjacent to the left
        self.right= None         #Peak index of peak adjacent to the right
        self.is_empty = True     #whether Peak object has any values within
        self.width = 0
    
    def c(self):                #sets up the Peak centroid value
        cent = 0
        for i in range(len(self.mz)): 
           cent += self.mz[i] * self.intensity[i]
        cent = cent / sum(self.intensity)
        self.centroid = cent
   
    def find_shoulders(self):     #returns a list of each peak's range of indices in spectrum terms. also removes smaller peaks from peak mz
        global threshold
        global adj_factor
        temp_indices = []
        local_max = []         #local_max and local_min hold the indices in spectrum terms
        local_min = []         # of the local maxima and minima, respectively
        adjthresh = threshold
        shift = self.indices[0]
        for i in range(len(self.mz)): 
            if i==0 or i == len(self.mz)-1:
                is_max, is_min = False, True
            else:
                is_max = self.intensity[i] > self.intensity[i-1] and self.intensity[i] > self.intensity[i+1]
                is_min = self.intensity[i] < self.intensity[i-1] and self.intensity[i] < self.intensity[i+1]
            if is_max:
                local_max.append(i+shift)
            if is_min:
                local_min.append(i+shift)

        local_max.sort()
        local_min.sort()
        j=0 #j iterates through the local_max and local_min lists
        if len(local_min) <= 2: return([self.indices])
        while 1:
            if j >= len(local_max) or j >= len(local_min): break;
            if self.intensity[local_max[j]-shift] > adjthresh:
                adjthresh = adj_factor * self.intensity[local_max[j]-shift]
                temp_indices.append(range(local_min[j], local_min[j+1]+1))
            j+=1

        return(temp_indices)
    
    def h(self):
        if self.intensity == []:
            print "ERROR: Peak is Empty", self.simulated, self.is_empty, self.indices, self.mz, self.intensity
        self.height = max(self.intensity)
    
    def w(self):
        for i in range(len(self.mz)):
                if self.intensity[i] > 0.5 * self.height: break;
        for j in range(i, len(self.mz)):
                if self.intensity[j] < 0.5 * self.height: break;
        self.width = self.mz[j] - self.mz[i]
    
    def PeakFit(self, UseDiff, SimNumber): #this is fitted over m/z, not charge. it doesn't use self.intensity so it can accomodate subtracted spectra
        if UseDiff:
            intensity = DifferenceSpectrum(self.mz, SimNumber)
        else:
            intensity = self.intensity
        #print intensity
        if len(intensity) >=3:
            try:
                self.fit_std = 1 #in case the peak was used before
                popt, pcov= curve_fit(norm, self.mz, intensity, [self.height, self.centroid, self.fit_std])
                self.fit_height = popt[0]
                self.fit_center = popt[1]
                self.fit_std = popt[2]

            except RuntimeError:
                self.fit_height = self.height
                self.fit_center= self.centroid
                self.fit_std = self.fit_std

        else: #if the peak is only 1-2 points then you can't fit a gaussian to it.
            self.fit_height = self.height
            self.fit_center= self.centroid
            self.fit_std = self.fit_std
    
    def PeakPlot(self):
        inputs = np.linspace(min(self.mz), max(self.mz), 40)
        outputs = norm(inputs, self.fit_height, self.fit_center, self.fit_std)
        #print "Inputs are", inputs, "Outputs are", outputs
        return(inputs, outputs)
        
class Envelope:
    """The Envelope object is used to find the correct charge state envelope for
        the chosen central peak. It stores the peak heights, charges, and
        locations, as well as the parameters of the Gaussian function that is
        fitted to the peak tops."""
    def _init_(self):
        global StdDevStartGuess
        self.number = 0 #the number of the envelope
        self.central_peak_index = 0
        self.center = 0
        self.height = 0
        self.central_fit_height = 0
        self.central_fit_std = 0
        self.central_fit_mz = 0
        self.mass = 0
        self.charge = 0
        self.peaks = [] #peak number in list Peaks
        self.peak_charges = []
        self.peak_mz = []
        self.peak_heights = []
        self.peak_std = []
        self.curve_scale= 0
        self.curve_center = 50 #on the z axis, not on the m/z axis
        self.curve_std = StdDevStartGuess #random beginning guess (also on z axis)
        self.abundance = 0
    
    def ChargeState(self, Peak, wrong_charges, UseDiff, TryCharge = []): 
        """Here, the central peak charge, the mass of the envelope, peak
        centroids, heights, charges, and indices are found. It also marks off the
        explained peaks as simulated so they won't be chosen for the central peak
        of the next envelope."""
        
        global max_charge  
        global min_charge
        global max_mass #max mass used to narrow down the charge states to be tested
        global min_mass
        global PeakNumber
        global MassTolerance
        global ScoreLimit
        global min_peak_number
        global StdDevStartGuess
        
        #the following block of code is to reset values in case we're doing an interative fit
        self.peaks = []
        self.peak_charges = []
        self.peak_heights = []
        self.peak_mz = [] 
        self.peak_std = []
        max_match_right = 0
        max_match_left = 0
        self.charge = 0
        self.curve_std = StdDevStartGuess
        
        
        for z in range(min_charge, max_charge+1):
            if TryCharge == [] or z in TryCharge:
                if not z in wrong_charges:
                    mass = z*(self.central_fit_mz - 1.00794)
                    if min_mass< mass <max_mass:
                        right_match = 0 #how many peaks to the right simulation gets right
                        left_match = 0 #how many peaks to the left simulation gets right
                        charge_state=z+1
                        
                        while 1:
                            if charge_state > max_charge: break;
                            mz = mass/charge_state + 1.00794
                            
                            for i in range(PeakNumber):
                                in_peak = abs(mz - Peak[i].centroid) < MassTolerance
                                if in_peak:
                                    left_match += 1
                                    break;
                            if not in_peak:#stops simulation when it doesn't explain a peak
                                break;
                            charge_state+=1
                        charge_state = z-1
                        while 1:
                            if charge_state < min_charge: break;
                            mz = mass/charge_state + 1.00794
                            for i in range(PeakNumber):
                                in_peak = abs(mz - Peak[i].centroid) < MassTolerance
                                if in_peak:
                                    right_match += 1
                                    break;
                            if not in_peak:
                                break;
                            charge_state-=1
                        if right_match + left_match > max_match_right + max_match_left:
                            max_match_right= right_match
                            max_match_left = left_match
                            self.charge = z

    
        self.mass = self.charge * (self.central_fit_mz - 1.00794) #here mass is calculated from supposed charge state
        self.peak_charges = range(self.charge - max_match_right, self.charge + max_match_left+1)
        
        if len(self.peak_charges)<min_peak_number:
            return(False, ScoreLimit)
        
        peak_intensities = []
        for z in self.peak_charges:
            indiv_height = []
            mz = self.mass/z + 1.00794
            for i in range(PeakNumber):
                if abs(mz - Peak[i].centroid) < MassTolerance:
                    Peak[i].PeakFit(UseDiff, self.number)
                    #print "Hello", Peak[i].centroid, Peak[i].fit_height
                    indiv_height.append(Peak[i].fit_height)
                    #print "Hi", indiv_height
            #print "yup, charge is", z, "height is", indiv_height, "and mz is", mz
            peak_intensities.append(sum(indiv_height)/len(indiv_height))
        popt, pcov= curve_fit(norm, np.array(self.peak_charges), np.array(peak_intensities), [self.height, self.charge, self.curve_std])
        
        
        resid = []         #list of the deviations between sim peak and real peak
        correct_peaks = []
        for z in self.peak_charges: #THIS PIECE OF CODE CALCULATES ALL DEVIATIONS BETWEEN SIMULATED PEAK AND ACTUAL PEAK
            peakind = []
            peakstd = []
            sim_height = norm(z, popt[0], popt[1], popt[2])
            sim_center = self.mass/z +1.00974
            sim_std = self.central_fit_std * self.charge / z
            for i in range(PeakNumber):
                if abs(sim_center - Peak[i].centroid) < MassTolerance:
                    Peak[i].PeakFit(UseDiff, self.number)
                    peakind.append(i)
            for j in range(len(peakind)):
                peakstd.append(ScoreFunction(sim_center, Peak[peakind[j]].fit_center, sim_height, Peak[peakind[j]].fit_height, sim_std, Peak[peakind[j]].fit_std, False))

            resid.append(min(peakstd))
            correct_peaks.append(peakind[peakstd.index(min(peakstd))])

        for l in correct_peaks:
            self.peaks.append(l)
            self.peak_heights.append(Peak[l].fit_height)
            self.peak_mz.append(Peak[l].fit_center)
            self.peak_std.append(Peak[l].fit_std)

            if TryCharge != []: 
                Peak[l].simulated = True
                
        self.mass = 0
        for i in range(len(self.peaks)): #this allows each m/z peak to inform the estimate of mass's value
            self.mass += self.peak_charges[i]* (self.peak_mz[i]-1.00794)
        self.mass /= len(self.peaks)

        return (True, sum(resid)/len(resid))
        


                        
    def CurveFit(self, Peak):
        global ScoreLimit
        global StdDevStartGuess
        self.curve_scale = self.height
        self.curve_std = StdDevStartGuess
        if len(self.peak_heights) >=3:

            popt, pcov= curve_fit(norm, np.array(self.peak_charges), np.array(self.peak_heights), [self.curve_scale, self.charge, self.curve_std])
            self.curve_center = self.charge
            self.curve_scale = popt[0]
            self.curve_center = popt[1]
            self.curve_std = popt[2]

        else:
            self.curve_scale = self.height
            self.curve_center =  sum(self.peak_charges)/len(self.peak_charges) #Averages peak charges if you can't fit a gaussian
            self.curve_std = StdDevStartGuess
        total_std = []
        for j in range(len(self.peaks)):
            sim_height = norm(self.peak_charges[j], self.curve_scale, self.curve_center, self.curve_std)
            sim_center = self.mass/self.peak_charges[j] +1.00974
            sim_std = self.central_fit_std * self.charge / self.peak_charges[j]
            total_std.append(ScoreFunction(sim_center, self.peak_mz[j], sim_height, self.peak_heights[j], sim_std, self.peak_std[j], False))
        try:
            returner = sum(total_std)/len(total_std)
        except ZeroDivisionError:
            returner = ScoreLimit
        return(returner)
    
    
    def Function(self, x): #outputs the y-value of the peak spectrum given the mz-value
        output = 0

        for z in self.peak_charges:
            peaksim_height = norm(z, self.curve_scale, self.curve_center, self.curve_std) #gets peak height on the gaussian fitted to envelope distribution
            output += norm(x, peaksim_height, (self.mass/z) + 1.00794, self.central_fit_std * self.charge / z) #height on this peak

        return(output)
    
    def SimulationSetUp(self, Peak, number):
        self.number = number
        index = Max_Unsimulated_Index(Peak)
        if index == 0.5 or index == 0:
            return False
        else:
            self.central_peak_index = index #the central index is unsimulated but other peaks in the envelope can have been simulated
            self.center = Peak[self.central_peak_index].centroid
            self.height = Peak[self.central_peak_index].height
            return True
    
    def CentralPeakFit(self, Peak, UseDiff):
        Peak[self.central_peak_index].PeakFit(UseDiff, self.number)
        self.central_fit_height = Peak[self.central_peak_index].fit_height
        self.central_fit_std = Peak[self.central_peak_index].fit_std
        self.central_fit_mz = Peak[self.central_peak_index].fit_center
        print "Considering the peak at ", round(self.center), ": \n"
    
    def CompleteFit(self, Peak, UseDiff):
        global BigDisplay
        global mzarray
        global intensityarray
        global ScoreLimit
        global Automatic
        global ChargeDisplayTime
        global LetterDict
        global AttemptLimit
        global RuntimeErrorLimit
        
        wrong_charges = []
        scores = []
        tried_charges = []
        charge_works = True #confirms that the supposed charge state explains enough peaks to be fitted to a gaussian
        possible_peak_mz  = []
        possible_peak_heights = []

        RuntimeErrorCounter = 0
        while True:
            try:
                charge_works, stddev = self.ChargeState(Peak, wrong_charges, UseDiff) #just use this stddev or refit?
                print "\npossible charge is", self.charge
                print "corresponding mass is", round(self.mass)
                print "corresponding peak mzs are", [round(m, 1) for m in self.peak_mz]
                print "corresponding peak charges are", self.peak_charges
                print "score is", round(stddev, 4)
		if self.peak_mz == [] or self.charge == 0:
                    print "Error: no more possible charges found."
                    break;
                else:
                    if Automatic:
                        for z in tried_charges:
                            if z % self.charge ==0:
                                print "charge", self.charge, "is a divisor of earlier charge", z
                                charge_works = False
                                break;
                if charge_works:
                    stddev = self.CurveFit(Peak)

                else:
                    stddev = ScoreLimit
                    print "this charge state fails."

                scores.append(stddev)

                wrong_charges.append(self.charge)
                tried_charges.append(self.charge)
                sim_peak_mz = [self.mass/z +1.00794 for z in self.peak_charges]
                possible_peak_mz.append(sim_peak_mz)
                possible_peak_heights.append([self.Function(i) for i in sim_peak_mz])

                
                if len(tried_charges) >= AttemptLimit: #THIS CAN BE CHANGED IF NECESSARY
                    break;

            except RuntimeError:
                print "RuntimeError for fitting peak", self.center, "to charge", self.charge
                wrong_charges.append(self.charge)
                RuntimeErrorCounter +=1
                if RuntimeErrorCounter > RuntimeErrorLimit:
                    print "No more possible charges found."
                    break;
        print "\n"
            #error means doesn't converge when fitting to gaussian. go back and switch the charge state
        if Automatic:
            print "for peak", round(self.center), "possible charges are", tried_charges
            print "corresponding masses are", [round((self.center -1.00794)* z) for z in tried_charges]
            print "corresponding scores are", [round(m) for m in scores]
        if len(tried_charges) == 0:
            return False
        if sum(np.array(scores)>=ScoreLimit) >=len(scores): #checks if no correct charge states were found (all would be greater than/equal to ScoreLimit then)
            return False
        else:
            if Automatic:
                self.charge = tried_charges[scores.index(min(scores))] 
            else:
                col = [tried_charges, scores, possible_peak_mz, possible_peak_heights]
                row = zip(*col)
                sort_row = sorted(row, key = operator.itemgetter(1))
                sort_col = list(zip(*sort_row))
                print "for peak", round(self.center), "possible charges are", sort_col[0]
                print "corresponding masses are", [round((self.center -1.00794)* z) for z in sort_col[0]]
                print "list of scores is", [round(m, 4) for m in sort_col[1]]
                while 1:
                    plt.figure()
                    ax = plt.subplot(1,1,1)
                    ax.axis([0.9*min([min(i) for i in possible_peak_mz]), 1.1* max([max(i) for i in possible_peak_mz]), 0, 1.2* max([max(j) for j in possible_peak_heights])])
                    ax.plot(mzarray, intensityarray, 'k')
                    ax.text(self.center, self.height, "Central Peak")
                    d = range(len(tried_charges))
                    d.reverse()
                    for i in d:
                        ax.plot(sort_col[2][i], sort_col[3][i], color = ColorDict[i], marker = 'o', linestyle = 'None', label = "%d" % (sort_col[0][i]))
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[::-1], labels[::-1], loc = 'best')
                    if BigDisplay:
                        figManager = plt.get_current_fig_manager()
                        figManager.window.showMaximized()
                        figManager.window.raise_()
                    plt.show()
                    plt.pause(ChargeDisplayTime)
                    while 1:
                        print "Enter / if you would like to continue viewing. Enter n to cancel this simulation. Or, enter one of the possible charges:", sort_col[0]
                        entry = raw_input()
                        if entry == '/' or entry == 'n' or entry in [str(z) for z in sort_col[0]]:
                            break;
                        print "Sorry, that is not a valid entry. Please reenter."

                    plt.close()
                    if entry != '/':
                        break;
                if entry == 'n':
                    return False
                self.charge = int(entry)
            
            self.ChargeState(Peak, [], UseDiff, [self.charge])
            self.CurveFit(Peak)

            self.abundance = sum([math.sqrt(2*math.pi) * norm(self.peak_charges[i], self.curve_scale, self.curve_center, self.curve_std) * self.peak_std[i]  for i in range(len(self.peak_std))])
            print "mass is", self.mass
            print "m/z center is", self.center
            print "peak charges are", self.peak_charges
            print "peaks are", self.peak_mz
            print "abundance is", self.abundance
            print "\n"
            #print Peak[self.central_peak_index].simulated
            return True            

                
def ScoreFunction(m0, m1, h0, h1, sd0, sd1, printplease): #returns the value of the scoring function that weights each of the errors
    global MassErrWeight
    global HeightErrWeight
    global StdDevErrWeight

    if printplease:
        print "masses", m0, m1
        print "mass score is", MassErrWeight*abs((m0-m1)/m0)
        print "heights", h0, h1
        print "height score is", HeightErrWeight*abs((h0-h1)/h0)
        print "stddev", sd0, sd1
        print "stddev score is", StdDevErrWeight*abs((abs(sd0)-abs(sd1))/sd0)
        print "score is", MassErrWeight*abs((m0-m1)/m0) + HeightErrWeight*abs((h0-h1)/h0) + StdDevErrWeight*abs((abs(sd0)-abs(sd1))/sd0)
    return (MassErrWeight*abs((m0-m1)/m0) + HeightErrWeight*abs((h0-h1)/h0) + StdDevErrWeight*abs((abs(sd0)-abs(sd1))/sd0))

def Max_Unsimulated_Index(Peak): #returns the index of the maximum value not simulated
    global BigDisplay
    global PeakNumber
    global Automatic
    global PeakDisplayTime
    if Automatic:
        index = 0                    #no check on if you've run out of peaks. it would just return 0
        height=0
        for i in range(PeakNumber):
            if not Peak[i].simulated:
                if Peak[i].height > height:
                    height = Peak[i].height
                    index=i
        return index
    else:
        indices = []
        while 1:
            index = 0                    #no check on if you've run out of peaks. it would just return 0
            height=0
            for i in range(PeakNumber):
                if not Peak[i].simulated and i not in indices:
                    if Peak[i].height > height:
                        #print "highest peak so far is", Peak[i].mz, Peak[i].height
                        height = Peak[i].height
                        index=i
            indices.append(index)
            if len(indices)>=5: break;
        while 1:
            plt.figure()
            plt.axis([0.9*min([Peak[j].centroid for j in indices]), 1.1* max([Peak[j].centroid for j in indices]), 0, 1.2* max([Peak[j].height for j in indices])])
            plt.plot(mzarray, intensityarray, 'k')
            for i in range(len(indices)):
                plt.plot(Peak[indices[i]].centroid, Peak[indices[i]].height, color = ColorDict[i], marker = 'o', linestyle = 'None', label = "Peak %d: %.0f" % (i+1, Peak[indices[i]].centroid))
            plt.legend(loc = 'best')
            if BigDisplay:
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()
                figManager.window.raise_()
            plt.show()
            plt.pause(PeakDisplayTime)
            while 1:
                print "Enter / if you would like to continue viewing. Enter n to stop simulations. Or, enter one of the peaks:", range(1, len(indices)+1)
                entry = raw_input()
                if entry == '/' or entry == 'n' or entry in [str(i+1) for i in range(len(indices))]:
                    break;
                print "Sorry, that is not a valid entry. Please reenter."
            plt.close()
            if entry != '/':
                break;
        if entry == 'n':
            return(0.5)
        else:
            return(indices[int(entry)-1])

def norm(x, height, center, std):
  return(height*np.exp(-(x - center)**2/(2*std**2)))

def Simulate(x, a, b, c, d, e, f, g, h, i, j):#returns y-value of an mz for the combined simulated spectrum
    # global SimulationNumber (is simulation number needed? since all start with curve_scale at 0, it should be fine)
    global Simulation
    y = 0
    y += a * Simulation[0].Function(x)
    y += b * Simulation[1].Function(x)
    y += c * Simulation[2].Function(x)
    y += d * Simulation[3].Function(x)
    y += e * Simulation[4].Function(x)
    y += f * Simulation[5].Function(x)
    y += g * Simulation[6].Function(x)
    y += h * Simulation[7].Function(x)
    y += i * Simulation[8].Function(x)
    y += j * Simulation[9].Function(x)
    return(y)
    
def SimulationFit(): #A least-squares fitting of a linear combination of the simulated envelopes to the spectrum.
                    #It finds the scaling factor for each envelope
    global mzarray
    global intensityarray
    global SimulationNumber
    global Simulation
    popt, pcov = curve_fit(Simulate, mzarray, intensityarray, [1,1,1,1,1,1,1,1,1,1])
    for i in range(SimulationNumber):
        Simulation[i].curve_scale *= popt[i]
        Simulation[i].abundance *= popt[i]
        #print Simulation[i].curve_scale
    
def DifferenceSpectrum(mzval, EnvelopeNumber): #mzval is the m/z region of question
    global mzarray
    global intensityarray
    b = [i != EnvelopeNumber for i in range(10)]
    Spectrum = Simulate(mzval, b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9])
    Spectrum = intensityarray[np.where(mzarray == min(mzval))[0][0] : np.where(mzarray == max(mzval))[0][0] + 1] - Spectrum
    return(Spectrum) 

def find_overlaps(mz, intensity, indices, use_half_left, use_half_right, ysg = None):
    r"""This function uses the second derivative to find the number of underlying
    peaks in the peak range, with a threshold value. Then, it fits up to 6 gaussians
    at a time to the peak range. If more than 6 underlying peaks are found, it
    splits the 6th peak in half and fits the two ranges separately. Then, it
    combines the parameters of the 6th peak. Finally, it returns the peak m/z,
    intensity, and index values in the returnmz, returnintensity, and returnindices
    tuples."""
    
    global PeakNumber
    global DerivThresh
    global NeedZeroCrossing
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    global WeightWidths
    thresh = DerivThresh # for second derivative. multiply to the median of negative derivatives
    height_fraction = 0.001 # for the peak m/z range. if a data point is below height_fraction * the peak height, then it's not in the peak m/z range

    if ysg == None:
        ysg = savitzky_golay(intensity, window_size=5, order=2, deriv = 2)
    thresh *= np.median(np.array([i for i in ysg if i <0]))
    tempind = []
    peakind = []
    ysg_temp = []
    inside_peak= False
    
    for i in range(len(ysg)):
        if ysg[i]<0:
            if 0<i<len(ysg)-1:
                if ysg[i]<thresh and ysg[i] < ysg[i+1] and ysg[i] < ysg[i-1]:
                    tempind.append(i)
                    ysg_temp.append(ysg[i])
                    inside_peak = True
            elif i ==0:
                if ysg[i]<thresh and ysg[i] < ysg[i+1]:
                    tempind.append(i)
                    ysg_temp.append(ysg[i])
                    inside_peak = True
            elif i == len(ysg)-1:
                if ysg[i]<thresh and ysg[i] < ysg[i-1]:
                    tempind.append(i)
                    ysg_temp.append(ysg[i])
                    inside_peak = True
        else:
            if inside_peak:
                if NeedZeroCrossing:
                    peakind.append(tempind[ysg_temp.index(min(ysg_temp))])
                else:
                    peakind.append(j for j in tempind)
                tempind = []
                ysg_temp = []
                inside_peak = False
    if inside_peak: #IN CASE THE LAST ONE IS STIL NEGATIVE, THEN THE LAST INDICES WON'T BE INCLUDED
        peakind.append(tempind[ysg_temp.index(min(ysg_temp))])
        tempind = []
        ysg_temp = []
        inside_peak = False
    

    
    res_use_half_left = use_half_left
    res_use_half_right = use_half_right
        
    if use_half_left:
        if 0 not in peakind:
            peakind.insert(0, 0)
        half_gaussian_center_left = mz[peakind[0]]
    else:
        half_gaussian_center_left = None
    
    if use_half_right:
        if use_half_right not in peakind:
            peakind.append(use_half_right)
        #print "Use half right is", use_half_right
        #print "last mz is", mz[use_half_right]
        #print "length of mz is", len(mz)
        #
        half_gaussian_center_right = mz[peakind[-1]]
    else:
        half_gaussian_center_right = None


    b = range(len(mz))
    step1= mz[0]
    step2 = mz[len(mz)-1]
    for i in b:
        if i != 0 and i != len(b)-1:
            if intensity[i] < max(intensity)/2 < intensity[i+1]:
                step1 = mz[i]
                break;
    b.reverse()
    for j in b:
        if j != 0 and j != len(b)-1:
            if intensity[j] < max(intensity)/2 < intensity[j-1]:
                step2 = mz[j]
                break;
    sd_guess = (step2-step1)/(2*1.1774) #The full width at half height of a Gaussian is 1.1774 standard deviations from the center, so we use this as a guess



    if len(peakind)  == 1:
        return([mz], [intensity], [indices])
    elif len(peakind) == 2:
        sd_guess /=2
        if WeightWidths:
            h1, h2, m1, m2 = [intensity[peakind[0]], intensity[peakind[1]], mz[peakind[0]], mz[peakind[1]]] #initial guesses
            sd1, sd2 = np.array([h1, h2]) * sd_guess / sum(np.array([h1, h2]))
        else:
            h1, h2, m1, m2, sd1, sd2 = [intensity[peakind[0]], intensity[peakind[1]], mz[peakind[0]], mz[peakind[1]], sd_guess, sd_guess] #initial guesses
        try:
            if not (res_use_half_left or res_use_half_right):
                p = [h1, m1, sd1, h2, m2, sd2]
                plsq = leastsq(res2, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5]
            elif res_use_half_left and res_use_half_right:
                p = [h1, sd1, h2, sd2]
                plsq = leastsq(res2, p, args = (intensity, mz))
                h1, sd1, h2, sd2 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3]
            elif res_use_half_left:
                p = [h1, sd1, h2, m2, sd2]
                plsq = leastsq(res2, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4]
            elif res_use_half_right:
                p = [h1, m1, sd1, h2, sd2]
                plsq = leastsq(res2, p, args = (intensity, mz))
                h1, m1, sd1, h2, sd2 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4]      
        except TypeError:
            pass
        
        returnmz = [] #each element is a list of three elements: list of m/zs, list of intensities, list of indices. remove 0s at end.
        returnintensity = []
        returnindices = []            
        for i in range(2):
            returnmz.append([])
            returnintensity.append([])
            returnindices.append([])
        for i in range(len(mz)):
            indiv_intensity = np.zeros(2)
            if norm(mz[i], h1, m1, sd1) > height_fraction*h1:
                indiv_intensity[0] += norm(mz[i], h1, m1, sd1)
            if norm(mz[i], h2, m2, sd2) > height_fraction * h2:
                indiv_intensity[1] += norm(mz[i], h2, m2, sd2)
            for j in range(len(indiv_intensity)):
                if indiv_intensity[j] > 0:
                    returnmz[j].append(mz[i])
                    returnintensity[j].append(intensity[i]*indiv_intensity[j]/sum(indiv_intensity)) #multiply the actual intensity by the ratio
                    returnindices[j].append(i+indices[0])
      
    elif len(peakind) == 3:
        sd_guess/= 3
        if WeightWidths:
            h1, h2, h3, m1, m2, m3 = [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]]] #initial guesses
            sd1, sd2, sd3 = np.array([h1, h2, h3]) * sd_guess / sum(np.array([h1, h2, h3]))
        else:
            h1, h2, h3, m1, m2, m3, sd1, sd2, sd3 = [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], sd_guess, sd_guess, sd_guess] #initial guesses
        try:
            if not (res_use_half_left or res_use_half_right):
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3]
                plsq = leastsq(res3, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8]
            elif res_use_half_left and res_use_half_right:
                p = [h1, sd1, h2, m2, sd2, h3, sd3]
                plsq = leastsq(res3, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, sd3 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6]
            elif res_use_half_left:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3]
                plsq = leastsq(res3, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7]
            elif res_use_half_right:
                p = [h1, m1, sd1, h2, m2, sd2, h3, sd3]
                plsq = leastsq(res3, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, sd3 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7]
        except TypeError:
            pass
        returnmz = [] #each element is a list of three elements: list of m/zs, list of intensities, list of indices. remove 0s at end.
        returnintensity = []
        returnindices = []            
        for i in range(3):
            returnmz.append([])
            returnintensity.append([])
            returnindices.append([])
        for i in range(len(mz)):
            indiv_intensity = np.zeros(3)
            if norm(mz[i], h1, m1, sd1) > height_fraction*h1:
                indiv_intensity[0] += norm(mz[i], h1, m1, sd1)
            if norm(mz[i], h2, m2, sd2) > height_fraction * h2:
                indiv_intensity[1] += norm(mz[i], h2, m2, sd2)
            if norm(mz[i], h3, m3, sd3) > height_fraction * h3:
                indiv_intensity[2] += norm(mz[i], h3, m3, sd3)
            for j in range(len(indiv_intensity)):
                if indiv_intensity[j] > 0:
                    returnmz[j].append(mz[i])
                    returnintensity[j].append(intensity[i]*indiv_intensity[j]/sum(indiv_intensity)) #multiply the actual intensity by the ratio
                    returnindices[j].append(i+indices[0])
    elif len(peakind) == 4:
        sd_guess /= 4
        if WeightWidths:
            h1, h2, h3, h4, m1, m2, m3, m4= [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]],intensity[peakind[3]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], mz[peakind[3]] ] #initial guesses
            sd1, sd2, sd3, sd4 = np.array([h1, h2, h3, h4]) * sd_guess / sum(np.array([h1, h2, h3, h4]))
        else:
            h1, h2, h3, h4, m1, m2, m3, m4, sd1, sd2, sd3, sd4 = [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]],intensity[peakind[3]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], mz[peakind[3]], sd_guess, sd_guess, sd_guess, sd_guess] #initial guesses
        try:
            if not (res_use_half_left or res_use_half_right):
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4]
                plsq = leastsq(res4, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11]
            elif res_use_half_left and res_use_half_right:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, sd4]
                plsq = leastsq(res4, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, sd4 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9]
            elif res_use_half_left:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4]
                plsq = leastsq(res4, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10]
            elif res_use_half_right:
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, sd4]
                plsq = leastsq(res4, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, sd4 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10]
        except TypeError:
            pass
        
        
        returnmz = [] #each element is a list of three elements: list of m/zs, list of intensities, list of indices. remove 0s at end.
        returnintensity = []
        returnindices = []            
        for i in range(4):
            returnmz.append([])
            returnintensity.append([])
            returnindices.append([])
        for i in range(len(mz)):
            indiv_intensity = np.zeros(4)
            if norm(mz[i], h1, m1, sd1) > height_fraction*h1:
                indiv_intensity[0] += norm(mz[i], h1, m1, sd1)
            if norm(mz[i], h2, m2, sd2) > height_fraction * h2:
                indiv_intensity[1] += norm(mz[i], h2, m2, sd2)
            if norm(mz[i], h3, m3, sd3) > height_fraction * h3:
                indiv_intensity[2] += norm(mz[i], h3, m3, sd3)
            if norm(mz[i], h4, m4, sd4) > height_fraction * h4:
                indiv_intensity[3] += norm(mz[i], h4, m4, sd4)
            for j in range(len(indiv_intensity)):
                if indiv_intensity[j] > 0:
                    returnmz[j].append(mz[i])
                    returnintensity[j].append(intensity[i]*indiv_intensity[j]/sum(indiv_intensity)) #multiply the actual intensity by the ratio
                    returnindices[j].append(i+indices[0])
        
    elif len(peakind) == 5:
        sd_guess /= 5
        if WeightWidths:
            h1, h2, h3, h4, h5, m1, m2, m3, m4, m5= [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]],intensity[peakind[3]], intensity[peakind[4]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], mz[peakind[3]], mz[peakind[4]]] #initial guesses
            sd1, sd2, sd3, sd4, sd5 = np.array([h1, h2, h3, h4, h5]) * sd_guess / sum(np.array([h1, h2, h3, h4, h5]))
        else:
            h1, h2, h3, h4, h5, m1, m2, m3, m4, m5, sd1, sd2, sd3, sd4, sd5 = [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]],intensity[peakind[3]], intensity[peakind[4]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], mz[peakind[3]], mz[peakind[4]], sd_guess, sd_guess, sd_guess, sd_guess, sd_guess] #initial guesses
        try:
            if not (res_use_half_left or res_use_half_right):
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5]
                plsq = leastsq(res5, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13], plsq[0][14]
            elif res_use_half_left and res_use_half_right:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, sd5]
                plsq = leastsq(res5, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, sd5 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12]
            elif res_use_half_left:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5]
                plsq = leastsq(res5, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13]
            elif res_use_half_right:
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, sd5]
                plsq = leastsq(res5, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, sd5 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13]
        except TypeError:
            pass        

        
        returnmz = [] #each element is a list of three elements: list of m/zs, list of intensities, list of indices. remove 0s at end.
        returnintensity = []
        returnindices = []            
        for i in range(5):
            returnmz.append([])
            returnintensity.append([])
            returnindices.append([])
        for i in range(len(mz)):
            indiv_intensity = np.zeros(5)
            if norm(mz[i], h1, m1, sd1) > height_fraction*h1:
                indiv_intensity[0] += norm(mz[i], h1, m1, sd1)
            if norm(mz[i], h2, m2, sd2) > height_fraction * h2:
                indiv_intensity[1] += norm(mz[i], h2, m2, sd2)
            if norm(mz[i], h3, m3, sd3) > height_fraction * h3:
                indiv_intensity[2] += norm(mz[i], h3, m3, sd3)
            if norm(mz[i], h4, m4, sd4) > height_fraction * h4:
                indiv_intensity[3] += norm(mz[i], h4, m4, sd4)
            if norm(mz[i], h5, m5, sd5) > height_fraction * h5:
                indiv_intensity[4] += norm(mz[i], h5, m5, sd5)
            for j in range(len(indiv_intensity)):
                if indiv_intensity[j] > 0:
                    returnmz[j].append(mz[i])
                    returnintensity[j].append(intensity[i]*indiv_intensity[j]/sum(indiv_intensity)) #multiply the actual intensity by the ratio
                    returnindices[j].append(i+indices[0])
        
    elif len(peakind) == 6:
        sd_guess /= 6

        if WeightWidths:
            h1, h2, h3, h4, h5, h6, m1, m2, m3, m4, m5, m6= [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]],intensity[peakind[3]], intensity[peakind[4]], intensity[peakind[5]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], mz[peakind[3]], mz[peakind[4]], mz[peakind[5]]] #initial guesses
            sd1, sd2, sd3, sd4, sd5, sd6 = np.array([h1, h2, h3, h4, h5, h6]) * sd_guess / sum(np.array([h1, h2, h3, h4, h5, h6]))
        else:
            h1, h2, h3, h4, h5, h6, m1, m2, m3, m4, m5, m6, sd1, sd2, sd3, sd4, sd5, sd6 = [intensity[peakind[0]], intensity[peakind[1]], intensity[peakind[2]],intensity[peakind[3]], intensity[peakind[4]], intensity[peakind[5]], mz[peakind[0]], mz[peakind[1]], mz[peakind[2]], mz[peakind[3]], mz[peakind[4]], mz[peakind[5]], sd_guess, sd_guess, sd_guess, sd_guess, sd_guess, sd_guess] #initial guesses
        try:
            if not (res_use_half_left or res_use_half_right):
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, m6, sd6]
                plsq = leastsq(res6, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, m6, sd6 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13], plsq[0][14], plsq[0][15], plsq[0][16], plsq[0][17]
            elif res_use_half_left and res_use_half_right:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, sd6]
                plsq = leastsq(res6, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, sd6 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13], plsq[0][14], plsq[0][15]
            elif res_use_half_left:
                p = [h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, m6, sd6]
                plsq = leastsq(res6, p, args = (intensity, mz))
                h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, m6, sd6 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13], plsq[0][14], plsq[0][15], plsq[0][16]
            elif res_use_half_right:
                p = [h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, sd6]
                plsq = leastsq(res6, p, args = (intensity, mz))
                h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, sd6 = plsq[0][0], plsq[0][1], plsq[0][2], plsq[0][3], plsq[0][4], plsq[0][5], plsq[0][6], plsq[0][7], plsq[0][8], plsq[0][9], plsq[0][10], plsq[0][11], plsq[0][12], plsq[0][13], plsq[0][14], plsq[0][15], plsq[0][16]
        except TypeError:
            pass

        
        returnmz = [] #each element is a list of three elements: list of m/zs, list of intensities, list of indices. remove 0s at end.
        returnintensity = []
        returnindices = []            
        for i in range(6):
            returnmz.append([])
            returnintensity.append([])
            returnindices.append([])
        for i in range(len(mz)):
            indiv_intensity = np.zeros(6)
            if norm(mz[i], h1, m1, sd1) > height_fraction*h1:
                indiv_intensity[0] += norm(mz[i], h1, m1, sd1)
            if norm(mz[i], h2, m2, sd2) > height_fraction * h2:
                indiv_intensity[1] += norm(mz[i], h2, m2, sd2)
            if norm(mz[i], h3, m3, sd3) > height_fraction * h3:
                indiv_intensity[2] += norm(mz[i], h3, m3, sd3)
            if norm(mz[i], h4, m4, sd4) > height_fraction * h4:
                indiv_intensity[3] += norm(mz[i], h4, m4, sd4)
            if norm(mz[i], h5, m5, sd5) > height_fraction * h5:
                indiv_intensity[4] += norm(mz[i], h5, m5, sd5)
            if norm(mz[i], h6, m6, sd6) > height_fraction * h6:
                indiv_intensity[5] += norm(mz[i], h6, m6, sd6)
            for j in range(len(indiv_intensity)):
                if indiv_intensity[j] > 0:
                    returnmz[j].append(mz[i])
                    returnintensity[j].append(intensity[i]*indiv_intensity[j]/sum(indiv_intensity)) #multiply the actual intensity by the ratio
                    returnindices[j].append(i+indices[0])        
        
    elif len(peakind)> 6:
        print "\nMore than 6 overlapping peaks; refitting.\n"
        mz1 = mz[: peakind[5] + 1]
        intensity1 = intensity[: peakind[5] + 1]
        indices1 = indices[: peakind[5] + 1]
        #LEFT OFF HERE 8_25_14
        ysg1 = ysg[: peakind[5] + 1]
        
        mz2 = mz[peakind[5] : ]
        intensity2 = intensity[peakind[5] : ]
        indices2 = indices[peakind[5] : ]
        ysg2 = ysg[peakind[5] : ]
        
        split1 = find_overlaps(mz1, intensity1, indices1, use_half_left, peakind[5], ysg = ysg1)
        
        if use_half_right:
            split2 = find_overlaps(mz2, intensity2, indices2, True, use_half_right-peakind[5], ysg = ysg2)
        else:
            split2 = find_overlaps(mz2, intensity2, indices2, True, use_half_right, ysg = ysg2)
        
        #print "the first split:", split1[0]
        #print "end of first split:", split1[0][-1]
        #print "begining of second split:", split2[0][0][1:]
        returnmz = list(split1[0][:-1]) + [split1[0][-1] +split2[0][0][1:]]+ list(split2[0][1:]) #YOU LEFT OFF HERE 8_17_13. Issues to be resolved: clean the middle half one. Recheck if residual functions work. also, turn tuple into list.
        returnintensity = list(split1[1][:-1]) + [split1[1][-1] +split2[1][0][1:]]+ list(split2[1][1:]) #REPLACE THE INTENSITY AT THE HALF CENTER CALCULATED BY THE RIGHT HALF
        returnindices = list(split1[2][:-1]) + [split1[2][-1] +split2[2][0][1:]]+ list(split2[2][1:])
    

    
    return(returnmz, returnintensity, returnindices)
    
    

def PeakSetup(Peak, mz):
    global mzarray
    global intensityarray
    global PeakNumber
    PeakNumber = 0
    mz.sort()
    adjacent_number = 0 # limits the number of adjacent
    for i in range(len(mz)):
        center = None
        for k in range(len(mzarray)):
            if mzarray[k]<=mz[i]<mzarray[k+1]:
                center = k
                break;
        counter = [center, center]
        if center == None:
            print "uhoh, center is None in PeakSetup"
        if intensityarray[center] <intensityarray[center + 1] and intensityarray[center]<intensityarray[center-1]:
            print 'PEAK SETUP ERROR: FOUND PEAK M/Z', mzarray[center], 'IS LOCAL MINIMUM'
            center1 = center - 1 #If Local Minimum, set up both the left peak and the right peak
            center2 = center + 1
            counter1 = [center1, center1]
            counter2 = [center2, center2]
            
            #LEFT PEAK
            while 1:
                counter1[0] -= 1
                if counter1[0] <= 0:
                    break;
                if intensityarray[counter1[0]] <= 0.01:
                    break;
                if intensityarray[counter1[0]] <intensityarray[counter1[0] +1] and intensityarray[counter1[0]] <intensityarray[counter1[0] -1]:
                    break;
            while 1:
                counter1[1] += 1
                if counter1[1] >= len(intensityarray) -1:
                    break;
                if intensityarray[counter1[1]] <= 0.01:
                    break;
                if intensityarray[counter1[1]] <intensityarray[counter1[1] +1] and intensityarray[counter1[1]] <intensityarray[counter1[1] -1]:
                    break;
            for j in range(counter1[0], counter1[1]+1):
                Peak[PeakNumber].mz.append(mzarray[j])
                Peak[PeakNumber].intensity.append(intensityarray[j])
                Peak[PeakNumber].indices.append(j)
            Peak[PeakNumber].is_empty = False
            if PeakNumber >0:
                if Peak[PeakNumber].intensity[0] >=0.001 and min(Peak[PeakNumber].indices) == max(Peak[PeakNumber -1].indices): #find adjacent peaks
                    Peak[PeakNumber].left = PeakNumber -1
                    Peak[PeakNumber-1].right = PeakNumber
                    adjacent_number += 1
                else:
                    adjacent_number = 0
            PeakNumber +=1
            
            #RIGHT PEAK
            while 1:
                counter2[0] -= 1
                if counter2[0] <= 0:
                    break;
                if intensityarray[counter2[0]] <= 0.01:
                    break;
                if intensityarray[counter2[0]] <intensityarray[counter2[0] +1] and intensityarray[counter2[0]] <intensityarray[counter2[0] -1]:
                    break;
            while 1:
                counter2[1] += 1
                if counter2[1] >= len(intensityarray) -1:
                    break;
                if intensityarray[counter2[1]] <= 0.01:
                    break;
                if intensityarray[counter2[1]] <intensityarray[counter2[1] +1] and intensityarray[counter2[1]] <intensityarray[counter2[1] -1]:
                    break;
            for j in range(counter2[0], counter2[1]+1):
                Peak[PeakNumber].mz.append(mzarray[j])
                Peak[PeakNumber].intensity.append(intensityarray[j])
                Peak[PeakNumber].indices.append(j)
            Peak[PeakNumber].is_empty = False
            if PeakNumber >0:
                if Peak[PeakNumber].intensity[0] >=0.001 and min(Peak[PeakNumber].indices) == max(Peak[PeakNumber -1].indices): #find adjacent peaks
                    Peak[PeakNumber].left = PeakNumber -1
                    Peak[PeakNumber-1].right = PeakNumber
                    adjacent_number += 1
                else:
                    adjacent_number = 0
            PeakNumber +=1
            
            
            
            #------------------- 
        else:
            while 1:
                counter[0] -= 1
                if counter[0] <= 0:
                    break;
                if intensityarray[counter[0]] <= 0.01:
                    break;
                if intensityarray[counter[0]] <intensityarray[counter[0] +1] and intensityarray[counter[0]] <intensityarray[counter[0] -1]:
                    break;
            while 1:
                counter[1] += 1
                if counter[1] >= len(intensityarray) -1:
                    break;
                if intensityarray[counter[1]] <= 0.01:
                    break;
                if intensityarray[counter[1]] <intensityarray[counter[1] +1] and intensityarray[counter[1]] <intensityarray[counter[1] -1]:
                    break;
            for j in range(counter[0], counter[1]+1):
                Peak[PeakNumber].mz.append(mzarray[j])
                Peak[PeakNumber].intensity.append(intensityarray[j])
                Peak[PeakNumber].indices.append(j)
            Peak[PeakNumber].is_empty = False
            if PeakNumber >0:
                if Peak[PeakNumber].intensity[0] >=0.001 and min(Peak[PeakNumber].indices) == max(Peak[PeakNumber -1].indices): #find adjacent peaks
                    Peak[PeakNumber].left = PeakNumber -1
                    Peak[PeakNumber-1].right = PeakNumber
                    adjacent_number += 1
                else:
                    adjacent_number = 0
            PeakNumber +=1

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """


    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def res2(p, y, x): #This is the residual function for fitting two gaussians to a peak range.
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    if not (res_use_half_left or res_use_half_right):
        h1, m1, sd1, h2, m2, sd2 = p
    elif res_use_half_left and res_use_half_right:
        h1, sd1, h2, sd2 = p
        m1 = half_gaussian_center_left
        m2 = half_gaussian_center_right
    elif res_use_half_left:
        h1, sd1, h2, m2, sd2 = p
        m1 = half_gaussian_center_left
    elif res_use_half_right:
        h1, m1, sd1, h2, sd2 = p
        m2 = half_gaussian_center_right
    y_fit = norm(x, h1, m1, sd1) + norm(x, h2, m2, sd2)
    err = y - y_fit
    return err

def res3(p, y, x):
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    if not (res_use_half_left or res_use_half_right):
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3 = p
    elif res_use_half_left and res_use_half_right:
        h1, sd1, h2, m2, sd2, h3, sd3 = p
        m1 = half_gaussian_center_left
        m3 = half_gaussian_center_right
    elif res_use_half_left:
        h1, sd1, h2, m2, sd2, h3, m3, sd3 = p
        m1 = half_gaussian_center_left
    elif res_use_half_right:
        h1, m1, sd1, h2, m2, sd2, h3, sd3 = p
        m3 = half_gaussian_center_right
    
    y_fit = norm(x, h1, m1, sd1) + norm(x, h2, m2, sd2) + norm(x, h3, m3, sd3)
    err = y - y_fit
    return err
  
def res4(p, y, x):
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    if not (res_use_half_left or res_use_half_right):
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4 = p
    elif res_use_half_left and res_use_half_right:
        h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, sd4 = p
        m1 = half_gaussian_center_left
        m4 = half_gaussian_center_right
    elif res_use_half_left:
        h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4 = p
        m1 = half_gaussian_center_left
    elif res_use_half_right:
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, sd4 = p
        m4 = half_gaussian_center_right
    y_fit = norm(x, h1, m1, sd1) + norm(x, h2, m2, sd2) + norm(x, h3, m3, sd3) + norm(x, h4, m4, sd4)
    err = y - y_fit
    return err

def res5(p, y, x):
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    if not (res_use_half_left or res_use_half_right):
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5 = p
    elif res_use_half_left and res_use_half_right:
        h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, sd5 = p
        m1 = half_gaussian_center_left
        m5 = half_gaussian_center_right
    elif res_use_half_left:
        h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5 = p
        m1 = half_gaussian_center_left
    elif res_use_half_right:
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, sd5 = p
        m5 = half_gaussian_center_right
    y_fit = norm(x, h1, m1, sd1) + norm(x, h2, m2, sd2) + norm(x, h3, m3, sd3) + norm(x, h4, m4, sd4) + norm(x, h5, m5, sd5)
    err = y - y_fit
    return err
  
def res6(p, y, x):
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    if not (res_use_half_left or res_use_half_right):
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, m6, sd6 = p
    elif res_use_half_left and res_use_half_right:
        h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, sd6 = p
        m1 = half_gaussian_center_left
        m6 = half_gaussian_center_right
    elif res_use_half_left:
        h1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, m6, sd6 = p
        m1 = half_gaussian_center_left
    elif res_use_half_right:
        h1, m1, sd1, h2, m2, sd2, h3, m3, sd3, h4, m4, sd4, h5, m5, sd5, h6, sd6 = p
        m6 = half_gaussian_center_right
    y_fit = norm(x, h1, m1, sd1) + norm(x, h2, m2, sd2) + norm(x, h3, m3, sd3) + norm(x, h4, m4, sd4) + norm(x, h5, m5, sd5) + norm(x, h6, m6, sd6)
    err = y - y_fit
    return err

        
def MakeFile(file_name, revise_name, beginning, mzarray, DiffSpec): #Writes the subtracted spectrum into another text file.
	temp_path = file_name.partition('.txt')[0] + '_' + revise_name + '.txt'
	file = open(temp_path, 'w')
	file.write(beginning)
	for i in range(len(mzarray)):
	    file.write(str(mzarray[i]))
	    file.write('\t')
	    file.write(str(DiffSpec[i]))
	    file.write('\n')
	file.close()
	print 'Subtracted spectrum saved to directory ', temp_path   


def main():
    
    #set up constants
    t = time.time()
    global threshold
    global adj_factor
    global max_mass
    global min_mass
    global max_charge
    global min_charge
    global PeakDisplayTime
    global ChargeDisplayTime
    global BigDisplay
    global Automatic
    MaxWidth = None
    MinWidth = None

    
    global PeakNumber
    global MaxSimulations
    global SimulationNumber
    global Simulation
    

    
    global DerivThresh
    global NeedZeroCrossing
    global WeightWidths
    
    global ColorDict
    ColorDict = dict([(0, 'r'), (1 , 'y'), (2, 'g'), (3, 'c'), (4,'m'), (5, 'b'), (6, 'purple'), (7, 'SaddleBrown'), (8, 'chartreuse'), (9, 'orange') ])
    global LetterDict
    LetterDict = dict([(0, 'A'), (1 , 'B'), (2, 'C'), (3, 'D'), (4,'E'), (5, 'F'), (6, 'G'), (7, 'H'), (8,'I'), (9, 'I')])
    
    global half_gaussian_center_left
    global half_gaussian_center_right
    global res_use_half_left
    global res_use_half_right
    
    global PeakDetectionType
    global MassTolerance
    global ScoreLimit
    global MassErrWeight
    global HeightErrWeight
    global StdDevErrWeight
    global AttemptLimit
    global RuntimeErrorLimit
    global min_peak_number
    
    global StdDevStartGuess
    
    opener = open('/Users/jlu96/Desktop/Python27/OPTIONS.txt') #THIS MUST BE THE DIRECTORY OF THE OPTIONS FILE THAT CAME WITH THE PROGRAM.
    index = 0
    for line in opener:
        value = line.partition('\t')[2].partition('\t')[0]
        #print "For index", index, "value is", value
        if index == 0:
            name = line.partition('\t')[2].partition('\n')[0]
        elif index == 1:
            Automatic = eval(value)
        elif index == 2:
            PeakDisplayTime = eval(value)
        elif index == 3:
            ChargeDisplayTime = eval(value)
        elif index == 4:
            BigDisplay = eval(value)
        elif index == 5:
            MaxMz = eval(value)
        elif index == 6:
            MinMz = eval(value)
        elif index == 7:
            max_mass = eval(value)
        elif index == 8:
            min_mass = eval(value)
        elif index == 9:
            max_charge = eval(value)
        elif index == 10:
            min_charge = eval(value)
        elif index == 11:
            MaxWidth = eval(value)
        elif index == 12:
            MinWidth = eval(value)

        #PROCESSING OPTIONS BEING READ IN
        
        elif index == 16:
            BkgrdRemove = eval(value)
        elif index == 17:
            background_remove_window = eval(value)
        elif index == 18:
            background_smooth_window = eval(value)
        elif index == 19:
            DoSmooth = eval(value)
        elif index == 20:
            smooth_window = eval(value)
        elif index == 21:
            DoSavitzky = eval(value)
        elif index == 22:
            Savitz_Window = eval(value)
        elif index == 23:
            Savitz_Order = eval(value)
        elif index == 24:
            Savitz_Times = eval(value)
        elif index == 25:
            UseProcessed = eval(value)
        
        #PEAK DETECTION OPTIONS BEING READ IN
            
        elif index == 29:
            PeakDetectionType = value
        elif index == 30:
            threshold = eval(value)
        elif index == 31:
            JustPeaks = eval(value)
        elif index == 32:
            FindOverlaps = eval(value)
        elif index == 33:
            DerivThresh = eval(value)
        elif index == 34:
            NeedZeroCrossing = eval(value)
        elif index == 35:
            WeightWidths = eval(value)
        elif index == 36:
            adj_factor = eval(value)
        elif index == 37:
            cwt_widths = np.arange(eval(value.partition('-')[0]), eval(value.partition('-')[2].partition(',')[0]), eval(value.partition(',')[2]))
        elif index == 38:
            cwt_wavelet = eval(value) 
        elif index == 39:
            cwt_max_distances = eval(value) 
        elif index == 40:
            cwt_gap_thresh = eval(value)  
        elif index == 41:
            cwt_min_length = eval(value)
        elif index == 42:
            cwt_min_snr = eval(value)
        elif index == 43:
            cwt_noise_perc = eval(value)
        
        #CHARGE STATE ASSIGNMENT OPTIONS BEING READ IN
        
        elif index == 47:
            MassTolerance = eval(value)
        elif index == 48:
            ScoreLimit = eval(value)
        elif index == 49:
            MassErrWeight = eval(value)
        elif index == 50:
            HeightErrWeight = eval(value)
        elif index == 51:
            StdDevErrWeight = eval(value)
        elif index == 52:
            AttemptLimit = eval(value)
        elif index == 53:
            RuntimeErrorLimit = eval(value)
        elif index == 54:
            min_peak_number = eval(value)
        elif index == 55:
            MaxSimulations = eval(value)
        elif index == 56:
            RepeatSearch = eval(value)
        elif index == 57:
            UseSubtract = eval(value)
        elif index == 58:
            SubtractAuto = eval(value)
        elif index == 59:
            NumberRefit = eval(value)
        elif index == 60:
            StdDevStartGuess = eval(value)    
            
        #SAVING SUBTRACTED SPECTRUM BEING READ IN
        
        elif index == 64:
            SaveSubtract = eval(value)
        elif index == 65:
            NegToZeros = eval(value)
        elif index == 66:
            SubtractName = value
        
        #SAVING MASS INFORMATION BEING READ IN
        
        elif index == 69:
            SaveMasses = eval(value)
        elif index == 70:
            MassName = value
        index += 1
    opener.close()

    global mzarray
    global intensityarray
    mzarray = []
    intensityarray = []

    #READ IN RAW DATA INTO MZARRAY AND INTENSITYARRAY***************************
    index = 0
    global beginning
    beginning = ''
    opener = open(name)
    for line in opener:
        if index <8:
            beginning += line
        else:
            
            a = line.partition('\t')
            mzarray.append(float(a[0]))
            b= a[2].partition('\n')
            intensityarray.append(float(b[0]))
        index+=1;
    opener.close()
    
    #in case you want to zoom
    if MinMz != None:
        for i in range(len(mzarray)):
            if mzarray[i] >= MinMz:
                mzarray = mzarray[i:]
                intensityarray = intensityarray[i:]
                break;
        #if MinMz and MaxMz are too big it doesn't matter-- nothing happens to the arrays
    if MaxMz != None:
        for i in range(len(mzarray)):
            if mzarray[i] >= MaxMz:
                mzarray = mzarray[:i]
                intensityarray = intensityarray[:i]
                break;

    mzarray = np.array(mzarray)
    intensityarray=np.array(intensityarray)
    
    t_read = time.time()
    
    #SPECTRUM PROCESSING********************************************************
    processed_intensity = list(intensityarray[:]) #This is just used as the intensity value for peak detection, whether it is actually processed or not.
    if DoSmooth:
        processed_intensity = Smooth(processed_intensity, smooth_window)
    if DoSavitzky:
        for i in range(Savitz_Times):
            processed_intensity = savitzky_golay(processed_intensity, Savitz_Window*2 + 1, Savitz_Order)
    if BkgrdRemove:
        processed_intensity = BackgroundRemove(processed_intensity, background_remove_window, background_smooth_window)    
        #plt.figure()
        #plt.plot(mzarray, processed_intensity)
    
    t_process = time.time()
    
    #PEAK INITIALIZATION *******************************************************
    global Peak
    Peak = range(10000) #Peak is a list of Peak objects which are indexed.
    for i in range(10000):
   	Peak[i] = Peaks()
   	Peak[i]._init_()
    threshold *= max(processed_intensity) #gives the threshold a value on the intensity
    
    if UseProcessed:
        intensityarray = np.array(processed_intensity[:])
    
    #PEAK FINDING ************************************************************* SHOULD RETURN LIST OF M/Z VALUES AT WHICH THERE IS A PEAK.
    if PeakDetectionType == 'i': #Local Maxima Detection
        mz = []
        for i in range(len(mzarray)-1):
            if processed_intensity[i] > threshold and processed_intensity[i] > processed_intensity[i+1] and processed_intensity[i] > processed_intensity[i-1]:
                mz.append(mzarray[i])
    
    elif PeakDetectionType == 'c': #Continuous Wavelet Transform
        mz = []
        
        peakind = signal.find_peaks_cwt(processed_intensity, cwt_widths, cwt_wavelet, cwt_max_distances, cwt_gap_thresh, cwt_min_length, cwt_min_snr, cwt_noise_perc)
        i = 0
        #print "peak indices were", peakind
        #This piece of code removes extra ones that are right next to each other.
        newpeakind = []
        for k, g in groupby(enumerate(peakind), lambda (i,x):i-x):
            b = map(operator.itemgetter(1), g)
            newpeakind.append(b[len(b)/2])
        
        peakind = newpeakind
        
        for i in range(len(peakind)):
            if processed_intensity[peakind[i]] > threshold:
                mz.append(mzarray[peakind[i]])
    
    elif PeakDetectionType == 'm': #Massign Peak Detection: Fixed and Adjusted Threshold
        inside_peak = False
        PeakNumber = 0
        for i in range(len(mzarray)):
            if processed_intensity[i] > threshold:
                inside_peak = True
                Peak[PeakNumber].mz.append(mzarray[i])
                Peak[PeakNumber].intensity.append(processed_intensity[i])
                Peak[PeakNumber].indices.append(i)
            else:
                if inside_peak:
                    PeakNumber += 1
                    inside_peak = False
    
        for i in range(PeakNumber):
            new_indices = Peak[i].find_shoulders()
            if len(new_indices)>1: #checks if any shoulder peaks were found
                Peak[i].indices = new_indices[0]
                Peak[i].mz = []
                Peak[i].intensity = []
                for j in Peak[i].indices:
                    Peak[i].mz.append(mzarray[j])
                    Peak[i].intensity.append(processed_intensity[j])
                for k in range(1, len(new_indices)):
                    Peak[PeakNumber].indices = new_indices[k]
                    for l in new_indices[k]:
                        Peak[PeakNumber].mz.append(mzarray[l])
                        Peak[PeakNumber].intensity.append(processed_intensity[l])
                    PeakNumber +=1
            elif len(new_indices)==1:
                Peak[i].indices = new_indices[0]
                Peak[i].mz = []
                Peak[i].intensity = []
                for j in Peak[i].indices:
                    Peak[i].mz.append(mzarray[j])
                    Peak[i].intensity.append(processed_intensity[j])
            else:
                pass
    else:
        print "ERROR: INVALID PEAK DETECTION TYPE."

    #PEAK SETUP*****************************************************************
    if PeakDetectionType == 'i' or PeakDetectionType == 'c':
        PeakSetup(Peak,mz) #sorts as well
        print "Number of Peaks initially found is", PeakNumber
    
    t_peak = time.time()
    
    #PEAK OVERLAP FINDING ******************************************************
    
    if FindOverlaps: #"  and (PeakDetectionType != 'm'):
        print "Now, looking for overlaps."
        i = 0
        while 1:
            if i >= PeakNumber:
                break;
            left_index = i
            while 1:
                if Peak[i].right is None:
                    break;
                i+=1
            right_index = i
            mz = []
            intensity = []
            indices = []
            half_gaussian_center_left = None
            half_gaussian_center_right = None
            res_use_half_left = False
            res_use_half_right = False
            
            
            if Peak[left_index].is_empty != True:
                for j in range(left_index, right_index+1):
                    mz += Peak[j].mz[:-1]
                    #print "Peak", j, "mz is", Peak[j].mz
                    intensity += Peak[j].intensity[:-1]
                    indices += Peak[j].indices[:-1]
                    if j == right_index:
                        mz.append(Peak[j].mz[-1])
                        intensity.append(Peak[j].intensity[-1])
                        indices.append(Peak[j].indices[-1])
                newmz, newintensity, newindices = find_overlaps(mz, intensity, indices, None, None)

                k= 0
                while 1:
                    if k>=len(newmz):
                        break; #THIS PIECE OF CODE REMOVES THE PEAKS FROM CONSIDERATION IF NOTHING IS FOUND
                    if newmz[k] == []:
                        #print "Empty peak removed from consideration." #THIS PIECE OF CODE REMOVES THE PEAKS FROM CONSIDERATION IF NOTHING IS FOUND
                        newmz.remove([])
                        newintensity.remove([])
                        newindices.remove([])
                        newmz.append([0])
                        newintensity.append([1])
                        newindices.append([0])
                    else:
                        k+=1
    
                for l in range(left_index, right_index+1):
                    if l - left_index >= len(newmz): #this is the case where it finds less peaks than we had when we inputted. mark off remaining peaks as empty.
                        Peak[l].mz = [0]
                        Peak[l].intensity = [1]
                        Peak[l].indices = [0]
                        Peak[l].is_empty = True
                    else:
                        Peak[l].mz = newmz[l-left_index]
                        Peak[l].intensity = newintensity[l-left_index]
                        Peak[l].indices = newindices[l-left_index]
                
                if len(newmz) > right_index-left_index + 1: #THIS PIECE OF CODE AGAIN REMOVES THE PEAKS FROM CONSIDERATION IF NOTHING IS FOUND
                    for p in range(right_index-left_index+1, len(newmz)):
                        Peak[PeakNumber].mz = newmz[p]
                        Peak[PeakNumber].intensity = newintensity[p]
                        Peak[PeakNumber].indices = newindices[p]
                        if newmz[p] == [0] or newmz[p] == []:
                            Peak[PeakNumber].is_empty = True
                            Peak[PeakNumber].mz = [0]
                            Peak[PeakNumber].intensity = [1]
                            Peak[PeakNumber].indices = [0]
                            print "Empty peak removed from consideration."
                            print "Updated number of peaks is", PeakNumber
                        PeakNumber+=1
            i += 1
        
        i = 0
        while 1:
            if i >=PeakNumber:
                break;
            if Peak[i].mz == [0] or Peak[i].mz == []:
                Peak.remove(Peak[i])
                PeakNumber -=1
                print "Empty peak removed from consideration."
                print "Updated number of peaks is", PeakNumber
            else:
                i+=1
        

    #PEAK HEIGHT, CENTROID, AND WIDTH CALCULATION*******************************
    for i in range(PeakNumber):
        Peak[i].h()
        Peak[i].c()
        Peak[i].w()
    
    i = 0
    while 1:
        if i >=PeakNumber:
            break;
        if (MinWidth != None and Peak[i].width < MinWidth)  or (MaxWidth != None and Peak[i].width > MaxWidth):
            Peak.remove(Peak[i])
            PeakNumber -=1
            print "Peak of width outside range removed."
            print "Updated number of peaks is", PeakNumber
        else:
            i+=1

        
    #CHECK PEAK FITTING*********************************************************
    for i in range(PeakNumber):
        Peak[i].PeakFit(False, 0)
    
    t_overlap = time.time()
    
    #JUST PEAK DISPLAYS*********************************************************    

    if JustPeaks:
        print "The threshold intensity is ", threshold
        print "The number of peaks found is ", PeakNumber
        print "\n"
        if FindOverlaps: #Displays peaks found and second derivative
            plt.figure()
            
            ax = plt.subplot(311)
            plt.title('Original Peak')
            plt.axis( [min(mzarray), max(mzarray), min(intensityarray)/max(intensityarray), 1.2])
            plt.plot(mzarray, intensityarray/max(intensityarray), 'k')
            plt.ylabel('Relative Intensity')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            
            
            ax = plt.subplot(312)
            plt.title('Second Derivative')
            secderiv = savitzky_golay(intensityarray, window_size=5, order=2, deriv = 2)
            plt.axis( [min(mzarray), max(mzarray), -1.2, 1.2])
            plt.plot(mzarray, secderiv/max(abs(secderiv)), 'm')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            plt.ylabel('Relative Intensity')            
            
            ax = plt.subplot(313)
            plt.title('Fitted Gaussians')
            plt.axis( [min(mzarray), max(mzarray), min(intensityarray)/max(intensityarray), 1.2])
            plt.plot(mzarray, intensityarray/max(intensityarray), 'k')
            for i in range(PeakNumber):
                Peak[i].h()
                Peak[i].c()
                Peak[i].PeakFit(False, 0)
                inputs, outputs = Peak[i].PeakPlot()
                plt.plot(inputs, outputs/max(intensityarray), 'c')
                plt.plot(Peak[i].fit_center, Peak[i].fit_height/max(intensityarray), 'r.')
            plt.xlabel('m/z')
            plt.ylabel('Relative Intensity')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            if BigDisplay:
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()
                figManager.window.raise_()
            plt.show()
        
        else: #Just displays found peaks
            plt.figure()
            ax = plt.subplot(1,1,1)
            plt.axis( [min(mzarray), max(mzarray), min(intensityarray)/max(intensityarray), 1.2])
            plt.plot(mzarray, intensityarray/max(intensityarray), 'k')
            for i in range(PeakNumber):
                plt.plot(Peak[i].centroid, Peak[i].height/max(intensityarray), 'ro')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            plt.title("Found Peaks")
            plt.xlabel('m/z')
            plt.ylabel('Relative Intensity')
            if BigDisplay:
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()
                figManager.window.raise_()
            plt.show()
        
        t_plot = time.time()
        print "\nTotal time used: ", t_plot-t
        print "Time to read in: ", t_read - t
        print "Time to process: ", t_process - t_read
        print "Time to detect peaks: ", t_peak - t_process
        print "Time to find overlaps: ", t_overlap - t_peak
        print "Time to plot: ", t_plot - t_overlap

    else:
        #ENVELOPE SIMULATION SETUP**************************************************
        TruePeakNumber = PeakNumber
        
        Simulation = range(20)
        for i in range(20):
           Simulation[i] = Envelope()
           Simulation[i]._init_()


        ExtraSimulations = []
        
        SearchNumber = 0
        while 1:
        
            SimulationNumber = 0
            masses = [] #list of masses found, makes sure there are no repeats within the mass tolerance
            continue_counter = 0
            while True:
                if not Simulation[SimulationNumber].SimulationSetUp(Peak, SimulationNumber):
                        Simulation[SimulationNumber]._init_()  
                        break;
                Simulation[SimulationNumber].CentralPeakFit(Peak, False)
                good = Simulation[SimulationNumber].CompleteFit(Peak, False)
                if not good:
                    Peak[Simulation[SimulationNumber].central_peak_index].simulated = True # DO for the other ones too? and add to the iterative subtracting?
                    Simulation[SimulationNumber]._init_()                                  # Should make these checks uniform
                    continue_counter +=1
                elif len(Simulation[SimulationNumber].peaks) < min_peak_number:
                    print "ERROR: The number of peaks found by the simulation above is too small. It will not be considered."
                    Simulation[SimulationNumber]._init_()
                    continue_counter += 1
                else:
                    masses.append(Simulation[SimulationNumber].mass) 
                    SimulationNumber +=1
                if continue_counter >2:
                    break;
                if SimulationNumber >= MaxSimulations: 
                    break;
        
            print "BEFORE SUBTRACTING AND REFITTING TO ACCOUNT FOR OVERLAPS:\n"
            for i in range(SimulationNumber):
                print "Mass", i+1, "is", Simulation[i].mass, "and abundance is", Simulation[i].abundance
                print "at center", Simulation[i].center
                print "for peaks", Simulation[i].peak_mz
                print "at charges", Simulation[i].peak_charges
            print "\n"
            #Choose to fit charge envelopes using levenberg-marquardt OR iteratively subtract and refit
            
            if not UseSubtract:
                SimulationFit()
            else:
                Automatic = SubtractAuto
                for j in range(NumberRefit):
                    b = range(SimulationNumber)
                    b.reverse()  # fit the smallest, most insignificant ones first
                    for i in b:
                        Simulation[SimulationNumber].CentralPeakFit(Peak, True)
                        Simulation[i].CompleteFit(Peak, True)
            
            i = 0
            while 1:
                if Simulation[i].curve_scale < threshold:
                    Simulation.remove(Simulation[i])
                    SimulationNumber -= 1
                else:
                    i += 1
                if i >= SimulationNumber - 1:
                    break;
            
            SearchNumber += 1
            if SearchNumber >= RepeatSearch + 1:
                break;
            else:
                for i in range(SimulationNumber):
                    ExtraSimulations.append(Simulation[i])
                for i in range(SimulationNumber):
                    Simulation.remove(Simulation[0])
                SimulationNumber = 0
                
                i = 0
                while 1:
                    if Peak[i].simulated == True:
                        Peak.remove(Peak[i])
                        PeakNumber -= 1
                    else:
                        i +=1
                    if i >= PeakNumber:
                        break;
                
        for i in range(len(ExtraSimulations)):
            Simulation.insert(0, ExtraSimulations[len(ExtraSimulations) - 1 - i])
            SimulationNumber +=1
        
                
        t_charge = time.time()          
                
        if SimulationNumber >=1:
            #PLOTTING*******************************************************************
            outputs = np.zeros(len(intensityarray))
            for i in range(SimulationNumber):
                outputs += Simulation[i].Function(mzarray)
            DiffSpec = intensityarray - outputs
            
            #FIX THE MAX AND MIN. ONLY REPORT RELATIVE ABUNDANCES.
            plt.figure()
            ax = plt.subplot(311)
            plt.title('Original Spectrum')
            plt.axis( [min(mzarray), max(mzarray), min(intensityarray)/max(intensityarray), 1.2])
            plt.plot(mzarray, intensityarray/max(intensityarray), 'k')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
    
            ax = plt.subplot(312)
            plt.title('Simulated Spectrum')
            plt.axis( [min(mzarray), max(mzarray), min(intensityarray)/max(intensityarray), 1.2])
            plt.plot(mzarray, outputs/max(intensityarray), 'b')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)    
            
            ax = plt.subplot(313)
            plt.title('Subtracted Spectrum')
            plt.axis( [min(mzarray), max(mzarray), -0.2, 1.2])
            plt.plot(mzarray, (DiffSpec)/max(intensityarray), 'k')
            plt.xlabel('m/z')
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            
            if BigDisplay:
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()
                figManager.window.raise_()
            plt.figure()
            
            for i in range(SimulationNumber):
                ax = plt.subplot(SimulationNumber, 2, 2*(i+1) -1)
                plt.axis( [min(mzarray), max(mzarray), min(intensityarray)/max(intensityarray), 1.2])
                ax.yaxis.set_minor_locator(MultipleLocator(0.25))
                plt.locator_params(axis = 'y', tight = True, nbins = 3)
                plt.plot(mzarray, intensityarray/max(intensityarray), 'k')
                plt.plot(mzarray, Simulation[i].Function(mzarray)/max(intensityarray), color = ColorDict[i], label = LetterDict[i] + str(round(Simulation[i].mass)) )
                for j in range(len(Simulation[i].peak_charges)):
                    plt.text(Simulation[i].mass/Simulation[i].peak_charges[j] + 1.00794, norm(Simulation[i].peak_charges[j], Simulation[i].curve_scale, Simulation[i].curve_center, Simulation[i].curve_std), LetterDict[i]  + str(Simulation[i].peak_charges[j]), ha = 'center', va = 'bottom')
                plt.legend(loc = 'best', frameon = False)        
            plt.xlabel('m/z')
            
            ax = plt.subplot(2, 2, 2)
            x = np.linspace(min([Simulation[i].mass for i in range(SimulationNumber)])-20000, max([Simulation[i].mass for i in range(SimulationNumber)]) + 20000, 1000)
            mass_intensity = np.zeros(len(x))
            for i in range(SimulationNumber):
                mass_intensity += norm(x, Simulation[i].abundance, Simulation[i].mass, Simulation[i].central_fit_std * Simulation[i].charge)
            plt.plot(x, mass_intensity/max(mass_intensity))
            for i in range(SimulationNumber):
                plt.text(Simulation[i].mass, Simulation[i].abundance/max([Simulation[j].abundance for j in range(SimulationNumber)]), LetterDict[i]  + str(round(Simulation[i].mass)))
            plt.title("Deconvolved masses")
            plt.xlabel('Mass')
            plt.ylabel('Relative Abundance')
            plt.axis([min(x), max(x), min(intensityarray)/max(intensityarray), 1.1])
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'x', tight = True, nbins = 7)
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            
            
            ax = plt.subplot(2,2,4)
            plt.axis( [min(mzarray), max(mzarray), 0, 1.2])
            plt.plot(mzarray, intensityarray/max(intensityarray), 'k', label = "Raw")
            plt.plot(mzarray, outputs/max(intensityarray), 'b', label = "Fitted")
            plt.xlabel('m/z')
            plt.legend(loc = 'best', frameon = False)
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
            plt.locator_params(axis = 'y', tight = True, nbins = 3)
            
            if BigDisplay:
                figManager = plt.get_current_fig_manager()
                figManager.window.showMaximized()        
                figManager.window.raise_()
            plt.show()
            
            t_plot = time.time()
            
            print "FINAL RESULTS:\n"
            for i in range(SimulationNumber):
                print "Mass", i+1, "is", round(Simulation[i].mass), "and abundance is", round(Simulation[i].abundance)
                print "at m/z center", round(Simulation[i].center, 1)
                print "for peaks", [round(Simulation[i].peak_mz[j], 1) for j in range(len(Simulation[i].peak_mz))] 
                print "at charges", Simulation[i].peak_charges
                print "Full width at half maximum of mass peak is", round(Simulation[i].central_fit_std * Simulation[i].charge * 2 * 1.1774, 2)
                print "Charge center of simulation: ", round(Simulation[i].curve_center, 2)
                print "Height at this charge center: ", round(Simulation[i].curve_scale, 1)
                print "Full width at half maximum over charge domain is", round(Simulation[i].curve_std * 2*1.1774, 2)
                print "\n"
            print "Total peaks detected:\t", TruePeakNumber        
            
            AllPeaks = []
            for i in range(SimulationNumber):
                AllPeaks += Simulation[i].peak_mz #the index of AllPeaks should give the corresponding simulation number
            
            Overlaps = [o for o in set(AllPeaks) if AllPeaks.count(o) > 1]
            OverlapPeakMassDict = {Overlaps[i]: [] for i in range(len(Overlaps))}
            
            for i in range(SimulationNumber):
                for j in range(len(Overlaps)):
                    if Overlaps[j] in Simulation[i].peak_mz:
                        OverlapPeakMassDict[Overlaps[j]].append(i)
            
            for i in range(len(Overlaps)):
                print "Peak at", round(Overlaps[i],1) , "is shared by masses ", [round(Simulation[k].mass) for k in OverlapPeakMassDict[Overlaps[i]]]
            
            print "\n"
            
            if SaveSubtract:
                if NegToZeros:
                    for i in range(len(DiffSpec)):
                        if DiffSpec[i] <0:
                            DiffSpec[i] = 0.0
                MakeFile(name, SubtractName, beginning, mzarray, DiffSpec)
            
    
            if SaveMasses:
                mass_path = name.partition('.txt')[0] + '_' + MassName + '.txt'
                file = open(mass_path, 'w')
                for i in range(SimulationNumber):
                    file.write("Mass " + str(i+1) + " is " + str(round(Simulation[i].mass)) + " and abundance is " + str(round(Simulation[i].abundance)))
                    file.write(" and sum of heights is " + str(round(sum(Simulation[i].peak_heights))))
                    file.write('\n')
                    file.write("at m/z center:\t" + str(round(Simulation[i].center, 1)))
                    file.write('\n')
                    file.write("for peaks:\t" + str([round(Simulation[i].peak_mz[j], 1) for j in range(len(Simulation[i].peak_mz))]))
                    file.write('\n')
                    file.write("at charges:\t" + str(Simulation[i].peak_charges))
                    file.write('\n')
                    file.write("Full width at half maximum of mass peak:\t" + str(round(Simulation[i].central_fit_std * Simulation[i].charge * 2 * 1.1774, 2)))
                    file.write('\n')
                    file.write("Charge center of simulation:\t" + str(round(Simulation[i].curve_center, 2)))
                    file.write('\n')
                    file.write("Height of fitted Gaussian envelope:\t" + str(round(Simulation[i].curve_scale, 1)))
                    file.write('\n')
                    file.write("Full width at half maximum of fitted envelope over charge domain:\t" + str(round(Simulation[i].curve_std * 2 * 1.1774, 2)))
                    file.write('\n')
                    file.write('\n')
                file.write('Total Peaks Detected:\t' + str(TruePeakNumber))
                file.write('\n')
                for i in range(len(Overlaps)):
                    file.write("Peak at " + str(round(Overlaps[i], 1)) + " is shared by masses " + str([round(Simulation[k].mass) for k in OverlapPeakMassDict[Overlaps[i]]]))
                    file.write('\n')
                file.close()
                print 'Mass information saved to directory ', mass_path  
            
            t_write = time.time()
            
            print "Total time used: ", t_write - t
            print "Time to read in: ", t_read - t
            print "Time to process: ", t_process - t_read
            print "Time to detect peaks: ", t_peak - t_process
            print "Time to find overlaps: ", t_overlap - t_peak
            print "Time to assign charge states: ", t_charge - t_overlap
            print "Time to plot: ", t_plot - t_charge
            print "Time to write to file: ", t_write - t_plot
    
    
    
if __name__ == '__main__':
    main()
