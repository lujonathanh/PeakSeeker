import numpy as np
import matplotlib.pyplot as plt
import random
import math

class Envelope:
    def _init_(self):
        self.number=0
        self.mass_std = 0 #modified, used to be self.central_fit_std
        self.mass = 0
        self.charge = 1 #modified, used to be 0
        self.peaks = [] #peak number in list Peaks
        self.peak_charges = []
        self.peak_mz = []
        self.peak_heights = []
        self.peak_std = []
        self.curve_scale= 0
        self.curve_center = 50 #on the z axis, not on the m/z axis
        self.curve_std = 0 #random beginning guess (also on z axis)
        self.abundance = 0
    
    def Function(self, x): #outputs the y-value of the peak spectrum given the mz-value
        output = 0

        for z in self.peak_charges:
            peaksim_height = norm(z, self.curve_scale, self.curve_center, self.curve_std) #gets peak height on the gaussian fitted to envelope distribution
            output += norm(x, peaksim_height, (self.mass/z) + 1.00794, self.mass_std / z) #height on this peak

        return(output)
    
    def DisplayParameters(self):
        print "Mass is: ", self.mass
        print "Charge states: ", self.peak_charges
        print "Full width at half maximum of mass peak: ", self.mass_std * 2.3548
        print "Charge center of simulation: ", self.curve_center
        print "Height at this charge center: ", self.curve_scale
        print "Full width at half maximum of the charge envelope: ", self.curve_std * 2.3548

def SimulationSetup(number):
    global Simulation
    print "Enter, in one line separated by semicolons, the following information:\n"
    print "Mass; charge states of the mass as a list, e.g. [3,4,5,6]; full width at half maximum of mass peak;\n"
    print "charge center of simulation (can be decimal when central peak isn't at the center of the Gaussian);\n"
    print "height at this charge center; full width at half maximum of the charge envelope."
    info = raw_input()
    Simulation[number].mass = eval(info.partition(';')[0])
    #print "Enter the charge states of the mass as a list; e.g. [3,4,5,6]."
    Simulation[number].peak_charges = eval(info.partition(';')[2].partition(';')[0])
    #print "Enter standard deviation of mass peak."
    Simulation[number].mass_std = eval(info.partition(';')[2].partition(';')[2].partition(';')[0]) / 2.3548
    #print "Enter the charge center of simulation. This can be a decimal in cases where the central peak isn't at the center of the Gaussian."
    Simulation[number].curve_center = eval(info.partition(';')[2].partition(';')[2].partition(';')[2].partition(';')[0])
    #print "Enter the height at this charge center."
    Simulation[number].curve_scale = eval(info.partition(';')[2].partition(';')[2].partition(';')[2].partition(';')[2].partition(';')[0])
    #print "Enter the standard deviation of the charge envelope."
    Simulation[number].curve_std = eval(info.partition(';')[2].partition(';')[2].partition(';')[2].partition(';')[2].partition(';')[2].partition(';')[0]) / 2.3548
    Simulation[number].abundance = sum([math.sqrt(2*math.pi) * norm(Simulation[number].peak_charges[i], Simulation[number].curve_scale, Simulation[number].curve_center, Simulation[number].curve_std) * Simulation[number].mass_std / Simulation[number].peak_charges[i]  for i in range(len(Simulation[number].peak_charges))])


def norm(x, height, center, std):
  return(height*np.exp(-(x - center)**2/(2*std**2)))
  
def Simulate(x, a, b, c, d, e):#returns y-value of an mz for the combined simulated spectrum
    # global SimulationNumber (is simulation number needed? since all start with curve_scale at 0, it should be fine)
    global Simulation
    y = 0
    y += a * Simulation[0].Function(x)
    y += b * Simulation[1].Function(x)
    y += c * Simulation[2].Function(x)
    y += d * Simulation[3].Function(x)
    y += e * Simulation[4].Function(x)
    return(y)

def MakeFile(file_name, beginning, mzarray, DiffSpec):

	"""
		MakeFile(file_name): makes a file.
	"""

	temp_path = "C:\Python27\Simulated Data\\" + file_name.partition('.txt')[0] + '.txt'
	file = open(temp_path, 'w')
	file.write(beginning)
	for i in range(len(mzarray)):
	    file.write(str(mzarray[i]))
	    file.write('\t')
	    file.write(str(DiffSpec[i]))
	    file.write('\n')
	file.close()
	print 'Simulated spectrum saved to directory', temp_path   

def MakeEnvelopes(adjust = True):
    global Simulation
    global SimulationNumber
    global MaxSimulations
    
    if adjust:
        print "Edit/Remove Past Simulations? Enter 'y'."
        entry = raw_input()
        if entry == 'y':
            i = 0
            while 1:
                if i >= SimulationNumber:
                    break;
                print "Parameters of Simulation", i
                Simulation[i].DisplayParameters()
                print "Edit? Enter 'y'. Remove? Enter 'r'. Move to the next one? Enter 'n'. Stop editting altogether? Enter anything else."
                entry = raw_input()
                if entry == 'y':
                    SimulationSetup(i)
                    i +=1
                elif entry == 'r':
                    Simulation.remove(Simulation[i])
                    SimulationNumber -=1
                elif entry == 'n':
                    i +=1
                else:
                    break;
             
    while 1:
        print "Add Simulation", SimulationNumber + 1, "? Enter 'y'."
        entry = raw_input()
        if entry != 'y':
            break;
        else:
            SimulationSetup(SimulationNumber)
            SimulationNumber += 1
            if SimulationNumber >= MaxSimulations:
                break;
    
def main():
    global SimulationNumber
    SimulationNumber = 0
    global MaxSimulations
    MaxSimulations = 5
    global Simulation
    global minmz
    global maxmz
    global numdata

    
    global x
    global intensity
    mode = 'm' #d is deconvoluting spectrum. m is making spectrum.
    ChargeDisplayTime = 10
    
    if mode == 'm':
        minmz = 0
        maxmz = 12000
        numdata = 10000

        noise_std = 500
        x = np.linspace(minmz, maxmz, numdata)

    if mode == 'd':
        directory = 'C:\Python27\\Raw Data\\tric2_4.txt'
        opener = open(directory)
        mzarray = []
        intensityarray = []
        index = 0
        for line in opener:
            if index>= 8:
                a = line.partition('\t')
                mzarray.append(float(a[0]))
                b= a[2].partition('\n')
                intensityarray.append(float(b[0]))
            index+=1;
        opener.close()
        mzarray = np.array(mzarray)
        intensityarray=np.array(intensityarray)
        x = mzarray
    

    
    
    Simulation = range(MaxSimulations)
    for i in range(MaxSimulations):
   	Simulation[i] = Envelope()
   	Simulation[i]._init_()
   
    MakeEnvelopes(False)
    

    while 1:   
        intensity = np.zeros(len(x))
        for i in range(SimulationNumber):
            intensity += Simulation[i].Function(x)
        
        if mode == 'm':
            intensity += np.array([abs(random.gauss(0, noise_std)) for i in range(len(x))])
        
        while 1:
            plt.figure()
            if mode == 'd':
                plt.subplot(211)
                plt.title("Overlaid Spectrum")
                plt.plot(mzarray, intensityarray, 'k')
                plt.plot(x, intensity, 'c')
                for i in range(SimulationNumber):
                    outputs = Simulation[i].Function(x)
                    plt.plot(x, outputs, label = round(Simulation[i].mass))
                plt.xlabel('m/z')
                plt.ylabel('Intensity')
                plt.legend(loc = 'best')
                plt.subplot(212)
                plt.title("Subtracted Spectrum")
                plt.plot(mzarray, intensityarray-intensity)
                
                plt.show()
                
                plt.pause(ChargeDisplayTime)
                plt.close()
            
            if mode == 'm':
                plt.plot(x, intensity)
                
                
                plt.title("Simulated Spectrum")
                plt.xlabel('m/z')
                plt.ylabel('Intensity')
                plt.show()
                plt.pause(ChargeDisplayTime)
                plt.close()
            if (raw_input("Enter 'y' to display again.") != 'y'):
                break;
        
        if (raw_input("Edit/Remove/Add envelopes? Enter 'y'. ") == 'y'):
            MakeEnvelopes(True)
        else:
            break;
            
            
    plt.figure()
    for i in range(SimulationNumber):
        plt.subplot(SimulationNumber, 1, i)
        plt.axis( [min(x), max(x), min(intensity), 1.2*max(intensity)])
        plt.plot(x, Simulation[i].Function(x))
        plt.xlabel('m/z')
        plt.ylabel('intensity')
        #for j in range(len(Simulation[i].peak_charges)):
        #    plt.text(Simulation[i].peak_mz[j], Simulation[i].peak_heights[j], Simulation[i].peak_charges[j])
        plt.title('%.0f' % Simulation[i].mass)
    
    if mode == 'm':
        if (raw_input("Enter 'y' to save spectrum.") == 'y'):
            file_name = raw_input("Enter a file name.")
            beginning = 'SIMULATED SPECTRUM - MS\n' + file_name + '\nMasses:' + str([Simulation[i].mass for i in range(SimulationNumber)]) +  '\n Abundances:' + str([Simulation[i].abundance for i in range(SimulationNumber)]) + '\nCharges of masses:'+ str([Simulation[i].peak_charges for i in range(SimulationNumber)]) + ' Full Width at Half Maximum of Mass Peak: ' + str([2*1.1774 * Simulation[i].mass_std for i in range(SimulationNumber)])+ ' Charge Center of Simulation: '+ str([Simulation[i].curve_center for i in range(SimulationNumber)]) + ' Height at charge center: ' + str([Simulation[i].curve_scale for i in range(SimulationNumber)]) + ' Full Width at Half Maximum of charge envelope: '  + str([2*1.1774* Simulation[i].curve_std for i in range(SimulationNumber)]) + '\nNoise level: ' + str(noise_std) + '\nData points: ' + str(numdata)+ '\nMass\tIntensity\n'
            MakeFile(file_name, beginning, x, intensity)
    for i in range(SimulationNumber):
        Simulation[i].peak_mz = [Simulation[i].mass/j + 1.00794 for j in Simulation[i].peak_charges]
        print "Total abundance of mass " + str(Simulation[i].mass) + " is " + str(Simulation[i].abundance) + "\n"
    
if __name__ == '__main__':
    main()
    
