#!/home/kuppel/anaconda2/bin/python

##!/usr/bin/env python2

##Updated: Feb 27, 2024
##Python2 script
###/home/kuppel/anaconda2/bin/python
import sys
import numpy as np
#import struct
import matplotlib.pyplot as plt
#from astropy import units as u
import os
import time
import argparse
#from astropy.io import ascii
#from astropy.table import Table

# Set up argparse
parser = argparse.ArgumentParser(
    description="Calculates Power values [in dBm] for S262 Day1 Experiment.",
    usage="%(prog)s <filename> [--debug]"
)
parser.add_argument("filename", type=str, help="The gqrx filename to process")
parser.add_argument('--debug', action='store_true', help='Print Power values to terminal')

# Parse the command line arguments
args = parser.parse_args()

# read the gqrx file and data about it from its name
#filename = sys.argv[1]#'gqrx_20170602_134841_199987000_960000_fc.raw'
#infile = open(filename,'rb')
#f0 =float(filename.split("_")[3])#  u.Hz
#sample_rate = float(filename.split("_")[4]) # I believe also Hz

try:
    filename = args.filename  # Expecting the filename as the first command-line argument
    with open(filename, 'rb') as infile:
        # Extracting parts of the filename assuming a specific format
        parts = filename.split("_")
        f0 = float(parts[3])  # Assuming this part represents a frequency in Hz
        sample_rate = float(parts[4])  # Assuming this part represents the sample rate in Hz
except:
    print "Usage: script.py <filename> \n The file is either unavailable or the gqrx filename was renamed/altered"

filesize = os.path.getsize(filename)*8 # bits
fc = 1.0 # multiplicative factor

#we split in second sized chunks...
chunksize =int( sample_rate * fc )
nchunks = int(np.floor(filesize/chunksize/2/32)) 

#we use this timestamp for output file
timestr = time.strftime('%Y-%m-%d_%H-%M-%S')

print "reading file..."

#file is a byte array, buffer is 4 bytes each (32 bit floats)
#we get I and Q values in float, so each time it's 8 bytes 

with open(filename,'rb') as infile:
	byte_array = bytearray(infile.read()) #read file
	data_array = np.ndarray(shape = (len(byte_array)//8,), dtype = "f, f", buffer = byte_array) #put into array n,2 type float

print "processing..."
data_complex = (data_array['f0'] + data_array['f1']*1j)[:-(len(data_array) % chunksize)].reshape(nchunks,-1) #make sure one can get the split in chunks by discarding end, make complex 
P_array_db = 10*np.log10(np.sum(np.abs(data_complex)**2/chunksize,axis=1)) # each chunk is a second and becomes a total power reading in db
time = range(nchunks) # s we generate the time array

if args.debug:
    print "Time [s]\tPower [dBm]" 

    # Iterate and print each pair of time and power values
    for t, p_dB in zip(time, P_array_db):
        #print "{t}\t{p_dB}" #Python3
        print("{}\t{}".format(t, p_dB)) #Python2

# Open a file named 'output.txt' in write mode
write_filename = "Hot_cold_"+timestr

with open(write_filename+".txt", 'w') as file:
    # Write the header to the file
    file.write("Time [s]\tPower [dBm]\n")

    # Iterate and write each pair of time and power values to the file
    for t, p_dB in zip(time, P_array_db):
        #file.write(f"{t}\t{p_dB}\n") #Python3
        file.write("{}\t{}\n".format(t, p_dB)) #Python2

print "done!"

#plot it
plt.plot(time,P_array_db)
plt.xlabel("Time [s]")
plt.ylabel("Total Power in Band [dB]")
plt.title("Hot Cold Measurement Plot "+timestr)
plt.savefig("Hot_cold_"+timestr+".eps")	
plt.show()
plt.close()

#Toma Badescu, Ankur Dev 2024... 
