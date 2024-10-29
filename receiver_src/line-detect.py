#!/home/kuppel/anaconda2/bin/python
import sys
import numpy as np
import struct
import matplotlib.pyplot as plt
from astropy import units as u
import os
import time
from astropy.io import ascii
from astropy.table import Table

filename = sys.argv[1]#'gqrx_20170602_134841_199987000_960000_fc.raw'
infile = open(filename,'rb')
f0 =float(filename.split("_")[3])# 0937000. #* u.Hz
sample_rate = float(filename.split("_")[4])
#f0 = 199987000
filesize = os.path.getsize(filename)*8 # bits
fc = 1.0 # multiplicative factor
chunksize =int( sample_rate * fc )
nchunks = int(np.floor(filesize/chunksize/2/32))
freq_chunk =f0+np.fft.fftfreq(chunksize,d = 1./sample_rate)
d_freq = freq_chunk[1] - freq_chunk[0]
spec_power_t = []
spec_power_t_db = []


#TOTAL POWER CALCULATION!

timestr = time.strftime('%Y-%m-%d_%H-%M-%S')
for k in range(nchunks):
	try:
		sig_chunk = np.array([struct.unpack('f',infile.read(4))[0]+struct.unpack('f',infile.read(4))[0]*1j for i in range(chunksize)])
	except:
		break
	P_spec = abs(np.fft.fft(sig_chunk))**2/len(sig_chunk)**2
	P_spec_db = 10*np.log10(P_spec)
#	P_sig = np.sum(abs(sig_chunk)**2)/len(sig_chunk)
#	P_sig_db = 10*np.log10(P_sig)
#	print k,P_sig_db
	
	spec_power_t.append(P_spec)
	spec_power_t_db.append(P_spec_db)
stack_spec = np.mean(spec_power_t,axis=0)
stack_spec_db = 10*np.log10(np.mean(spec_power_t,axis=0))

data_tab = Table([freq_chunk,stack_spec],names=('Frequency [Hz]','Intensity [W/Hz]'),dtype=(np.float,np.float))
ascii.write(data_tab,output="Specral_line"+timestr+".txt",format='fixed_width',delimiter=' ')

plt.plot(freq_chunk,stack_spec)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Intensity [W/Hz]")
plt.title("Spectral Line  "+timestr)
plt.savefig("spectral_line_"+timestr+".pdf")	
plt.show()
#readblock = sample_rate.value/2
#lengthins =np.floor(filesize/8/sample_rate)
#meanvals = []
#mean_std = []
#slices = int(lengthins.value)
#
#for j in range(slices):
#	ii=[]
#	qq=[]
#	for i in np.arange(j*readblock,(j+1)*readblock):
#		ii.append(struct.unpack('f',infile.read(4))[0])
#		qq.append(struct.unpack('f',infile.read(4))[0])
#	t = (j*len(ii)+np.arange(len(ii)))*1/sample_rate
#	signal = ii*np.cos(2*np.pi*(f0*t).decompose()*u.radian)-qq*np.sin(2*np.pi*(f0*t).decompose()*u.radian)
#	fourier = np.fft.rfft(signal)
#	winlen = len(signal)
#	freq = np.fft.rfftfreq(winlen,d = 1./sample_rate.value)*u.Hz
##	plt.subplot(slices/4,4,j+1)
#	dbspectrum = 10*np.log10(abs(fourier)*(freq[1].value-freq[0].value))
##	plt.plot(freq+f0,dbspectrum)
#	tempmean = np.mean(dbspectrum)
#	tempstd =  np.std(dbspectrum)
#	meanvals.append(tempmean)
#	mean_std.append(tempstd)
#	print 100*float(j+1)/slices,"% done.",tempmean,tempstd
#plt.errorbar(np.arange(slices)+1,meanvals,yerr=mean_std)
#plt.xlabel("seconds")
#plt.ylabel("dB")
#
#plt.savefig("hotcold_test1.eps")

#fftamp = [np.fft.rfft(a) for a in amp]
#sample_rate = wavefile[0]
#winlen = len(ii[0])
#freq = np.fft.rfftfreq(winlen, d=1./sample_rate)
#
#for i,f import n zip(range(1,len(fftamp)),fftamp):
#	plt.subplot(5,6,i)
#	plt.plot(freq[1:],abs(f[1:]))
#
#plt.show()
