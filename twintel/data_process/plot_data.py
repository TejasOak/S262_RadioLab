#!/home/kuppel2/anaconda3/envs/py2env/bin/python
from datetime import datetime
from time import time
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pickle
import sys
datain = []


#read the datafile
try:
	datatab = ascii.read(sys.argv[1])
except Exception,e:
	print e

#assign the columns to arrays
#coords = coordinates of the pointings
#fluxa and fluxb are fluxes from ant 1 and 2
#Q and I are components of the visibility

coords = [SkyCoord(d[0]+" "+d[1],unit=(u.hourangle,u.deg))  for d in zip(datatab['col3'],datatab['col4'])]
fluxa = list(datatab['col5'])
fluxb = list(datatab['col8'])
I = list(datatab['col6'])
Q = list(datatab['col7'])
nr = list(datatab['col8'])
coordpair = zip(coords,nr)
loadpickle = 0

#create lists which hold the recorded values associated with an
#index number
pfluxa =[list(f) for f in zip(fluxa,range(len(fluxa)))]
pfluxb =[ list(f) for f in zip(fluxb,range(len(fluxb)))]
pQ =[list(f) for f in  zip(Q,range(len(Q)))]
pI =[list(f) for f in  zip(I,range(len(I)))]


#despike by subtracting from the original data the same data
#shifted by 1 index to the left or to the right
#this eliminates the slow variations due to the sun or other source
#being scanned and (hopefully) leaves only the spikes in the data
#which are then removed from the arrays
#iterate 6 times doing this, since the largest spikes are 5-6 positions wide
i=8
while i:
	ind = np.where(abs(fluxa - np.roll(fluxa,1**i)) > 60)[0]
	pfluxa = np.delete(pfluxa,ind,axis=0)
	fluxa = np.delete(fluxa,ind)
	print len(ind)
	ind = np.where(abs(fluxb - np.roll(fluxb,1**i)) > 60)[0]
	pfluxb = np.delete(pfluxb,ind,axis=0)
	fluxb = np.delete(fluxb,ind)

	ind = np.where(abs(Q - np.roll(Q,1**i)) > 30)[0]
	pQ = np.delete(pQ,ind,axis=0)
	Q = np.delete(Q,ind)
	
	ind = np.where(abs(I - np.roll(I,1**i)) > 30)[0]
	pI = np.delete(pI,ind,axis=0)
	I = np.delete(I,ind)


	i-=1



#plot the raw flux data so we can decide from where to where the data spans
plt.plot(fluxa)
plt.title("Decide on start and end of the data")
plt.xlabel("Time")
plt.ylabel("Amp")
plt.show()

#choose the data range
try:
	print "=== Choose data range after looking at the plot ==="
	data_start = int(raw_input("Starting position [int]:"))
	data_end =  int(raw_input("Ending position [int]:"))
	fluxa = fluxa[data_start:data_end] 
	fluxb = fluxb[data_start:data_end]
	Q     = Q[data_start:data_end]
	I     = I[data_start:data_end]
	pfluxa = pfluxa[data_start:data_end]
	pfluxb = pfluxb[data_start:data_end]
	pQ     = pQ[data_start:data_end]
	pI     = pI[data_start:data_end]
	coords = coords[data_start:data_end]
except Exception,e:
	print e


#set to zero mean
meanfluxa = np.mean(fluxa)
meanfluxb = np.mean(fluxb)
meanQ = np.mean(Q)
meanI = np.mean(I)
#using the value index pairs interpolate over the missing data 

ifluxa = np.interp(range(data_start,data_end),[p[1] for p in pfluxa],[p[0]-meanfluxa for p in pfluxa])
ifluxb = np.interp(range(data_start,data_end),[p[1] for p in pfluxb],[p[0]-meanfluxb for p in pfluxb])
iQ = np.interp(range(data_start,data_end),[p[1] for p in pQ],[p[0]-meanQ for p in pQ])
iI = np.interp(range(data_start,data_end),[p[1] for p in pI],[p[0]-meanI for p in pI])

#since we take multiple measurements at the same coordinate position 
# because we sample so fast that the coordinates don't change fast enough
# we create a dictionary where the keys are the coordinates and the values
# are the recorded data at each coordinate. We take the median value from these
# data and use those.

keys = np.unique(coords)
datdic = {}
for c,fa,fb,ii,qq in zip(coords,ifluxa,ifluxb,iI,iQ):
	try:
		datdic[(c.ra.value,c.dec.value)].append([fa,fb,ii,qq])
		#print "a"
	except:
		datdic[(c.ra.value,c.dec.value)]=[[fa,fb,ii,qq]]

mediandat = {}
dcoords = []
dvalues = []
for d in datdic:
	medi = np.median(datdic[d],axis=0)
	mediandat[d] = medi
	dcoords.append(d)
	dvalues.append(medi)

#load file
if loadpickle:
	with open("cleansig.txt",'r') as fin:
		print "Loading pickled file"
		dvalues = pickle.load(fin)
		dcoords = [coords[c[1]] for c in dvalues[0]]


#do some interf


#interpolate

#gridding: we create a 2d grid for our data
# here we define the extents of the grid in ra dec

ramax = max([c.ra.value for c  in coords])
ramin = min([c.ra.value for c  in coords])
decmax = max([c.dec.value for c  in coords])
decmin = min([c.dec.value for c  in coords])
#create the grid
grid_ra,grid_dec = np.mgrid[ramin:ramax:50j,decmin:decmax:50j]
#use a time stamp for newly created files so we don't overwrite:
tt = str(time())
#then we interpolate in 2d
grid1 = griddata(dcoords,[d[0] for d in dvalues],(grid_ra,grid_dec),method = 'cubic')
plt.imshow(grid1.T[2:-1,2:-1],extent = [ramin,ramax,decmin,decmax])
plt.title("Antenna A Sun Scan")
plt.savefig("AntennaAscan"+tt+".png")
plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.show()

grid1 = griddata(dcoords,[d[1] for d in dvalues],(grid_ra,grid_dec),method = 'cubic')
plt.imshow(grid1.T[2:-1,2:-1],extent = [ramin,ramax,decmin,decmax])
plt.title("Antenna B Sun Scan")
plt.savefig("AntennaBscan"+tt+".png")
plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.show()

grid_ra,grid_dec = np.mgrid[ramin:ramax:150j,decmin:decmax:100j]
grid1 = griddata(dcoords,[abs(d[2]+d[3]*1j)**2 for d in dvalues],(grid_ra,grid_dec),method = 'cubic')
plt.imshow(grid1.T[2:-1,2:-1],extent = [ramin,ramax,decmin,decmax])
plt.tricontour([d[0] for d in dcoords],[d[1] for d in dcoords],[abs(d[2]+d[3]*1j)**2 for d in dvalues])
plt.title("Interferometer Scan")
plt.savefig("Interfscan"+tt+".png")
plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.show()




#plt.plot(fluxa)
#plt.show()

#plt.tricontour([d.ra.value for d in coords],[d.dec.value for d in coords],fluxa,20)
#plt.savefig("sunplot"+str(time())+".pdf")
#plt.scatter([d.ra.value for d in datatab["RA Dec"]],[d.dec.value for d in datatab["RA Dec"]],s=[d for d in datatab["TelAflux"]]/np.max([d for d in datatab["TelAflux"]]))
