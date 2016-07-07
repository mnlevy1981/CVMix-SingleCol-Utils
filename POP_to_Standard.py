#!/usr/bin/env python

# This script converts POP netcdf output to uniform format needed to compare
# CVMix tests to other GCMs. Two parts:
# 1) Change variable names (and add time dimension to some POP invariant vars)
# 2) Convert units from cgs to mks

import argparse
import os.path
import sys
from netCDF4 import Dataset
import numpy

def min_max(x):
  print numpy.amin(x), numpy.amax(x)

# Require POP history filename
parser = argparse.ArgumentParser(description='Covert POP history file to CVMix standard format')
parser.add_argument('filename_in',  metavar='INFILE',  help='POP history file')
parser.add_argument('filename_out', metavar='OUTFILE', help='Output netCDF file')
args = parser.parse_args()
POPhist = args.filename_in
OutFile = args.filename_out

try:
  data = Dataset(POPhist)
except:
  print "ERROR:", POPhist, "could not be found or is not a netCDF file"
  sys.exit(1)

print "Reading data from", POPhist
loc_time     = (data.variables['time'][:] - 366.0) * 86400.0 # days -> sec
loc_zt       = data.variables['z_t'][:]*0.01 # cm -> m
loc_zw       = numpy.concatenate(([0],data.variables['z_w_bot'][:]))*0.01 # cm -> m
loc_blot     = data.variables['HBLT'][:,0,0]*0.01 # cm -> m
loc_thetao   = data.variables['TEMP'][:,:,0,0]
loc_so       = data.variables['SALT'][:,:,0,0]
loc_uo       = data.variables['UVEL'][:,:,0,0]*0.01 # cm/s -> m/s
loc_vo       = data.variables['VVEL'][:,:,0,0]*0.01 # cm/s -> m/s
loc_difvho   = data.variables['VDC_T'][:,:,0,0]*0.0001 # cm^2/s -> m^2/s
#loc_difvso   = data.variables['VDC_S'][:,:,0,0]*0.0001 # cm^2/s -> m^2/s
loc_difvmo   = data.variables['VVC'][:,:,0,0]*0.0001 # cm^2/s -> m^2/s
loc_difvhonl = data.variables['KPP_NONLOCAL_T'][:,:,0,0]
loc_difvsonl = data.variables['KPP_NONLOCAL_S'][:,:,0,0]
#loc_difvmonl = data.variables['KPP_NONLOCAL_UV'][:,:,0,0]

if False:
  min_max(loc_time)
  min_max(loc_zt)
  min_max(loc_zw)
  min_max(loc_blot)
  min_max(loc_thetao)
  min_max(loc_so)
  min_max(loc_uo)
  min_max(loc_vo)
  min_max(loc_difvho)
  #min_max(loc_difvso)
  min_max(loc_difvmo)
  min_max(loc_difvhonl)
  min_max(loc_difvsonl)
  #min_max(loc_difvmonl)

# (2) Write output file
if os.path.isfile(OutFile):
  os.remove(OutFile)

print "Writing data to", OutFile
ncfile = Dataset(OutFile,'w')
# (2a) Define dimensions
ncfile.createDimension('time',loc_time.size)
ncfile.createDimension('zt',loc_zt.size)
ncfile.createDimension('zw',loc_zw.size)

# (2b) Define / set variables
time           = ncfile.createVariable('time','f8',('time',))
time.long_name = 'Elapsed Time'
time.units     = 'seconds'
time[:]        = loc_time

zt           = ncfile.createVariable('zt','f8',('time','zt'))
zt.long_name = 'Cell Center Depth'
zt.units     = 'm'
for i in range(0, loc_time.size):
  zt[i,:]  = loc_zt

zw           = ncfile.createVariable('zw','f8',('time','zw'))
zw.long_name = 'Cell Interface Depth'
zw.units     = 'm'
for i in range(0, loc_time.size):
  zw[i,:]  = loc_zw

blot           = ncfile.createVariable('blot','f8',('time',))
blot.long_name = 'Boundary Layer Ocean Thickness'
blot.units     = 'm'
blot[:]  = loc_blot

# (2c) Close file to actually write
ncfile.close()

