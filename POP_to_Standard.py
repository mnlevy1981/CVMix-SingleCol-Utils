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

def min_max(x,name):
  print "Range of", name, "is", numpy.amin(x), "to", numpy.amax(x)

# Require POP history filename
parser = argparse.ArgumentParser(description='Covert POP history file to CVMix standard format')
parser.add_argument('filename_in',  metavar='INFILE',  help='POP history file')
parser.add_argument('filename_out', metavar='OUTFILE', help='Output netCDF file')
parser.add_argument('--verbose', dest='verbose', action='store_const', const=True, default=False, help='Verbose stdout')
args = parser.parse_args()
POPhist = args.filename_in
OutFile = args.filename_out
verbose = args.verbose

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
loc_difvmo   = data.variables['VVC'][:,:,0,0]*0.0001 # cm^2/s -> m^2/s
loc_difvhonl = data.variables['KPP_NONLOCAL_T'][:,:,0,0]
loc_difvsonl = data.variables['KPP_NONLOCAL_S'][:,:,0,0]
# Two variables that may not be in POP output
has_difvso = True
try:
  loc_difvso   = data.variables['VDC_S'][:,:,0,0]*0.0001 # cm^2/s -> m^2/s
except:
  print "Warning: difvso is not available in POP history file"
  has_difvso = False
has_difvmonl = True
try:
  loc_difvmonl = data.variables['KPP_NONLOCAL_UV'][:,:,0,0]
except:
  print "Warning: difvmonl is not available in POP history file"
  has_difvmonl = False

if verbose:
  min_max(loc_time, 'time')
  min_max(loc_zt, 'zt')
  min_max(loc_zw, 'zw')
  min_max(loc_blot, 'blot')
  min_max(loc_thetao, 'thetao')
  min_max(loc_so, 'so')
  min_max(loc_uo, 'uo')
  min_max(loc_vo, 'vo')
  min_max(loc_difvho, 'difvho')
  if has_difvso:
    min_max(loc_difvso, 'difvso')
  min_max(loc_difvmo, 'difvmo')
  min_max(loc_difvhonl, 'difvhonl')
  min_max(loc_difvsonl, 'difvsonl')
  if has_difvmonl:
    min_max(loc_difvmonl, 'difvsonl')

# (2) Write output file
if os.path.isfile(OutFile):
  os.remove(OutFile)

print "Writing data to", OutFile
ncfile = Dataset(OutFile,'w')
# (2a) Define dimensions
ncfile.createDimension('time',loc_time.size)
ncfile.createDimension('nlev',loc_zt.size)
ncfile.createDimension('nface',loc_zw.size)

# (2b) Define / set variables
time           = ncfile.createVariable('time','f8',('time',))
time.long_name = 'Elapsed Time'
time.units     = 'seconds'
time[:]        = loc_time

zt           = ncfile.createVariable('zt','f8',('time','nlev'))
zt.long_name = 'Cell Center Depth'
zt.units     = 'm'
for i in range(0, loc_time.size):
  zt[i,:]  = loc_zt

zw           = ncfile.createVariable('zw','f8',('time','nface'))
zw.long_name = 'Cell Interface Depth'
zw.units     = 'm'
for i in range(0, loc_time.size):
  zw[i,:]  = loc_zw

blot           = ncfile.createVariable('blot','f8',('time',))
blot.long_name = 'Boundary Layer Ocean Thickness'
blot.units     = 'm'
blot[:]        = loc_blot

thetao           = ncfile.createVariable('thetao','f8',('time','nlev'))
thetao.long_name = 'Sea Water Potential Temperature'
thetao.units     = 'deg C'
thetao[:,:]      = loc_thetao

so           = ncfile.createVariable('so','f8',('time','nlev'))
so.long_name = 'Sea Water Salinity'
so.units     = 'Absolute Salinity (g/kg)'
so[:,:]      = loc_so

uo           = ncfile.createVariable('uo','f8',('time','nlev'))
uo.long_name = 'Sea Water x-Velocity'
uo.units     = 'm/s'
uo[:,:]      = loc_uo

vo           = ncfile.createVariable('vo','f8',('time','nlev'))
vo.long_name = 'Sea Water y-Velocity'
vo.units     = 'm/s'
vo[:,:]      = loc_vo

difvho           = ncfile.createVariable('difvho','f8',('time','nface'))
difvho.long_name = 'Ocean Vertical Heat Diffusivity'
difvho.units     = 'm^2/s'
difvho[:,0]      = 0.0
difvho[:,1:]     = loc_difvho

if has_difvso:
  difvso           = ncfile.createVariable('difvso','f8',('time','nface'))
  difvso.long_name = 'Ocean Vertical Salt Diffusivity'
  difvso.units     = 'm^2/s'
  difvso[:,0]      = 0.0
  difvso[:,1:]     = loc_difvso

difvmo           = ncfile.createVariable('difvmo','f8',('time','nface'))
difvmo.long_name = 'Ocean Vertical Momentum Diffusivity'
difvmo.units     = 'm^2/s'
difvmo[:,0]      = 0.0
difvmo[:,1:]     = loc_difvmo

difvhonl           = ncfile.createVariable('difvhonl','f8',('time','nface'))
difvhonl.long_name = 'Non-Local Component of Ocean Vertical Heat Diffusivity'
difvhonl.units     = 'unitless'
difvhonl[:,0]      = 0.0
difvhonl[:,1:]     = loc_difvhonl

difvsonl           = ncfile.createVariable('difvsonl','f8',('time','nface'))
difvsonl.long_name = 'Non-Local Component of Ocean Vertical Salt Diffusivity'
difvsonl.units     = 'unitless'
difvsonl[:,0]      = 0.0
difvsonl[:,1:]     = loc_difvsonl

if has_difvmonl:
  difvmonl           = ncfile.createVariable('difvmonl','f8',('time','nface'))
  difvmonl.long_name = 'Non-Local Component of Ocean Vertical Momentum Diffusivity'
  difvmonl.units     = 'unitless'
  difvmonl[:,0]      = 0.0
  difvmonl[:,1:]     = loc_difvmonl

# (2c) Close file to actually write
ncfile.close()

