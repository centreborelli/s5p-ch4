#!/usr/bin/python3

# DISCLAIMER:
#
# The people behind the HDF5 format should be tarred and feathered.
# If I am willing to subject myself to the torture of dealing with this
# shitfuckery is only because somebody has to extract the damn numbers from
# the damn files.
# --eml


# extract input argument
import sys
filename_in  = sys.argv[1]
dataset_id   = sys.argv[2]
filename_out = sys.argv[3]

# load file
import h5py
f = h5py.File(filename_in, "r")

# select dataset with the lovely syntax of the h5py object
x = f[dataset_id][()]

# in case the array is four dimensional but the first dimension is trivial
# (this is actually the typical convention for HDF5 multispectral images)
# then reove the first dimension and keep the rest of them
if len(x.shape) == 4 and x.shape[0] == 1:
	x = x[0,:,:,:]

# same thing for three-dimensional array
if len(x.shape) == 3 and x.shape[0] == 1:
	x = x[0,:,:]

# save the array to the output file
import iio
iio.write(filename_out, x)
