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
par_x0       = int(sys.argv[3])
par_x1       = int(sys.argv[4])
par_y0       = int(sys.argv[5])
par_y1       = int(sys.argv[6])
filename_out = sys.argv[7]

# load file
import h5py
f = h5py.File(filename_in, "r")

# extract dataset section (only works for four-dimensional data)
x = f[dataset_id][0, par_x0:par_x1, par_y0:par_y1, :]

# save the array to the output file
import iio
iio.write(filename_out, x)
