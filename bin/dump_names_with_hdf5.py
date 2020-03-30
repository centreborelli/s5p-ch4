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
filename_in = sys.argv[1]

# load file
import h5py
f = h5py.File(filename_in, "r")

# recursively print keys
def recursively_print_keys(d, p):
	for k in d.keys():
		v = d[k]
		n = f"{p}/{k}"
		if hasattr(v, "keys"):
			recursively_print_keys(v, n)
		else:
			# The following conditionals are necessary because,
			# contrarily to what is explained on the tutorial,
			#  ( http://docs.h5py.org/en/stable/quick.html )
			# some Datasets do *not* have shape nor size.
			s = v.size  if hasattr(v, "size")  else 0
			z = v.shape if hasattr(v, "shape") else (0,)
			print(f"{s}\t{'x'.join(str(x) for x in z)}\t{n}")

#recursively_print_keys(f, filename_in)
recursively_print_keys(f, "")
