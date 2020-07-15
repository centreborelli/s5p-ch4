
# This is an example on how to download a L1 / L2 pair and crop them at the
# same localization

# Input: orbit number, longitude/latitude bounding box, grid size, frequencies
# Output: two images, with the L2 methane data and the L1 frequency slice
# output files:
#	out_ch4.npy
#	out_l1.npy

# Note: since tsd does not allow (yet) to search for orbit number, changing the
# values of this script needs some manual trickery, by looking at the scihub
# interface for the orbit number corresponding to a date/site.  Then, you need
# to edit the "CFG" variable below.

ORBIT=11550
LONLAT="-10 20 10 40"
GRID="1000 1000"
FREQUENCIES="2358.91 2359.17 2369.55 2370.4"


# STEP 0. silly config (I don't know how to call the tsd command line otherwise)
export PYTHONPATH=$HOME/src/tsd
tsd="python3 $HOME/src/tsd/tsd/search_scihub.py"


# STEP 1. Build list of URLs using TSD
# NOTE: change the values of the CFG variable if you need to change the orbit
# NOTE2: a future version of TSD will accept an orbit parameter directly
#
CFG="--lon 2 --lat 48 -s 2020-1-5 --satellite Sentinel-5P"
CRED="--user=s5pguest --password=s5pguest"
$tsd $CFG --product-type L1B_RA_BD8 > l1_b8.j
$tsd $CFG --product-type L2__CH4___ > l2_ch4.j
cat l*j                             |\
	tr , '\n'                   |\
	grep -E '(link|title)'      |\
	sed 's/.*\"://'             |\
	sed 's/\$/\\\\$/'           |\
	xargs -n 2                  |\
	grep _${ORBIT}_             |\
	while read x y; do
		echo wget $CRED \"$y\" -O $x.nc
done > wgets.sh

# STEP 2. Download the urls (takes about 40 minutes)
. wgets.sh


# STEP 3. give reasonalbe names to the downloaded files
ln -s S5P_*_L2__CH4_*_${ORBIT}_*.nc i${ORBIT}_CH4.nc
ln -s S5P_*_RA_BD8_*_${ORBIT}_*.nc i${ORBIT}_BD8.nc


# STEP 4. localize the L1 and L2 images to the requested grid
DSET=methane_mixing_ratio_bias_corrected
./localize $LONLAT $GRID i${ORBIT}_CH4.nc -s crop_ch4.npy -d $DSET
./localize $LONLAT $GRID i${ORBIT}_BD8.nc -s crop_l1.npy  -sf "$FREQUENCIES"

# optional: set nodata value to NAN
plambda crop_ch4.npy "x 1e30 > nan x if" -o crop_ch4.npy
