
# input filename
F=img/L1_5.nc

# output prefix for dumped tiff files (the suffix .tiff is always appended)
O=t2/o

python3 dump_names_with_hdf5.py $F | cut -f3 | while read x; do
	echo python3 dump_one_dataset.py $F $x $O`echo $x|tr / -`.tif
done | parallel -j 7
