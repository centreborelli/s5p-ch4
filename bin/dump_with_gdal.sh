
mkdir -p tmp/L1
mkdir -p tmp/L2_NO2

gdal_translate -sds img/L1.nc tmp/L1/i
gdal_translate -sds img/L2_NO2.nc tmp/L2_NO2/i
