cd /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/trend_corrected || exit
module load cdo

for filename in *.nc
do
  cdo timmean $filename annual/annual_$filename
done

rm annual/annual_CERA20C_correcteds100_1905_representative.nc  #do not contain a full year of data
rm annual/annual_CERA20C_correcteds100_2005_representative.nc

chmod 400 annual/*

# annual means of 20CRv3 input data
cd /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/20CRv3 || exit
module load nco
for filename in *.nc
do
  ncpdq -a time,latitude,longitude,number $filename annual/reordered_$filename
  cdo timmean annual/reordered_$filename annual/annual_$filename
  chmod 400 annual/annual_$filename
  rm annual/reordered*
done