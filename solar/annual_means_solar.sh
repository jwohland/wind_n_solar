
cd /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/solar_power || exit

module load cdo

# delete from here onwards
for filename in *number_0.nc
do
  cdo timmean $filename annual/annual_$filename
done
rm annual/*1905*.nc  #do not contain a full year of data
rm annual/*2005*.nc
