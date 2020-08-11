
cd /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/solar_power || exit

module load cdo

cd both_constant || exit
  for filename in *.nc
  do
    cdo timmean $filename annual/annual_$filename
  done
cd .. || exit

cd tilt_constant || exit
  for filename in *.nc
  do
    cdo timmean $filename annual/annual_$filename
  done
cd .. || exit

cd neither_constant || exit
  for filename in *.nc
  do
    cdo timmean $filename annual/annual_$filename
  done
cd .. || exit