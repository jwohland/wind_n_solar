#!/bin/sh
#BSUB -J calculate_low-pass[1-54]
#BSUB -n 1
#BSUB -R "rusage[mem=4096]"
#BSUB -W 04:00
#BSUB -r
#BSUB -o /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/docs/filter/calculate_low-pass.txt

if [ "${LSB_JOBINDEX}" -le "30" ]; then
  noaa_index=-1
  if [ "${LSB_JOBINDEX}" -le "10" ];then
    cera_index=$( expr ${LSB_JOBINDEX} - 1 )
  elif [ "${LSB_JOBINDEX}" -le "20" ];then
    cera_index=$( expr ${LSB_JOBINDEX} - 11 )
  else
    cera_index=$( expr ${LSB_JOBINDEX} - 21 )
  fi
else
  cera_index=-1
  if [ "${LSB_JOBINDEX}" -le "38" ];then
    noaa_index=$( expr ${LSB_JOBINDEX} - 31 )
  elif [ "${LSB_JOBINDEX}" -le "46" ];then
    noaa_index=$( expr ${LSB_JOBINDEX} - 39 )
  else
    noaa_index=$( expr ${LSB_JOBINDEX} - 47 )
  fi
fi

if [ "${LSB_JOBINDEX}" -le "10" ]; then
  lat_index=15
elif [ "${LSB_JOBINDEX}" -le "20" ]; then
  lat_index=25
elif [ "${LSB_JOBINDEX}" -le "30" ]; then
  lat_index=35
elif [ "${LSB_JOBINDEX}" -le "38" ]; then
  lat_index=15
elif [ "${LSB_JOBINDEX}" -le "46" ]; then
  lat_index=25
else
  lat_index=35
fi

python calculate_low-pass_filtered.py $cera_index $noaa_index $lat_index

