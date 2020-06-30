#!/bin/sh
#BSUB -J calculate_low-pass_rep[1-42]
#BSUB -n 1
#BSUB -R "rusage[mem=15000]"
#BSUB -W 04:00
#BSUB -r
#BSUB -o /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/docs/filter/calculate_low-pass_rep.txt

cera_index=7
noaa_index=5
lat_index=$((${LSB_JOBINDEX}-1))

python calculate_low-pass_filtered.py -1 $noaa_index $lat_index
python calculate_low-pass_filtered.py $cera_index -1 $lat_index