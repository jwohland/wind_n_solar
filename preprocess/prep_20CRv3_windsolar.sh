#!/bin/sh
#BSUB -J prep_20CRv3_windsolar
#BSUB -n 1
#BSUB -W 12:00
#BSUB -o ../logs/prep_20CRv3_windsolar.txt

# 1) restrict 20CR data to European domain (75/-15/30/35)
# 2) regrid 20CR to CERA20C grid


module load cdo

data_dir=/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/20CRv3
target_grid=/cluster/work/apatt/wojan/data/CERA20C/monthly/CERA20C_ssrc_1901.nc

cd /cluster/work/apatt/wojan/data/20CRv3/wind.100m || exit

for year in {1901..2009}
do
  echo $year
  cd $year || exit
  for filename in *.nc
  do
    if [ -f  ${data_dir}/${filename}_CERAgrid.nc ]; then
      echo "exists already"
    else
      # 1) select box that is 5 degrees larger than target grid to avoid boundary effects
      cdo sellonlatbox,-20,40,25,80 UGRD100m.1901_mem001.nc ${data_dir}/tmp.nc
      # 2) bilinear interpolation to CERA20C grid
      cdo remapbil,${target_grid} ${data_dir}/tmp.nc ${data_dir}/${filename}_CERAgrid.nc
      rm ${data_dir}/tmp.nc
    fi
  done
done