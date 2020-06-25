# This script has to run on a NERSC machine
# Please use the machines dtn[01 - 04].nersc.gov to issue these commands, as they are optimized for HPSS extraction.
# ssh jwohland@dtn01.nersc.gov

cd /global/cscratch1/sd/jwohland || exit
for year in {1903..2009}
do
  for mem in {001..008}
  do
    # U component
    if [ -f "${year}/UGRD100m.${year}_mem${mem}.nc" ]; then
      echo "${year}/UGRD100m.${year}_mem${mem}.nc exists."
    else
      htar -Hnostage -xvf /home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/subdaily/UGRD100m/UGRD100m_${year}.tar ${year}/UGRD100m.${year}_mem${mem}.nc
    fi
    # V component
    if [ -f "${year}/VGRD100m.${year}_mem${mem}.nc" ]; then
      echo "${year}/VGRD100m.${year}_mem${mem}.nc exists."
    else
      htar -Hnostage -xvf /home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/subdaily/VGRD100m/VGRD100m_${year}.tar ${year}/VGRD100m.${year}_mem${mem}.nc
    fi
  done
done
