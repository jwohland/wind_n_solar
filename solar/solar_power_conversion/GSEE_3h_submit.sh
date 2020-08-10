#!/bin/sh
#BSUB -J GSEE_3h[1-3270]
#BSUB -R "rusage[mem=2000]"
#BSUB -W 04:00
#BSUB -r
#BSUB -oo ../../../logs/GSEE_3h_%I.txt

i=0
for year in {1..109}
do
  for ensemble_member in {1..10}
  do
    for scenario in {1..3}
    do
      #echo ${scenario}
      if [[ i -eq ${LSB_JOBINDEX} ]] ; then
        echo ${year} ${ensemble_member} ${scenario}
        python GSEE_3h.py  ${year} ${ensemble_member} ${scenario}
      fi
      i=$(( $i + 1))
    done
  done
done
