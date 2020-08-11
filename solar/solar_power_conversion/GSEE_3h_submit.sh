#!/bin/sh
#BSUB -J GSEE_3h[1-109]
#BSUB -R "rusage[mem=2000]"
#BSUB -W 72:00
#BSUB -r
#BSUB -oo ../../../logs/GSEE_3h_%I.txt

for ensemble_member in {1..10}
do
  echo ${ensemble_member}
  for scenario in {1..3}
  do
    echo ${scenario}
    python GSEE_3h.py ${LSB_JOBINDEX} ${ensemble_member} ${scenario}
  done
done

