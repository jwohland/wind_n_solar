#!/bin/sh
#BSUB -J GSEE_3h[1-109]
#BSUB -R "rusage[mem=8000]"
#BSUB -W 24:00
#BSUB -r
#BSUB -oo ../../logs/GSEE_3h_%I.txt

for ensemble_member in {1..10}
do
  echo ${ensemble_member}
  python GSEE_3h.py  ${LSB_JOBINDEX} ${ensemble_member}
done

