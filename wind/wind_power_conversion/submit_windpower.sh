#!/bin/sh
#BSUB -J calculate_windpower_[1-101]
#BSUB -n 1
#BSUB -R "rusage[mem=7000]"
#BSUB -W 24:00
#BSUB -r
#BSUB -oo ../../../logs/calculate_windpower_%I.txt

python windpower.py $LSB_JOBINDEX