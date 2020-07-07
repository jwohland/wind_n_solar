#!/bin/sh
#BSUB -J calculate_windpower_test
#BSUB -n 1
#BSUB -R "rusage[mem=10000]"
#BSUB -W 24:00
#BSUB -r
#BSUB -o ../../../logs/calculate_windpower_test.txt

python windpower.py 10