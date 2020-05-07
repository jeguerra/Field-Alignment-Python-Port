#!/bin/bash
#PBS -A det
#PBS -N memFA
#PBS -l procs=2
#PBS -l walltime=8:00:00
#PBS -d .
/apps/matlab/R2014b/bin/matlab -nodisplay -nodesktop -nojvm -nosplash < foo.txt > out 2>&1
echo "qsub ok"
