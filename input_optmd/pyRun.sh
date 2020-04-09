#!/bin/sh
#$ -cwd
#$ -N pyvasp
#$ -o output
#$ -V
#$ -S /bin/bash
#$ -pe x40 80
source /opt/intel/impi/2018.2.199/intel64/bin/mpivars.sh
export I_MPI_FALLBACK_DEVICE=0
export I_MPI_PIN=0
export I_MPI_MPD_RSH=/usr/bin/rsh
export I_MPI_HYDRA_BOOTSTRAP=rsh

python ontheflymd_setup.py