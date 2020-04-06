#!/bin/sh
#$ -cwd
#$ -N JOBNAME
#$ -o output
#$ -V
#$ -S /bin/bash
#$ -pe x40 40
source /opt/intel/impi/2018.2.199/intel64/bin/mpivars.sh
export I_MPI_FALLBACK_DEVICE=0
export I_MPI_PIN=0
export I_MPI_MPD_RSH=/usr/bin/rsh
export I_MPI_HYDRA_BOOTSTRAP=rsh

mpirun -genv I_MPI_FABRICS=shm:dapl -genv I_MPI_DEBUG=5 -genv I_MPI_PIN=1 -genv I_MPI_PIN_MODE=pm -np 40 /usr/local/vasp.5.3.5_intelmpi/bin/vasp 