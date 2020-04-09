import os

# set up executable
label = 'C'
no_cpus = 32
npool = 32
from ase.calculators.vasp import Vasp2


os.environ['VASP_PP_PATH'] = "/home/f-gimbert/potvasp"
os.environ['ASE_VASP_COMMAND'] = "mpirun -genv I_MPI_FABRICS=shm:dapl -genv I_MPI_DEBUG=5 -genv I_MPI_PIN=1 -genv I_MPI_PIN_MODE=pm -np 40 /home/share/vasp.5.4.4/bin/vasp_std"


# create ASE calculator
mydir = 'Test_OptMD'    # Directory where we will do the calculations

# Make self-consistent ground state
dft_calc = Vasp2(kpts=(4, 4, 4), directory=mydir)