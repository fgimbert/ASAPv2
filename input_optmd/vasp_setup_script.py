import os

# set up executable
label = 'C'
no_cpus = 32
npool = 32
from ase.calculators.vasp import Vasp2


os.environ['VASP_PP_PATH'] = VASPPPPATH
os.environ['ASE_VASP_COMMAND'] = ASEVASPCOMMAND


# create ASE calculator
mydir = 'Test_OptMD'    # Directory where we will do the calculations

# Make self-consistent ground state
dft_calc = Vasp2(kpts=(4, 4, 4), directory=mydir)