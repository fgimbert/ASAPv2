from ase.build import bulk
from ase.calculators.vasp import Vasp2

import os

os.environ['VASP_PP_PATH'] = VASPPPPATH
os.environ['ASE_VASP_COMMAND'] = ASEVASPCOMMAND

print('Bonjour Ici asepy.py !')

si = bulk('Si')

print(os.environ['VASP_PP_PATH'])

mydir = 'Test'    # Directory where we will do the calculations

# Make self-consistent ground state
calc = Vasp2(kpts=(4, 4, 4), directory=mydir)

si.set_calculator(calc)
si.get_potential_energy()  # Run the calculation

# Non-SC calculation along band path
kpts = {'path': 'WGX',     # The BS path
        'npoints': 30}     # Number of points along the path

calc.set(isym=0,           # Turn off kpoint symmetry reduction
         icharg=11,        # Non-SC calculation
         kpts=kpts)

# Run the calculation
si.get_potential_energy()

