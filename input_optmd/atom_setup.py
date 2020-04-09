import numpy as np
from ase import units
from ase.spacegroup import crystal
from ase.build import bulk

np.random.seed(12345)

a = 3.52678
super_cell = bulk('C', 'diamond', a=a, cubic=True)