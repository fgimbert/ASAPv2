import numpy as np

from flare import otf, kernels
from flare.gp import GaussianProcess
from flare.mgp.mgp import MappedGaussianProcess
from flare.ase.calculator import FLARE_Calculator
import flare.kernels.mc_simple as mc_simple

# ---------- create gaussian process model -------------------
kernel = mc_simple.two_plus_three_body_mc
kernel_grad = mc_simple.two_plus_three_body_mc_grad
energy_force_kernel = mc_simple.two_plus_three_mc_force_en
energy_kernel = mc_simple.two_plus_three_mc_en

hyps = np.array([0.1, 1., 0.001, 1, 0.03])
two_cut = 4.0
three_cut = 4.0
cutoffs = np.array([two_cut, three_cut])
hyp_labels = ['sig2', 'ls2', 'sig3', 'ls3', 'noise']
opt_algorithm = 'BFGS'

gp_model = GaussianProcess(kernel, kernel_grad, hyps, cutoffs,
                           energy_kernel=energy_kernel,
                           energy_force_kernel=energy_force_kernel,
                           hyp_labels=hyp_labels,
                           opt_algorithm=opt_algorithm, par=True)



# ------------ create ASE's flare calculator -----------------------
flare_calc = FLARE_Calculator(gp_model, par=True, use_mapping=False)
