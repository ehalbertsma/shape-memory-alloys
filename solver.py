import numpy as np

import math_functions
import plot_solution
import model_parameters

# initialize all parameters
model = Model_Parameters()
# parameters form the comparam.m file TODO
model.compute_parameters()

model.integrate()

# run as a script


#
#
#plot_solution()