import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

energy = np.logspace(-4, 3, 1000)
temperature = np.linspace(1e4, 1e8, 4)
colors = mpl.cm.inferno(np.linspace(0, 1, len(temperature)))

# cgs units
boltz = 1.380658e-16 # erg K^-1
speed_of_light = 2.99792458e10 # cm s^-1
planck = 6.62607015e-27 # erg s













