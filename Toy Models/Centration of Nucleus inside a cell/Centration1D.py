#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Centration1D.py
Created on Sat Sep 10, 2016
@author: akkimura (Python 3.5)
Description: simple 1D simulation for nuclear centration 
"""

import numpy as np
from numpy import pi, exp
import matplotlib.pyplot as plt

# Initialize parameters
L = 25  # half length of the cell [um]
N = 100  # number of microtubules (half)
R = 5  # the Stokes radius of the nucleus [um]
eta = 1  # the viscosity of the cytoplasm [pN s/um^2]
f = 1  # pulling force of a single motor [pN]
c = 0.01  # the density of the motors on microtubules [/um]
Vg = 1  # growth speed of microtubule [um/s]
Fcat = 0.2  # the catastrophe frequency of microtubules [/s]
kappa = 10  # flexural rigidity of microtubules [pN um^2]
K = 6 * pi * eta * R  # coefficient for nucleus drag [pN s/um]
dt = 2  # sec per simulation step [s]
totSTEP = 120  # the number of simulation steps

X0 = 20  # initial position of the nucleus [um]

# Initialize variables
Xcyto = np.zeros(totSTEP + 1)
Xcort = np.zeros(totSTEP + 1)
Xpush = np.zeros(totSTEP + 1)
Xcyto[0] = X0
Xcort[0] = X0
Xpush[0] = X0

# Calculate
for st in range(totSTEP):  # i starts from 0 to totSTEP-1
    
    # cytoplasmic pulling model
    Fcyto = -2 * f * c * N * Xcyto[st]
    Xcyto[st + 1] = Xcyto[st] + dt * Fcyto / K 
    
    # cortical pulling model
    if Xcort[st] > L:
        Fcort_plus = 0 
    else:
        Fcort_plus = f * N * exp(-Fcat * (L - Xcort[st]) / Vg)
    if Xcort[st] < -L:
        Fcort_minus = 0
    else:
        Fcort_minus = -f * N * exp(-Fcat * (L + Xcort[st]) / Vg)
    Fcort = Fcort_plus + Fcort_minus
    Xcort[st + 1] = Xcort[st] + dt * Fcort / K
    
    # pushing model
    Fpush = N * pi * pi * kappa * (1 / ((L + Xpush[st]) ** 2) - 1 / ((L - Xpush[st]) ** 2))
    Xpush[st + 1] = Xpush[st] + dt * Fpush / K

# Plot results
t_values = np.linspace(0, totSTEP * dt, totSTEP + 1)
plt.plot(t_values, Xcyto, color="b", label="Cyto Pull")
plt.plot(t_values, Xcort, color="r", label="Cort Pull")
plt.plot(t_values, Xpush, color="g", label="Push")

# Graph modifications
ax = plt.gca()  # get current axis
ax.set_title("Centration Simulation", size=24, weight='bold')
ax.set_xlabel("Time [s]", size=18)
ax.set_ylabel("Position of the Nucleus [um]", size=18)
plt.legend(loc='center right')  # show legends

# Save the figure
fig = plt.gcf()  # get current figure
plt.savefig("centration_test.png")
