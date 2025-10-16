#purpose: 

#import parameters
from defineParameters import params
#import packages
from numpy.random import rand
import numpy as np
from scipy.integrate import solve_ivp
from time import time # allows us to check run times

I_saved = []
t_saved = []

#------------------------
# Unpack parameters
#------------------------
n = params['n']
To = params['T_0']
tspan = params['tspan']

Lx = params['L_x'] # Size of material
Ly = params['L_y'] # Size of material
electrode_width = params['electrode_width']
electrode_height = params['electrode_height']

#------------------------
# Randomize initial particle positions
#------------------------

xo = 10. * rand(n);  # initial x-positions
yo = electrode_height * rand(n) - 0.5*electrode_height  # initial y-positions

#------------------------
# Zero all initial particle velocities
#------------------------
vox = np.zeros(n)
voy = np.zeros(n)

#------------------------
# Initial config vector
#------------------------

initial = np.concatenate((xo, vox, yo, voy)) # initial state vector

#------------------------
# Run ode simulation of the differentiable force function
#------------------------

# Define the function to be solved
from calculateForces import calculateForces

print('starting the solver')

V = params['V']

start_time = time() # start clock

t, states = solve_ivp(calculateForces,tspan,initial,'RK45')

end_time = time() # end clock

elapsed_time = end_time - start_time # simulation run time

# TODO: Go over this
#------------------------

X_pos = states[0:n:1]      # First n entries are x positions
X_vel = states[n:2*n:1]  # Second group of n entries are x-velocities
Y_pos = states[2*n:3*n:1]  # Third group of n entries are y-positions
Y_vel = states[3*n:4*n:1]  # Last group of n entries are y-velocities

#produce 2D plot of particle positions, save as gif

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots()

# Set axis limits and labels once

ax.set_xlim(-electrode_width - 5, Lx+electrode_width + 5)
ax.set_ylim(-electrode_height/2 - 5, electrode_height/2 + 5)
ax.set_xlabel('X Position (units to be determined)')
ax.set_ylabel('Y Position (units to be determined)')

def animate(): #define animation function
    return

