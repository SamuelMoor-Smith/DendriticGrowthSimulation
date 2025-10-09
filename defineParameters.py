# Ported into python from matlab by Ben Campbell
# Original Author: Samuel Moore-Smith and Jason Yuan
# Further updates by Ben Campbell
'''
Purpose: Define Parameters and Setup Simulation for Gold Memristor Simulation
'''

import numpy as np
from numpy.random import rand
from scipy.constants import Boltzmann, elementary_charge

params = { # create a dictionary to hold parameters
    #--------------------------------
    # Basic Simulation Parameters
    #--------------------------------

    'i' : 1, # Simulation number
    "fin" : 1, #If 1, stops movement of particles once they hit the other end

    #--------------------------------
    # Well-understood Physical Constants
    #--------------------------------

    'k_B' : Boltzmann, # Boltzmann's constant in J/K
    "T_coeff" : 3000, # !!! Temperature coefficient [K^-1] -- this is a guess, not sure what to use here

    #--------------------------------
    # Simulation Length and Particle Count
    # --------------------------------
    
    "n" : 500, # Number of particles
    "simulation_length" : 6, # in seconds
    'dt' : 1e-6, # Time step for simulation

    #--------------------------------
    # Initial position box and material size
    #--------------------------------
    'L_x' : 100, # Size of material in x-direction
    'L_y' : 100, # Size of material in y-direction

    #--------------------------------
    # Heat related parameters
    #--------------------------------
    'T_0' : 300, # Initial temperature in K
    'Q' : 1, # Thermal Energy [J]
    'alpha' : 1e-5, # Thermal diffusivity in m^2/s
    # Jasons code?
    'CT' : 0.24*5e-9, # !!! Heat capacity = specific heat capacity*mass [J/K] http://www2.ucdsb.on.ca/tiss/stretton/database/Specific_Heat_Capacity_Table.html -- check for gold

    #--------------------------------
    # Pinning Force Parameters
    #--------------------------------
    'm' : 100, # number of pinning sites
    # Pinning Site Locations and Forces
    "wpa_repulse" :  2000,
    "wpa_attract" : 500,

    #--------------------------------
    # Applied Electric Field Parameters
    #--------------------------------

    #--------------------------------
    #Residual Stress Parameters
    #--------------------------------
    'E' : 323.5e9, # Young's modulus of AlN in Pa = N/m^2
    'nu' : 0.23, # Poisson's ratio of AlN, not sure if this should be used or be constant
    "k_C" : 8.99e9, # Coulomb's constant in N*m^2/C^2
    "q" : elementary_charge, # Elementary charge in C
    "r_Ag" : 126e-12 # Ionic radius of Ag in m
}

# add to params

# simulation parameters
params["tspan"] = np.arange(0, params["simulation_length"], params["dt"]), # Time span for simulation output

# material and electrode dimensions
params["electrode_width"] = params["L_x"] / 2 # Width of electrodes on material
params["electrode_height"] = params["L_y"] # Height of electrodes on material

#pinning parameters
params['wp_attract']=(-(rand(1,0.5*params['m']))*params['wpa_attract'] - params['wpa_attract']/2) # update the value from a scalar to an array
params['wp_repulse']=((rand(1,0.5*params['m']))*params['wpa_repulse'] + params['wpa_repulse']/2)
params['w_pin'] = np.concatenate((params['wp_attract'], params['wp_repulse']), axis=None) # Combine attractive and repulsive pinning forces into one array
params['x_pin'] = rand(params['m'],1) * params['L_x'] # Randomly distribute pinning sites in x-direction
params['y_pin'] = rand(params['m'],1) * params['L_y'] - 0.5*params['L_y'] # Randomly distribute pinning sites in y-direction
#params['w_pin'] = 100 # Pinning potential amplitude
params['R_pin'] = rand(params['m'],1) * 10 # Pinning potential distance

# Jasons code?
params["k"] = 420/params['L_x']; # !!! Heat transfer coefficient [W/m^2K] from https://www.spiraxsarco.com/learn-about-steam/steam-engineering-principles-and-heat-transfer/heat-transfer

#--------------------------------
#
#--------------------------------
