#####################################################################
##################### RUN DENDRITIC SIMULATION ######################
#####################################################################

from defineParameters import params
from runSimulation import run_simulation

if __name__ == "__main__": # Only run the following code if this file is being run directly, not when it’s imported by another script
    t, states = run_simulation(params)
    print("Simulation finished successfully!")
    print(f"Final time: {t[-1]:.3f}")
    print(f"State vector shape: {states.shape}")