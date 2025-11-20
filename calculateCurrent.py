from defineParameters import params
import numpy as np
from getNextIndex import get_next_index

def calculate_current(x, y, L, V, lambda_, Rt, steps, num_e):
    """
    Translated from MATLAB calculateCurrent.m
    Simulates electron hopping and returns current I.
    """
    n = len(x)
    threshold_dist = 50

    # Skip simulation if device too long
    min_distance_to_end = np.min(L - x) # minimum distance of a particle to the end electrode
    if min_distance_to_end > threshold_dist:
        return 0

    LEFT_ELECTRODE = n + 1
    RIGHT_ELECTRODE = n + 2

    # Start all electrons at left electrode
    curr_indexes = np.full(num_e, LEFT_ELECTRODE) # array of length of the number of simulated electrons filled with n+1
    currx = np.zeros(num_e)
    curry = np.full(num_e, np.nan) # array of length of the number of simulated electrons filled with nan

    for _ in range(steps):
        dx = np.tile(x, (num_e, 1)) - currx[:, np.newaxis]
        dy = np.tile(y, (num_e, 1)) - curry[:, np.newaxis]
        dy[np.isnan(dy)] = 0
        d = np.sqrt(dx ** 2 + dy ** 2) # not sure why this is here

        positive = dx >= 0
        resist = Rt * np.exp(d / lambda_) # resistance 
        dist_weights = positive * (1 / resist)

        sigmoid_end = np.ones(num_e)
        resist_end = Rt * np.exp((L - currx) / lambda_)
        to_end = sigmoid_end / resist_end
        to_end[curr_indexes == RIGHT_ELECTRODE] = 1000  # stay at end
        to_begin = np.zeros(num_e)

        tot_weights = np.hstack([dist_weights, to_begin[:, None], to_end[:, None]])

        curr_indexes, currx, curry = get_next_index(tot_weights, x, y, currx, curry, curr_indexes, L)

    I = np.sum(currx == L)
    return I