import numpy as np

import numpy as np

def get_next_index(tot_weights, x, y, currx, curry, curr_indexes, L):
    """
    Translated from MATLAB getNextIndex.m
    Determines the next index and corresponding (x,y) position for each electron.
    """
    num_e, cols = tot_weights.shape
    n = cols - 2  # number of particle targets (minus 2 electrodes)

    LEFT_ELECTRODE = n + 1
    RIGHT_ELECTRODE = n + 2

    # Normalize hopping probabilities
    row_sums = np.sum(tot_weights, axis=1, keepdims=True)
    W_normalized = tot_weights / row_sums
    W_cdf = np.cumsum(W_normalized, axis=1)

    # Random choice via CDF
    rands = np.random.rand(num_e, 1)
    C = rands < W_cdf
    next_indices = np.argmax(C, axis=1) + 1  # MATLAB is 1-based, Python is 0-based; adjust later if needed

    # Initialize next positions
    nextx = np.zeros(num_e)
    nexty = np.zeros(num_e)

    # Beginning electrode (x=0, y=nan)
    mask_left = next_indices == LEFT_ELECTRODE
    nextx[mask_left] = 0
    nexty[mask_left] = np.nan

    # Right electrode (x=L, y=nan)
    mask_right = next_indices == RIGHT_ELECTRODE
    nextx[mask_right] = L
    nexty[mask_right] = np.nan

    # Particles in between
    mask_particles = (~mask_left) & (~mask_right)
    nextx[mask_particles] = x[next_indices[mask_particles] - 1]
    nexty[mask_particles] = y[next_indices[mask_particles] - 1]

    return next_indices, nextx, nexty


"""def getNextIndex(tot_weights, x, y, currx, curry, curr_indexes, L):
    '''
    Python version of the MATLAB function getNextIndex.
    Computes the next hopping indices and positions for each electron
    based on transition probabilities.
    '''

    num_e = tot_weights.shape[0]              # total number of electrons
    n = tot_weights.shape[1] - 2              # number of possible particle targets (minus 2 electrodes)
    LEFT_ELECTRODE = n                        # n+1 in MATLAB (zero-based index)
    RIGHT_ELECTRODE = n + 1                   # n+2 in MATLAB (zero-based index)

    # Normalize hopping probabilities per electron
    W_normalized = tot_weights / np.sum(tot_weights, axis=1, keepdims=True)
    W_cdf = np.cumsum(W_normalized, axis=1)   # cumulative distribution function for transitions
    rands = np.random.rand(num_e, 1)          # random numbers for each electron

    # Find first index where rands < CDF
    C = rands < W_cdf
    next_indices = np.argmax(C, axis=1)       # MATLAB's max(...,[],2) → argmax along axis 1

    # Initialize nextx and nexty
    nextx = np.zeros(num_e)
    nexty = np.zeros(num_e)

    # Set beginning positions to (0, 0) [y = nan]
    mask_left = next_indices == LEFT_ELECTRODE
    nextx[mask_left] = 0
    nexty[mask_left] = np.nan

    # Set end positions to (L, 0) [y = nan]
    mask_right = next_indices == RIGHT_ELECTRODE
    nextx[mask_right] = L
    nexty[mask_right] = np.nan

    # Set nextx and nexty for electrons hopping to actual particles
    not_begin_or_end = ~(mask_left | mask_right)
    nextx[not_begin_or_end] = x[next_indices[not_begin_or_end]]
    nexty[not_begin_or_end] = y[next_indices[not_begin_or_end]]

    return next_indices, nextx, nexty
"""