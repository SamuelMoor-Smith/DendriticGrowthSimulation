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
    min_distance_to_end = np.min(L - x)
    if min_distance_to_end > threshold_dist:
        return 0

    LEFT_ELECTRODE = n + 1
    RIGHT_ELECTRODE = n + 2

    # Start all electrons at left electrode
    curr_indexes = np.full(num_e, LEFT_ELECTRODE)
    currx = np.zeros(num_e)
    curry = np.full(num_e, np.nan)

    for _ in range(steps):
        dx = np.tile(x, (num_e, 1)) - currx[:, np.newaxis]
        dy = np.tile(y, (num_e, 1)) - curry[:, np.newaxis]
        dy[np.isnan(dy)] = 0
        d = np.sqrt(dx ** 2 + dy ** 2)

        positive = dx >= 0
        resist = Rt * np.exp(d / lambda_)
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


"""def calculateCurrent(x, y, L, V, lam, Rt, steps, num_e):
    
    n = len(x) # number of particles
    threshold_dist = 50

    # If the minimum distance to the end is greater than the threshold distance,
    # then set current to 0 and return to avoid unnecessary calculations
    # probably makes sense due to exponential decay of quantum tunneling
 
    min_distance_to_end = np.min(L - x)
    if min_distance_to_end > threshold_dist:
        I = 0
        return

    # Start with electrons at random particles in device
    #curr_indexes = np.random.choice(np.arange(100), size=p, replace=False) # Depreciated

    # start with num_e electrons at left electrode
    LEFT_ELECTRODE = n
    curr_indexes = LEFT_ELECTRODE * np.ones(num_e) # left electrode is index n

    # Get the current positions of the electrons
    currx = np.zeros(num_e, dtype=float) # start at x=0
    curry = np.zeros((num_e,), dtype=float) # start at y = nan' (y doesn't matter at electrode)


    # in case num_e changes - not happening anymore
    # num_e = length(currx);

    for step in range(1, steps + 1):
        dx = np.tile(np.transpose(x),(num_e,1)) - currx
        dy = np.tile(np.transpose(y),(num_e,1)) - curry

        d = np.sqrt(dx**2 + dy**2) # distance of each electron
        
        '''if type(dx) == list:
            dx = np.array(dx) #allows the use of the bollean array gernation below''' # not needed, above generates list
        positive = (dx >= 0) # bollean array
        
        # voltage = V * dx / L;
        # numerator = 1 + 2000*max(voltage, 0);
        # sigmoid = 1./(1+exp(-voltage));
        # sigmoid = V * ones(size(dx));    % Default to V
        # sigmoid(dx == 0) = 1;            % Set to 1 where dx == 0
        resist = Rt*np.exp(d/lam)
        
        dist_weights = dx[positive] / resist
        # dist_weights(isnan(dist_weights)) = 0;
        
        # voltage_end = V.*(L-currx')./L;
        # sigmoid_end = 1./(1+exp(-voltage_end));
        sigmoid_end = np.ones(num_e)    # Default to V
        # sigmoid_end = V;
        resist_end = Rt*np.exp((L-currx)/lam)
        
        to_end = (sigmoid_end/resist_end) * np.ones(num_e,1)
        to_end(curr_indexes==n+1) = 1000 # if already at end, always stay there
        
        to_begin = np.zeros(num_e)

        tot_weights = np.array([dist_weights, to_begin, to_end])
        
        curr_indexes,currx,curry = getNextIndex(tot_weights,x,y,currx,curry,curr_indexes, L)

    I = sum(currx==L)

    return [I]"""