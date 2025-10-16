import warp as wp
from numpy import arange, array, len, append, float32



# RK4 integration algorithm
def rk4(f:function, initial_conditions:list, t_final:float, dt:float): #integrates with RK4. initial conditions must be given with time first - ODE must take time as first argument
  
    # f should output an n-list of floats. Each entry of this list gives the derivative of the ith 1st order ODE.
    tpoints = arange(initial_conditions[0],t_final,dt,dtype=np.float32) # time points to evaluate at
    tpoints = wp.from_numpy(tpoints, dtype=wp.float32) # move to a warp array
    # we want this to generalize to a system of n ODEs
    n = len(initial_conditions) - 1 # subtract one because we don't want to count t_0

    @wp.kernel # warp kernel to do the rk4 integration step t_i at once on the GPU
    def rk4_step(t_i: wp.float32, xpoints: dict, dt: float, k1: wp.vec(n, dtype=wp.float32), k2: wp.vec(n, dtype=wp.float32), k3: wp.vec(n, dtype=wp.float32), k4: wp.vec(n, dtype=wp.float32)):
        i = wp.tid() # index of the ODE we are integrating

        xpoints["x_%i"%i].append(current_conditions[i])
        xpoints["x_%i"%i][-1] += (k1[i]+2*k2+2*k3+k4)/6 # i have reservation about the efficiency of all this list generation and unpacking but I can't think of a better way to do this in complete generality


    xpoints = { "x_%i"%i : [initial_conditions[i+1]] for i in range(n)} # initialize the list of coordinates with the initial conditions
    # use a dictionary so that it is variable length
    for t_i in range(len(tpoints)-1):
        current_conditions = [ xpoints["x_%i"%j][t_i] for j in range(n) ]   #take the last of each coordinate; use t_i index so that on loop t_i we only use the conditions from t_i-1 not a mix (index of -1 would cause a mix).
        k1 = dt*f(tpoints[t_i],*current_conditions) # since f may be a function of arbitray many coordinates we unpack the list
        k2 = dt*f(tpoints[t_i]+0.5*dt,*[current_conditions[j] + 0.5*k1 for j in range(n)])
        k3 = dt*f(tpoints[t_i]+0.5*dt,*[current_conditions[j] + 0.5*k2 for j in range(n)])
        k4 = dt*f(tpoints[t_i]+dt,*[current_conditions[j] + k3 for j in range(n)])
        # each of the ki's should be numpy arrays of length n at this point - need to convert to warp vectors
        n_vec = wp.vec(n, dtype=wp.float32) # create a warp vector type of length n
        k1 = n_vec(k1) # convert each ki to a warp vector
        k2 = n_vec(k2)
        k3 = n_vec(k3)
        k4 = n_vec(k4)
        #launch the rk4 step kernel on the GPU
        wp.launch(rk4_step, dim=n,inputs=[tpoints[t_i],xpoints,dt, k1, k2, k3, k4]) # launch the rk4 step kernel on the GPU
    return tpoints, [xpoints["x_%i"%i] for i in range(n)] # times and list of lists