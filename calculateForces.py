import numpy as np
from defineParameters import params
from calculateCurrent import calculateCurrent

# Global simulation state
tindex = 0
V = 0
I_saved, t_saved, I_last = [], [], 0
rand_dirs_global_x, rand_dirs_global_y = None, None


def calculate_forces(t, states, params):
    """
    Python translation of MATLAB calculateForces.m
    Computes the time derivative dxdt for all particle states.
    """
    global tindex, V, I_saved, t_saved, I_last
    global rand_dirs_global_x, rand_dirs_global_y

    n = params["n"]
    eta = params["eta"]

    # Initialize global direction vectors if empty
    if rand_dirs_global_x is None:
        rand_dirs_global_x = 2 * np.random.randint(0, 2, n) - 1
        rand_dirs_global_y = 2 * np.random.randint(0, 2, n) - 1

    # ------------------------
    # Unpack state variables
    # ------------------------
    x_p = states[0:n]
    x_v = states[n:2 * n]
    y_p = states[2 * n:3 * n]
    y_v = states[3 * n:4 * n]
    T = states[4 * n]

    # ------------------------
    # Compute current occasionally
    # ------------------------
    if tindex < len(params["tspan"]) and t > params["tspan"][tindex]:
        rand_dirs_global_x = 2 * np.random.randint(0, 2, n) - 1
        rand_dirs_global_y = 2 * np.random.randint(0, 2, n) - 1

        if tindex % 2 == 0:
            I_last = calculateCurrent(
                x_p, y_p, params["L_x"], V,
                params["lambda"], params["Rt"], params["steps"], params["num_e"]
            )

        tindex += 1

        if tindex > 1000 and V != 0:
            V = 0
            print("Voltage set to zero after 1000 iterations")

        print(f"Time index {tindex} / {len(params['tspan'])}")
        print(f"Computed I at t = {tindex}, I = {I_last}")

        I_saved.append(I_last)
        t_saved.append(t)

    # ------------------------
    # Forces
    # ------------------------
    Fa_x = applied_force(n, x_p, params["alpha"], V, params["L_x"])
    Fd_x, Fd_y = drag_force(n, x_v, y_v, params["eta"], params["Cd"])
    FI_x = interfacial_force(n, x_p, y_p, params["wI"], params["RI"], params["L_x"])
    Fp_x, Fp_y = pinning_force(
        n, x_p, y_p, params["w_pin"], params["x_pin"], params["y_pin"], params["R_pin"], params["L_x"]
    )
    Ft_x, Ft_y = temperature_fluctuations(
        n, eta, params["T_coeff"], T, rand_dirs_global_x, rand_dirs_global_y
    )
    Fr_x, Fr_y = residual_force(n, x_p, y_p, params["L_x"], params["L_y"])

    Fc_x = np.zeros(n)
    Fc_y = np.zeros(n)

    fin_array = finishing_array(x_p, params["L_x"], params["fin"])

    # Total forces
    forces_x = Fa_x + Fd_x + FI_x + Fp_x + Ft_x + Fr_x + Fc_x
    forces_y = Fd_y + Fp_y + Ft_y + Fr_y + Fc_y

    # ------------------------
    # Solve for dx/dt
    # ------------------------
    dxdt = np.zeros(4 * n + 1)
    dxdt[0:n] = x_v * fin_array
    dxdt[n:2 * n] = (forces_x / eta) * fin_array
    dxdt[2 * n:3 * n] = y_v * fin_array
    dxdt[3 * n:4 * n] = (forces_y / eta) * fin_array
    dxdt[4 * n] = (params["CT"] * params["Q"]) - params["k"] * (T - params["T_0"])

    return dxdt


# -----------------------------------------------------
# Subfunctions
# -----------------------------------------------------

def distances(x1, x2, y1, y2):
    dx = x1 - x2
    dy = y1 - y2
    d = np.sqrt(dx ** 2 + dy ** 2)
    return d, dx, dy


def applied_force(n, x_p, alpha, V, Lx):
    Fa_x = np.zeros(n)
    inside = (x_p >= 0) & (x_p <= Lx)
    Fa_x[inside] = alpha[inside].flatten() * V / (12 * Lx)
    return Fa_x


def drag_force(n, x_v, y_v, eta, Cd):
    Fd_x = -eta * Cd * x_v
    Fd_y = -eta * Cd * y_v
    return Fd_x, Fd_y


def interfacial_force(n, x_p, y_p, wI, RI, Lx):
    FI_x = np.zeros(n)
    inside = (x_p >= 0) & (x_p <= Lx)
    FI_x[inside] = -(2 * wI * RI ** 2) * (
        x_p[inside] * np.exp(-(x_p[inside] ** 2) / RI ** 2)
        + (x_p[inside] - Lx) * np.exp(-((x_p[inside] - Lx) ** 2) / RI ** 2)
    )
    return FI_x


def pinning_force(n, x_p, y_p, w_pin, x_pin, y_pin, R_pin, Lx):
    Fp_x = np.zeros(n)
    Fp_y = np.zeros(n)
    for i in range(len(w_pin)):
        d, dx, dy = distances(x_p, x_pin[i], y_p, y_pin[i])
        F = (2 * w_pin[i] * R_pin[i] ** 2) * np.exp(-(d ** 2) / R_pin[i] ** 2)
        Fp_x += F * dx
        Fp_y += F * dy
    return Fp_x, Fp_y


def temperature_fluctuations(n, eta, T_coeff, T, rand_x, rand_y):
    noise_scale = np.sqrt(eta * T_coeff * T)
    Ft_x = rand_x * noise_scale
    Ft_y = rand_y * noise_scale
    return Ft_x, Ft_y


def residual_force(n, x_p, y_p, Lx, Ly):
    cell_size = 5.0
    w_resid = 1.0
    Nx = max(1, int(np.ceil(Lx / cell_size)))
    Ny = max(1, int(np.ceil(Ly / cell_size)))

    y_shifted = y_p - np.min(y_p)
    ix = np.clip(np.floor(x_p / cell_size).astype(int), 0, Nx - 1)
    iy = np.clip(np.floor(y_shifted / cell_size).astype(int), 0, Ny - 1)

    count = np.zeros((Ny, Nx))
    for i in range(n):
        count[iy[i], ix[i]] += 1

    count[tuple((iy, ix))] = np.maximum(count[tuple((iy, ix))] - 1, 0)

    def shift(arr, dx, dy):
        return np.roll(np.roll(arr, dy, axis=0), dx, axis=1)

    count_left = np.pad(count[:, :Nx - 1], ((0, 0), (1, 0)), mode="edge")
    count_right = np.pad(count[:, 1:], ((0, 0), (0, 1)), mode="edge")
    count_up = np.pad(count[1:, :], ((0, 1), (0, 0)), mode="edge")
    count_down = np.pad(count[:-1, :], ((1, 0), (0, 0)), mode="edge")

    lin_idx = (iy, ix)
    rho_left = count_left[lin_idx]
    rho_right = count_right[lin_idx]
    rho_up = count_up[lin_idx]
    rho_down = count_down[lin_idx]

    Fr_x = w_resid * (rho_left - rho_right)
    Fr_y = w_resid * (rho_down - rho_up)
    return Fr_x, Fr_y


def finishing_array(x_p, L, fin):
    if fin == 1:
        return np.where(x_p >= 2 * L - 5e-6, 0, 1)
    else:
        return np.ones_like(x_p)


"""# ------------------------
# Global variables
# ------------------------
tindex = 0
V = 0.0

t_saved = []
I_saved = []
t_saved = []
I_last = 0.0

n = params['n']

rand_dirs_global_x = 2 * np.random.randint(0, 2, size=(n, 1)) - 1 # chooses +-1 with uniform dist.
rand_dirs_global_y = 2 * np.random.randint(0, 2, size=(n, 1)) - 1


eta = params['eta'] # viscosity

# ------------------------
# Main ODE RHS function
# ------------------------
def calculateForces(t, states, params):

    # Initialize globals if needed
    if I_saved is None or len(I_saved) == 0:
        I_saved = []
        t_saved = []
        I_last = 0.0

    if rand_dirs_global_x is None or rand_dirs_global_y is None:
        rand_dirs_global_x = 2 * np.random.randint(0, 2, params["n"]) - 1
        rand_dirs_global_y = 2 * np.random.randint(0, 2, params["n"]) - 1

    # ------------------------
    # Unpack state vector
    # ------------------------
    x_p = states[0:n]
    x_v = states[n:2 * n]
    y_p = states[2 * n:3 * n]
    y_v = states[3 * n:4 * n]
    T = states[4 * n] # temperature

    # ------------------------
    # Compute current periodically
    # ------------------------
    if tindex <= len(params["tspan"]) and t > params["tspan"][tindex]: # what is the point of this if statement?
        rand_dirs_global_x = 2 * np.random.randint(0, 2, n) - 1
        rand_dirs_global_y = 2 * np.random.randint(0, 2, n) - 1

        if tindex % 2 == 0:
            I_last = calculateCurrent(
                x_p, y_p,
                params["Lx"], V,
                params["lambda"], params["Rt"],
                params["steps"], params["num_e"]
            )

        tindex += 1

        if tindex > 1000 and V > 0:
            V = 0
            print("Voltage set to zero after 1000 iterations")

        print(f"Time index: {tindex} / {len(params['tspan'])}")
        print(f"Computed I at t = {tindex} | I = {I_last}")

        I_saved.append(I_last)
        t_saved.append(t)

    # ------------------------
    # Compute forces
    # ------------------------
    Fa_x = applied_force(n, x_p, params["alpha"], V, params["Lx"])
    Fd_x, Fd_y = drag_force(n, x_v, y_v, params["eta"], params["Cd"])
    FI_x = interfacial_force(n, x_p, y_p, params["wI"], params["RI"], params["Lx"])
    Fp_x, Fp_y = pinning_force(
        n, x_p, y_p, params["w_pin"],
        params["x_pin"], params["y_pin"],
        params["R_pin"], params["Lx"]
    )
    Ft_x, Ft_y = temperature_fluctuations(
        n, eta, params["T_coeff"], T,
        rand_dirs_global_x, rand_dirs_global_y
    )

    Fr_x, Fr_y = residual_force(n, x_p, y_p, params["Lx"], params["Ly"])
    Fc_x = np.zeros(n)
    Fc_y = np.zeros(n)

    fin_array = finishing_array(x_p, params["Lx"], params["fin"])

    forces_x = Fa_x + Fd_x + FI_x + Fp_x + Ft_x + Fr_x + Fc_x
    forces_y = Fd_y + Fp_y + Ft_y + Fr_y + Fc_y

    dxdt = np.zeros_like(states)

    # x dynamics
    dxdt[0:n] = x_v * fin_array
    dxdt[n:2 * n] = (forces_x / eta) * fin_array

    # y dynamics
    dxdt[2 * n:3 * n] = y_v * fin_array
    dxdt[3 * n:4 * n] = (forces_y / eta) * fin_array

    # temperature
    dxdt[4 * n] = params["CT"] * params["Q"] - params["k"] * (T - params["To"])

    return dxdt


# ------------------------
# Distance calculation
# ------------------------
def distances(x1, x2, y1, y2):
    assert np.all(np.isfinite(x1)) and np.all(np.isfinite(x2)), "x input contains non-finite values"
    assert np.all(np.isfinite(y1)) and np.all(np.isfinite(y2)), "y input contains non-finite values"

    dx = x1 - x2
    dy = y1 - y2
    d = np.sqrt(dx ** 2 + dy ** 2)
    return d, dx, dy


# ------------------------
# Applied Force
# ------------------------
def applied_force(n, x_p, alpha, V, Lx):
    Fa_x = np.zeros(n)
    inside = x_p < Lx
    Fa_x[inside] = alpha[inside] * V / (0.5 * Lx)
    if np.any(Fa_x != 0):
        assert V > 0, "Applied force should be zero if voltage is zero"
    return Fa_x


# ------------------------
# Drag Force
# ------------------------
def drag_force(n, x_v, y_v, eta, Cd):
    Fd_x = -eta * Cd * x_v
    Fd_y = -eta * Cd * y_v
    return Fd_x, Fd_y


# ------------------------
# Interfacial Force
# ------------------------
def interfacial_force(n, x_p, y_p, wI, RI, Lx):
    FI_x = np.zeros(n)
    inside = (x_p > 0) & (x_p < Lx)
    FI_x[inside] = -(2 * wI / RI ** 2) * (
        x_p[inside] * np.exp(-x_p[inside] ** 2 / RI ** 2)
        + (x_p[inside] - Lx) * np.exp(-(x_p[inside] - Lx) ** 2 / RI ** 2)
    )
    return FI_x


# ------------------------
# Pinning Force
# ------------------------
def pinning_force(n, x_p, y_p, w_pin, x_pin, y_pin, R_pin, Lx):
    Fp_x = np.zeros(n)
    Fp_y = np.zeros(n)

    for i in range(len(w_pin)):
        d, dx, dy = distances(x_p, x_pin[i], y_p, y_pin[i])
        F = (2 * w_pin[i] / (R_pin[i] ** 2)) * np.exp(-d ** 2 / (R_pin[i] ** 2))
        Fp_x += F * dx
        Fp_y += F * dy
    return Fp_x, Fp_y


# ------------------------
# Temperature Fluctuations
# ------------------------
def temperature_fluctuations(n, eta, T_coeff, T, rand_x, rand_y):
    noise_scale = np.sqrt(eta * T_coeff * T)
    Ft_x = rand_x * noise_scale
    Ft_y = rand_y * noise_scale
    return Ft_x, Ft_y


# ------------------------
# Residual Force
# ------------------------
def residual_force(n, x_p, y_p, Lx, Ly):
    cell_size = 5
    w_resid = 1

    Nx = max(1, int(np.ceil(Lx / cell_size)))
    Ny = max(1, int(np.ceil(Ly / cell_size)))

    y0 = np.min(y_p)
    y_shifted = y_p - y0

    ix = np.floor(x_p / cell_size).astype(int)
    iy = np.floor(y_shifted / cell_size).astype(int)

    ix = np.clip(ix, 0, Nx - 1)
    iy = np.clip(iy, 0, Ny - 1)

    count = np.zeros((Ny, Nx))
    for i in range(n):
        count[iy[i], ix[i]] += 1
    count[iy, ix] = np.maximum(count[iy, ix] - 1, 0)

    count_left = np.hstack((count[:, [0]], count[:, :-1]))
    count_right = np.hstack((count[:, 1:], count[:, [-1]]))
    count_up = np.vstack((count[1:, :], count[[-1], :]))
    count_down = np.vstack((count[[0], :], count[:-1, :]))

    lin_idx = iy * Nx + ix
    rho_left = count_left.flat[lin_idx]
    rho_right = count_right.flat[lin_idx]
    rho_up = count_up.flat[lin_idx]
    rho_down = count_down.flat[lin_idx]

    Fr_x = w_resid * (rho_left - rho_right)
    Fr_y = w_resid * (rho_down - rho_up)
    return Fr_x, Fr_y


# ------------------------
# Finishing Array
# ------------------------
def finishing_array(x_p, L, fin):
    if fin == 1:
        fin_array = (x_p < 2 * L - 5e-6).astype(float)
    else:
        fin_array = np.ones_like(x_p)
    return fin_array"""


    