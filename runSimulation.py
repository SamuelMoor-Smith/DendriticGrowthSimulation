import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import PillowWriter, FuncAnimation
from scipy.integrate import solve_ivp
import imageio

from calculateForces import calculate_forces


def run_simulation(params):
    """
    Translated from MATLAB runSimulation.m
    Runs the full dendritic growth simulation with visualization and GIF output.
    """

    # Initialize globals (same as MATLAB)
    from calculateForces import I_saved, t_saved
    I_saved.clear()
    t_saved.clear()

    np.random.seed(1)

    # ------------------------
    # Unpack parameters
    # ------------------------
    n = params["n"]
    To = params["T_0"]
    tspan = params["tspan"]

    Lx = params["L_x"]
    Ly = params["L_y"]
    electrode_width = params["electrode_width"]
    electrode_height = params["electrode_height"]

    # ------------------------
    # Randomize initial particle positions
    # ------------------------
    xo = 10 * np.random.rand(n)
    yo = electrode_height * np.random.rand(n) - 0.5 * electrode_height

    # ------------------------
    # Zero all initial particle velocities
    # ------------------------
    vox = np.zeros(n)
    voy = np.zeros(n)

    # ------------------------
    # Initial state vector
    # ------------------------
    initial = np.concatenate([xo, vox, yo, voy, [To]])

    # ------------------------
    # Run ODE simulation
    # ------------------------
    print("Starting solver...")
    from calculateForces import tindex, V
    tindex = 1
    V = params["V"]

    sol = solve_ivp(
        fun=lambda t, y: calculate_forces(t, y, params),
        t_span=(tspan[0], tspan[-1]),
        y0=initial,
        t_eval=tspan,
        method="RK45",
        vectorized=False,
        rtol=1e-6,
        atol=1e-9,
    )

    t = sol.t
    states = sol.y.T

    X_pos = states[:, 0:n]
    Y_pos = states[:, 2 * n:3 * n]

    # ------------------------
    # Plot setup
    # ------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    plt.tight_layout(pad=3)

    # ---- Top plot: particle evolution ----
    ax1.set_xlim([-electrode_width - 5, Lx + electrode_width + 5])
    ax1.set_ylim([-electrode_height / 2 - 5, electrode_height / 2 + 5])
    ax1.set_xlabel("x-position")
    ax1.set_ylabel("y-position")
    ax1.set_title("Dendritic Growth Simulation")

    # Draw electrodes
    ax1.fill_betweenx(
        [-electrode_height / 2, electrode_height / 2],
        -electrode_width,
        0,
        color=(1.0, 0.8431, 0.0),
        alpha=0.3,
    )
    ax1.fill_betweenx(
        [-electrode_height / 2, electrode_height / 2],
        Lx,
        Lx + electrode_width,
        color=(1.0, 0.8431, 0.0),
        alpha=0.3,
    )

    # ---- Pinning potential background ----
    x = np.linspace(0, Lx, 100)
    y = np.linspace(-Ly / 2, Ly / 2, 100)
    X, Y = np.meshgrid(x, y)
    U = np.zeros_like(X)
    for k in range(len(params["w_pin"])):
        dx = X - params["x_pin"][k]
        dy = Y - params["y_pin"][k]
        U += params["w_pin"][k] * np.exp(-(dx**2 + dy**2) / (params["R_pin"][k] ** 2))

    pcm = ax1.pcolormesh(X, Y, U, cmap="pink", shading="auto", alpha=0.3)
    cb = plt.colorbar(pcm, ax=ax1)
    cb.set_label("Pinning Potential U")

    # ---- Initialize particle positions ----
    particles, = ax1.plot(X_pos[0, :], Y_pos[0, :], "b.", markersize=4)

    # ---- Bottom plot: current vs time ----
    ax2.set_xlim([t[0], t[-1]])
    ax2.set_ylim([0, params["num_e"]])
    ax2.set_xlabel("time")
    ax2.set_ylabel("Current I")
    ax2.grid(True)
    current_line, = ax2.plot([], [], "r-", lw=1.5)

    # ------------------------
    # Animation update function
    # ------------------------
    def update(frame):
        particles.set_data(X_pos[frame, :], Y_pos[frame, :])

        from calculateForces import I_saved, t_saved
        if len(t_saved) > 0:
            mask = np.array(t_saved) <= t[frame]
            if np.any(mask):
                current_line.set_data(np.array(t_saved)[mask], np.array(I_saved)[mask])

        ax1.set_title(f"t = {t[frame]:.2f}")
        return particles, current_line

    # ------------------------
    # Create and save animation
    # ------------------------
    anim = FuncAnimation(fig, update, frames=len(t), interval=100, blit=False)

    filename = f"dendrite_growth_simulation.gif"
    print(f"Saving GIF to {filename}...")
    writer = PillowWriter(fps=10)
    anim.save(filename, writer=writer)

    plt.close(fig)
    print("Simulation complete. GIF saved.")

    return t, states
