# Solve 1d advection with radiation as heat source using Finite Difference method
import numpy as np
import matplotlib.pyplot as plt

def plotstate(f, x, T):
    plt.clf()
    plt.plot(x, T, 'o-', linewidth=2, color='black', label='numerical')
    plt.xlabel(r'position, $x$', fontsize=20)
    plt.ylabel(r'state, $T$', fontsize=20)
    plt.axis([0, 1, 150, 400])
    plt.figure(f.number);
    plt.tick_params(axis='both', labelsize=16)
    f.tight_layout(); 
    plt.draw(); plt.pause(0.04)

def FOU(x_grids, dx, dt, T, cfl, n_iter, C, T_space, f):
    """
    First-Order Upwind
    """
    for n in range(n_iter):
        T_old = T.copy()

        # source term
        source = C * (T_old[1:]**4 - T_space**4)

        T[1:] = T_old[1:] - cfl * (T_old[1:] - T_old[:-1]) + dt * source
        # T[0] stays constant (boundary condition)

        plotstate(f, x_grids, T)

    return T

def solve(a, L, N, C, T_in, T_space, walltime):
    """
    N is number of intervals.
    """
    cfl = 0.5
    
    dx = L/N
    dt = cfl * dx / a

    # adjust dt and cfl
    n_iter = int(np.ceil(walltime/dt))
    dt = walltime/n_iter
    cfl = a * dt / dx
    print(f"N = {N}   dt = {dt}   cfl = {cfl}")

    x_grid = np.linspace(0., L, N+1)
    T_init = np.ones(N+1) * T_in  # initial temperature field
    T = T_init.copy()

    fig = plt.figure(figsize=(10,5))


    # solutions at t=T using three different methods
    sol_FOU = FOU(x_grid, dx, dt, T, cfl, n_iter, C, T_space, fig)


def build_problem():
    # constants
    sigma = 5.67e-8  # Stefan-Boltzmann constant, W/(m^2*K^4)
    T_space = 4.0  # surrounding space temperature, K

    # inlet temperature
    T_inlet = 298.0  # inlet temperature of fluid, K

    # fluid properties
    rho = 1000.0  # density of fluid, kg/m^3
    cp = 1500.0  # specific heat capacity of fluid, J/(kg*K)
    u = 0.03368  # advection velocity, m/s
    
    # radiator channel properties
    wet_p = np.pi*0.012  # wet perimeter of fluid channel, m
    A = 0.25 * np.pi * 0.012**2  # cross-sectional area of fluid channel, m^2
    L = 1.0  # length of fluid channel, m
    emissivity = 0.8  # emissivity of radiating surface

    # simulation parameters
    walltime = 500  # total simulation time, s
    N = 64  # number of spatial intervals
    
    # compute radiation coefficient
    C = -(emissivity * sigma * wet_p) / (rho * cp * A)

    solve(u, L, N, C, T_inlet, T_space, walltime)


def plot():
    fig, ax = plt.subplots()
    
    ax.plot(x_grid, sol_FOU, label="FOU")
    ax.plot(x_grid, sol_LW,  label="LW")
    ax.plot(x_grid, sol_BW,  label="BW")
    ax.plot(x_grid, sol_exact, '--', label="exact solution")
    
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.set_title(fr"Solutions at  t={T} for three methods       $N_x$={N}")
    ax.legend()
    plt.savefig("Q1_uA_100.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    plt.show()


if __name__ == "__main__":
    build_problem()
