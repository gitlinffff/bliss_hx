# Solve 1d advection with radiation as heat source using Finite Volume method
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



def Roe_flux(UL, UR, a):
    """
    Roe flux for linear advection equation.
    """
    a_hat = a  # Roe average speed (for linear advection, it's just a)
    return 0.5 * a * (UL + UR) - 0.5 * abs(a_hat) * (UR - UL)

def solve_FV(a, xc, T0, C, T_space, T_inlet, walltime, CFL, flux_func):
    dx = xc[1]-xc[0]  # assume dx is same for all cells
    dt = CFL*dx/a; Nt = int(np.ceil(walltime/dt)); dt = walltime/Nt
    Ne = len(T0)  # number of elements/cells
    

    T = T0.copy()
    R = T0.copy()  # Flux residual
    f = plt.figure(figsize=(10,5))
    for n in range(Nt):  # loop over time steps
        plotstate(f, xc, T)
        R *= 0. # zero out residual vector, over all cells
        for j in range(Ne-1):  # loop over interior interfaces
            TL = T[j]   # left state
            TR = T[j+1] # right state
            Fhat = flux_func(TL, TR, a)  # numerical flux
            R[j] += Fhat    # increment left cell residual
            R[j+1] -= Fhat  # and right cell residual
        
        # add boundary fluxes to flux residual
        R[0] -= flux_func(T_inlet, T[0], a)
        R[-1] += flux_func(T[-1], T[-1], a)

        # add source term contribution to residual
        R -= C * (T**4 - T_space**4)
        
        # for nonlinear problems, set dt here using wave speeds calculated from the flux evaluations
        T -= dt/dx * R  # update state to next time step

        # print total residual for monitoring convergence
        tot_residual = np.sum(np.abs(R))
        print(f"n = {n}, T = {(n+1)*dt}s, total residual = {tot_residual:.3e}")
        if tot_residual < 1e-6:
            print(f"Converged at time step {n}, T = {(n+1)*dt}s, total residual = {tot_residual:.3e}")
            break

    plt.pause(3)
    plt.close(f)
    return T


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
    Ne = 64  # number of elements
    cfl = 0.5
    
    # compute radiation coefficient
    C = -(emissivity * sigma * wet_p) / (rho * cp * A)

    
    x_grid = np.linspace(0., L, Ne+1)
    xc = 0.5*(x_grid[0:Ne] + x_grid[1:Ne+1]) # cell centers
    T_init = np.ones(Ne) * T_inlet  # initial temperature field
    

    solve_FV(u, xc, T_init, C, T_space, T_inlet, walltime, cfl, Roe_flux)


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
