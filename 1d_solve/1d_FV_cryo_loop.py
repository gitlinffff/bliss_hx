# Solve 1d advection with radiation as heat source using Finite Volume method
import numpy as np
import matplotlib.pyplot as plt

def plotstate(f, x, T, N_ptr, savefig=False):
    dx = x[1]-x[0]
    L_ptr = N_ptr*dx

    plt.clf()
    plt.plot(x, T, 'o-', ms=1, linewidth=0.5, color='black', label='numerical')
    plt.axvline(x=L_ptr, color='red', linestyle='--')
    plt.text(x[0], T[0], f"{T[0]:.2f}", color='red', fontsize=16)
    plt.text(L_ptr+0.02, T[N_ptr-1], f"{T[N_ptr-1]:.2f}", color='red', fontsize=16)
    plt.text(x[-1]+0.02, T[-1], f"{T[-1]:.2f}", color='red', fontsize=16)
    plt.xlabel(r'position, $x$', fontsize=20)
    plt.ylabel(r'state, $T$', fontsize=20)
    plt.axis([0, 3, 150, 400])
    plt.figure(f.number);
    plt.tick_params(axis='both', labelsize=16)
    f.tight_layout(); 
    if savefig:
        plt.savefig(f"temperature.png", dpi=200, bbox_inches='tight', pad_inches=0.1)
    else:
        plt.draw(); plt.pause(0.04)


def Roe_flux(UL, UR, a):
    """
    Roe flux for linear advection equation.
    """
    a_hat = a  # Roe average speed (for linear advection, it's just a)
    return 0.5 * a * (UL + UR) - 0.5 * abs(a_hat) * (UR - UL)

def solve_FV(a, xc, T0, C_rad, C_ptr, T_space, T_ptr, N_ptr, walltime, CFL, flux_func):
    dx = xc[1]-xc[0]  # assume dx is same for all cells
    dt = CFL*dx/a; Nt = int(np.ceil(walltime/dt)); dt = walltime/Nt
    Ne = len(T0)  # number of elements
    

    T = T0.copy()
    R = T0.copy()  # Flux residual
    f = plt.figure(figsize=(10,5))
    for n in range(Nt):  # loop over time steps
        plotstate(f, xc, T, N_ptr)
        R *= 0. # zero out residual vector, over all cells
        for j in range(Ne-1):  # loop over interior interfaces
            TL = T[j]   # left state
            TR = T[j+1] # right state
            Fhat = flux_func(TL, TR, a)  # numerical flux
            R[j] += Fhat    # increment left cell residual
            R[j+1] -= Fhat  # and right cell residual
        
        # periodic boundary
        R[0] -= flux_func(T[-1], T[0], a)
        R[-1] += flux_func(T[-1], T[0], a)

        # add source term contribution to residual
        R[:N_ptr] -= C_ptr * (T_ptr - T[:N_ptr])     # convection at PTR
        R[N_ptr:] -= C_rad * (T[N_ptr:]**4 - T_space**4) # radiation at radiator
        
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

def get_Nu(Re, Pr):
    if Re < 2300:
        # Laminar flow, use Sieder-Tate correlation
        Nu = 3.66
    else:
        # error if flow is turbulent, as we are assuming laminar flow in this problem
        raise ValueError("Flow is turbulent, Nusselt number not defined for this regime.")
    return Nu



def build_problem():
    # =======================================
    # ====== Define problem parameters ======
    # =======================================
    # constants
    sigma = 5.67e-8  # Stefan-Boltzmann constant, W/(m^2*K^4)
    T_space = 240.15  # surrounding space temperature, K

    # fluid properties
    rho = 1000.0  # density of fluid, kg/m^3
    cp = 1500.0  # specific heat capacity of fluid, J/(kg*K)
    u = 0.03368  # advection velocity, m/s
    k = 0.12     # thermal conductivity of fluid, W/(m*K)
    mu = 23.9 * 1e-4   # dynamic viscosity of fluid, Pa*s

    # radiator channel properties
    d_rad = 0.012  # diameter of fluid channel, m
    wp_rad = np.pi*d_rad  # wet perimeter of fluid channel, m
    Ac_rad = 0.25 * np.pi * d_rad**2  # cross-sectional area of fluid channel, m^2
    L_rad = 1.0  # length of fluid channel, m
    emissivity = 0.8  # emissivity of radiating surface

    # PTR parameters
    d_ptr = 0.012  # diameter of channel at PTR, m
    wp_ptr = np.pi*d_ptr  # wet perimeter of channel at PTR, m
    Ac_ptr = 0.25 * np.pi * d_ptr**2  # cross-sectional area of channel at PTR, m^2
    L_ptr = 1.5  # length of channel at PTR, m
    T_ptr = 298.0  # assume constant temperature on inner surface, K

    # simulation parameters
    walltime = 1000  # total simulation time, s
    Ne_list = [256, 512, 1024, 2048, 4096]  # number of elements
    cfl = 0.5
    
    # ========================================
    # ====== Compute derived parameters ======
    # ========================================
    # determine convection heat transfer coefficient at PTR channel
    Re = (rho * u * d_ptr) / mu  # Reynolds number
    Pr = (cp * mu) / k  # Prandtl number
    Nu = get_Nu(Re, Pr)  # Nusselt number, function of Re and Pr
    h = (Nu * k) / d_ptr  # convection heat transfer coefficient, W/(m^2*K)

    # compute coefficients of heat source term from radiation and convection
    C_rad = -(emissivity * sigma * wp_rad) / (rho * cp * Ac_rad)*0.7
    C_ptr = (h * wp_ptr) / (rho * cp * Ac_ptr)

    # ========================================
    # =========== build grid and IC ==========
    # ========================================
    solutions = {}
    for Ne in Ne_list:

        L = L_rad + L_ptr
        N_ptr = int(np.floor(Ne * (L_ptr/L))) # set number of cells in PTR region
        x_grid = np.linspace(0., L, Ne+1)
        xc = 0.5*(x_grid[0:Ne] + x_grid[1:Ne+1]) # cell centers
        T_init = np.ones(Ne) * 250  # initial temperature field
        

        T_sol = solve_FV(u, xc, T_init, C_rad, C_ptr, T_space, T_ptr, N_ptr, walltime, cfl, Roe_flux)

        # compute heat transfer rates at PTR and radiator
        Q_rad = cp * rho * Ac_rad * u * (T_sol[-1] - T_sol[N_ptr-1])
        
        solutions[Ne] = {
            "x_grid": xc,
            "T_sol": T_sol,
            "Q_rad": Q_rad
        }

    # plot solution for all Ne in one figure
    f = plt.figure(figsize=(10,5))
    for Ne, sol in solutions.items():
        xc = sol["x_grid"]
        T_sol = sol["T_sol"]

        plt.plot(xc, T_sol, label=f"Ne={Ne}")
    plt.xlabel(r'position, $x$', fontsize=20)
    plt.ylabel(r'temperature, $T$', fontsize=20)
    plt.savefig("Q1_T_solutions.png", dpi=200, bbox_inches='tight', pad_inches=0.1)




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
