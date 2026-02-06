import numpy as np
import matplotlib.pyplot as plt

def Nu(Re, Pr):
    # laminar flow inside tubes


    # Dittus-Boelter equation for turbulent flow inside tubes
    return 0.023 * Re**0.8 * Pr**0.4

def simulate_radiator_cooling(
    T,           # fluid temperature (K)
    mt, # Mass flow rate of fluid (kg/s)
    Q_in,
    cp_fluid,       # Specific heat capacity of fluid (J/(kg*K))
    x_grids, # Total length of the radiator flow path (m)
    radiator_width, # Radiating surface width (m)
    emissivity,     # Radiator surface emissivity (dimensionless)
    num_elements=50 # Number of discrete elements (for 1D discretization)
):
    
    sigma = 5.670374e-8  # Stefan-Boltzmann constant (W/(m^2*K^4))
    T_sink = 4           # Deep space sink temperature (K) - negligible but included
    
    tolerance = 1e-5

    # The radiating area for an element is (L_element) * (Width)
    dx = x_grids[1] - x_grids[0]
    area_per_element = dx * radiator_width
   
    f = plt.figure(figsize=(10,5))

    iteration_or_not = True
    while iteration_or_not:
        plotstate(f, x_grids, T)

        T_elem_ave = 0.5*(T[1:] + T[:-1])
        Q_rad_elem = emissivity * sigma * area_per_element * (T_elem_ave**4 - T_sink**4)
        dT = Q_rad_elem / (mt * cp_fluid)
        T[1:] = T[:-1] - dT
        
        # calculate current total radiation flux
        Q_rad_total = np.sum(Q_rad_elem)
        
        totQ_res = np.abs(Q_rad_total - Q_in) 
        print(f"{T[0]:.2f}   {T[-1]:.2f}   {totQ_res:.4f}", flush=True)
        if totQ_res <= tolerance: 
            iteration_or_not = False
        else:
            # update inlet temperature due to input Q
            T[0] = T[-1] + Q_in / (mt * cp_fluid)

    plt.pause(3)
    #plt.close(f)

    return T, Q_rad_total

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


# --- USER INPUT PARAMETERS ---
def problem():
    # Parameters based on the user's calculations:
    N = 50   # number of nodes
    L = 1    # length of radiator (m)
    width = 0.35  # width of radiator (m)

    Q_hot_watt = 40  # Heat to reject (W)
    mass_flow_rate = 0.768 / 3600 # kg/s
    cp = 1500  # J/(kg*K) (Approximate value for EG/Water)
    emissivity = 0.5 # emissivity of the coating
    T = np.ones(N) * 300.0  # initial temperature
    x = np.linspace(0, L, N)

    simulate_radiator_cooling(
        T,           # Inlet fluid temperature (K)
        mass_flow_rate, # Mass flow rate of fluid (kg/s)
        Q_hot_watt,
        cp,       # Specific heat capacity of fluid (J/(kg*K))
        x, #
        width, # Radiating surface width (m)
        emissivity
    )


if __name__ == "__main__":
    problem()


