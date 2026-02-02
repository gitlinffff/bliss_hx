import numpy as np
import matplotlib.pyplot as plt
from fluid_property import fluids_dict

def calc_pump_power(fluid, Q, delT, L, D, N_elbow):
    """
    calculate fluid power delivered. Pump efficiency not considered.
    Inputs:
        fluid : dictionary of fluid properties
        Q :     Heat load from PTR (W)
        delT :  temperature change of fluid (K)
        L :     total pipe length (m)
        D :     Inner Pipe Diameter (m) 
        N_elbow : number of 90° elbows
    """
    # fluid parameters
    cp    = fluid['cp']     # heat capacity  J/(kg K)
    rho_f = fluid['rho_f']  # fluid density  kg/m3
    mu    = fluid['mu']     # dynamic viscosity (Pa s)
    
    # flow calculations
    mt = Q / (cp*delT)       # mass flow rate
    Vt = mt / rho_f          # volumetric flow rate
    Ac = 0.25*np.pi*D*D      # pipe inner cross section area
    v = mt / (rho_f*Ac)      # flow speed
    Re = rho_f * v * D / mu  # Reynolds Number

    # find friction factor
    if Re < 2300:
        print(f"laminar, Re = {Re}")
        f = 64/Re
    elif Re >= 4000:
        print(f"Turbulent, Re = {Re}")
        f = (0.79*np.log(Re)-1.64)**(-2)
    else:
        print("undefined Re.")
    
    # =========================
    # ===== Pressure Loss =====
    # =========================
    # pressure loss from friction in the straight pipe 
    delP_pipe = f * (L / D) * (0.5*rho_f*v*v)
    # pressure loss from components
    L_D_elbow = 30 # equivalent L/D of 90° elbows  /*https://tameson.com/pages/pressure-drop*/
    delP_comp = f * (N_elbow * L_D_elbow) * (0.5*rho_f*v*v) # (gamma * 0.5*rho_f*v*v, gamma is resistance coefficient by test or vendor specification)
    # total pressure loss
    delP_tot = delP_pipe + delP_comp
    
    # power needed to be delivered to fluid
    P_f = delP_tot * Vt
    
    # fluid mass
    mass_f = Ac * L * rho_f
    
    return Vt, mt, delP_comp, delP_tot, P_f, mass_f

def main():
    # system parameters
    fluid = fluids_dict['Ethylene_Glycol']
    Q = 40.                         # Heat load from PTR (W)
    delT = 25.                      # temperature change of fluid (K)
    L = np.linspace(2., 10., 3)     # total pipe length (m)
    D = 0.00635                     # Inner Pipe Diameter (m)
    N_elbow = np.linspace(2, 15, 3)  # number of 90° elbows
    #N_elbow = 10                    # number of 90° elbows

    # Calculate pressure drops over the meshgrid
    L_mesh, N_elbow_mesh = np.meshgrid(L, N_elbow)
    Vt, mt, delP_comp, delP_tot, P_f, mass_f = calc_pump_power(fluid, Q, delT, L_mesh,
                                                               D, N_elbow_mesh)
    
    # Create the filled contour plot
    plt.figure(figsize=(8, 6))

    # Levels can be adjusted or removed to let matplotlib decide the scale automatically
    cp = plt.contourf(L_mesh, N_elbow_mesh, delP_tot, levels=20, cmap='viridis')

    # Add a colorbar to show the pressure scale
    cbar = plt.colorbar(cp)
    cbar.set_label('Total Pressure Drop $\\Delta P_{tot}$ (Pa)')

    # Labeling the axes
    plt.title('Pressure Drop vs. Pipe Length and Number of Elbows')
    plt.xlabel('Total Pipe Length $L$ (m)')
    plt.ylabel('Number of 90° Elbows $N_{elbow}$')

    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

if __name__ == "__main__":
    main()
