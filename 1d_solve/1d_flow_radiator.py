import numpy as np
import matplotlib.pyplot as plt

def simulate_radiator_cooling(
    T_in,           # Inlet fluid temperature (K)
    Q_hot,          # Total heat to reject (W) - calculated load (39.75 W)
    mass_flow_rate, # Mass flow rate of fluid (kg/s)
    cp_fluid,       # Specific heat capacity of fluid (J/(kg*K))
    radiator_length, # Total length of the radiator flow path (m)
    radiator_width, # Radiating surface width (m)
    emissivity,     # Radiator surface emissivity (dimensionless)
    num_elements=50 # Number of discrete elements (for 1D discretization)
):
    """
    Simulates the 1D temperature profile of a fluid flowing through a radiator.
    
    The simulation calculates the required radiating area based on the heat load (Q_hot) 
    and then simulates the temperature drop along the length using that required area.
    
    Returns:
        tuple: (temperatures, distances, required_area)
    """
    
    # 1. Constants
    sigma = 5.670374e-8  # Stefan-Boltzmann constant (W/(m^2*K^4))
    T_sink = 4           # Deep space sink temperature (K) - negligible but included
    
    # 2. Calculate the REQUIRED total radiating area (A_req)
    # This calculation assumes an average temperature (T_avg) to find the size needed.
    # The actual temperature drop is then simulated over that area.
    
    # Estimate the final temperature assuming the radiator works:
    T_out_estimated = T_in - (Q_hot / (mass_flow_rate * cp_fluid))
    T_avg_design = (T_in + T_out_estimated) / 2
    
    # Radiative Heat Equation (simplified Q = sigma*eps*A*T^4 for deep space)
    # A = Q / (sigma * eps * (T_avg^4 - T_sink^4))
    required_area = Q_hot / (emissivity * sigma * (T_avg_design**4 - T_sink**4))
    
    # 3. Discretization Setup
    
    # Area per element (A_elem): assuming a single continuous plate (L x W)
    # The radiating area for an element is (L_element) * (Width)
    element_length = radiator_length / num_elements
    area_per_element = element_length * radiator_width
    
    # Check if the calculated area matches the geometry:
    # We must ensure the element area is sized according to the required area.
    # We will assume the geometry (L x W) is used to calculate the element area, 
    # but the total size is validated against A_req.
    
    # If the provided L x W geometry is the required area, use it:
    if np.isclose(radiator_length * radiator_width, required_area):
        pass # Geometry matches required area
    else:
        # Use the required area to define the element area for consistency
        # Assuming the total length 'radiator_length' is fixed, adjust the effective width
        effective_width = required_area / radiator_length
        area_per_element = required_area / num_elements
        
    
    # 4. Simulation Loop (Steady State)
    
    temperatures = [T_in]
    distances = [0.0]
    Q_rad_total = 0.0
    
    # Initial temperature is the inlet temp
    T_current = T_in
    
    for i in range(1, num_elements + 1):
        
        # Radiated power from the current element (Q_rad)
        # Assuming the element surface temperature is the fluid temperature (T_current)
        Q_rad_elem = emissivity * sigma * area_per_element * (T_current**4 - T_sink**4)
        
        # Check for convergence / stop condition
        # If Q_rad_elem is very small or temperature is too low, stop.
        if Q_rad_elem < 1e-6:
            Q_rad_elem = 0
        
        # Heat balance for the fluid element:
        # Q_rad_elem = dQ_fluid = mass_flow_rate * cp_fluid * dT
        
        # Calculate the temperature change (dT) across this element
        dT = Q_rad_elem / (mass_flow_rate * cp_fluid)
        
        # New temperature is lower (cooling)
        T_new = T_current - dT
        
        # Store results
        T_current = T_new
        temperatures.append(T_current)
        distances.append(i * element_length)
        Q_rad_total += Q_rad_elem

    return np.array(temperatures), np.array(distances), Q_rad_total, required_area

# --- USER INPUT PARAMETERS ---

# Parameters based on the user's calculations:
Q_hot_watt = 39.75  # Heat to reject (W)
mass_flow_rate_kg_s = 0.768 / 3600  # 0.768 kg/h converted to kg/s (~2.133e-4 kg/s)
cp_ethylene_glycol = 3000  # J/(kg*K) (Approximate value for EG/Water)
emissivity_radiator = 0.85 # High-emissivity coating (e.g., White Paint)
T_inlet_K = 298.0  # Inlet temperature from PTR hot end (K)

# Radiator Geometry (example, assuming a square radiator)
Radiator_L_m = 0.5   # 50 cm
Radiator_W_m = 0.5   # 50 cm
# Total geometric area is 0.25 m^2 (2500 cm^2)

# --- RUN SIMULATION ---

T_profile, X_distances, Q_total, A_required = simulate_radiator_cooling(
    T_in=T_inlet_K,
    Q_hot=Q_hot_watt,
    mass_flow_rate=mass_flow_rate_kg_s,
    cp_fluid=cp_ethylene_glycol,
    radiator_length=Radiator_L_m,
    radiator_width=Radiator_W_m,
    emissivity=emissivity_radiator
)

# --- OUTPUT AND PLOT RESULTS ---

# Calculate actual heat rejection
Q_actual = mass_flow_rate_kg_s * cp_ethylene_glycol * (T_profile[0] - T_profile[-1])

print(f"--- Radiator Simulation Results ---")
print(f"Heat Load (Q_hot): {Q_hot_watt:.2f} W")
print(f"Mass Flow Rate (m_dot): {mass_flow_rate_kg_s:.4e} kg/s")
print(f"Fluid Cp: {cp_ethylene_glycol} J/(kg*K)")
print("-" * 35)
print(f"Inlet Temperature (T_in): {T_profile[0]:.2f} K")
print(f"Outlet Temperature (T_out): {T_profile[-1]:.2f} K")
print(f"Total Temperature Drop (Delta T): {(T_profile[0] - T_profile[-1]):.2f} K")
print(f"Actual Heat Rejected: {Q_actual:.2f} W")
print(f"Actual Heat Rejected 2: {Q_total:.2f} W")

print(f"Required Design Area (A_req): {A_required:.4f} m^2 ({A_required * 10000:.1f} cm^2)")


# Plotting the temperature profile
plt.figure(figsize=(10, 6))
plt.plot(X_distances, T_profile - 273.15, marker='o', markersize=4, linestyle='-')
plt.title('Fluid Temperature Profile Along Radiator Length')
plt.xlabel('Distance Along Radiator (m)')
plt.ylabel('Fluid Temperature (°C)')
plt.grid(True)
plt.show()
