# Thermodynamics properties of fluids
fluids_dict = {
    
    # imaginary fluid
    'imag_fluid': {
        'cp': 1500,              # J/(kg K)
        'rho_f': 1000,         # kg/m3       # NIST https://webbook.nist.gov/cgi/fluid.cgi?ID=C75694&Action=Page
        'mu': 23.9 * 1e-4      # (Pa s)      # NIST https://webbook.nist.gov/cgi/fluid.cgi?ID=C75694&Action=Page
    },
    
    # CFC 11  (Reference temperature: 242.5 K)
    'CFC11': {
        'cp': 840.38,            # J/(kg K)    # Osborne, The Heat Capacity, Entropy, Heats of Fusion 
                                         # and Vaporization and Vapor Pressure of Fluorotrichloromethane
        'rho_f': 1601.9,         # kg/m3       # NIST https://webbook.nist.gov/cgi/fluid.cgi?ID=C75694&Action=Page
        'mu': 7.8308 * 1e-4      # (Pa s)      # NIST https://webbook.nist.gov/cgi/fluid.cgi?ID=C75694&Action=Page
    },
    
    # Therminol®59  (Reference temperature: 223.15 K)
    'Therminol59': {
        'cp': 1460,              # J/(kg K)    # https://domxoloda.ru/oils/docs/HTF-59.PDF
        'rho_f': 1025,           # kg/m3       # https://domxoloda.ru/oils/docs/HTF-59.PDF
        'mu': 25043.1e-4     # (Pa s)      # https://domxoloda.ru/oils/docs/HTF-59.PDF
    },
    
    # Therminol®LT  (Reference temperature: 223.15 K)
    'TherminolLT': {
        'cp': 1530,              # J/(kg K)    # https://www.sintelub.com/wp-content/uploads/PDS/29.Therminol_LT.pdf
        'rho_f': 921,            # kg/m3       # https://www.sintelub.com/wp-content/uploads/PDS/29.Therminol_LT.pdf
        'mu': 39.9e-4,           # (Pa s)      # https://www.sintelub.com/wp-content/uploads/PDS/29.Therminol_LT.pdf
        'nu': 4.33,              # mm2/s       # https://www.sintelub.com/wp-content/uploads/PDS/29.Therminol_LT.pdf
        'k': 0.1384              # W/(m K)     # https://www.sintelub.com/wp-content/uploads/PDS/29.Therminol_LT.pdf
    },
    
    # Syltherm®XLT  (Reference temperature: 242.5 K)
    'SylthermXLT': {
        'cp': 1666,              # J/(kg K)    # https://www.npl.washington.edu/TRIMS/sites/sand.npl.washington.edu.TRIMS/files/manuals-documentation/syltherm-xlt-technical-data-sheet.pdf
        'rho_f': 907.2,          # kg/m3       # https://www.npl.washington.edu/TRIMS/sites/sand.npl.washington.edu.TRIMS/files/manuals-documentation/syltherm-xlt-technical-data-sheet.pdf
        'mu': 130 * 1e-4,        # (Pa s)      # https://www.npl.washington.edu/TRIMS/sites/sand.npl.washington.edu.TRIMS/files/manuals-documentation/syltherm-xlt-technical-data-sheet.pdf
        'k': 0.1212              # W/(m K)     # https://www.npl.washington.edu/TRIMS/sites/sand.npl.washington.edu.TRIMS/files/manuals-documentation/syltherm-xlt-technical-data-sheet.pdf
    },
    
    # Therminol®VLT  (Reference temperature: 222.15 K)
    'TherminolVLT': {
        'cp': 1650,              # J/(kg K)    # https://www.sintelub.com/wp-content/uploads/PDS/30.Therminol_VLT.pdf
        'rho_f': 809,            # kg/m3       # https://www.sintelub.com/wp-content/uploads/PDS/30.Therminol_VLT.pdf
        'mu': 23.9 * 1e-4,        # (Pa s)      # https://www.sintelub.com/wp-content/uploads/PDS/30.Therminol_VLT.pdf
        'k': 0.1181              # W/(m K)     # https://www.sintelub.com/wp-content/uploads/PDS/30.Therminol_VLT.pdf
    },
    
    # R-236fa  (Reference temperature: 242.5 K)
    'R_236fa': {
        'cp': 1166.1,            # J/(kg K)    # NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C690391&Units=SI
        'rho_f': 1533.9,         # kg/m3       # NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C690391&Units=SI
        'mu': 6.1904e-4          # (Pa s)      # NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C690391&Units=SI
    },

    # NH3  (Reference temperature: 223.15 K  pressure: 0.101MPa)
    'NH3_1e-1MPa': {
        'cp':    4404.,          # J/(kg K)    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        'rho_f': 702.,           # kg/m3       # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        'mu':    3.2532e-4,      # (Pa s)      # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        'k':     0.59694         # (W/m K)     # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
    },

    # NH3  (Reference temperature: 223.15 K  pressure: 1MPa)
    'NH3_1MPa': {
        'cp':    4401.,          # J/(kg K)    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        'rho_f': 702.,           # kg/m3       # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        'mu':    3.2664e-4,      # (Pa s)      # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
        'k':     0.59828         # (W/m K)     # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417
    },

    # Ethylene Glycol based Water Solutions (223.15 K, 50%.Vol)
    'Ethylene_Glycol': {
        'cp': 3040,            # J/(kg K)    # https://www.engineeringtoolbox.com/ethylene-glycol-d_146.html
        'rho_f': 1127,         # kg/m3       # https://www.engineeringtoolbox.com/ethylene-glycol-d_146.html
        'mu': 644e-4           # (Pa s)      # https://www.mokon.com/products/fluids/glycol-solutions/pdf/Ethylene-Glycol-Technical-Data.pdf
    },
}
