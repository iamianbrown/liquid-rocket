import numpy as np
import matplotlib.pyplot as plt
import TCA
import Chamber_Size_Test


class Injector:
    # -------------------------------------------Initial Data-------------------------------------------------#

    '''
    _r designates radial propellant
    _z designates axial propellant
    '''

    def __init__(self, TCAobj, CHAMBERobj, injector_props):  
        self.TCAobj = TCAobj
        self.Chamberobj = CHAMBERobj
        self.injector_props = injector_props

        
        
      
        
    # def pintleParams(cea_data, gap_size, d_o, P_c):
    #     # cea_data is a list containing all relevant CEARUN parameters (T, OF, V_exit)
    #     # calculate the following parameters
    #     return theta, delta_p, BF, mdot_t

    P_c = 300 * 6894.76  # CHAMBER PRESSURE (Pascals)
    mdot_t = 1.989416305  # TOTAL MASS FLOW RATE (kg/s)
    OF = 3  # OXIDIZER TO FUEL RATIO
    d_c = .13  # CHAMBER DIAMETER (m)
    rho_r = 1141  # RADIAL PROPELLANT DENSITY (kg/m^3)
    rho_z = 810  # AXIAL PROPELLANT DENSITY (kg/m^3)
    d1 = 1.3  # ORIFICE ROW 1 DIAMETER (mm)
    d2 = 1.3  # ORIFICE ROW 2 DIAMETER (mm)
    C_d = 0.75  # DISCHARGE COEFFICIENT
    delta_P = 0.25 * P_c  # FUEL PRESSURE DROP
    delta_P_o = 0.15 * P_c  # OXIDIZER PRESSURE DROP

    # Calculate radial and axial propellant mass flow rates
    mdot_r = mdot_t / (1 + OF) * OF
    mdot_z = mdot_t / (1 + OF)

    # Calculate pintle diameter and radius using chamber-to-pintle ratio
    d_p = d_c / 4
    r_p = d_p / 2

    # Gap tolerance as provided by machine shop
    gap_tolerance = 0.05 / 1000

    # Constants needed for LMR and Spray Angle relation
    alpha = 0.7
    beta = 2.0

    # ---------------------------------------ORIFICE SIZING-------------------------------------------#

    # Orifice propellant discharge area

    '''
    _m designates value in meters
    _mm designates value in millimeters
    '''

    A_r_m = mdot_r / (C_d * np.sqrt(2 * rho_r * delta_P_o))  # m^2
    A_r_mm = A_r_m * 1000000  # mm^2

    a1 = np.pi * ((d1 * 0.5) ** 2)  # area of each orifice in row 1
    a2 = np.pi * ((d2 * 0.5) ** 2)  # area of "" row 2

    orifice_pairs = a1 + a2  # area per pair of orifices

    n_orifice_pairs = round(A_r_mm / orifice_pairs)

    percent_error = 100 * (n_orifice_pairs - (A_r_mm / orifice_pairs)) / (A_r_mm / orifice_pairs)

    BF = (n_orifice_pairs * (d1 + d2)) / (2 * np.pi * ((d_p * 1000) * 0.5))

    # -------------------------------------ANNULAR GAP SIZING-----------------------------------------#

    A_z_m = mdot_z / (C_d * np.sqrt(2 * rho_z * delta_P))
    A_z_mm = A_z_m * 1000000

    gap = np.sqrt((A_z_mm / np.pi) + ((r_p * 1000) ** 2)) - (r_p * 1000)
    gap_low = gap - (gap_tolerance * 1000)
    gap_high = gap + (gap_tolerance * 1000)

    # ---------------------------------------CALCULATE LMR AND THETA---------------------------------------#

    A_lr_c = np.pi * ((d1 / 2000) ** 2 + (d2 / 2000) ** 2)  # cross-sectional area for one orifice pair
    A_lz_c = (gap / 1000) * (
                d1 / 1000 + d2 / 1000)  # cross-sectional area of annular stream that impinges on the orifice pair

    U_r_c = mdot_r / (rho_r * A_r_m)
    U_z_c = mdot_z / (rho_z * A_z_m)

    LMR_c = (rho_r * (U_r_c ** 2) * A_lr_c) / (rho_z * (U_z_c ** 2) * A_lz_c)
    theta_c = alpha * np.arctan(beta * LMR_c) * (180 / np.pi) + 20 #theta_c = converging angle of nozzle

    # -------------------------------------------SETTING UP PLOTS-------------------------------------------------#

    # Array of gap sizes between 0.05 mm and 0.8 mm
    a = np.linspace(0.05 / 1000, 0.8 / 1000, 100000)

    A_lr = np.pi * ((d1 / 2 / 1000) ** 2)
    A_lz = (gap / 1000) * (d1 / 1000)
    U_r = mdot_r / (rho_r * A_r_m)
    U_z = mdot_z / (rho_z * np.pi * ((r_p + a) ** 2 - (r_p) ** 2))

    LMR = (rho_r * (U_r ** 2) * A_lr) / (rho_z * (U_z ** 2) * A_lz)
    theta = (alpha * np.arctan(beta * LMR) * (180 / np.pi)) + 20
