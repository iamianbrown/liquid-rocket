import numpy as np
import matplotlib.pyplot as plt
from TCA import *
from Chamber_Size_Test import *


class Injector:
    # -------------------------------------------Initial Data-------------------------------------------------#

    '''
    _r designates radial propellant
    _z designates axial propellant
    '''

    def __init__(self, TCAobj, CHAMBERobj, injector_props):  
        
        self.TCAobj = TCAobj
        self.CHAMBERobj = CHAMBERobj
        self.injector_props = injector_props #user will pass this in as a dictionary with the values for rho_r, rho_z, d1, d2, C_d, delta_P, delta_P_o

        PINTLEvalues = self.pintle_system_calc()
        INJvalues = self.sizingcalc()
        self.injmerge = {
            'mdot_r': PINTLEvalues['mdot_r'],
            'mdot_z': PINTLEvalues['mdot_z'],
            'D_p': PINTLEvalues['D_p'],
            'R_p': PINTLEvalues['R_p'],
            'L_pintle': PINTLEvalues['L_pintle'],
            'A_r_mm': INJvalues['A_r_mm'],
            'A_z_mm': INJvalues['A_z_mm'],
            'n_oriface_pairs': INJvalues['n_oriface_pairs'],
            'BF': INJvalues['BF'],
            'gap': INJvalues['gap'],
            'LMR': INJvalues['LMR'],
            'theta_c': INJvalues['theta_c']
        }
        self.injector_props.update(self.injmerge)

    def pintle_system_calc(self):
        
        mdot = self.TCAobj.mdot
        OF_Ratio = self.TCAobj.OF_ratio
        R_c = self.CHAMBERobj.geometric_props['R_c']
        L_c = self.CHAMBERobj.geometric_props['L_c']
        
        # Calculate radial and axial propellant mass flow rates
        mdot_r = mdot / (1+OF_Ratio) * OF_Ratio 
        mdot_z = mdot / (1+OF_Ratio) 

        #Calculate pintle diameter and radius using chamber-to-pintle ratio
        D_p = R_c / 4 #Pintle diameter
        R_p = D_p / 2 #Pintle radius
        L_pintle = L_c/3  # Pintle length

        pintle_system_dict = {
            'mdot_r': mdot_r,
            'mdot_z': mdot_z,
            'D_p': D_p,
            'R_p': R_p,
            'L_pintle': L_pintle,
        }
        return(pintle_system_dict)

    def sizingcalc(self):
        
        PINTLEvalues = self.pintle_system_calc()
        mdot_r = PINTLEvalues['mdot_r']
        mdot_z = PINTLEvalues['mdot_z']
        D_p = PINTLEvalues['D_p']
        R_p = PINTLEvalues['R_p']
        gap_tolerance = 0.05 / 1000 # Gap tolerance as provided by machine shop
        alpha = 0.7 # Constant needed for LMR and Spray Angle relation
        beta = 2.0 # Constant needed for LMR and Spray Angle relation

        # ---------------------------------------ORIFICE SIZING-------------------------------------------#

        # Orifice propellant discharge area

        '''
        _m designates value in meters
        _mm designates value in millimeters
        '''

        A_r_m = mdot_r / (self.injector_props['C_d'] * np.sqrt(2 * self.injector_props['rho_r'] * self.injector_props['delta_P_o']))  # m^2
        A_r_mm = A_r_m * 1000000  # mm^2

        a1 = np.pi * ((self.injector_props['d1'] * 0.5) ** 2)  # area of each orifice in row 1
        a2 = np.pi * ((self.injector_props['d2'] * 0.5) ** 2)  # area of "" row 2

        orifice_pairs = a1 + a2  # area per pair of orifices

        n_orifice_pairs = round(A_r_mm / orifice_pairs)

        percent_error = 100 * (n_orifice_pairs - (A_r_mm / orifice_pairs)) / (A_r_mm / orifice_pairs)

        BF = (n_orifice_pairs * (self.injector_props['d1'] + self.injector_props['d2'])) / (2 * np.pi * ((D_p * 1000) * 0.5))

        # -------------------------------------ANNULAR GAP SIZING-----------------------------------------#

        A_z_m = mdot_z / (self.injector_props['C_d'] * np.sqrt(2 * self.injector_props['rho_z'] * self.injector_props['delta_P']))
        A_z_mm = A_z_m * 1000000

        gap = np.sqrt((A_z_mm / np.pi) + ((R_p * 1000) ** 2)) - (R_p * 1000)
        gap_low = gap - (gap_tolerance * 1000)
        gap_high = gap + (gap_tolerance * 1000)

        # ---------------------------------------CALCULATE LMR AND THETA---------------------------------------#

        A_lr_c = np.pi * ((self.injector_props['d1'] / 2000) ** 2 + (self.injector_props['d2'] / 2000) ** 2)  # cross-sectional area for one orifice pair
        A_lz_c = (gap / 1000) * (
                    self.injector_props['d1'] / 1000 + self.injector_props['d2'] / 1000)  # cross-sectional area of annular stream that impinges on the orifice pair

        U_r_c = mdot_r / (self.injector_props['rho_r'] * A_r_m)
        U_z_c = mdot_z / (self.injector_props['rho_z'] * A_z_m)

        LMR_c = (self.injector_props['rho_r'] * (U_r_c ** 2) * A_lr_c) / (self.injector_props['rho_z'] * (U_z_c ** 2) * A_lz_c)
        theta_c = alpha * np.arctan(beta * LMR_c) * (180 / np.pi) + 20 #theta_c = converging angle of nozzle

    # -------------------------------------------SETTING UP PLOTS-------------------------------------------------#

        # Array of gap sizes between 0.05 mm and 0.8 mm 
        a = np.linspace(0.05 / 1000, 0.8 / 1000, 100000) # ALL THE MATH IN THIS 'SETTING UP PLOTS' SECTION NEEDS TO BE LOOKED AT ---> THIS ARRAY MESSES UP MY STUFF 

        A_lr = np.pi * ((self.injector_props['d1'] / 2 / 1000) ** 2)
        A_lz = (gap / 1000) * (self.injector_props['d1'] / 1000)
        U_r = mdot_r / (self.injector_props['rho_r'] * A_r_m)
        U_z = mdot_z / (self.injector_props['rho_z'] * np.pi * ((R_p + a) ** 2 - (R_p) ** 2)) #THIS IS WHERE THE CODE BREAKS BECAUSE YOU'RE ADDING THAT 'a' ARRAY SO THE SIZES NEED TO MATCH WITH 'R_p'

        

        LMR = (self.injector_props['rho_r'] * (U_r ** 2) * A_lr) / (self.injector_props['rho_z'] * (U_z ** 2) * A_lz)
        theta = (alpha * np.arctan(beta * LMR) * (180 / np.pi)) + 20

        sizingdict = {
            'A_r_mm': A_r_mm,
            'A_z_mm': A_z_mm,
            'n_oriface_pairs': n_orifice_pairs,
            'BF': BF,
            'gap': gap,
            'LMR': LMR_c,
            'theta_c': theta_c
        }
        return(sizingdict)
