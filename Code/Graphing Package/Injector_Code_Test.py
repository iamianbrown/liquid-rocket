import numpy as np
import matplotlib.pyplot as plt
from TCA import *
from Chamber_Size_Test import *


class Injector:
    # -------------------------------------------Initial Data-------------------------------------------------#
    def __init__(self, COMBUSTIONobj, CHAMBERobj, injector_props):  

        '''
        _r designates radial propellant
        _z designates axial propellant

        This class is considered to be last in the general math flow of this program
        You must pass in a 'Combustion' object, a 'Chamber' object, and one dictionary: injector_props with the keys of...
            rho_r, rho_z, d1, d2, C_d, delta_P, delta_P_o
        This Injector object holds a dictionary with all of the geometric paramters of your pintle injector
        '''

        self.COMBUSTIONobj = COMBUSTIONobj
        self.CHAMBERobj = CHAMBERobj
        self.injector_props = injector_props #user will pass this in as a dictionary with the values for rho_r, rho_z, d1, d2, C_d, delta_P, delta_P_o

        PINTLEvalues = self.pintle_system_calc()
        INJvalues = self.sizingcalc()
        self.injmerge = {
            'D_p': PINTLEvalues['D_p'],
            'R_p': PINTLEvalues['R_p'],
            'L_pintle': PINTLEvalues['L_pintle'],
            'A_r_mm': INJvalues['A_r_mm'],
            'A_z_mm': INJvalues['A_z_mm'],
            'orifice_pairs': INJvalues['orifice_pairs'],
            'n_orifice_pairs': INJvalues['n_orifice_pairs'],
            'BF': INJvalues['BF'],
            'gap': INJvalues['gap'],
            'LMR': INJvalues['LMR'],
            'LMR_c': INJvalues['LMR_c'],
            'theta_c': INJvalues['theta_c'],
            #'theta': INJvalues['theta']
        }
        self.injector_props.update(self.injmerge)

    def pintle_system_calc(self):
        
        '''
        This function calculates the big picture pintle dimensions
        Input: strictly utilizes the CHAMBERobj's geometric_props dictionary
        Outputs: returns a dictionary containing...
            D_p, R_p, L_pintle
        '''
        
        R_c = self.CHAMBERobj.geometric_props['R_c']
        L_c = self.CHAMBERobj.geometric_props['L_c']

        #Calculate pintle diameter and radius using chamber-to-pintle ratio
        D_p = R_c / 3 #Pintle diameter 
        R_p = D_p / 2 #Pintle radius
        L_pintle = L_c/3  # Pintle length
        skip_length = D_p * self.injector_props['skip_distance']

        pintle_system_dict = {
            'D_p': D_p,
            'R_p': R_p,
            'L_pintle': L_pintle,
        }
        return(pintle_system_dict)

    def sizingcalc(self):
        
        PINTLEvalues = self.pintle_system_calc()
        D_p = PINTLEvalues['D_p'] 
        R_p = PINTLEvalues['R_p']
        alpha = 0.7 # Constant needed for LMR and Spray Angle relation
        beta = 2.0 # Constant needed for LMR and Spray Angle relation

        # ---------------------------------------ORIFICE SIZING-------------------------------------------#

        '''
        _m designates value in meters
        _mm designates value in millimeters
        '''

        A_r_m = self.COMBUSTIONobj.mdot_r / (self.injector_props['C_d'] * np.sqrt(2 * self.injector_props['rho_r'] * self.injector_props['delta_P_o']))  # Radial Discharge area [m^2]
        A_r_mm = A_r_m * 1000000  # mm^2

        a1 = np.pi * ((self.injector_props['d1'] * 0.5) ** 2)  # area of each orifice in row 1
        a2 = np.pi * ((self.injector_props['d2'] * 0.5) ** 2)  # area of "" row 2

        orifice_pairs = a1 + a2  # area per pair of orifices
        n_orifice_pairs = round(A_r_mm / orifice_pairs)
        BF = (n_orifice_pairs * (self.injector_props['d1'] + self.injector_props['d2'])) / (2 * np.pi * ((D_p * 1000) * 0.5))

        # -------------------------------------ANNULAR GAP SIZING-----------------------------------------#

        A_z_m = self.COMBUSTIONobj.mdot_z / (self.injector_props['C_d'] * np.sqrt(2 * self.injector_props['rho_z'] * self.injector_props['delta_P']))
        A_z_mm = A_z_m * 1000000

        gap = np.sqrt((A_z_mm / np.pi) + ((R_p * 1000) ** 2)) - (R_p * 1000) #Annular gap size (.5 addition correction factor to increase gap size thus decreasing spray angle)

        # ---------------------------------------CALCULATE LMR AND THETA---------------------------------------#

        A_lr_c = np.pi * ((self.injector_props['d1'] / 2000) ** 2 + (self.injector_props['d2'] / 2000) ** 2)  # cross-sectional area for one orifice pair

        A_lz_c = (gap / 1000) * (
                    self.injector_props['d1'] / 1000 + self.injector_props['d2'] / 1000)  # cross-sectional area of annular stream that impinges on the orifice pair


        U_r_c = self.COMBUSTIONobj.mdot_r / (self.injector_props['rho_r'] * A_r_m) #Velocity of radial propellant
        U_z_c = self.COMBUSTIONobj.mdot_z / (self.injector_props['rho_z'] * A_z_m)

        LMR_c = (self.injector_props['rho_r'] * (U_r_c ** 2) * A_lr_c) / (self.injector_props['rho_z'] * (U_z_c ** 2) * A_lz_c)
        theta_c = alpha * np.arctan(beta * LMR_c) * (180 / np.pi) + 20 #theta_c = spray angle

    # -------------------------------------------SETTING UP PLOTS-------------------------------------------------#

        # Array of gap sizes between 0.05 mm and 0.8 mm 
        a = np.linspace(0.05 / 1000, 0.8 / 1000, 100000) # ALL THE MATH IN THIS 'SETTING UP PLOTS' SECTION NEEDS TO BE LOOKED AT ---> THIS ARRAY MESSES UP MY STUFF 

        A_lr = np.pi * ((self.injector_props['d1'] / 2 / 1000) ** 2)
        A_lz = (gap / 1000) * (self.injector_props['d1'] / 1000)
        U_r = self.COMBUSTIONobj.mdot_r / (self.injector_props['rho_r'] * A_r_m)
        U_z = self.COMBUSTIONobj.mdot_z / (self.injector_props['rho_z'] * np.pi * ((R_p + a) ** 2 - (R_p) ** 2)) #THIS IS WHERE THE CODE BREAKS BECAUSE YOU'RE ADDING THAT 'a' ARRAY SO THE SIZES NEED TO MATCH WITH 'R_p'

        LMR = (self.injector_props['rho_r'] * (U_r ** 2) * A_lr) / (self.injector_props['rho_z'] * (U_z ** 2) * A_lz)
        #theta = (alpha * np.arctan(beta * LMR) * (180 / np.pi)) + 20

        sizingdict = {
            'A_r_mm': A_r_mm,
            'A_z_mm': A_z_mm,
            'orifice_pairs': orifice_pairs,
            'n_orifice_pairs': n_orifice_pairs,
            'BF': BF,
            'gap': gap,
            'LMR_c': LMR_c,
            'LMR': LMR,
            'theta_c': theta_c, #the spray angle
        }
        return(sizingdict)

    def pintle_graph(self):

        a = np.linspace(0.05/1000, 0.8/1000, 100000)
        alpha = 0.7 # Constant needed for LMR and Spray Angle relation
        beta = 2.0 # Constant needed for LMR and Spray Angle relation
        gap_tolerance = 0.05/1000
        vals = self.sizingcalc()
        gap_low = vals['gap'] - (gap_tolerance * 1000)
        gap_high = vals['gap'] + (gap_tolerance * 1000)
        percent_error = 100*(vals['n_orifice_pairs'] - (vals['A_r_mm'] / vals['orifice_pairs'])) / (vals['A_r_mm'] / vals['orifice_pairs'])
        

        fig, axs = plt.subplots(2, 1, figsize=(12,12))
        fig.suptitle('LMR and Theta vs annular gap')

        axs[0].plot(a, vals['LMR'])
        axs[1].plot(a, vals['theta'])

        axs[0].set_xlabel('annular gap size')
        axs[0].set_ylabel('LMR')
        axs[1].set_xlabel('annular gap size')
        axs[1].set_ylabel('spray angle')
        axs[0].locator_params(axis="x", nbins=15)
        axs[0].locator_params(axis="y", nbins=15)
        axs[1].locator_params(axis="x", nbins=15)
        axs[1].locator_params(axis="y", nbins=20)
        axs[0].set_xlim(a[0])
        axs[1].set_xlim(a[0])
        axs[0].set_ylim(0)
        axs[1].set_ylim(0)

        #--------------------------------LMR AND THETA TOLERANCED VALES---------------------------------------#

        LMR_low = np.interp(((vals['gap']/1000) - gap_tolerance), a, vals['LMR'])
        LMR_high = np.interp(((vals['gap']/1000) + gap_tolerance), a, vals['LMR'])

        vals['theta']

        theta_low = alpha*np.arctan(beta*LMR_low)*(180/np.pi) + 20
        theta_high = alpha*np.arctan(beta*LMR_high)*(180/np.pi) + 20

        axs[0].vlines(x = vals['gap']/1000, label = 'With gap determined by DA eqn', linewidth = 1, linestyle = 'dashed', color = 'red', ymin = 0, ymax = vals['LMR_c'])
        axs[0].vlines(x = (vals['gap']/1000) - gap_tolerance, linewidth = 1,linestyle = 'dashed', color = 'red', ymin = 0, ymax = LMR_low)
        axs[0].vlines(x = (vals['gap']/1000) + gap_tolerance, linewidth = 1,linestyle = 'dashed', color = 'red', ymin = 0, ymax = LMR_high)

        axs[1].vlines(x = vals['gap']/1000, label = 'With gap determined by DA eqn', linewidth = 1,linestyle = 'dashed', color = 'red', ymin = 0, ymax = vals['theta_c'])
        axs[1].vlines(x = (vals['gap']/1000) - gap_tolerance, linewidth = 1,linestyle = 'dashed', color = 'red', ymin = 0, ymax = theta_low)
        axs[1].vlines(x = (vals['gap']/1000) + gap_tolerance, linewidth = 1,linestyle = 'dashed', color = 'red', ymin = 0, ymax = theta_high)

        axs[0].hlines(y = LMR_low, xmin = 0, xmax = (vals['gap']/1000) - gap_tolerance, linewidth = 1,linestyle = 'dashed', color = 'red')
        axs[0].hlines(y = LMR_high, xmin = 0, xmax = (vals['gap']/1000) + gap_tolerance,linewidth = 1,linestyle = 'dashed', color = 'red')
        axs[1].hlines(y = theta_low, xmin = 0, xmax = (vals['gap']/1000) - gap_tolerance,linewidth = 1,linestyle = 'dashed', color = 'red')
        axs[1].hlines(y = theta_high, xmin = 0, xmax = (vals['gap']/1000) + gap_tolerance,linewidth = 1,linestyle = 'dashed', color = 'red')
        axs[0].hlines(y =vals['LMR_c'], xmin = 0, xmax = (vals['gap']/1000),linewidth = 1,linestyle = 'dashed', color = 'red')
        axs[1].hlines(y = vals['theta_c'], xmin = 0, xmax = (vals['gap']/1000),linewidth = 1,linestyle = 'dashed', color = 'red')

        plt.show()

        #------------------------------------------OUTPUTS------------------------------------------------#

        print('Radial discharge area = {} mm^2'.format(round(vals['A_r_mm'], 2)))
        print('Axial discharge area = {} mm^2'.format(round(vals['A_z_mm'], 2)))
        print('Number of orifice pairs = {}'.format(vals['n_orifice_pairs']))
        print('Orifice area percent error = {}%'.format(round(percent_error, 2)))
        print('BF = {}'.format(round(vals['BF'], 2)))
        print('------------------------------------------')
        print('USING THE DISCHARGE AREA EQUATION:')
        print('Gap size = {}mm < {}mm < {}mm'.format(round(gap_low, 3), round(vals['gap'], 3), round(gap_high, 3)))
        print('LMR = {} < {} < {}'.format(round(LMR_low, 3), round(vals['LMR_c'], 3), round(LMR_high, 3)))
        print('Spray Angle = {}deg < {}deg < {}deg'.format(round(theta_low, 2), round(vals['theta_c'], 2), round(theta_high, 2)))
