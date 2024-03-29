import numpy as np
import matplotlib.pyplot as plt
from Injector_Code_Test import *
from rocketcea.cea_obj import CEA_Obj
from TCA import *

class Chamber:
    
    def __init__(self, COMBUSTIONobj, geometric_props): 
        
        '''
            This class is considered to be second in the general math flow of this program.
            You must pass in a 'Combustion' object as well as two dictionaries:
                1) geometric_props with: 'R_c' and 'L_characteristic'
                2) bartz_props with: 'mu' and 'm' and 'w' and 'M'
            This object holds two dictionaries for the geometric calculations and the bartz calculations
        '''

        self.COMBUSTIONobj = COMBUSTIONobj 

        self.geometric_props = geometric_props
        GEOvalues = self.geocalc()
        self.geomerge = {
            'A_t': GEOvalues['A_t'],    # Throat area
            'eps_c': GEOvalues['eps_c'],# compression ratio
            'V_c': GEOvalues['V_c'],    # Volume of combustion chamber
            'L_c': GEOvalues['L_c'],    # Length of combustion chamber
            'R_t': GEOvalues['R_t'],    # Throat radius
            'R_e': GEOvalues['R_e'],    # Exit radius
            'R_b': GEOvalues['R_b'],    # Radius of bigger arc, usually 1.5 times the throat radius
            'R_s': GEOvalues['R_s'],    # Radius of smaller arc, usually 0.4 times the throat radius
            }
        self.geometric_props.update(self.geomerge) 

    def geocalc(self):
        
        '''
            This function does the general geometry calculations for the combustion chamber.
            It returns a dictionary containing some geometric properties of the combustion chamber.
            Input: Utilizes COMBUSTIONobj and dictionary geometric_props that user passes in
            Outputs: Returns a dictionary containing...
                A_t, eps_c, V_c, L_c, R_t, R_e, R_b, R_s
        '''

        mdot = self.COMBUSTIONobj.tca_props['mdot']
        p_c = self.COMBUSTIONobj.tca_props['p_c']
        T_c = self.COMBUSTIONobj.tca_props['T_c']
        gamma = self.COMBUSTIONobj.tca_props['gamma']
        R = self.COMBUSTIONobj.tca_props['R']
        eps = self.COMBUSTIONobj.tca_props['eps']
        c_star = self.COMBUSTIONobj.tca_props['cstar'] / 3.281
        
        A_t = mdot*c_star/(p_c*6894.76) #Throat Area

        # The compression ratio will be
        eps_c = (np.pi*self.geometric_props['R_c']**2)/A_t

        V_c = self.geometric_props['L_characteristic']*A_t # Volume of combustion chamber [m^3]

        # If we want the combustion chamber radius to be R_c , then we need a
        # length:
        L_c = V_c/(np.pi*self.geometric_props['R_c']**2) # Length of combustion chamber

        ## STEP 3: GEOMETRICAL PARAMETERS FOR THE NOZZLE

        R_t = np.sqrt(A_t/np.pi)  # Throat radius
        R_e = np.sqrt(eps)*R_t  # Exit radius

        R_b = 1.5*R_t     # Radius of bigger arc, usually 1.5 times the throat radius
        R_s = 0.382*R_t     # Radius of smaller arc, usually 0.4 times the throat radius

        geodict = {
            'A_t': A_t, #throat area 
            'eps_c': eps_c, #expansion ratio
            'V_c': V_c, #volume of chamber
            'L_c': L_c, #length of chamber
            'R_t': R_t, #throat radius
            'R_e': R_e, #exit radius
            'R_b': R_b, # Radius of bigger arc, usually 1.5 times the throat radius
            'R_s': R_s, # Radius of smaller arc, usually 0.4 times the throat radius
            'c_star': c_star, # c star [ft/s]
            'mdot' : mdot,
            'T_c': T_c,
            'gamma': gamma,
            'p_c': p_c,
            'R': R,
        }

        return(geodict)
