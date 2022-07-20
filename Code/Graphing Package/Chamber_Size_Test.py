import numpy as np
import matplotlib.pyplot as plt
import Injector_Code_Test
from rocketcea.cea_obj import CEA_Obj
import TCA


class Chamber:
    
    
    def __init__(self, TCAobj, geometric_props, bartz_props): 
        
        self.TCAobj = TCAobj 

        self.geometric_props = geometric_props
        # geometric_props should be an input dictionary of R_c and L_characteristic
        GEOvalues = self.geocalc()
        geomerge = {
            'A_t': GEOvalues['A_t'],    # Throat area
            'eps_c': GEOvalues['eps_c'],# compression ratio
            'V_c': GEOvalues['V_c'],    # Volume of combustion chamber
            'L_c': GEOvalues['L_c'],    # Length of combustion chamber
            'R_t': GEOvalues['R_t'],    # Throat radius
            'R_e': GEOvalues['R_e'],    # Exit radius
            'R_b': GEOvalues['R_b'],    # Radius of bigger arc, usually 1.5 times the throat radius
            'R_s': GEOvalues['R_s'],    # Radius of smaller arc, usually 0.4 times the throat radius
            'R_pintle': GEOvalues['R_pintle'],  #Pintle radius
            'L_pintle': GEOvalues['L_pintle'],  #Pintle length
            }
        geomerge.update(self.geometric_props) 
        
        self.bartz_props = bartz_props
        # bartz_props should be an input dictionary containing mu, m, w, M
        BARTZvalues = self.Bartzcalc()
        bartzmerge = {
            'D_t': BARTZvalues['D_t'],
            'Pr': BARTZvalues['Pr'],
            'c_t': BARTZvalues['c_t'],
            'T_w': BARTZvalues['T_w'],
            'sigma': BARTZvalues['sigma'],
            }
        bartzmerge.update(bartz_props)
        
    
    def geocalc(self):
        A_t = TCA(self.TCAobj).mdot/(TCA(self.TCAobj).p_c/np.sqrt(TCA(self.TCAobj).T_c)*np.sqrt(TCA(self.TCAobj).gamma/R*(2/(TCA(self.TCAobj).gamma+1))**((TCA(self.TCAobj).gamma+1)/(TCA(self.TCAobj).gamma-1))))

        # The compression ratio will be
        eps_c = (np.pi*self.geometric_props['R_c']**2)/A_t

        V_c = self.geometric_props['L_characteristic']*A_t # Volume of combustion chamber [m^3]

        # If we want the combustion chamber radius to be R_c , then we need a
        # length:
        L_c = V_c/(np.pi*self.geometric_props['R_c']**2) # Length of combustion chamber

        ## STEP 3: GEOMETRICAL PARAMETERS FOR THE NOZZLE

        R_t = np.sqrt(A_t/np.pi)  # Throat radius
        R_e = np.sqrt(TCA(self.TCAobj).eps)*R_t  # Exit radius

        R_b = 1.5*R_t     # Radius of bigger arc, usually 1.5 times the throat radius
        R_s = 0.4*R_t     # Radius of smaller arc, usually 0.4 times the throat radius

        R_pintle = self.geometric_props['R_c']/4  # Pintle radius
        L_pintle = L_c/3  # Pintle length

        geodict = {
            'A_t': A_t,
            'eps_c': eps_c,
            'V_c': V_c,
            'L_c': L_c,
            'R_t': R_t,
            'R_e': R_e,
            'R_b': R_b,
            'R_s': R_s,
            'R_pintle': R_pintle,
            'L_pintle': L_pintle,
        }

        return(geodict)

    def Bartzcalc(self):
        GEOvaleus = self.geocalc()
        D_t = 2*GEOvaleus['R_t']
        Pr = 4*TCA(self.TCAobj).gamma/(9*TCA(self.TCAobj).gamma - 5) # prandtl number given in Bartz paper
        c_t = np.sqrt((1/TCA(self.TCAobj).gamma)*((TCA(self.TCAobj).gamma+1)/2)**((TCA(self.TCAobj).gamma+1)/(TCA(self.TCAobj).gamma-1))*R*TCA(self.TCAobj).T_c) # characteristic velocity
        T_w = TCA(self.TCAobj).T_c                # assume for now
        sigma = 1/((1/2)*(T_w/TCA(self.TCAobj).T_c)*(1+((TCA(self.TCAobj).gamma-1)/2)*self.bartz_props['M']**2)+1/2)**(0.8-(self.bartz_props['w']/5))*(1+((TCA(self.TCAobj).gamma-1)/2)*self.bartz_props['M']**2)**(self.bartz_props['w']/5) # dimensionless factor

        bartzdict = {
            'D_t': D_t,
            'Pr': Pr,
            'c_t': c_t,
            'T_w': T_w,
            'sigma': sigma,
        }

        return(bartzdict)
    
        
