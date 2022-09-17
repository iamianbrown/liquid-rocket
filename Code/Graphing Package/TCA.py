import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from Chamber_Size_Test import *
from Injector_Code_Test import *

#Thrust Chamber Assembly Class
class Combustion:
    def __init__(self, tca_props, bartz_props): 
        
        ''' 
            This class contains the everything related to the combustion parameters for the overall engine.
            You must pass in two dictionaries relating to tca_props (overall engine properties) and bartz_props (bartz calculations)
            INPUTS - 
                tca_props: top level engine/combustion properties
                bartz_props: parameters needed for bartz calculations
        '''

        #Construction of your objects
        self.tca_props = tca_props 
        self.bartz_props = bartz_props
        
        self.CEAvalues = self.CEArun() #Runs the CEArun() function to calculate CEArun parameters
        self.mdot = self.CEAvalues['mdot']
        self.mdot_r = self.CEAvalues['mdot_r']
        self.mdot_z = self.CEAvalues['mdot_z']
        self.eps = self.CEAvalues['eps']
        self.v_e = self.CEAvalues['v_e']
        self.T_c = self.CEAvalues['T_c']
        self.cstar = self.CEAvalues['cstar']
        self.gamma = self.CEAvalues['gamma']
        self.Cp = self.CEAvalues['Cp']
        self.Cv = self.CEAvalues['Cv']
        self.R = self.CEAvalues['R']
        self.tca_props.update(self.CEAvalues) #Merges the CEArun dictionary output into tca_props

        BARTZvalues = self.Bartzcalc() #Runs the Bartzcalc() function to calculate Bartz calculation variables
        self.bartzmerge = {
            'Pr': BARTZvalues['Pr'],
            'T_w': BARTZvalues['T_w'],
            'sigma': BARTZvalues['sigma'],
            }
        self.bartz_props.update(self.bartzmerge) #Merges the Bartzcalc dictionary outputs into bartz_props

#   ------------------------------------------- CEArun Stuff -------------------------------------------

    def CEArun(self):
            
                '''
                    This is the implementation of the rocketcea module which calculates the CEA outputs within python. 
                    It returns a dictionary of the CEArun outputs
                    Input: Utilizes the tca_props dictionary the user passes into the TCA class
                    Output: Dictionary containing...
                        cstar, T_c, v_e, eps, mdot, mdot_r, mdot_z, gamma, Cp, Cv, R
                '''

                # C = CEA_Obj( oxName='LOX', fuelName='JetA') ---> How to create the CEA object!

                C = CEA_Obj(oxName=self.tca_props['oxidizer'], fuelName=self.tca_props['fuel'])

                # print("O/F, Isp(m/s), Cstar(m/s), Temp(K)") ---> Not sure this works yet.
                
                
                pressure_ratio = self.tca_props['p_c']/14.7 #Dimensionless value
                opt_eps = C.get_eps_at_PcOvPe(Pc=self.tca_props['p_c'], MR = self.tca_props['OF_Ratio'], PcOvPe = pressure_ratio) #optimal expansion ratio
                Isp, mode = C.estimate_Ambient_Isp(Pc = self.tca_props['p_c'], MR = self.tca_props['OF_Ratio'], eps = opt_eps, Pamb = 14.7) #Specific Impulse
                

                cstar = C.get_Cstar(Pc=self.tca_props['p_c'], MR = self.tca_props['OF_Ratio'])
                Temps = C.get_Temperatures(Pc=self.tca_props['p_c'],MR=self.tca_props['OF_Ratio'],eps=opt_eps)
                mw_gamma_chamber = C.get_Chamber_MolWt_gamma(Pc=self.tca_props['p_c'], MR=self.tca_props['OF_Ratio'], eps = opt_eps)
                gamma = mw_gamma_chamber[1]
                T_c = Temps[0] / 1.8
                v_e = Isp * 9.8066
                mdot = self.tca_props['F']/v_e
                mdot_r = mdot / (1 + self.tca_props['OF_Ratio']) * self.tca_props['OF_Ratio'] 
                mdot_z = mdot / (1 + self.tca_props['OF_Ratio']) 
                Cp = C.get_Chamber_Cp(Pc=self.tca_props['p_c'], MR=self.tca_props['OF_Ratio'], eps=opt_eps)
                Cv = Cp/gamma
                R = Cp-Cv

                ceadict = {
                    'cstar': cstar,
                    'T_c': T_c,
                    'v_e': v_e,
                    'eps': opt_eps,
                    'mdot': mdot,
                    'mdot_r': mdot_r,
                    'mdot_z': mdot_z,
                    'gamma': gamma,
                    'Cp': Cp,
                    'Cv': Cv,
                    'R': R
                }
                return (ceadict)

    def Bartzcalc(self):

        '''
            This function contains the calculations using the bartz formulas.
            It returns a dictionaries with the answers to these bartz equations
            Input: Utilizes the bartz_props dictionary the user passes into the TCA class
            Output: Dictionary containing...
                Pr, T_w, sigma
        '''

        gamma = self.gamma
        T_c = self.T_c

        Pr = 4*gamma/(9*gamma - 5) # prandtl number given in Bartz paper
        T_w = T_c # assume for now
        sigma = 1/((1/2)*(T_w/T_c)*(1+((gamma-1)/2)*self.bartz_props['M']**2)+1/2)**(0.8-(self.bartz_props['w']/5))*(1+((gamma-1)/2)*self.bartz_props['M']**2)**(self.bartz_props['w']/5) # dimensionless factor

        bartzdict = {
            'Pr': Pr,
            'T_w': T_w,
            'sigma': sigma,
        }

        return(bartzdict)