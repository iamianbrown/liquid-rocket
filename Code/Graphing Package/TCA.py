import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from Chamber_Size_Test import *
from Injector_Code_Test import *

#Thrust Chamber Assembly Class
class TCA:
    def __init__(self, tca_props, CHAMBERobj=None, INJobj=None): #add INJobj and CHAMBERobj as paramaters
        
        self.tca_props = tca_props
        self.INJobj = INJobj
        self.CHAMBERobj = CHAMBERobj

        self.CEAvalues = self.CEArun()
        self.mdot = self.CEAvalues['mdot']
        self.eps = self.CEAvalues['eps']
        self.v_e = self.CEAvalues['v_e']
        self.T_c = self.CEAvalues['T_c']
        self.cstar = self.CEAvalues['cstar']
        self.gamma = self.CEAvalues['gamma']
        self.Cp = self.CEAvalues['Cp']
        self.Cv = self.CEAvalues['Cv']
        self.R = self.CEAvalues['R']

        self.tca_props.update(self.CEAvalues)

#   ------------------------------------------- CEArun Stuff -------------------------------------------

    def CEArun(self):
            
                # C = CEA_Obj( oxName='LOX', fuelName='JetA', isp_units='m/s', cstar_units='m/s', temperature_units='K')

                C = CEA_Obj(oxName=self.tca_props['oxidizer'], fuelName=self.tca_props['fuel'])

                # print("O/F, Isp(m/s), Cstar(m/s), Temp(K)")

                max_Isp = 0
                opt_eps = 0

                for eps in np.arange(2, 7, 0.1):
                    # Isp_th = C.get_Throat_Isp(Pc=Pc, MR=of_ratio)
                    # Isp, Cstar, Tc = C.get_IvacCstrTc(Pc=Pc, MR=of_ratio, eps=eps)

                    # print(f"{of_ratio:.2f}     {Isp_th:.0f}      {Cstar:.0f}        {Tc:.0f}   {C.get_PcOvPe(Pc=Pc, MR=2, eps=eps):.0f}")

                    Isp, mode = C.estimate_Ambient_Isp(Pc=self.tca_props['p_c'], MR=self.tca_props['OF_Ratio'], eps=eps, Pamb=14.7)
                    if Isp > max_Isp:
                        max_Isp = Isp
                        opt_eps = eps

                    #print(f"{eps:.2f}     {Isp:.0f}")

                cstar = C.get_Cstar(Pc=self.tca_props['p_c'])
                Temps = C.get_Temperatures(Pc=self.tca_props['p_c'],MR=self.tca_props['OF_Ratio'],eps=opt_eps)
                mw_gamma_chamber = C.get_Chamber_MolWt_gamma(Pc=self.tca_props['p_c'], MR=self.tca_props['OF_Ratio'], eps = opt_eps)
                gamma = mw_gamma_chamber[1]
                T_c = Temps[0]
                v_e = max_Isp
                mdot = self.tca_props['F']/v_e
                Cp = C.get_Chamber_Cp(Pc=self.tca_props['p_c'], MR=self.tca_props['OF_Ratio'], eps=opt_eps)
                Cv = Cp/gamma
                R = Cp-Cv

                ceadict = {
                    'cstar': cstar,
                    'T_c': T_c,
                    'v_e': v_e,
                    'eps': opt_eps,
                    'mdot': mdot,
                    'gamma': gamma,
                    'Cp': Cp,
                    'Cv': Cv,
                    'R': R
                }

                return (ceadict)

    def getind(self, independent):
        #return(self.tca_props[independent])
        if independent in self.tca_props.keys():
            return(self.tca_props[independent])
        elif independent in self.CHAMBERobj.geometric_props.keys():
            return(self.CHAMBERobj.geometric_props[independent])
        elif independent in self.CHAMBERobj.bartz_props.keys():
            return(self.CHAMBERobj.bartz_props[independent])
        elif independent in self.INJobj.injector_props.keys():
            return(self.INJobj.injector_props[independent])

    def getdep(self, dependent):
        # Dictvalues = self.CEArun()
        # return(Dictvalues[dependent])
        if dependent in self.tca_props.keys():
            Dictvalues = self.CEArun()
            return(Dictvalues[dependent])
        elif dependent in self.CHAMBERobj.geometric_props.keys():
            Dictvalues = self.CHAMBERobj.geocalc()
            return(Dictvalues[dependent])
        elif dependent in self.CHAMBERobj.bartz_props.keys():
            Dictvalues = self.CHAMBERobj.bartzcalc()
            return(Dictvalues[dependent])
        elif dependent in self.INJobj.injector_props.keys():
            Dictvalues = self.INJobj.sizingcalc()
            return(Dictvalues[dependent])

    def changeVal(self, independent, value):
        #self.tca_props[independent] = value
        if independent in self.tca_props.keys():
            self.tca_props[independent] = value
        elif independent in self.CHAMBERobj.geometric_props.keys():
            self.CHAMBERobj.geometric_props[independent] = value
        elif independent in self.CHAMBERobj.bartz_props.keys():
            self.CHAMBERobj.bartz_props[independent] = value
        elif independent in self.INJobj.injector_props.keys():
            self.INJobj.injector_props[independent] = value
        
        

    def plotparams(self, ind, dep, minVal, maxVal):

            list_ind = []
            list_dep = []

            for j in range(minVal, maxVal):
                self.changeVal(ind,j)
                list_dep.append(self.getdep(dep))
                list_ind.append(self.getind(ind))
            
            plt.plot(list_ind, list_dep, '-g')
            plt.xlabel(ind)
            plt.ylabel(dep)
            plt.show()