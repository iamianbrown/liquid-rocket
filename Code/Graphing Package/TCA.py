import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj

#Thrust Chamber Assembly Class
class TCA:
    def __init__(self, fuel, oxidizer, F, OF_ratio, p_c, chamber, injector):
        
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.F = F
        self.OF_ratio = OF_ratio
        self.p_c = p_c

        
        CEAvalues = self.CEArun()
        self.system_props = {

        }

#   ------------------------------------------- CEArun Stuff -------------------------------------------

    def CEArun(self):
            
                # C = CEA_Obj( oxName='LOX', fuelName='JetA', isp_units='m/s', cstar_units='m/s', temperature_units='K')

                C = CEA_Obj(oxName=self.geometric_props['oxidizer'], fuelName=self.geometric_props['fuel'])

                # print("O/F, Isp(m/s), Cstar(m/s), Temp(K)")

                max_Isp = 0
                opt_eps = 0

                for eps in np.arange(2, 7, 0.1):
                    # Isp_th = C.get_Throat_Isp(Pc=Pc, MR=of_ratio)
                    # Isp, Cstar, Tc = C.get_IvacCstrTc(Pc=Pc, MR=of_ratio, eps=eps)

                    # print(f"{of_ratio:.2f}     {Isp_th:.0f}      {Cstar:.0f}        {Tc:.0f}   {C.get_PcOvPe(Pc=Pc, MR=2, eps=eps):.0f}")

                    Isp, mode = C.estimate_Ambient_Isp(Pc=self.combustion_props['p_c'], MR=self.combustion_props['OF_ratio'], eps=eps, Pamb=14.7)
                    if Isp > max_Isp:
                        max_Isp = Isp
                        opt_eps = eps

                    #print(f"{eps:.2f}     {Isp:.0f}")

                cstar = C.get_Cstar(Pc=p_c)
                Temps = C.get_Temperatures(Pc=self.p_c,MR=self.OF_ratio,eps=opt_eps)
                T_c = Temps[0]
                v_e = max_Isp  
                mdot = self.F/v_e

                ceadict = {
                    'cstar': cstar,
                    'T_c': T_c,
                    'v_e': v_e,
                    'eps': opt_eps,
                    'mdot': mdot
                }

                return (ceadict)