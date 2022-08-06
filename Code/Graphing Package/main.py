from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *
from Graphing_Code import *

# listTCA = []

# for i in range(300,320):
#     list.append(TCA('JetA', 'LOX', 400, 2.1, i))

# for obj in list:
#     print(obj.mdot)

# listCHAMBER = []
# listINJECTOR = []



tca1 = TCA('JetA', 'LOX', 400, 2.2, 300) #TCA Obj
#print(tca1.mdot)

geo = {
    'R_c': np.linspace(.1,.25,100),
    'L_characteristic': 1,
}
bartz = {
    'mu': .2,
    'm': .2,
    'w': .2,
    'M': .2
}



chamber1 = Chamber(tca1, geometric_props=geo, bartz_props=bartz) #Chamber obj
print(chamber1.geocalc())
#print(chamber1.geometric_props['A_t']) 

inj = {
    'rho_r': 800,
    'rho_z': 1100,
    'd1': 1.3,
    'd2': 1.3,
    'C_d': .75,
    'delta_P': 75,
    'delta_P_o': 50,
}
#injector1 = Injector(tca1, chamber1, inj) 
#print(injector1.sizingcalc())





