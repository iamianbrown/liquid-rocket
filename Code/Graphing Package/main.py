from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *
from Graphing_Code import *


tca1 = TCA('JetA', 'LOX', 400, 2.2, 300) #TCA Obj
tca1.plotparams('F', 'mdot', 400, 450)

geo = {
    'R_c': .25,
    'L_characteristic': 1,
}

bartz = {
    'mu': .2,
    'm': .2,
    'w': .2,
    'M': .2,
}

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

# tca1 = tca(a, b, c, d)
# tca1.plotparams(ind, dep, step/linspace, misc)



