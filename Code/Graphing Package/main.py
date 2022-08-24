from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *
from Graphing_Code import *

tcaprops = {
    'fuel': 'JetA',
    'oxidizer': 'LOX',
    'F': 400,
    'OF_Ratio': 2.2,
    'p_c': 300
}

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

tca1 = TCA(tcaprops) #TCA Obj
chamber1 = Chamber(tca1, geo, bartz)
inj1 = Injector(tca1, chamber1, inj)
tca1 = TCA(tcaprops, chamber1, inj1)

tca1.plotparams('mdot', 'LMR', 1, 3)
#tca1.plotparams('F', 'mdot', 400, 450)





