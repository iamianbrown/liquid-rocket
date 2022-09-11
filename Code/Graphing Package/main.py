from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *
from Graphing_Code import *

'''
    The general mathematical 'flow' of this program goes from TCA -> Chamber -> Injector
'''

tcaprops = {
    'fuel': 'JetA',
    'oxidizer': 'LOX',
    'F': 4900,  #4900 N is Grunt's Thrust
    'OF_Ratio': 2.25,
    'p_c': 300
}

geo = {
    'R_c': .12,
    'L_characteristic': .9,
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
    'C_d': .9,
    'delta_P': 75,
    'delta_P_o': 45, 
}

combustion1 = Combustion(tcaprops, bartz_props=bartz) #TCA Obj
chamber1 = Chamber(combustion1, geo)
inj1 = Injector(combustion1, chamber1, inj)

MkII = Engine(COMBUSTIONobj=combustion1, CHAMBERobj=chamber1, INJobj=inj1)

MkII.plotparams('F', 'A_t', 4900, 5700, 1000)






