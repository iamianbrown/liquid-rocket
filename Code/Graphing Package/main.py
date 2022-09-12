from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *
from Graphing_Code import *

'''
    The general mathematical 'flow' of this program goes from TCA -> Chamber -> Injector
'''

combustionprops = {
    'fuel': 'JetA',
    'oxidizer': 'LOX',
    'F': 4900,  #4900 N is Grunt's Thrust
    'OF_Ratio': 2.25,
    'p_c': 300
}

chamber_geometry = {
    'R_c': .12,
    'L_characteristic': .9,
}

bartzprops = {
    'mu': .2,
    'm': .2,
    'w': .2,
    'M': .2,
}

injectorprops = {
    'rho_r': 800,
    'rho_z': 1100,
    'd1': 1.3,
    'd2': 1.3,
    'C_d': .9,
    'delta_P': 75,
    'delta_P_o': 45, 
}

combustion1 = Combustion(tca_props=combustionprops, bartz_props=bartzprops) #TCA Obj
chamber1 = Chamber(combustion1, geometric_props=chamber_geometry)
inj1 = Injector(combustion1, chamber1, injectorprops)

MkII = Engine(COMBUSTIONobj=combustion1, CHAMBERobj=chamber1, INJobj=inj1)

MkII.plotparams('R_c', 'BF', .1, .2, 1000)






