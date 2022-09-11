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

tca1 = TCA(tcaprops, bartz_props=bartz) #TCA Obj
chamber1 = Chamber(tca1, geo)
inj1 = Injector(tca1, chamber1, inj)
Engine = TCA(tca_props=tcaprops, bartz_props=bartz, CHAMBERobj=chamber1, INJobj=inj1)

print(Engine.CHAMBERobj.geometric_props['A_t'])

Engine.plotparams('R_c', 'BF', .1, .2, 1000)
#tca1.plotparams('F', 'mdot', 400, 450)





