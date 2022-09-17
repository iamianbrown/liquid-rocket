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
    'F': 5000,  #5000 N is Grunt's Thrust
    'OF_Ratio': 2.25, 
    'p_c': 300 #psi
}

chamber_geometry = {
    'R_c': .13, #m
    'L_characteristic': .9, #m
}

bartzprops = {
    'mu': .2,
    'm': .2,
    'w': .2,
    'M': .2,
}

injectorprops = {
    'rho_r': 1141, #kg/m^3
    'rho_z': 810, #kg/m^3
    'd1': 1.3, #mm
    'd2': 1.3, #mm
    'C_d': .75,
    'delta_P': 517107, #Pa
    'delta_P_o': 310264, #Pa
    #need to add skip distance input into injector
    'skip_distance': 0.9 #m

}

combustion1 = Combustion(tca_props=combustionprops, bartz_props=bartzprops) #TCA Obj
chamber1 = Chamber(combustion1, geometric_props=chamber_geometry)
inj1 = Injector(combustion1, chamber1, injectorprops)
MkII = Engine(COMBUSTIONobj=combustion1, CHAMBERobj=chamber1, INJobj=inj1)

#print(chamber1.geometric_props['A_t'])

#print(inj1.injector_props['L_pintle'])

MkII.plotparams('A_r_mm', 'theta_c', 50, 80, 1000)
# MkII.plotparams('F', 'mdot', 4900, 5500, 1000)






