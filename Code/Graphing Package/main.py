from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *
from Graphing_Code import *

'''
    The general mathematical 'flow' of this program goes from TCA -> Chamber -> Injector
'''

combustionprops = { #goes into combustion class
    'fuel': 'JetA',
    'oxidizer': 'LOX',
    'F': 5000,  #5000 N is Grunt's Thrust
    'OF_Ratio': 2.25, 
    'p_c': 300 #psi - chamber pressure 
}

chamber_geometry = { #goes into chamber class
    'R_c': .065, #Chamber radius [m]
    'L_characteristic': .9, #Characteristic length [m]
}

bartzprops = { #goes into combustion class
    'mu': .2,
    'm': .2,
    'w': .2,
    'M': .2,
}

injectorprops = { #goes into injector class
    'rho_r': 1141, #kg/m^3 - density of radial propellant 
    'rho_z': 810, #kg/m^3 - density of annular propellant
    'd1': 1.3, #mm - row 1 orifice diameters for pintle
    'd2': 1.3, #mm - row 2 orifice diameters for pintle
    'C_d': .75, #Discharge coefficient
    'delta_P': 517000, #Pa (fuel pressure drop) - usually 25% of p_c  #15% is 310264 and 25% is 517000
    'delta_P_o': 310264, #Pa (oxidizer pressure drop) - usually 15% of p_c
    #need to add skip distance input into injector
    'skip_distance': 0.9 #m

}

combustion1 = Combustion(tca_props=combustionprops, bartz_props=bartzprops) #TCA Obj
chamber1 = Chamber(combustion1, geometric_props=chamber_geometry)
inj1 = Injector(combustion1, chamber1, injectorprops)
MkII = Engine(COMBUSTIONobj=combustion1, CHAMBERobj=chamber1, INJobj=inj1)


#print(combustion1.tca_props['eps'])
#print(chamber1.geometric_props['A_t'])
#print(inj1.injector_props['BF'])
#print(MkII.COMBUSTIONobj.tca_props['Cp'])

#inj1.pintle_graph()


#MkII.plotparams('R_c', 'theta_c', .07, .12, 1000)
#MkII.plotparams('F', 'mdot', 4900, 5500, 1000)
MkII.engineVisual()






