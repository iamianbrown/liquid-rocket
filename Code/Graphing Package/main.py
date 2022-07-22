from Injector_Code_Test import *
from TCA import *
from Chamber_Size_Test import *

A = TCA('JetA', 'LOX', 400, 2.2, 400)
#print(A.mdot)

geo = {
    'R_c': np.linspace(.1,.25,10),
    'L_characteristic': 1
}
bartz = {
    'mu': .2,
    'm': .2,
    'w': .2,
    'M': .2
}
B = Chamber(A, geometric_props=geo, bartz_props=bartz)
print(B.geocalc())
#print(B.geometric_props['A_t']) 

# inj = {
#     'rho_r': 800,
#     'rho_z': 1100,
#     'd1': 1.3,
#     'd2': 1.3,
#     'C_d': .75,
#     'delta_P': 75,
#     'delta_P_o': 50,
# }
# C = Injector(A, B, inj) 
# print(C.sizingcalc())



