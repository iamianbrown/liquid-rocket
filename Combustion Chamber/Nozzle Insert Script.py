# Assuming a submerged nozzle insert design for calculation 
# purposes but with an extension to make it flush, should  
# still make most calculations the same since the extended 
# area isn't affected by the support

import numpy as np

class NozzleInsert:
    
    def __init__(self, R_c, thickness, length, A_t, P_c, shear_str, compressive_str):
        self.R_c = R_c
        self.thick = thickness
        self.L = length
        self.A_t = A_t
        self.P_c = P_c
        self.shear_str = shear_str
        self.compressive_str = compressive_str

    def A_compressive(self):  # area of bearing ring
        area = np.pi * (self.R_c - self.thick) * self.thick
        return area

    def stress_compressive(self):  # compressive stress on nozzle
        stress = ((np.pi/4) * (self.R_c**2) * self.P_c)/self.A_compressive()
        return stress
        
    def SF_compressive(self):
        safety_factor = self.compressive_str/self.stress_compressive()
        return safety_factor
    
    def A_shear(self):
        area = np.pi * self.R_c * self.L
        return area
    
    def stress_shear(self):
        stress = ((np.pi/4) * (R_c**2) * P_c)/self.A_shear()
        return stress
    
    def SF_shear(self):
        safety_factor = self.shear_str/self.stress_shear()
        return safety_factor
    
    def is_sizing_good(self):  
        safety_compressive = self.SF_compressive()
        safety_shear = self.SF_shear()

        if safety_compressive > 1.5 and safety_shear > 1.5:
            print("Good sizing")
        elif safety_compressive < 1.5:
            print("Compressive safety is small")
        else:
            print("Shear safety is small")
                  
if __name__ == '__main__':
    R_c = 0.06 # combustion chamber radius [m]
    thick = 0.005   # thickness of compression area [m]
    L = 0.015   # length of shear line [m]

    A_t = 0.0016719730201810484 # throat area [m^2]
    P_c = 2e6 # combustion chamber pressure [Pa]

    ## Material properties:
    graphite = {
                "Shear strength": 3.37e6,
                "Compressive strength": 114.5e6,
                }
            

nozzle1 = NozzleInsert(R_c, thick, L, A_t, P_c, graphite["Shear strength"], graphite["Compressive strength"], )        

nozzle1.is_sizing_good()