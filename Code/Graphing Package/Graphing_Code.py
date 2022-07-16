#Don't be scared... 
import Chamber_Size_Test
import Injector_Code_Test
import numpy as np

ind = 0
dep = 0

class graphing:
    def __init__(self,Injector,ChamberSizing,independent,dependent):

        self.Injector = Injector
        self.ChamberSizing = ChamberSizing
        self.independent = independent
        self.dependent = dependent

    def prompt(self):
        
        oxidizer = input("Enter Oxidizer: ")
        fuel = input("Enter Fuel: ")
        F = float(input("Enter thrust force [N]: "))
        Pe = float(input("Enter exit pressure [psi]: "))
        Pc = float(input("Enter chamber pressure [psi]: "))
        OFratio = float(input("Enter OFratio: "))
        p_0 = float(input("Enter p_0 [Pa]: "))
        R_c = float(input("Enter radius of chamber [m]: "))
        L_characteristic = float(input("Enter characteristic length: "))
        Theta_c = float(input("Enter angle of converging section: "))
        Theta_e = float(input("Enter angle of diverging section: "))
        spray_angle = float(input("Enter spray angle: "))
        mu = float(input("Enter mu: "))
        m = float(input("Enter m: "))
        w = float(input("Enter w: "))
        M = float(input("Enter M: "))
        d_c = float(input("Enter chamber diameter [m]: "))
        rho_r = float(input("Enter radial propellant density [kg/m^3]: "))
        rho_z = float(input("Enter axial propellant density [kg/m^3]: "))
        d1 = float(input("Enter orifice row 1 diameter: "))
        d2 = float(input("Enter orifice row 2 diameter: "))
        C_d = float(input("Enter discharge coefficient: "))
        delta_P = float(input("Enter fuel pressure drop: "))
        delta_P_o = float(input("Enter oxidizer pressure drop: "))

        # test = Chamber_Size_Test.ChamberSizing.__init__()

        pedro = np.array([Pe,Pc,OFratio,p_0,R_c,L_characteristic,Theta_c,Theta_e,spray_angle,mu,m,w,M,d_c,rho_r,rho_z,d1,d2,C_d,delta_P,delta_P_o])
        # cham = np.array(F,OFratio,p_0,R_c,L_characteristic,Theta_c,Theta_e,spray_angle,mu,m,w,M)
        

        for i in pedro:
            if pedro[i] == 0.1:
                
                independentrange = float(input("Enter max value to test independent variable: "))
                ind = np.linspace(0,independentrange,1000)
                pedro[i] = ind
        
        step1 = Chamber_Size_Test.ChamberSizing(F,OFratio,p_0,R_c,L_characteristic,Theta_c,Theta_e,spray_angle,mu,m,w,M)


        
        #if user inputs 0.1 for any numerical input above ---> code will assume that variable to be the independent variable graphed
        #if user inputs 0 for any numerical input above ---> code will assume that variable to be the dependent variable graphed


