import numpy as np
import matplotlib.pyplot as plt
import Injector_Code_Test
from rocketcea.cea_obj import CEA_Obj


class Chamber:
    
    def __init__(self, oxidizer, fuel, F, OF_ratio, p_c, R_c, L_characteristic, Theta_c, Theta_e, spray_angle, mu, m, w, M):

        CEAvalues = self.CEArun()
        self.combustion_props = {
            'F': F,                     # Thrust force [N]
            'OF_ratio': OF_ratio,       # Oxidizer/fuel ratio
            'p_c':p_c,                  # Stagnation pressure at the combustion chamber [Pa]
            'T_c': CEAvalues['T_c'],
            'cstar': CEAvalues['cstar'],
            'v_e': CEAvalues['v_e'],
            'eps': CEAvalues['eps'],
            #calculate R, Cv, Cp, gamma, Exit velocity
            }

        self.geometric_props = {
            'R_c': R_c,                 # Radius of combustion chamber [m]
            'L_characteristic': L_characteristic, # Characteristic length
            'Theta_c': Theta_c,         # Angle of converging section
            'Theta_e': -Theta_e,        # Angle of diverging section
            'spray_angle': spray_angle, # Spray angle
            }
        
        self.bartz_props = {
            'mu': mu, 
            'm': m,
            'w': w,
            'M': M,
            }
        
        self.geocalc() #this line
        CEA_Output = self.CEArun()
        self.m_dot = m_dot
        self.injector = Injector_Code_Test.Injector(d_c,rho_r,rho_z,d1,d2,C_d,delta_P,delta_P_o,mdot) #injector object is being passed into the chamber class

        #add a function that calculates the mass flow 
    
    def CEArun(self):
        
            # C = CEA_Obj( oxName='LOX', fuelName='JetA', isp_units='m/s', cstar_units='m/s', temperature_units='K')

            C = CEA_Obj(oxName=self.oxidizer, fuelName=self.fuel)

            # print("O/F, Isp(m/s), Cstar(m/s), Temp(K)")

            max_Isp = 0
            opt_eps = 0

            for eps in np.arange(2, 7, 0.1):
                # Isp_th = C.get_Throat_Isp(Pc=Pc, MR=of_ratio)
                # Isp, Cstar, Tc = C.get_IvacCstrTc(Pc=Pc, MR=of_ratio, eps=eps)

                # print(f"{of_ratio:.2f}     {Isp_th:.0f}      {Cstar:.0f}        {Tc:.0f}   {C.get_PcOvPe(Pc=Pc, MR=2, eps=eps):.0f}")

                Isp, mode = C.estimate_Ambient_Isp(Pc=self.p_c, MR=self.OF_ratio, eps=eps, Pamb=14.7)
                if Isp > max_Isp:
                    max_Isp = Isp
                    opt_eps = eps

                #print(f"{eps:.2f}     {Isp:.0f}")

            cstar = C.get_Cstar(Pc=p_c)
            Temps = C.get_Temperatures(Pc=self.p_c,MR=self.OF_ratio,eps=opt_eps)
            T_c = Temps[0]
            v_e = max_Isp
            eps = opt_eps   

            ceadict = {
                'cstar': cstar,
                'T_c': T_c,
                'v_e': v_e,
                'eps': eps,

            }

            return (ceadict)
            
    
    
    def geocalc(self):
        m_dot = self.combustion_props['F']/self.CEArun()['v_e']  # Needed mass flow for the required thrust
        A_t = self.m_dot/(self.combustion_props['p_c']/np.sqrt(self.CEArun()['T_0'])*np.sqrt(self.CEArun()['gamma']/R*(2/(self.CEArun()['gamma']+1))**((self.CEArun()['gamma']+1)/(self.CEArun()['gamma']-1))))

        # The compression ratio will be
        eps_c = (np.pi*self.geometric_props['R_c']**2)/self.A_t

        V_c = self.geometric_props['L_characteristic']*self.A_t # Volume of combustion chamber [m^3]

        # If we want the combustion chamber radius to be R_c , then we need a
        # length:
        L_c = self.V_c/(np.pi*self.geometric_props['R_c']**2) # Length of combustion chamber

        ## STEP 3: GEOMETRICAL PARAMETERS FOR THE NOZZLE

        R_t = np.sqrt(self.A_t/np.pi)  # Throat radius
        R_e = np.sqrt(self.CEArun()['eps'])*self.R_t  # Exit radius

        R_b = 1.5*self.R_t     # Radius of bigger arc, usually 1.5 times the throat radius
        R_s = 0.4*self.R_t     # Radius of smaller arc, usually 0.4 times the throat radius

        R_pintle = self.geometric_props['R_c']/4  # Pintle radius
        L_pintle = self.L_c/3  # Pintle length

    def Bartzcalc(self):
        D_t = 2*self.R_t
        Pr = 4*self.CEArun()['gamma']/(9*self.CEArun()['gamma']-5) # prandtl number given in Bartz paper
        c_t = np.sqrt((1/self.CEArun()['gamma'])*((self.CEArun()['gamma']+1)/2)**((self.CEArun()['gamma']+1)/(self.CEArun()['gamma']-1))*R*self.CEArun()['T_0']) # characteristic velocity
        T_w = self.CEArun()['T_0']                # assume for now
        sigma = 1/((1/2)*(self.T_w/self.CEArun()['T_0'])*(1+((self.CEArun()['gamma']-1)/2)*self.bartz_props['M']**2)+1/2)**(0.8-(self.bartz_props['w']/5))*(1+((self.CEArun()['gamma']-1)/2)*self.bartz_props['M']**2)**(self.bartz_props['w']/5) # dimensionless factor
    
    def values(self):
        Cv = self.CEArun()['Cp']/self.CEArun()['gamma']
        R = self.CEArun()['Cp']-Cv              # Specific gas constant
        
        return Cv, R
        
    def plot(self):
        
        Cv, R, m_dot, A_t, eps_c, V_c, L_c, R_t, R_e, R_b, R_s, R_pintle, L_pintle, Pr, c_t, T_w, sigma = self.values()
        
        def atan(angle):
            return np.tan(angle*np.pi/180)
        def asin(angle):
            return np.sin(angle*np.pi/180)
        def acos(angle):
            return np.cos(angle*np.pi/180)

        def circle(x, R):    
            return (-np.sqrt(R**2-x**2)+(R+R_t))
        
        def pintle_end(x):
            return np.sqrt(R_pintle**2-(x-(chamber_end+L_pintle))**2)

        y0 = R_t
        x0 = 0

        y1 = R_e
        x1 = -(y1 - y0 + x0*atan(self.geometric_props['Theta_e']))/atan(self.geometric_props['Theta_e'])

        ym1 = self.geometric_props['R_c']
        xm1 = -(ym1 - y0 + x0*atan(self.geometric_props['Theta_c']))/atan(self.geometric_props['Theta_c'])

        ym2 = ym1
        xm2 = -L_c +(ym2 - ym1 + xm1*atan(self.geometric_props['Theta_c']))/atan(self.geometric_props['Theta_c'])

        xb = np.linspace(-R_b*asin(self.geometric_props['Theta_c']),0,50)  # Horizontal limits of big arc as function of Theta_c
        xs = np.linspace(0,R_s*asin(-self.geometric_props['Theta_e']),50)  # Horizontal limits of small arc as function of Theta_e

        # Combustion chamber walls
        linex0 = np.linspace(xm2-((ym1+R_b-(R_b*acos(self.geometric_props['Theta_c']))-ym2)/asin(self.geometric_props['Theta_c'])), 
                             xm1-((ym1+R_b-(R_b*acos(self.geometric_props['Theta_c']))-ym2)/asin(self.geometric_props['Theta_c'])), 10)
        liney0 = np.linspace(ym2, 
                             ym1, 10)
        chamber_end = xm2-((ym1+R_b-(R_b*acos(self.geometric_props['Theta_c']))-ym2)/asin(self.geometric_props['Theta_c']))

        # Converging nozzle
        linex1 = np.linspace(xm1-R_b*asin(self.geometric_props['Theta_c'])+((ym1+R_b-(R_b*acos(self.geometric_props['Theta_c']))-ym2)/atan(self.geometric_props['Theta_c'])), 
                             x0-R_b*asin(self.geometric_props['Theta_c']), 10)
        liney1 = np.linspace(ym2, 
                             y0+(R_b-R_b*acos(self.geometric_props['Theta_c'])), 10)

        # Diverging nozzle
        linex2 = np.linspace(x0+R_s*asin(-self.geometric_props['Theta_e']), 
                             x1+R_s*asin(-self.geometric_props['Theta_e']), 10)
        liney2 = np.linspace(y0+(R_s-R_s*acos(self.geometric_props['Theta_e'])),
                             y1+(R_s-R_s*acos(self.geometric_props['Theta_e'])), 10)

        # Spray angle
        spray_angle = abs(90-self.geometric_props['spray_angle']) # Correction
        sprayx = np.linspace(chamber_end+L_pintle,chamber_end+L_pintle+atan(spray_angle)*(self.geometric_props['R_c']-R_pintle),10)
        sprayy = np.linspace(R_pintle, self.geometric_props['R_c'],10)

        
        ###### Bartz
        
        D_t = 2*R_t
        def Bartz(y_values):
            h = ((0.026/(D_t**0.2))*((self.bartz_props['mu']**0.2*self.CEArun()['Cp'])/(Pr**0.6))*(self.combustion_props['p_c']/c_t)**0.8*(D_t/R_b)**0.1)*(A_t/(np.pi*y_values**2))**0.9*sigma
            return h
        
        
        ##### Shared x-axis plots
        
        fig, (ax_comb, ax_h) = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

        ax_comb.plot(linex0, liney0, 'k-')
        ax_comb.plot(linex1, liney1, 'k-')
        ax_comb.plot(linex2, liney2, 'k-')
        ax_comb.plot(linex0, -liney0, 'k-')
        ax_comb.plot(linex1, -liney1, 'k-')
        ax_comb.plot(linex2, -liney2, 'k-')

        ax_comb.plot(xb,circle(xb,R_b), 'k-')
        ax_comb.plot(xs,circle(xs,R_s), 'k-')
        ax_comb.plot(xb,-circle(xb,R_b), 'k-')
        ax_comb.plot(xs,-circle(xs,R_s), 'k-')

        ax_comb.plot(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50),
                 pintle_end(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50)), 'k-')
        ax_comb.plot(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50),
                 -pintle_end(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50)), 'k-')

        ax_comb.plot(sprayx,sprayy,'b--')
        ax_comb.plot(sprayx,-sprayy,'b--')

        ax_comb.hlines(y=R_pintle, xmin = chamber_end, xmax = chamber_end+L_pintle, color = 'black')
        ax_comb.hlines(y=-R_pintle, xmin = chamber_end, xmax = chamber_end+L_pintle, color = 'black')

        ax_comb.vlines(x=chamber_end, ymin=-ym2, ymax=ym2, color = 'black')
        ax_comb.vlines(x=chamber_end+L_pintle, ymin = -R_pintle, ymax = R_pintle, color = 'black')
        ax_comb.set_ylabel('Distance $[m]$', fontsize = 18)
        
        ax_h.plot(linex0, Bartz(liney0), 'k-')
        ax_h.plot(linex1, Bartz(liney1), 'k-')
        ax_h.plot(linex2, Bartz(liney2), 'k-')
        ax_h.plot(xb,Bartz(circle(xb,R_b)), 'k-')
        ax_h.plot(xs,Bartz(circle(xs,R_s)), 'k-')
        ax_h.set_xlabel('Distance $[m]$', fontsize = 18)
        ax_h.set_ylabel('Heat flux $[W/m^2]$', fontsize = 18)
        plt.show(block=False)
    
    def print_values(self):

        Cv, R, m_dot, A_t, eps_c, V_c, L_c, R_t, R_e, R_b, R_s, R_pintle, L_pintle, Pr, c_t, T_w, sigma =self.values()

        print(f'Cv                        = {Cv}')
        print(f'R                         = {R}')
        print(f'm_dot                     = {m_dot}')
        print(f'Throat area               = {A_t}')
        print(f'compression ratio         = {eps_c}')
        print(f'Combustion chamber volume = {V_c}')
        print(f'Combustion chamber length = {L_c}')
        print(f'Throat radius             = {R_t}')
        print(f'Exit radius               = {R_e}')
        print(f'Pintle radius             = {R_pintle}')
        print(f'Pintle length             = {L_pintle}')
        print()
        print(f'~~~~~ Bartz equation ~~~~~')
        print(f'Prandtl number            = {Pr}')
        print(f'Characteristic velocity   = {c_t}')
        print(f'Sigma                     = {sigma}')

    def dict_values(self):
        Cv, R, m_dot, A_t, eps_c, V_c, L_c, R_t, R_e, R_b, R_s, R_pintle, L_pintle, Pr, c_t, T_w, sigma = self.values()
        dictionary = {
            'Cv': Cv,
            'R': R,
            'm_dot': m_dot,
            'Throat area': A_t,
            'compression ratio': eps_c,
            'Combustion chamber volume': V_c,
            'Combustion chamber length': L_c,
            'Throat radius': R_t,
            'Exit radius': R_e,
            'Pintle radius': R_pintle,
            'Pintle length': L_pintle,
            'Prandtl number': Pr,
            'Characteristic velocity': c_t,
            'Sigma': sigma
        }
        return dictionary