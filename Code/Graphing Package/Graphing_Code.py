#Don't be scared... 
from TCA import *
from Chamber_Size_Test import *
from Injector_Code_Test import *
import matplotlib.pyplot as plt
import numpy as np
np.seterr(invalid='ignore')

class Engine:
    def __init__(self, COMBUSTIONobj, CHAMBERobj, INJobj):
        
        '''
        This class compiles all of the objects into one to create an overall Engine.
            It houses the combustion, chamber, and injector paramters all into one.
        This class was made to make the graphing process easier.
        '''

        self.COMBUSTIONobj = COMBUSTIONobj
        self.CHAMBERobj = CHAMBERobj
        self.INJobj = INJobj

    def getind(self, independent):
 
        '''
        This function matches the independent variable to the string passed in by the user
        '''

        if independent in self.COMBUSTIONobj.tca_props.keys():
            return(self.COMBUSTIONobj.tca_props[independent])
        
        elif independent in self.CHAMBERobj.geometric_props.keys():
            return(self.CHAMBERobj.geometric_props[independent])
        
        elif independent in self.COMBUSTIONobj.bartz_props.keys():
            return(self.COMBUSTIONobj.bartz_props[independent])
        
        elif independent in self.INJobj.injector_props.keys():
            return(self.INJobj.injector_props[independent])

    def getdep(self, dependent):

        '''
        This function just reruns the calculations for the dependent variable that the user inputs
        '''

        if dependent in self.COMBUSTIONobj.tca_props.keys():
            Dictvalues = self.COMBUSTIONobj.CEArun()
            return(Dictvalues[dependent])
        
        elif dependent in self.CHAMBERobj.geometric_props.keys():
            step1 = Combustion(tca_props=self.COMBUSTIONobj.tca_props,bartz_props=self.COMBUSTIONobj.bartz_props)
            step2 = Chamber(COMBUSTIONobj=step1, geometric_props=self.CHAMBERobj.geometric_props)
            Dictvalues = step2.geocalc()
            return(Dictvalues[dependent])
        
        elif dependent in self.COMBUSTIONobj.bartz_props.keys():
            Dictvalues = self.COMBUSTIONobj.bartzcalc()
            return(Dictvalues[dependent])
        
        elif dependent in self.INJobj.injector_props.keys():
            step1 = Combustion(tca_props=self.COMBUSTIONobj.tca_props, bartz_props=self.COMBUSTIONobj.bartz_props)
            step2 = Chamber(COMBUSTIONobj=step1, geometric_props=self.CHAMBERobj.geometric_props)
            step3 = Injector(COMBUSTIONobj=step1, CHAMBERobj=step2, injector_props=self.INJobj.injector_props)
            Dictvalues = step3.sizingcalc()
            return(Dictvalues[dependent])

    def changeVal(self, independent, value):

        '''
        This function iteratively changes the value of the independent variable using the linspace in plotparams()
        '''

        if independent in self.COMBUSTIONobj.tca_props.keys():
            self.COMBUSTIONobj.tca_props[independent] = value
        
        elif independent in self.CHAMBERobj.geometric_props.keys():
            self.CHAMBERobj.geometric_props[independent] = value
        
        elif independent in self.COMBUSTIONobj.bartz_props.keys():
            self.COMBUSTIONobj.bartz_props[independent] = value
        
        elif independent in self.INJobj.injector_props.keys():
            self.INJobj.injector_props[independent] = value
        
    def plotparams(self, ind, dep, minVal, maxVal, step):

            '''
                This is the method you call in the main.py file when you want to graph two variables against each other.
                
                You must pass in an independent variable to graph along with its dependent variable in string format, as well as
                a minumum to maximum value with a defined step to graph the independent variable between
            '''

            list_ind = []
            list_dep = []
            lin = np.linspace(minVal, maxVal, step)

            for j in lin:
                self.changeVal(ind,j)
                list_dep.append(self.getdep(dep))
                list_ind.append(self.getind(ind))
                
            plt.plot(list_ind, list_dep, '-g')
            plt.xlabel(ind)
            plt.ylabel(dep)
            plt.title(dep + ' as a function of ' + ind)
            plt.grid()
            plt.show()

    def engineVisual(self):
        
        ## FUNCTIONS
        def atan(angle):
            return np.tan(angle*np.pi/180)
        def asin(angle):
            return np.sin(angle*np.pi/180)
        def acos(angle):
            return np.cos(angle*np.pi/180)

        def circle(x,R):    # Function for throat arcs with positions (0, R+R_t)
            return (-np.sqrt(R**2-x**2)+(R+R_t))

        def pintle_end(x):
            return np.sqrt(R_pintle**2-(x-(chamber_end+L_pintle))**2) 

        #### GEOMETRICAL VARIABLES
        R_c = self.CHAMBERobj.geometric_props['R_c']
        L_c = self.CHAMBERobj.geometric_props['L_c'] / 100
        R_t = self.CHAMBERobj.geometric_props['R_t'] / 10
        R_e = self.CHAMBERobj.geometric_props['R_e'] / 10
        
        Theta_c = 35      # Angle for converging nozzle [deg]
        Theta_n = 33     # Angle for diverging nozzle [deg]
        Theta_e = 7      # Angle for exit of nozzle [deg]
        e_ratio = self.COMBUSTIONobj.tca_props['eps']

        R_b = self.CHAMBERobj.geometric_props['R_b'] / 10
        R_s = self.CHAMBERobj.geometric_props['R_s'] / 10

        R_pintle = self.INJobj.injector_props['R_p']
        L_pintle = self.INJobj.injector_props['L_pintle'] / 100

        spray_angle = self.INJobj.injector_props['theta_c']  # [deg]

        def parabolic_nozzle(t):

            Nx = 0.382 * acos(Theta_n-90)
            Ny = 0.382 * asin(Theta_n-90) + 1.382

            Ex = self.CHAMBERobj.geometric_props['L_characteristic']*((np.sqrt(e_ratio)-1)/atan(15))
            Ey = np.sqrt(e_ratio)

            m1 = atan(Theta_n)
            m2 = atan(Theta_e)

            C1 = Ny - m1*Nx
            C2 = Ey - m2*Ex

            Qx = (C2-C1)/(m1-m2)
            Qy = (m1*C2 - m2*C1)/(m1-m2)

            x_par = R_t*((1-t)**2*Nx + 2*(1-t)*t*Qx + t**2*Ex)
            y_par = R_t*((1-t)**2*Ny + 2*(1-t)*t*Qy + t**2*Ey)

            return x_par, y_par

        def entrance(theta):
            x_par = 1.5 * R_t * acos(theta)
            y_par = 1.5 * R_t * asin(theta) + 2.5 * R_t
            return x_par, y_par

        def exit(theta):
            x_par = 0.382 * R_t * acos(theta)
            y_par = 0.382 * R_t * asin(theta) + 1.382 * R_t
            return x_par, y_par

        x_entrance_vals = []
        y_entrance_vals = []
        entrance_steps = np.linspace(-Theta_c-90, -90, 50)

        for i in entrance_steps: 
            x_par, y_par = entrance(i)[0],entrance(i)[1]
            x_entrance_vals.append(x_par)
            y_entrance_vals.append(y_par)

        x_exit_vals = []
        y_exit_vals = []
        exit_steps = np.linspace(-90, Theta_n - 90, 50)

        for i in exit_steps: 
            x_par, y_par = exit(i)[0],exit(i)[1]
            x_exit_vals.append(x_par)
            y_exit_vals.append(y_par)

        x_noz_vals = []
        y_noz_vals = []

        t_steps = np.linspace(0,1,50)

        for i in t_steps: 
            x_par, y_par = parabolic_nozzle(i)[0],parabolic_nozzle(i)[1]
            x_noz_vals.append(x_par)
            y_noz_vals.append(y_par)

        y0 = R_t
        x0 = 0

        y1 = R_e
        x1 = -(y1 - y0 + x0*atan(Theta_n))/atan(Theta_n)

        ym1 = R_c
        xm1 = -(ym1 - y0 + x0*atan(Theta_c))/atan(Theta_c)

        ym2 = ym1
        xm2 = -L_c +(ym2 - ym1 + xm1*atan(Theta_c))/atan(Theta_c)

        xb = np.linspace(-R_b*asin(Theta_c),0,50)  # Horizontal limits of big arc as function of Theta_c
        xs = np.linspace(0,R_s*asin(-Theta_n),50)  # Horizontal limits of small arc as function of Theta_n

        # Combustion chamber walls
        linex0 = np.linspace(xm2-((ym1+R_b-(R_b*acos(Theta_c))-ym2)/asin(Theta_c)), 
                            xm1-((ym1+R_b-(R_b*acos(Theta_c))-ym2)/asin(Theta_c)), 10)
        liney0 = np.linspace(ym2, 
                            ym1, 10)
        chamber_end = xm2-((ym1+R_b-(R_b*acos(Theta_c))-ym2)/asin(Theta_c))

        # Converging nozzle
        linex1 = np.linspace(xm1-R_b*asin(Theta_c)+((ym1+R_b-(R_b*acos(Theta_c))-ym2)/atan(Theta_c)), 
                            x0-R_b*asin(Theta_c), 10)
        liney1 = np.linspace(ym2, 
                            y0+(R_b-R_b*acos(Theta_c)), 10)

        # Spray angle
        spray_angle = abs(90-spray_angle) # Correction
        sprayx = np.linspace(chamber_end+L_pintle,chamber_end+L_pintle+atan(spray_angle)*(R_c-R_pintle),10)
        sprayy = np.linspace(R_pintle, R_c, 10)

        # Plot
        plt.figure(figsize=(14, 12))

        plt.plot(x_entrance_vals, y_entrance_vals, 'k-')
        plt.plot(x_exit_vals, y_exit_vals, 'k-')
        plt.plot(x_noz_vals,y_noz_vals, 'k-')
        plt.plot(x_entrance_vals, np.negative(y_entrance_vals), 'k-')
        plt.plot(x_exit_vals, np.negative(y_exit_vals), 'k-')
        plt.plot(x_noz_vals,np.negative(y_noz_vals), 'k-')

        plt.plot(linex0, liney0, 'k-')
        plt.plot(linex1, liney1, 'k-')
        plt.plot(linex0, -liney0, 'k-')
        plt.plot(linex1, -liney1, 'k-')

        plt.plot(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50),
                pintle_end(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50)), 'k-')
        plt.plot(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50),
                -pintle_end(np.linspace(chamber_end+L_pintle, chamber_end+L_pintle+R_pintle, 50)), 'k-')
        plt.plot(sprayx,sprayy,'b--')
        plt.plot(sprayx,-sprayy,'b--')

        plt.hlines(y=R_pintle, xmin = chamber_end, xmax = chamber_end+L_pintle, color = 'black')
        plt.hlines(y=-R_pintle, xmin = chamber_end, xmax = chamber_end+L_pintle, color = 'black')

        plt.vlines(x=chamber_end, ymin=-ym2, ymax=ym2, color = 'black')
        plt.vlines(x=chamber_end+L_pintle, ymin = -R_pintle, ymax = R_pintle, color = 'black')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
