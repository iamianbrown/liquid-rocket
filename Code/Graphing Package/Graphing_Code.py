#Don't be scared... 
from TCA import *
from Chamber_Size_Test import *
from Injector_Code_Test import *
import matplotlib.pyplot as plt
import numpy as np

class Engine:
    def __init__(self, COMBUSTIONobj=None, CHAMBERobj=None, INJobj=None):
        
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
            step2 = Chamber(TCAobj=step1, geometric_props=self.CHAMBERobj.geometric_props)
            Dictvalues = step2.geocalc()
            return(Dictvalues[dependent])
        
        elif dependent in self.COMBUSTIONobj.bartz_props.keys():
            Dictvalues = self.COMBUSTIONobj.bartzcalc()
            return(Dictvalues[dependent])
        
        elif dependent in self.INJobj.injector_props.keys():
            step1 = Combustion(tca_props=self.COMBUSTIONobj.tca_props, bartz_props=self.COMBUSTIONobj.bartz_props)
            step2 = Chamber(TCAobj=step1, geometric_props=self.CHAMBERobj.geometric_props)
            step3 = Injector(TCAobj=step1, CHAMBERobj=step2, injector_props=self.INJobj.injector_props)
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

        