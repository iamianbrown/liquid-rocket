#Don't be scared... 
from TCA import *
from Chamber_Size_Test import *
from Injector_Code_Test import *
import numpy as np

class graph:
    def __init__(self, thing):
        self.thing = thing
    
    def graphing(self):
        list1 = []
        list2 = []
        geo = {
        'R_c': .25,
        'L_characteristic': 1
        }
        bartz = {
        'mu': .2,
        'm': .2,
        'w': .2,
        'M': .2
        }

        for i in range(300,320):
            list1.append(TCA('JetA', 'LOX', 400, 2.1, i))
        
        for j in list1:
            list2.append(Chamber(list1[j],))
        
        return(list2)
           
        
    






