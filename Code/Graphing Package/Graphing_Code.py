#Don't be scared... 
from pickle import NONE
from TCA import *
from Chamber_Size_Test import *
from Injector_Code_Test import *
import matplotlib.pyplot as plt
import numpy as np

class graph:
    def __init__(self, ENGINEobj=None, CHAMBERobj=None, INJobj=None):
        
        self.ENGINEobj = ENGINEobj
        self.CHAMBERobj = CHAMBERobj
        self.INJobj = INJobj

        