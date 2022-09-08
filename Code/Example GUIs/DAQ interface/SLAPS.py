from telnetlib import GA
import tkinter as tk
from tkinter import ttk
from GUIBaseClasses import *
import time

class MainApp(tk.Tk):
    def __init__(self):
        super().__init__()

        self.starttime = time.time()
        self.style = DarkStyle()
        
        self.grid_rowconfigure(0, w=1)
        self.grid_rowconfigure(1, w=1)
        self.grid_columnconfigure(1, w=1)

        self.temp = Info(self, label='Temperature', plottable=True)
        self.press = Info(self, label='Pressure', plottable=True)
        self.graphNotebook = ttk.Notebook(self)

        self.temp.grid(row=0, column=0, sticky='news')
        self.graphNotebook.grid(row=0, column=1, rowspan=2, sticky='news')
        self.press.grid(row=1, column=0, sticky='news')

        # add TCs to temp
        self.TC1 = TC(self, {'Pin': 'AIN0'})
        self.TC2 = TC(self, {'Pin': 'AIN1'})
        self.temp.addrow('TC1', self.TC1)
        self.temp.addrow('TC2', self.TC2)

        # add gauges to press
        self.pirani = Gauge(self, {'Relay': 'F101', 'Output':'AIN2'}, 'Lesker Pirani')
        self.press.addrow('Lesker Pirani', self.pirani)
    
        # make canvassed plot
        self.tempPlot = CanvasedPlot(self, [self.TC1, self.TC2], ['TC1', 'TC2'], ylabel='Temperature', notebook=self.graphNotebook, tablabel='Temperature')
        self.pressPlot = CanvasedPlot(self, [self.pirani], ['Pirani Gauge'], ylabel='Pressure', notebook=self.graphNotebook, tablabel='Pressure')

if __name__ == '__main__':
    m = MainApp()
    m.mainloop()