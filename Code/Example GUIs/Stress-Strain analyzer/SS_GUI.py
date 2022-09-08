import tkinter as tk
from tkinter import WORD, ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showerror
import os
from  matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.pyplot as plt
import pandas as pd
import SS_calculation


class fileSelector(tk.Tk):
    def __init__(self):
        super().__init__()
        self.geometry('1000x800')
        self.title('Stress-Strain Analysis')
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.data = None
        self.material_props = {'Young\'s Modulus (MPa)' : None, 'Yield Strength (MPa)' : None, 'Ultimate Strength (MPa)' : None, 'Fracture Stress (MPa)' : None}
        self.fig = None

        self.frame1 = ttk.Frame(self, height=200)
        self.frame2 = ttk.Frame(self)
        self.frame3 = ttk.Frame(self)

        #configure window
        self.grid_columnconfigure(0,w=4)
        self.grid_rowconfigure(0,w=1)
        self.grid_rowconfigure(1, w=1)
        self.grid_rowconfigure(2, w=1)

        # configure frame1
        self.frame1.grid(row=0, column=0, sticky='news', columnspan=2)
        self.frame1.grid_columnconfigure(0, w=1)
        self.frame1.grid_rowconfigure(0, w=1)
        self.frame1.grid_rowconfigure(1, w=1)
        self.frame1.grid_rowconfigure(2, w=1)


        # configure frame2
        self.frame2.grid(row=1, column=0, sticky='news')

        # configure frame3
        self.frame3.grid(row=2, column=0, sticky='news')

        # add buttons
        self.fileSelect = ttk.Button(self.frame1, text='Select Stress-Strain Data File', command=self.select_file)
        self.save_text = ttk.Button(self.frame1, text='Save Properties', command=self.write_params_to_csv)
        self.save_graph = ttk.Button(self.frame1, text='Save Image', command=self.save_image)

        # position buttons
        self.fileSelect.grid(pady=5, row=0, column=0, sticky='news')
        self.save_text.grid(pady=5, row=1, column=0, sticky='news')
        self.save_graph.grid(pady=5, row=2, column=0, sticky='news')

        self.filename = None

    def on_closing(self):
        plt.close('all')
        self.after(250, self.destroy)

    def select_file(self):
        filetypes = (('Text Files', '*.txt'), ('All Files', '*.*'))
        dialog = fd.askopenfilename(title='Select a file', initialdir=os.getcwd(), filetypes=filetypes)
        if self.filename != '':
            self.filename = dialog
            self.plot_data()
            self.show_params()

    def plot_data(self):
        if self.filename != None:
            self.data = SS_calculation.read_data(self.filename)
            self.fig, ax_1 = SS_calculation.plot_SS(self.data, figsize=(5,5))
            fig, ax_2 = SS_calculation.plot_SS(self.data)
            canvas = FigureCanvasTkAgg(fig, master=self.frame2)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        else:
            showerror(title='Error', message='Please input a valid filename')
    
    def show_params(self):
        columns = ('Parameter', 'Value')
        self.table = ttk.Treeview(self.frame3, columns=columns, show='headings')

        self.table.heading('Parameter', text='Parameter')
        self.table.heading('Value', text='Value')
        self.table.insert('', 'end', values=('Young\'s Modulus', str(SS_calculation.calculate_E(self.data)) + ' MPa'))
        self.table.insert('', 'end', values=('Ultimate Strength', str(SS_calculation.calculate_ultimate(self.data)) + ' MPa'))
        self.table.insert('', 'end', values=('Yield Strength', str(SS_calculation.calc_yield(self.data)) + 'MPa'))
        self.table.insert('', 'end', values=('Fracture Stress', str(SS_calculation.calculate_fracture(self.data)) + 'MPa'))

        self.table.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    def write_params_to_csv(self):
        self.material_props = {
        'Young\'s Modulus (MPa)' : SS_calculation.calculate_E(self.data), 
        'Yield Strength (MPa)' : SS_calculation.calc_yield(self.data), 
        'Ultimate Strength (MPa)' : SS_calculation.calculate_ultimate(self.data), 
        'Fracture Stress (MPa)' : SS_calculation.calculate_fracture(self.data)
        }

        f = fd.asksaveasfilename(initialfile = 'Untitled.txt',
        defaultextension=".txt",filetypes=[("All Files","*.*"),("Text Documents","*.txt")])
        with open(f, 'w') as file:
            for item in self.material_props.items():
                file.write(item[0] + '   ' + str(item[1]) + '\n')
    
    def save_image(self):
        f = fd.asksaveasfilename(initialfile = 'Untitled.png',
        defaultextension=".png",filetypes=[("All Files","*.*"),("PNG File","*.png"),("JPEG File","*.jpg")])
        self.fig.savefig(f)

if __name__ == '__main__':
    root = fileSelector()
    root.mainloop()