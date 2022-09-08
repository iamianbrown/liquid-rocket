import tkinter as tk
from tkinter import DISABLED, NORMAL, ttk
from tkinter import messagebox
import clickRecorder
import pyautogui
import pandas as pd
from pynput import mouse

class exporterGUI(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title('Button Sequencer')
        
        self.grid_columnconfigure(0,w=1)
        self.grid_rowconfigure(0,w=1)
        frame = tk.Frame(self)
        frame.grid(row=0, column=0, sticky='news')
        frame.grid_columnconfigure(0,w=1)
        frame.grid_columnconfigure(1,w=1)
        frame.grid_columnconfigure(2,w=1)
        frame.grid_rowconfigure(1,w=1)
        frame.grid_rowconfigure(0,w=1)
        self.Record = ttk.Button(frame, text='Record', command=self.record, state=NORMAL)
        self.Record.grid(row=1, column=0, sticky='ew')
        self.Stop = ttk.Button(frame, text='Stop', command=self.stop, state=DISABLED)
        self.Stop.grid(row=1, column=1, sticky='ew')
        self.Play = ttk.Button(frame, text='Play', command=self.play, state=DISABLED)
        self.Play.grid(row=1,column=2,sticky='ew')
        self.inputLabel = tk.Label(frame, text='Number of Loops:')
        self.inputLabel.grid(row=0, column=0,sticky='e')
        self.entry = tk.Entry(frame)
        self.entry.grid(row=0, column=1, sticky='ew')

    def record(self):
        clickRecorder.recording = True
        self.Stop['state'] = NORMAL
    
    def stop(self):
        if clickRecorder.recording == True:
            clickRecorder.write_clicks(clickRecorder.clicks) #write clicks to a .csv file
            clickRecorder.clicks.clear()
            clickRecorder.recording = False
        self.Stop['state'] = DISABLED
        self.Play['state'] = NORMAL

    def play(self):
        clicks = pd.read_csv('mouseClicks.csv', names=['x','y'])
        try:
            for i in range(int(self.entry.get())):
                for index, row in clicks.iterrows():
                    pyautogui.click(row['x'], row['y'])
        except:
            messagebox.showerror("Please input a number!")

if __name__ == '__main__':
    #start listening for clicks
    with mouse.Listener(on_click=clickRecorder.on_click) as listener:
        root = exporterGUI()
        root.mainloop()