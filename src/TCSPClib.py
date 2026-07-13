import os
import numpy as np
import pandas as pd
import tkinter as tk
try:
    from . import deflib1 as deflib
except ImportError:
    import deflib1 as deflib
try:
    from . import export1 as exportlib
except ImportError:
    import export1 as exportlib
try:
    from . import mathlib3 as mathlib
except ImportError:
    import mathlib3 as mathlib
import matplotlib.pyplot as plt

class TCSPCprocessor:
    def __init__(self, parentframe, filepath, savedir):
        self.parentframe = parentframe
        self.filepath_str = filepath
        self.filepath = self.filepath_str.get()
        self.savedir_str = savedir
        if self.savedir_str is None:
            self.savedir = f'{self.filepath}_processed'
        else:
            self.savedir = self.savedir_str.get()
        self.data = None
        self.figs = []
        self.axis = []
        self.timeaxis = [0, 1]
        self.wavelengthaxis = [0, 1]
        self.tresdata = np.array([[0, 0], [0, 0]])
    
    def build_frame(self):
        self.plotframe = tk.Frame(self.parentframe, width=400, height=200, borderwidth=5, relief="ridge")
        self.plotframe.pack(fill='both', expand=True)
        self.load_button = tk.Button(self.plotframe, text="Load TCSPC Data", command=self.load_tcspc)
        self.load_button.pack(pady=20)
        self.plot_tres_linear_button = tk.Button(self.plotframe, text="Plot TRES Linear", command=self.plot_tres_linear)
        self.plot_tres_linear_button.pack(pady=10)
        self.plot_tres_log_button = tk.Button(self.plotframe, text="Plot TRES Log", command=self.plot_tres_log)
        self.plot_tres_log_button.pack(pady=10)
    
    def load_tcspc(self):
        # check if maindir is valid and contains ..TIME.., ..TRES.., ..WAVE.. files
        # check if filepath is a directory
        if not os.path.isdir(self.filepath_str.get()):
            print('Error: TCSPC main directory is not valid.')
            return
        parentdirfilenames = os.listdir(self.filepath_str.get())
        self.filepath = self.filepath_str.get()
        self.savedir = self.savedir_str.get()
        if not os.path.exists(self.savedir):
            # create savedir as subdir of parent directory
            self.savedir = os.path.join(self.filepath, 'processed_TCSPC')
        if not os.path.exists(self.savedir):
            os.makedirs(self.savedir)
        timefile = None
        wavefile = None
        tresfile = None
        for filename in parentdirfilenames:
            if 'TIME' in filename.upper():
                timefile = os.path.join(self.filepath, filename)
            elif 'WAVE' in filename.upper():
                wavefile = os.path.join(self.filepath, filename)
            elif 'TRES' in filename.upper():
                tresfile = os.path.join(self.filepath, filename)
        if timefile and wavefile and tresfile:
            self.timeaxis = loadtimefile(timefile)
            self.wavelengthaxis = loadwavefile(wavefile)
            self.tresdata = loadtresfile(tresfile)
            try: 
                self.tresdata = np.rot90(self.tresdata)  # rotate if possible
            except:
                pass
            # further processing can be done here
            print('TCSPC data processed successfully.')
        else: 
            print('Error: Required TCSPC files not found in the specified directory.')

    def plot_tres_linear(self):
        self.fig, self.ax = plt.subplots()
        # imshow tresdata as image with wavelength on y-axis and time on x-axis
        self.im = self.ax.imshow(self.tresdata, aspect='auto', extent=[self.wavelengthaxis[0], self.wavelengthaxis[-1], self.timeaxis[-1], self.timeaxis[0]], cmap='viridis')
        self.ax.set_xlabel('Wavelength (nm)')
        self.ax.set_ylabel('Time (ps)')
        self.fig.colorbar(self.im, ax=self.ax, label='Intensitiy (counts)')
        self.ax.set_title('TRES Data (Linear Scale)')
        self.fig.tight_layout()
        # save figure to savedir
        figpath = os.path.join(self.savedir, 'TRES_Linear.png')
        print('Saving TRES linear plot to:', figpath)
        self.fig.savefig(figpath, dpi=600)

        self.fig.show()
    
    def plot_tres_log(self):
        self.fig, self.ax = plt.subplots()
        # imshow tresdata as image with wavelength on y-axis and time on x-axis, log scale
        self.im = self.ax.imshow(np.log(self.tresdata+1), aspect='auto', extent=[self.wavelengthaxis[0], self.wavelengthaxis[-1], self.timeaxis[-1], self.timeaxis[0]], cmap='viridis')
        self.ax.set_xlabel('Wavelength (nm)')
        self.ax.set_ylabel('Time (ps)')
        self.fig.colorbar(self.im, ax=self.ax, label='Log Intensity (counts)')
        self.ax.set_title('TRES Data (Log Scale)')
        self.fig.tight_layout()
        
        # save figure to savedir
        figpath = os.path.join(self.savedir, 'TRES_Log.png')
        print('Saving TRES log plot to:', figpath)

        self.fig.savefig(figpath, dpi=600)

        self.fig.show()

#load time axis in ps from timefile
def loadtimefile(timefile):
	timeaxis = np.loadtxt(timefile, delimiter='\t', dtype=float)
	return timeaxis

# load wavelength axis in nm from wavefile
def loadwavefile(wavefile):
	wavelengths = np.loadtxt(wavefile, delimiter='\t', dtype=float)
	return wavelengths

# load TCSPC data from tresfile
def loadtresfile(tresfile):
	tres = np.loadtxt(tresfile, delimiter='\t', dtype=float)
	return tres