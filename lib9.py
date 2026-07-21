import numpy as np
from PIL import Image
from matplotlib.figure import Figure
import sys, pickle, copy
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
#from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.widgets import Cursor
from matplotlib.widgets import Button
from matplotlib.ticker import AutoLocator
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from scipy.special import wofz
from tkinter import filedialog as tkfd
import threading as thre
from concurrent.futures import ThreadPoolExecutor, as_completed
import mathlib3 as matl # type: ignore
import deflib1 as deflib # type: ignore
import PMclasslib1 as PMlib # type: ignore
import os, gc
import traceback
import error_handler as error_handler  # Centralized error handling and logging
import hsi_normalization as hsi_normalization  # HSI normalization module
import memory_tracker as memory_tracker  # Memory tracking and logging
import roihandler as roihandler 
import copy

SpectDataFloats = ['Slit Width (µm)', 'Central Wavelength (nm)',
                   'Cooling Temperature (°C)',
                   'Exposure Time (s)',
                   #'Wavelength First Pixel (nm)', # is obtained from WL
                   #'Wavelength Last Pixel (nm)',
                   'Delta Wavelength (nm)',
                   'x-position',
                   'y-position',
                   'z-position',
                   'Short Wavelength (nm)',
                   'Long Wavelength (nm)',
                   'magnification']

# todo: make the load data function in the SpectrumData class more robust for statistical cosmic removal

class SpectrumData:
    def __init__(self, filename, WL, BG, loadeachbg = False, linearbg=False, removecosmics=False, cosmicthreshold=20, cosmicpixels=3, removecosmicmethod=list(deflib.cosmicfuncts.keys())[0], WL_eV=None, default_dataset='Spectrum (PL-BG)'):
        self.removecosmicsmethod = removecosmicmethod
        self.loadeachbg = loadeachbg
        self.linearbg = linearbg
        self.removecosmics = removecosmics
        self.cosmicthreshold = cosmicthreshold
        self.cosmicpixels = cosmicpixels
        self.WL = WL
        self.WL_eV = WL_eV
        self.dofit = False # set True if fit is done
        if self.loadeachbg == True:
            self.BG = []
        else:
            self.BG = BG
        self.filename = filename
        self.readinkeys = ['BG', 'PL'] # WL is defined in XYMap
        self.default_dataset = default_dataset
        self.openFstate = [] # open floats state from metadata
        self.openDstate = [] # open Data from spectrometer
        self.dataokay = False    # set True if everything openend properly
        self.data = {}
        self.roistore = {}
        self.PL = []
        self.PLB = [] # PL-BG
        self.Specdiff1 = None # first derivative dPLB/dWL (None until calculated)
        self.Specdiff2 = None # second derivative d2PLB/dWL2 (None until calculated)
        self._read_file()

        # init fit data
        self.fwhm = np.nan
        self.fitmaxX = np.nan
        self.fitmaxY = np.nan
        self.fitdata = [None] # only fit parameter are stored, not the fit itself
        # what is inside fitdata after the fit is done? 
        # a sub array for each fit function (look at matl.fitkeys)
        # each array contains the fit output parameters (fitkeys[function][4], which is the number of parameters), followed by:
        # ['ss_res', 'ss_tot', 'r_squared', 'fwhm', 'pixstart', 'pixend', 'wlstart', 'wlend', 'max_x', 'max_y']
        #  fitmaxX[-1] fitmaxY[-2] aqpixstart[-3] aqpixend[-4] r_squared[-7] 'ss_res'[-8] 'ss_tot'[-9] 
        self.fitparams = matl.buildfitparas() # store all fit parameters
        self.fitparamunits = matl.buildfitparas()

    def _read_file(self):
        """
        Read spectrum data from file.
        
        Uses centralized error handling via ErrorEngine when available.
        Uses get_default_error_engine() to access shared error handler instance.
        Maintains dataokay flag to signal file read success/failure.
        """
        # Get default error engine (shared singleton instance)
        # This is safe to use from any thread/context
        error_engine = error_handler.get_default_error_engine()
        #default_decoding = 'utf-8'  # Default encoding for reading files

        try:
            with open(self.filename, 'r') as file:
                lines = file.readlines()
        except Exception as e:
            # Use ErrorEngine for file open errors
            error_engine.error(
                exception=e,
                context="Opening spectrum file",
                filename=self.filename,
                action="file_open"
            )
            # Maintain existing dataokay behavior
            self.dataokay = False
            return

        # Process lines to store variables
        startreaddata = False
        for line in lines:
            if ':' in line:
                key, value = map(str.strip, line.split(':', 1))
                if key in SpectDataFloats:
                    try:
                        self.data[key] = float(value)
                    except Exception as e:
                        # Use ErrorEngine for parsing errors
                        error_engine.warning(
                            message=f"Could not parse float value for key '{key}'",
                            context="Parsing spectrum metadata",
                            filename=self.filename,
                            key=key,
                            value=value
                        )
                        self.data[key] = value
                        self.openFstate.append(False)
                else:
                    self.data[key] = value
            elif '\t' in line:  # Data lines with tabs
                parts = line.split()
                if startreaddata == False:
                    count = 0
                    # start reading if at least two keys in the line
                    for i in parts:
                        if i in self.readinkeys:
                            count += 1
                    if count > 1:
                        startreaddata = True
                elif startreaddata == True:
                    try:
                        if self.loadeachbg == True:
                            #self.BG.append(int(parts[1]))
                            print("Warning: loadeachbg is True, but BG data is not being loaded. Check implementation.")
                        #self.WL.append(float(parts[0]))  WL is only read once by XYMap since each SpectrumData has the same WL-axis
                        self.PL.append(int(parts[2])) # type: ignore
                    except Exception as e:
                        # Use ErrorEngine for data parsing errors
                        error_engine.error(
                            exception=e,
                            context="Parsing spectrum data line",
                            filename=self.filename,
                            line_content=line.strip()
                        )
        if self.loadeachbg == True and self.linearbg == True:
            av = np.mean(self.BG)
            for i in range(len(self.BG)):
                self.BG[i] = av

        # Convert PL (and BG for loadeachbg) to float32 numpy arrays to reduce RAM usage
        self.PL = np.array(self.PL, dtype=np.float32)
        if self.loadeachbg == True:
            self.BG = np.array(self.BG, dtype=np.float32)

        try:
            self.PLB = np.subtract(self.PL, self.BG).astype(np.float32) # add PLB = PL-BG
        except Exception as e:
            # Use ErrorEngine for background subtraction errors
            error_engine.error(
                exception=e,
                context="Background subtraction",
                filename=self.filename,
                PL_length=len(self.PL),
                BG_length=len(self.BG)
            )
        # write openstate list
        for i in SpectDataFloats:
            if i not in list(self.data.keys()):
                self.openFstate.append(False)
        if len(self.WL) == 0:
            self.openDstate.append(None)
        if len(self.BG) == 0:
            self.openDstate.append(None)
        if len(self.PL) == 0:
            self.openDstate.append(None)
        self.setOK()
        # remove cosmic rays todo: add nearest neighbor method
        if self.removecosmics == True:
            try:
                self.PLB = np.asarray(deflib.cosmicfuncts[self.removecosmicsmethod](self.PLB, self.cosmicthreshold, self.cosmicpixels), dtype=np.float32)
            except Exception as e:
                # Use ErrorEngine for cosmic removal errors
                error_engine.error(
                    exception=e,
                    context="Cosmic ray removal",
                    filename=self.filename,
                    method=self.removecosmicsmethod,
                    threshold=self.cosmicthreshold
                )
        
        # if 
        
        # important: clean up! delete everything that is not needed anymore
        del lines

    def setOK(self):
        if False in self.openFstate:
            pass
        elif None in self.openDstate:
            pass
        else:
            self.dataokay = True

    # return a specific attribute of this class
    def get_attribute(self, attr_name:str):
        try:
            return getattr(self, attr_name)
        except AttributeError as e:
            print("Attribute {} not found in class SpectrumData.".format(attr_name))
    
    def __del__(self):
        """Destructor to clean up resources"""
        # Clear large data arrays to help garbage collector
        if hasattr(self, 'PL'):
            self.PL = []
        if hasattr(self, 'BG'):
            self.BG = []
        if hasattr(self, 'PLB'):
            self.PLB = []
        # Note: WL is shared reference, don't clear it


# create XY Map that contains the Pixels 
class XYMap:
    def __init__(self, fnames, cmapframe, specframe, loadbg=False, linearbg=False, removecosmics=False, cosmicthreshold=20, cosmicpixels=3, cosmicmethod=list(deflib.cosmicfuncts.keys())[0], defentries=deflib.defaults, derivative_polynomarray=[None, None, None, None], calc_norm_and_derive=None, calc_norm_on_intensity=None, skip_gui_build=False, default_dataset='Spectrum (PL-BG)'):
        self.defentries = defentries
        self.remcosmicfunc = cosmicmethod
        self.removecosmics = removecosmics
        self.linearbg = linearbg
        self.cosmicthreshold = cosmicthreshold
        self.cosmicpixels = cosmicpixels
        self.fnames = fnames
        self.readinkeys = ['WL']
        self.loadeachbg = loadbg
        self.specs = []                                                     # Spectral Objects
        self.fontsize = self.defentries['fontsize']                         # Default Plot Font Size
        self.default_dataset = default_dataset
        # init tk variables
        self.HSI_fit_useROI = tk.BooleanVar()                               # use ROI for fit
        self.HSI_fit_useROI.set(self.defentries['Fit_use_ROI_mask'])        # set default value for HSI_fit_useROI
        self.colormap = tk.StringVar()                                      # Colormap
        self.colormap.set(self.defentries['default_colormap'])              # set default colormap
        self.WL_selection = tk.StringVar()                                  # Wavelength Selection
        self.WL_values= deflib.WL_values
        self.WL_selection.set(self.WL_values[self.WL_values.index(defentries['selected_WL_axis'])]) # set default WL selection
        self.HSI_from_fitparam_useROI = tk.BooleanVar()          # use ROI for fit from fit parameters
        self.HSI_from_fitparam_useROI.set(self.defentries['HSI_from_fitparam_useROI'])
        # derivative settings: [0] first derivative bool, [1] second derivative bool, [2] polynomial order, [3] N fit points
        self.derivative_polynomarray = derivative_polynomarray               # derivative settings from main
        # normalization settings for derivatives
        self.calc_norm_and_derive = calc_norm_and_derive                    # normalize, then calc derivatives checkbox
        self.calc_norm_on_intensity = calc_norm_on_intensity                # normalize on total intensity checkbox
        
        self.fitkeys = matl.fitkeys
        self.countthreshv = self.defentries['colormap_threshold']
        
        # Store frames for GUI building
        self.cmapframe = cmapframe                                          # Colormap Frame
        self.specframe = specframe                                          # Spectrum Frame
        
        # Initialize data structures
        self.PMdict = {}                                                    # Pixel Matrix Dictionary
        self._hsi_counter = 0                                               # Global counter for unique HSI naming
        self.disspecs = {}                                                  # displayed spectra Objects contains [spec][wl]{metadata}
        self.allfpnames = matl.getlistofallFitparameters()
        self.allfpnamesinone = matl.getlistofallFitparaminone()
        self.speckeys = {'Wavelength axis': 'WL', 'Background (BG)': 'BG',
                         'Counts (PL)': 'PL', 'Spectrum (PL-BG)': 'PLB', # PLB is BG corrected (=PL-BG)
                         'first derivative': 'Specdiff1', # first derivative of PLB: d(PLB)/dWL
                         'second derivative': 'Specdiff2', # second derivative of PLB: d2(PLB)/dWL2
                         'first derivative (normalized)': 'Specdiff1_norm', # first derivative normalized by signal
                         'second derivative (normalized)': 'Specdiff2_norm', # second derivative normalized by signal
                         'first derivative (norm on intensity, then derive)': 'Specdiff1_norm_intensity', # normalize on total intensity, then derive
                         'second derivative (norm on intensity, then derive)': 'Specdiff2_norm_intensity', # normalize on total intensity, then derive
                         'first derivative (norm on counts, then derive)': 'Specdiff1_norm_counts', # normalize on counts, then derive
                         'second derivative (norm on counts, then derive)': 'Specdiff2_norm_counts' # normalize on counts, then derive
                         } 
        self.windowfunctions = list(matl.fitkeys.keys())                    # Window Functions from matl.fitkeys
        
        # Load data if files provided
        if len(fnames) > 0:
            self.load_data()
        else:
            # Initialize empty structures for loading from saved state
            self.WL = []
            self.WL_eV = []
            self.BG = []
            self.DataSpecMin = 0
            self.DataSpecMax = 1000
            self.DataSpecdL = 0.1
            self.DataPixSt = 0
            self.wlstart = self.defentries['lowest_wavelength']
            self.wlend = self.defentries['highest_wavelength']
            self.wlstart_eV = matl.arr_nm2ev(np.array([self.wlstart]))[0]  # Convert to eV
            self.wlend_eV = matl.arr_nm2ev(np.array([self.wlend]))[0]  # Convert to eV
            self.SpecDataMatrix = []
            self.gdx = 0
            self.gdy = 0
            self.DataPixDX = 0
            self.DataPixDY = 0
            self.aqpixstart = 0
            self.aqpixend = -1
            self.aqpixstart_eV = 0
            self.aqpixstart_eV = -1
            self.mxcoords = []
            self.mycoords = []
            self.PixAxX = []
            self.PixAxY = []
            self.PMmetadata = {}
        
        # Build GUI only if not skipped
        if not skip_gui_build:
            self.build_gui()
    
    def load_data(self):
        """Load spectral data from files"""
        self.loadfiles()
        
        # Initialize data ranges
        if len(self.WL) > 0 and len(self.specs) > 0:
            self.DataSpecMin = round(float(np.amin(self.WL)), 1)                             # Spectrum Start
            self.DataSpecMax = round(float(np.amax(self.WL[-1])), 1)                         # Spectrum End
            self.DataSpecdL = self.specs[0].data['Delta Wavelength (nm)']   # delta Lambda
        else:
            self.DataSpecMin = 0
            self.DataSpecMax = 1000
            self.DataSpecdL = 0.1
        
        self.wlstart = self.defentries['lowest_wavelength']                 # lowest wavelength
        self.wlend = self.defentries['highest_wavelength']                  # highest wavelength
        
        self.autogenmatrix()                                                # generate emty grid and fill Data obj into Matrix
        self.calculate_derivatives()                                        # calculate derivatives if requested
        
        # Initialize pixel data
        if len(self.specs) > 0:
            self.DataPixSt = len(self.specs[0].WL)                          # Number of WL Pixels
        else:
            self.DataPixSt = 0
        
        self.DataPixDX = self.gdx if hasattr(self, 'gdx') else 0           # Plot Pixel dx
        self.DataPixDY = self.gdy if hasattr(self, 'gdy') else 0           # Plot Pixel dy
        self.aqpixstart = 0                                                 # evaluate window start dL
        self.aqpixend = -1                                                  # evaluate window end dL
        
        # Set metadata if we have data
        if len(self.PMdict) > 0:
            self.PMmetadata['HSI0'] = {'wlstart': self.wlstart, 'wlend': self.wlend, 'countthresh': self.countthreshv, 'aqpixstart': self.aqpixstart, 'aqpixend': self.aqpixend, 'unit': '', 'quantatity': '', 'fitmodel':'', 'fitparameter': ''}
    
    def build_gui(self):
        """Build the GUI elements - separated from data loading"""
        self.SpecButtons = self.build_button_frame(self.specframe)          # build Spec GUI
        self.buildMinMaxSpec(self.cmapframe)                                # build CMAP grid GUI
        self.build_PixMatrix_frame(self.cmapframe)                          # build Pixel Matrix GUI
        self.buildselectboxes(self.cmapframe, list(self.speckeys.keys()))
        self.build_plot_optionsframe(self.cmapframe)
        
        self.updatewl()
        self.UpdateHSIselect()
    
    def build_plot_optionsframe(self, parframe):
        #self.plot_optionsframe = tk.Frame(parframe, border=5, relief="raised")
        #tk.Label(self.plot_optionsframe, text="Plot Options").grid(row=0, column=0)
        self.plotoptions_frame = tk.Frame(parframe, border=5, relief="sunken")
        self.plotoptions_frame.grid(row=0, column=5, rowspan=6, sticky=tk.NW)
        # test: add text
        tk.Label(self.plotoptions_frame, text="Plot Options").grid(row=0, column=0)

        # add normalizatize derivatives on signal: 
        # normalizes the derivative intensity to the PLB intensity at each wavelength to enhance small features

        self.normalize_HSI_var = tk.BooleanVar()                    # normalize derivatives to selected criteria
        self.normalize_HSI_var.set(self.defentries['do_normalize_HSI_var'])  # set default value for normalize_derivatives_on_PLB
        tk.Checkbutton(self.plotoptions_frame, text="Normalize HSI", variable=self.normalize_HSI_var).grid(row=1, column=0)

        # selct what to normalize to
        tk.Label(self.plotoptions_frame, text="Normalize to:").grid(row=2, column=0)
        
        # Normalization method selector
        self.normalization_methods = ['None', 'Integrated Counts', 'Max Intensity', 'Counts at Wavelength', 'Normalize intern']
        self.normalize_method_var = tk.StringVar()
        self.normalize_method_var.set(self.defentries.get('normalize_method', 'None'))
        self.normalize_method_box = ttk.Combobox(self.plotoptions_frame, textvariable=self.normalize_method_var, values=self.normalization_methods, width=18)
        self.normalize_method_box.grid(row=3, column=0)
        
        # Normalization wavelength parameter (for methods that need it)
        tk.Label(self.plotoptions_frame, text="Norm. WL (nm):").grid(row=4, column=0, sticky=tk.W)
        self.normalize_wl_var = tk.StringVar()
        self.normalize_wl_var.set(self.defentries.get('normalize_wavelength', '600'))
        self.normalize_wl_entry = tk.Entry(self.plotoptions_frame, textvariable=self.normalize_wl_var, width=10)
        self.normalize_wl_entry.grid(row=5, column=0, sticky=tk.W)
        
        # Normalization dataset selector - allows normalizing based on a different dataset
        # e.g., normalize derivatives by counts from the original Spectrum (PL-BG)
        tk.Label(self.plotoptions_frame, text="Norm. Dataset:").grid(row=6, column=0, sticky=tk.W)
        self.normalize_dataset_var = tk.StringVar()
        self.normalize_dataset_var.set(self.defentries.get('normalize_dataset', 'Spectrum (PL-BG)'))
        self.normalize_dataset_box = ttk.Combobox(self.plotoptions_frame, textvariable=self.normalize_dataset_var, 
                                                   values=list(self.speckeys.keys()), width=18)
        self.normalize_dataset_box.grid(row=7, column=0, sticky=tk.W)
        

    def buildselectboxes(self, frame, values):
        tk.Label(frame, text="Select Data Set".format(self.DataSpecMax)).grid(row=0, column=1)
        self.selectspecbox = ttk.Combobox(frame, values=values, width=40)
        # Use the default from defentries if available, otherwise fall back to 'Spectrum (PL-BG)'
        DEFAULT_DATA_SET = 'Spectrum (PL-BG)'
        default_dataset = self.defentries.get('data_set', self.default_dataset)
        if default_dataset in values:
            self.selectspecbox.set(default_dataset)
        elif self.default_dataset in values:
            self.selectspecbox.set(DEFAULT_DATA_SET)  # Hardcoded fallback if defentries value not found
        self.selectspecbox.grid(row=1, column=1)
        tk.Label(frame, text="Select Colormap".format(self.DataSpecMax)).grid(row=0, column=2)
        self.selectcolmapbox = ttk.Combobox(frame, values=plt.colormaps(), textvariable=self.colormap)
        self.selectcolmapbox.set(self.defentries['default_colormap'])  # set default colormap from defentries
        self.selectcolmapbox.grid(row=1, column=2) 

    def updatecountthresh(self):
        try:
            self.countthreshv = int(self.countthresh.get())
        except Exception as e:
            print("Error", '{} Insert valid threshold of type int.'.format(str(e)))

    # Spectral Plot Input Update
    def updatewl(self):
        try:
            self.wlstart = round(float(self.proc_spec_min.get()), 1)
            self.wlend = round(float(self.proc_spec_max.get()), 1)
        except Exception as e:
            print("Error", '{} Insert valid spectral Borders of type float.'.format(str(e)))
        passt = False
        if self.wlstart > self.wlend:
            print('ERROR', 'Lowest Wavelength must be smaller than Highest Wavelength! Reconsider Input!')
            self.wlstart = round(float(self.DataSpecMin), 1)
            self.proc_spec_min.delete(0, tk.END)
            self.proc_spec_min.insert(0, str(self.wlstart))
            self.wlend = round(float(self.DataSpecMax), 1)
            self.proc_spec_max.delete(0, tk.END)
            self.proc_spec_max.insert(0, str(self.wlend))
        elif self.wlstart < self.DataSpecMin:
            passt = True
            print('ERROR', 'Lowest Wavelength is below data WL. Set WL to lowest datapoint.')
            self.wlstart = round(float(self.DataSpecMin), 1)
            self.proc_spec_min.delete(0, tk.END)
            self.proc_spec_min.insert(0, str(self.wlstart))
            if self.wlend > self.DataSpecMax:
                print('ERROR', 'Lowest Wavelength is below data WL. Set WL to highest datapoint.')
                self.wlend = round(float(self.DataSpecMax), 1)
                self.proc_spec_max.delete(0, tk.END)
                self.proc_spec_max.insert(0, str(self.wlend))
        elif self.wlend > self.DataSpecMax:
            passt = True
            print('ERROR', 'Lowest Wavelength is below data WL. Set WL to highest datapoint.')
            self.wlend = round(float(self.DataSpecMax), 1)
            self.proc_spec_max.delete(0, tk.END)
            self.proc_spec_max.insert(0, str(self.wlend))
            if self.wlstart < self.DataSpecMin:
                print('ERROR', 'Lowest Wavelength is below data WL. Set WL to lowest datapoint.')
                self.wlstart = round(float(self.DataSpecMin), 1)
                self.proc_spec_min.delete(0, tk.END)
                self.proc_spec_min.insert(0, str(self.wlstart))
        else:
            passt = True
        if passt == True:
            # convert lambdamin and lambdamax into pixels
            self.aqpixstart = int((self.wlstart-self.DataSpecMin)/self.DataSpecdL)
            self.aqpixend = int((self.wlend-self.DataSpecMin)/self.DataSpecdL)
            # update wl boundaries in eV
            self.wlstart_eV = matl.arr_nm2ev(np.array([self.wlstart]))[0]
            self.wlend_eV = matl.arr_nm2ev(np.array([self.wlend]))[0]
            # convert evmin and evmax into pixels
            self.aqpixstart_eV = int((self.wlstart_eV-matl.arr_nm2ev(np.array([self.DataSpecMax]))[0])/matl.arr_nm2ev(np.array([self.DataSpecdL]))[0])
            self.aqpixend_eV = int((self.wlend_eV-matl.arr_nm2ev(np.array([self.DataSpecMax]))[0])/matl.arr_nm2ev(np.array([self.DataSpecdL]))[0])
            # for testing: print updated boundaries    
    def updateproc_spec_max(self):
        try:
            self.proc_spec_min.delete(0, tk.END)
            self.proc_spec_min.insert(0, str(self.wlstart))
            self.proc_spec_max.delete(0, tk.END)
            self.proc_spec_max.insert(0, str(self.wlend))
        except Exception as e:
            print("Error", '{} Insert valid spectral Borders of type float.'.format(str(e)))
        self.updatewl()

    # min and max wl can be inserted here for preceed window
    def buildMinMaxSpec(self, frame):
        # ProcSpecMin input field
        tk.Label(frame, text="Lowest wavelength \\ nm\nMinimum: {:.1f} nm".format(float(self.DataSpecMin))).grid(row=0, column=0)
        self.proc_spec_min = tk.Entry(frame)
        self.proc_spec_min.grid(row=1, column=0)
        self.proc_spec_min.insert(0, self.defentries['lowest_wavelength'])

        # ProcSpecMax input field
        tk.Label(frame, text="Highest wavelength \\ nm\nMaximum: {:.1f} nm".format(float(self.DataSpecMax))).grid(row=2, column=0)
        self.proc_spec_max = tk.Entry(frame)
        self.proc_spec_max.grid(row=3, column=0)
        self.proc_spec_max.insert(0, self.defentries['highest_wavelength'])

        # Button to create colormap
        tk.Label(frame, text="Press button below \nto update data matrix").grid(row=2, column=1)
        b1 = tk.Button(frame, text="Create intensity colormap", command= lambda: self.buildandPlotIntCmap(datatype=self.selectspecbox.get()))
        b1.grid(row=3, column=1)
        b2 = tk.Button(frame, text="Create spectral maximum colormap", command= lambda: self.buildandPlotSpecCmap())
        b2.grid(row=4, column=1)

        b3 = tk.Checkbutton(frame, text="Fit: use Selected ROI mask", variable=self.HSI_fit_useROI).grid(row=5, column=1)

        b4 = tk.Button(frame, text="Update spectral fit maxima", command= lambda: self.updateandPlotSpecCmap('fitmaxX'))
        b4.grid(row=6, column=1)

        # Plot Font size
        tk.Label(frame, text="Plot font size".format(self.DataSpecMin)).grid(row=2, column=2)
        self.CMFont = tk.Entry(frame)
        self.CMFont.grid(row=3, column=2)
        self.CMFont.insert(0, str(self.fontsize))

        # scipy fit maxiter
        tk.Label(frame, text="Max tries for fit").grid(row=4, column=2)
        self.fitmaxiter = tk.Entry(frame)
        self.fitmaxiter.grid(row=5, column=2)
        self.fitmaxiter.insert(0, self.defentries['maxfev'])

        # threshold for colormap
        tk.Label(frame, text="Colormap threshold \\ Counts").grid(row=4, column=0)
        self.countthresh = tk.Entry(frame)
        self.countthresh.grid(row=5, column=0)
        self.countthresh.insert(0, self.defentries['colormap_threshold'])

    def build_PixMatrix_frame(self, parframe):
        frame = tk.Frame(parframe, border=5, relief="sunken")
        frame.grid(row=0, column=4, rowspan=6, sticky=tk.NW)
        tk.Label(frame, text="Process Pixel Matrix".format(self.DataSpecMin)).grid(row=0, column=0)
        b1= tk.Button(frame, text="fit rotated 2D Gaussian to Matrix", command= lambda: self.fit2drotgausstopixmatrix())
        b1.grid(row=1, column=0)
        b2= tk.Button(frame, text="fit 2D Gaussian to Matrix", command= lambda: self.fit2dgausstopixmatrix())
        b2.grid(row=2, column=0)

    def build_roi_frame(self, parframe, PMname='HSI0'):
        if not hasattr(self, 'roihandler'):
            self.roihandler = roihandler.Roihandler(cmap=self.colormap.get())#self)
        self.roisel = tk.StringVar()

        frame = tk.Frame(parframe, border=5, relief="raised")
        frame.grid(row=0, column=5, rowspan=6, sticky=tk.NW)
        tk.Label(frame, text="Select ROI").grid(row=0, column=0)
        self.roiselgui = ttk.Combobox(frame, textvariable=self.roisel)
        self.roiselgui.grid(row=1, column=0)
        
        # assign roiselgui to roihandler
        self.roihandler.construct_roiselgui(self.roiselgui)
        
        # Populate ROI list from existing data if available - fixes empty combobox after loading
        if hasattr(self.roihandler, 'roilist'): #and len(self.roihandler.roilist) > 0: 
            # the and condition is not needed, since an empty list is also fine
            self.roiselgui['values'] = list(self.roihandler.roilist.keys())
            # but if the list is not empty, select the last one
            if len(self.roihandler.roilist) > 0:
                self.roiselgui.set(list(self.roihandler.roilist.keys())[-1])

        try:
            b1 = tk.Button(frame, text="ROI Editing last Selection", command= lambda: self.roihandler.construct(self.PMdict[self.hsiselect.get()].PixMatrix, self.roiselgui))
        except: # select first HSI
            print('HSI selection failed. Selecting first HSI.')
            self.hsiselect.set(str(list(self.PMdict.keys())[0]))
            b1 = tk.Button(frame, text="ROI Editing last Selection", command= lambda: self.roihandler.construct(self.PMdict[self.hsiselect.get()].PixMatrix, self.roiselgui))
        b1.grid(row=2, column=0)
        b2 = tk.Button(frame, text="plot ROI", command= lambda: self.roihandler.plotroi())
        b2.grid(row=3, column=0)
        # select a color for the roiplot
        roicolor = tk.StringVar()
        tk.Label(frame, text="ROI Plot Color").grid(row=4, column=0)
        roicolorselect = ttk.Combobox(frame, textvariable=roicolor, values=['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta', 'yellow', 'teal', 'navy', 'maroon', 'olive', 'lime', 'coral', 'indigo', 'gold', 'silver'])
        roicolorselect.grid(row=5, column=0)
        roicolorselect.set(self.defentries['roi_plot_color']) # set default ROI plot color
        # button plot roi with selected color on selected HSI
        b_roi_hsi = tk.Button(frame, text="Plot ROI on HSI", command= lambda: self.roihandler.plot_roi_on_pixmatrix(
            self.PMdict[self.hsiselect.get()].PixMatrix, 
            self.roisel.get(), 
            vis_type=self.roivistypeselect.get(), 
            color=roicolor.get(),
            title=self._get_hsi_title(self.hsiselect.get()), 
            unit=self.PMdict[self.hsiselect.get()].metadata['unit'] 
                if hasattr(self.PMdict[self.hsiselect.get()].metadata, 'unit')
                else ''
            )
        )
        b_roi_hsi.grid(row=6, column=0)
        # select roi plot type
        tk.Label(frame, text="ROI Plot Type").grid(row=7, column=0)
        linemodes = ['cornerlines', 'overlay']
        self.roivistypeselect = ttk.Combobox(frame, values=linemodes, width=15)
        selmode = self.defentries['roi_vis_type'] # remove , if selmode not in linemodes else linemodes[0]
        if selmode not in linemodes:
            selmode = linemodes[0]
        self.roivistypeselect.set(selmode) # set default ROI plot type
        self.roivistypeselect.grid(row=8, column=0)

        # entry to enter, what roi (indicees 1,2,3 should be plotted at once, use colors 1,2,3 for that)
        tk.Label(frame, text="Plot ROI indices (e.g. 1,2,3)").grid(row=9, column=0)
        self.roi_plot_indices_entry = tk.Entry(frame)
        self.roi_plot_indices_entry.grid(row=10, column=0)
        # button to plot multiple rois with different colors on selected HSI
        b_plot_multiple_rois = tk.Button(frame, text="Plot multiple ROIs on HSI", command= lambda: self.plot_multiple_rois_on_hsi())
        b_plot_multiple_rois.grid(row=11, column=0)

        b3 = tk.Button(frame, text="delete ROI", command= lambda: self.roihandler.delete_roi())
        b3.grid(row=13, column=0)     
    
        # build second clumn for ROI
        tk.Label(frame, text="Select HSI Image").grid(row=0, column=1)
        self.hsiselect = ttk.Combobox(frame, width=40)
        self.hsiselect.grid(row=1, column=1)
        b4 = tk.Button(frame, text="Plot HSI", command= lambda: self.plotPixelMatrix(self.hsiselect.get())) # original command
                       # test
                       #command = lambda: [self.plotPixelMatrix(i) for i in range(len(self.PMdict.keys()))]) # test command
        b4.grid(row=2, column=1)
        b5 = tk.Button(frame, text="Multiply HSI with ROI", command= lambda: self.multiroitopixmatrix())
        b5.grid(row=3, column=1)
        b6 = tk.Button(frame, text="Delete selected HSI", command= lambda: self.delHSI())
        b6.grid(row=4, column=1)

        b7 = tk.Button(frame, text="Save selected HSI", command= lambda: self.saveHSI())
        b7.grid(row=5, column=1)

        b8 = tk.Button(frame, text="Load HSI", command= lambda: self.loadHSI())
        b8.grid(row=6, column=1)

        b9 = tk.Button(frame, text="Export HSI to .csv", command= lambda: self.exportHSI())
        b9.grid(row=7, column=1)

        b10 = tk.Button(frame, text="Export HSI to .itx", command= lambda: self.exportHSI_itx())
        b10.grid(row=8, column=1)

        # build third column for Spectral Data
        tk.Label(frame, text="Select Spectral Data").grid(row=0, column=2)
        self.specselect = ttk.Combobox(frame, width=40)
        self.specselect.grid(row=1, column=2)
        
        # Populate Spectrum list from existing data if available - fixes empty combobox after loading
        if hasattr(self, 'disspecs') and len(self.disspecs) > 0:
            self.specselect['values'] = list(self.disspecs.keys())
            self.specselect.set(list(self.disspecs.keys())[-1])
            
        b7 = tk.Button(frame, text="Plot Spectrum", command= lambda: self.plotSpectral(self.specselect.get()))
        b7.grid(row=2, column=2)
        b8 = tk.Button(frame, text="Save Spectrum to .txt", command= lambda: self.saveSpectrum(self.specselect.get()))
        b8.grid(row=3, column=2)
        b9 = tk.Button(frame, text="average hsi to spectrum", command= lambda: self.averageHSItoSpecData())
        b9.grid(row=5, column=2)
        b10 = tk.Button(frame, text="Delete selected Spectral Data", command= lambda:
                          self.delSpecData(self.specselect.get()))
        b10.grid(row=4, column=2)
        
        # Add button for exporting all averaged spectra
        b_export_all = tk.Button(frame, text="Export All Averaged Spectra", 
                                 command=lambda: self.exportAllAveragedSpectra())
        b_export_all.grid(row=6, column=2)
        
        # add entry for a correction spectrum
        tk.Label(frame, text="Correction Spectrum").grid(row=7, column=2)

        # Button to select correction spectrum file
        b11 = tk.Button(frame, text="Select Correction Spec File", command=self.select_correction_spectrum_file)
        b11.grid(row=8, column=2)
        b12 = tk.Button(frame, text="Correct Spectrum", command= lambda: self.correctSpectrum(self.specselect.get()))
        b12.grid(row=9, column=2)
    
    def build_plot_options_frame(self, parframe):
        """Build GUI frame for HSI plot options with scale bar and advanced formatting."""
        frame = tk.Frame(parframe, border=5, relief="ridge")
        frame.grid(row=0, column=6, rowspan=6, sticky=tk.NW)
        
        # Title
        tk.Label(frame, text="HSI Plot Options", font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)
        
        # Colormap selection
        tk.Label(frame, text="Colormap:").grid(row=1, column=0, sticky=tk.W, padx=5)
        self.hsi_cmap_var = tk.StringVar(value=self.defentries.get('hsi_cmap', 'hot'))
        self.hsi_cmap_combo = ttk.Combobox(frame, textvariable=self.hsi_cmap_var, 
                                           values=['hot', 'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'gray', 'jet'],
                                           width=15)
        self.hsi_cmap_combo.grid(row=1, column=1, padx=5, pady=2)
        
        # vmin/vmax for color scale
        tk.Label(frame, text="vmin (leave empty for auto):").grid(row=2, column=0, sticky=tk.W, padx=5)
        self.hsi_vmin_var = tk.StringVar(value=str(self.defentries.get('hsi_vmin', '')))
        self.hsi_vmin_entry = tk.Entry(frame, textvariable=self.hsi_vmin_var, width=18)
        self.hsi_vmin_entry.grid(row=2, column=1, padx=5, pady=2)
        
        tk.Label(frame, text="vmax (leave empty for auto):").grid(row=3, column=0, sticky=tk.W, padx=5)
        self.hsi_vmax_var = tk.StringVar(value=str(self.defentries.get('hsi_vmax', '')))
        self.hsi_vmax_entry = tk.Entry(frame, textvariable=self.hsi_vmax_var, width=18)
        self.hsi_vmax_entry.grid(row=3, column=1, padx=5, pady=2)
        
        # Scale bar options
        tk.Label(frame, text="Scale bar length (μm):").grid(row=4, column=0, sticky=tk.W, padx=5)
        self.hsi_scalebar_len_var = tk.DoubleVar(value=self.defentries.get('hsi_scalebar_length', 20.0))
        self.hsi_scalebar_len_entry = tk.Entry(frame, textvariable=self.hsi_scalebar_len_var, width=18)
        self.hsi_scalebar_len_entry.grid(row=4, column=1, padx=5, pady=2)
        
        tk.Label(frame, text="Scale bar width (μm):").grid(row=5, column=0, sticky=tk.W, padx=5)
        self.hsi_scalebar_width_var = tk.DoubleVar(value=self.defentries.get('hsi_scalebar_width', 2.0))
        self.hsi_scalebar_width_entry = tk.Entry(frame, textvariable=self.hsi_scalebar_width_var, width=18)
        self.hsi_scalebar_width_entry.grid(row=5, column=1, padx=5, pady=2)
        
        tk.Label(frame, text="Scale bar position X (μm):").grid(row=6, column=0, sticky=tk.W, padx=5)
        self.hsi_scalebar_pos_x_var = tk.DoubleVar(value=self.defentries.get('hsi_scalebar_pos_x', 2.0))
        self.hsi_scalebar_pos_x_entry = tk.Entry(frame, textvariable=self.hsi_scalebar_pos_x_var, width=18)
        self.hsi_scalebar_pos_x_entry.grid(row=6, column=1, padx=5, pady=2)
        
        tk.Label(frame, text="Scale bar position Y (μm):").grid(row=7, column=0, sticky=tk.W, padx=5)
        self.hsi_scalebar_pos_y_var = tk.DoubleVar(value=self.defentries.get('hsi_scalebar_pos_y', 2.0))
        self.hsi_scalebar_pos_y_entry = tk.Entry(frame, textvariable=self.hsi_scalebar_pos_y_var, width=18)
        self.hsi_scalebar_pos_y_entry.grid(row=7, column=1, padx=5, pady=2)
        
        # Figure size
        tk.Label(frame, text="Figure width (inches):").grid(row=8, column=0, sticky=tk.W, padx=5)
        self.hsi_figsize_width_var = tk.DoubleVar(value=self.defentries.get('hsi_figsize_width', 7.0))
        self.hsi_figsize_width_entry = tk.Entry(frame, textvariable=self.hsi_figsize_width_var, width=18)
        self.hsi_figsize_width_entry.grid(row=8, column=1, padx=5, pady=2)
        
        tk.Label(frame, text="Figure height (inches):").grid(row=9, column=0, sticky=tk.W, padx=5)
        self.hsi_figsize_height_var = tk.DoubleVar(value=self.defentries.get('hsi_figsize_height', 6.0))
        self.hsi_figsize_height_entry = tk.Entry(frame, textvariable=self.hsi_figsize_height_var, width=18)
        self.hsi_figsize_height_entry.grid(row=9, column=1, padx=5, pady=2)
        
        # Title
        tk.Label(frame, text="Plot title:").grid(row=10, column=0, sticky=tk.W, padx=5)
        self.hsi_title_var = tk.StringVar(value=self.defentries.get('hsi_title', ''))
        self.hsi_title_entry = tk.Entry(frame, textvariable=self.hsi_title_var, width=18)
        self.hsi_title_entry.grid(row=10, column=1, padx=5, pady=2)
        
        # Colorbar options
        tk.Label(frame, text="Colorbar unit:").grid(row=11, column=0, sticky=tk.W, padx=5)
        self.hsi_cbar_unit_var = tk.StringVar(value=self.defentries.get('hsi_cbar_unit', 'Kilo counts'))
        self.hsi_cbar_unit_entry = tk.Entry(frame, textvariable=self.hsi_cbar_unit_var, width=18)
        self.hsi_cbar_unit_entry.grid(row=11, column=1, padx=5, pady=2)
        
        self.hsi_show_colorbar_var = tk.BooleanVar(value=self.defentries.get('hsi_show_colorbar', True))
        tk.Checkbutton(frame, text="Show colorbar", variable=self.hsi_show_colorbar_var).grid(row=12, column=0, columnspan=2, pady=2)
        
        # Font size
        tk.Label(frame, text="Scale bar font size:").grid(row=13, column=0, sticky=tk.W, padx=5)
        self.hsi_scalebar_fontsize_var = tk.IntVar(value=self.defentries.get('hsi_scalebar_fontsize', 12))
        self.hsi_scalebar_fontsize_entry = tk.Entry(frame, textvariable=self.hsi_scalebar_fontsize_var, width=18)
        self.hsi_scalebar_fontsize_entry.grid(row=13, column=1, padx=5, pady=2)
        
        # Unit
        tk.Label(frame, text="Unit (e.g., $\\mu m$):").grid(row=14, column=0, sticky=tk.W, padx=5)
        self.hsi_unit_var = tk.StringVar(value=self.defentries.get('hsi_unit', '$\\mu m$'))
        self.hsi_unit_entry = tk.Entry(frame, textvariable=self.hsi_unit_var, width=18)
        self.hsi_unit_entry.grid(row=14, column=1, padx=5, pady=2)
        
        # Plot button
        tk.Button(frame, text="Plot HSI with Options", command=self.plot_hsi_with_options, 
                 bg='lightblue').grid(row=15, column=0, columnspan=2, pady=10, padx=5, sticky=tk.EW)
        
        # Save to file button
        tk.Button(frame, text="Save Plot to File", command=self.save_hsi_plot_to_file,
                 bg='lightgreen').grid(row=16, column=0, columnspan=2, pady=5, padx=5, sticky=tk.EW)
    
    def plot_multiple_rois_on_hsi(self, plotmode='cornerlines'):
        """Plot multiple ROIs on the selected HSI with different colors."""
        try:
            # Get selected HSI
            hsi_name = self.hsiselect.get()
            if hsi_name not in self.PMdict.keys():
                print('Error: No valid HSI selected.')
                return
            
            # Get the PMclass object
            pm_obj = self.PMdict[hsi_name]
            pixmatrix = pm_obj.PixMatrix

            if self.roi_plot_indices_entry.get().strip() == '': # if entry is empty, plot all rois
                roibynames = list(self.roihandler.roilist.keys()) # if no indices entered, plot all rois
            else:
                # Get ROI indices to plot
                roi_indices_str = self.roi_plot_indices_entry.get().split(',')
                if not roi_indices_str:
                    print('Error: No ROI indices entered. Please enter indices separated by commas (e.g., 1,2,3).')
                    return
                
                roi_indices = []
                for i in range(len(roi_indices_str)):
                    try: 
                        roi_indices.append(int(roi_indices_str[i].strip())) # split by comma and take first part, then convert to int
                    except Exception as e:
                        print(f'Error parsing ROI index "{roi_indices[i]}": {e}. Skipping this index.')
                roibynames = roihandler.roiindicees2roinames(self.roihandler, roi_indices)
            
            # Get selected color for each ROI (cycle through colors if more ROIs than colors)
            color_options = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta', 'yellow', 'teal', 'navy', 'maroon', 'olive', 'lime', 'coral', 'indigo', 'gold', 'silver']
            
            #for i, roi_idx in enumerate(roibynames):
            #    color = color_options[i % len(color_options)]
            #    self.roihandler.plot_roi_on_pixmatrix(pixmatrix, roi_idx, vis_type=self.defentries.get('roi_vis_type', 'cornerlines'), color=color)
            #plot_multiple_rois_on_pixmatrix(self, handler, pixmatrix, roinames, plotmodes, colors, fontsize=14):
            # get selected plotmode from ROI Plot Type combobox
            selplotmode = self.roivistypeselect.get()
            if selplotmode not in ['cornerlines', 'overlay']:
                print(f'Warning: Invalid ROI plot mode "{selplotmode}". Defaulting to "cornerlines".')
                selplotmode = defaultplotmode = self.defentries.get('roi_vis_type', 'cornerlines')

            plotmodes = [selplotmode] * len(roibynames)
            colors = [color_options[i % len(color_options)] for i in range(len(roibynames))]
            unit = 'Intensity'
            if hasattr(self.PMdict[hsi_name].metadata, 'unit'):
                unit = self.PMdict[hsi_name].metadata['unit']

            self.roihandler.plot_multiple_rois_on_pixmatrix(
                self.roihandler,
                pixmatrix, 
                roibynames, 
                plotmodes,
                colors,
                fontsize=self.defentries.get('roi_plot_fontsize', 14),
                title=self._get_hsi_title(hsi_name), 
                unit = unit
            )
        
        except Exception as e:
            print(f'Error plotting multiple ROIs on HSI: {e}')
            traceback.print_exc()

    def plot_hsi_with_options(self):
        """Plot the selected HSI using the options from the GUI."""
        try:
            # Get selected HSI
            hsi_name = self.hsiselect.get()
            if hsi_name not in self.PMdict.keys():
                print('Error: No valid HSI selected.')
                return
            
            # Get the PMclass object
            pm_obj = self.PMdict[hsi_name]
            data = pm_obj.PixMatrix
            
            # Parse vmin/vmax (empty string means None)
            vmin_str = self.hsi_vmin_var.get().strip()
            vmax_str = self.hsi_vmax_var.get().strip()
            vmin = float(vmin_str) if vmin_str else None
            vmax = float(vmax_str) if vmax_str else None
            
            # Prepare metadata
            metadata = {
                'dx': pm_obj.gdx,
                'dy': pm_obj.gdy,
                'unit': self.hsi_unit_var.get()
            }
            
            # Prepare parameters
            params = {
                'cmap': self.hsi_cmap_var.get(),
                'vmin': vmin,
                'vmax': vmax,
                'scalebarlength': self.hsi_scalebar_len_var.get(),
                'scalebarwidth': self.hsi_scalebar_width_var.get(),
                'scalebarpos': (self.hsi_scalebar_pos_x_var.get(), self.hsi_scalebar_pos_y_var.get()),
                'figsize': (self.hsi_figsize_width_var.get(), self.hsi_figsize_height_var.get()),
                'title': self.hsi_title_var.get(),
                'show_colorbar': self.hsi_show_colorbar_var.get(),
                'cbar_unit': self.hsi_cbar_unit_var.get(),
                'scalebarfontsize': self.hsi_scalebar_fontsize_var.get(),
                'enable_drag': False,
                'dx': pm_obj.gdx,
                'unit': self.hsi_unit_var.get()
            }
            
            # Initialize HSIPlotManager if not exists
            if not hasattr(self, 'hsi_plot_manager'):
                self.hsi_plot_manager = deflib.HSIPlotManager(None, params)
            else:
                self.hsi_plot_manager.params.update(params)
            
            # Create the plot
            fig, ax = self.hsi_plot_manager.construct_plot(data, metadata=metadata)
            plt.show()
            
        except Exception as e:
            print(f'Error plotting HSI with options: {e}')
            import traceback
            traceback.print_exc()
    
    def save_hsi_plot_to_file(self):
        """Save the HSI plot to a file."""
        try:
            # Get selected HSI
            hsi_name = self.hsiselect.get()
            if hsi_name not in self.PMdict.keys():
                print('Error: No valid HSI selected.')
                return
            
            # Ask for filename
            filename = tkfd.asksaveasfilename(
                defaultextension='.png',
                filetypes=[('PNG files', '*.png'), ('PDF files', '*.pdf'), ('SVG files', '*.svg'), ('All files', '*.*')]
            )
            
            if not filename:
                return
            
            # Get the PMclass object
            pm_obj = self.PMdict[hsi_name]
            data = pm_obj.PixMatrix
            
            # Parse vmin/vmax
            vmin_str = self.hsi_vmin_var.get().strip()
            vmax_str = self.hsi_vmax_var.get().strip()
            vmin = float(vmin_str) if vmin_str else None
            vmax = float(vmax_str) if vmax_str else None
            
            # Prepare metadata
            metadata = {
                'dx': pm_obj.gdx,
                'dy': pm_obj.gdy,
                'unit': self.hsi_unit_var.get()
            }
            
            # Prepare parameters
            params = {
                'cmap': self.hsi_cmap_var.get(),
                'vmin': vmin,
                'vmax': vmax,
                'scalebarlength': self.hsi_scalebar_len_var.get(),
                'scalebarwidth': self.hsi_scalebar_width_var.get(),
                'scalebarpos': (self.hsi_scalebar_pos_x_var.get(), self.hsi_scalebar_pos_y_var.get()),
                'figsize': (self.hsi_figsize_width_var.get(), self.hsi_figsize_height_var.get()),
                'title': self.hsi_title_var.get(),
                'show_colorbar': self.hsi_show_colorbar_var.get(),
                'cbar_unit': self.hsi_cbar_unit_var.get(),
                'scalebarfontsize': self.hsi_scalebar_fontsize_var.get(),
                'enable_drag': False,
                'dx': pm_obj.gdx,
                'unit': self.hsi_unit_var.get()
            }
            
            # Create plot
            fig, ax = deflib.plot_HSI(data, metadata=metadata, **params)
            
            # Save to file
            fig.savefig(filename, dpi=600, bbox_inches='tight')
            plt.close(fig)
            
            print(f'HSI plot saved to: {filename}')
            
        except Exception as e:
            print(f'Error saving HSI plot: {e}')
            traceback.print_exc()
    
    def select_correction_spectrum_file(self):
        self.correctionspecname = tkfd.askopenfilename(filetypes=[("Correction spectrum", "*")])
        self.loadcorrectionSpectrum()
    
    def saveHSI(self):
        # pickle the selected HSI to a file
        filename = tkfd.asksaveasfilename(defaultextension='.pkl', filetypes=[('Pickle files', '*.pkl')])
        if filename:
            # pickle the selected HSI to filename
            selhsi = self.hsiselect.get()
            if selhsi in self.PMdict.keys():
                with open(filename, 'wb') as f:
                    pickle.dump(self.PMdict[selhsi], f)
    
    def loadHSI(self):
        # load a pickle file and add it to the HSI list
        filename = tkfd.askopenfilename(filetypes=[('Pickle files', '*.pkl')])
        if filename:
            with open(filename, 'rb') as f:
                hsi = pickle.load(f)
            # Use counter for unique naming
            hsiname = f'HSI{self._hsi_counter}'
            self._hsi_counter += 1
            # add the loaded HSI to the HSI list
            self.PMdict[hsiname] = hsi
            self.hsiselect['values'] = list(self.PMdict.keys())
            self.hsiselect.set(hsiname)
    
    def exportHSI(self): # export the selected HSI to a .csv file
        filename = tkfd.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV files', '*.csv')])
        if filename:
            selhsi = self.hsiselect.get()
            # create a string from the metadata 'key': 'value' pairs
            metadata = '\n'.join([f'{key}: {value}' for key, value in self.PMdict[selhsi].metadata.items()])
            # write the metadata to the file
            with open(filename, 'w') as f:
                f.write(metadata + '\n')
                # write the HSI data to the file
                for row in self.PMdict[selhsi].PixMatrix:
                    # Convert np.nan to 'Nan', otherwise str(value)
                    formatted_row = [
                        'Nan' if (isinstance(val, float) and np.isnan(val)) else str(val)
                        for val in row
                    ]
                    f.write(';'.join(formatted_row) + '\n')
    
    def exportHSI_itx(self): # export the selected HSI to a .itx file
        '''
        Export the selected HSI to a .itx file.
        .itx file format: 
        IGOR
        WAVES/N=(x, y) <name>
        BEGIN
        <data>
        END
        X SetScale/P x 0,1,"" dataname
        X SetScale/P y 0,1,"" dataname
        X Note dataname, "metadata"
        metadata as key-value pairs: "key1: value1 ; key2: value2 ; ..."
        data are simply tab-seperated table, with NaN values replaced by NaN, 
        '''
        # get metadata ready:
        hsiname = self.hsiselect.get() # replace spaces with underscores for the name
        hsiname_igor = hsiname.replace(' ', '_')
        # this is essential, since Igor pro recognizes na me as wave 2 waves na and me
        metadata = ' ; '.join([f'{key}: {value} ' for key, value in self.PMdict[hsiname].metadata.items()])
        #data = np.transpose(self.PMdict[hsiname].PixMatrix)
        data = self.PMdict[hsiname].PixMatrix # testing, we must "rotate" the data by -90° clockwise, since Igor pro expects the first dimension to be x and the second to be y
        # data = np.flipud(np.transpose(data)) # rotate +90° clockwise <- this one was the wrong direcation lol
        data = np.rot90(data, k=-1) # rotate -90° clockwise
        xscaling = self.PMdict[hsiname].gdx
        yscaling = self.PMdict[hsiname].gdy
        # write to .itx file
        filename = tkfd.asksaveasfilename(defaultextension='.itx', filetypes=[('IGOR Pro files', '*.itx')])
        if filename:
            with open(filename, 'w') as f:
                f.write('IGOR\n')
                f.write(f'WAVES/N=({data.shape[0]}, {data.shape[1]}) {hsiname_igor}\n') # do not change the [0] and [1] - I figured this out!
                f.write('BEGIN\n')
                for row in data:
                    formatted_row = [
                        'NaN' if (isinstance(val, float) and np.isnan(val)) else str(val)
                        for val in row
                    ]
                    f.write('\t'.join(formatted_row) + '\n')
                f.write('END\n\n')
                # write the x ticks into a wave 
                f.write(f'WAVES/D xticks_{hsiname_igor}\nBEGIN\n')
                for i in range(data.shape[0]+1):
                    f.write(f'{i * xscaling}\n')
                f.write('END\n\n')
                # write the y ticks into a wave
                f.write(f'WAVES/D yticks_{hsiname_igor}\nBEGIN\n')
                for i in range(data.shape[1]+1):
                    f.write(f'{i * yscaling}\n')
                f.write('END\n\n')
                f.write(f'X SetScale/P x 0,{xscaling},"" {hsiname_igor}\n')
                f.write(f'X SetScale/P y 0,{yscaling},"" {hsiname_igor}\n')
                f.write(f'X Note {hsiname_igor}, "{metadata}"\n')
    
    def send_to_HSIPlot(self, hsiplotinstance):
        selhsi = self.hsiselect.get()
        if selhsi in self.PMdict.keys():
            hsiplotinstance.set_hsi_data(self.PMdict[selhsi].PixMatrix, self.PMdict[selhsi].metadata)
        # tk: switch to HSI Plot tab
        self.tab_control.select(self.tab_hsiplot) # type: ignore
    
    def exportHSIWithSpectra(self):
        """
        Export the selected HSI along with associated averaged spectra to CSV files.
        Creates multiple files:
        - <basename>_matrix.csv: HSI pixel matrix
        - <basename>_avg_spectra.csv: All averaged spectra from this HSI
        """
        filename = tkfd.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV files', '*.csv')])
        if not filename:
            return
        
        # Get base filename without extension
        base_name = filename.rsplit('.', 1)[0]
        selhsi = self.hsiselect.get()
        
        # Export HSI matrix
        matrix_filename = f"{base_name}_matrix.csv"
        metadata = '\n'.join([f'{key}: {value}' for key, value in self.PMdict[selhsi].metadata.items()])
        with open(matrix_filename, 'w') as f:
            f.write(metadata + '\n')
            for row in self.PMdict[selhsi].PixMatrix:
                formatted_row = [
                    'Nan' if (isinstance(val, float) and np.isnan(val)) else str(val)
                    for val in row
                ]
                f.write(';'.join(formatted_row) + '\n')
        print(f"Exported HSI matrix to: {matrix_filename}")
        
        # Export averaged spectra from this HSI
        spectra_filename = f"{base_name}_avg_spectra.csv"
        
        # Find all averaged spectra from this HSI
        hsi_spectra = {name: spec for name, spec in self.disspecs.items() 
                       if spec.parenthsi == selhsi}
        
        if len(hsi_spectra) > 0:
            with open(spectra_filename, 'w') as f:
                # Write header with metadata
                f.write(f"# Averaged spectra from HSI: {selhsi}\n")
                f.write(f"# Wavelength range: {self.wlstart} - {self.wlend} nm\n")
                f.write(f"# Number of spectra: {len(hsi_spectra)}\n")
                f.write("#\n")
                
                # Write column headers
                headers = ['Wavelength (nm)'] + list(hsi_spectra.keys())
                f.write(','.join(headers) + '\n')
                
                # Get wavelength array from first spectrum
                first_spec = list(hsi_spectra.values())[0]
                wavelengths = first_spec.WL
                
                # Write data rows
                for i, wl in enumerate(wavelengths):
                    row = [str(wl)]
                    for spec in hsi_spectra.values():
                        if i < len(spec.Spec):
                            row.append(str(spec.Spec[i]))
                        else:
                            row.append('NaN')
                    f.write(','.join(row) + '\n')
            
            print(f"Exported {len(hsi_spectra)} averaged spectra to: {spectra_filename}")
        else:
            print(f"No averaged spectra found for HSI: {selhsi}")
        
        print("Export complete!")
    
    def exportAllAveragedSpectra(self):
        """
        Export all averaged spectra (all data types) for the selected HSI.
        First generates averaged spectra for all data types if they don't exist,
        then exports them to a CSV file.
        """
        filename = tkfd.asksaveasfilename(defaultextension='.csv', 
                                          filetypes=[('CSV files', '*.csv')])
        if not filename:
            return
        
        selhsi = self.hsiselect.get()
        
        # Generate averaged spectra for all data types
        print(f"Generating averaged spectra for all data types from {selhsi}...")
        data_types = ['PLB', 'Specdiff1', 'Specdiff2', 'Specdiff1_norm', 'Specdiff2_norm', 'Specdiff1_norm_intensity', 'Specdiff2_norm_intensity', 'Specdiff1_norm_counts', 'Specdiff2_norm_counts']
        generated_spectra = self.averageHSItoSpecDataMultiple(data_types)
        
        if len(generated_spectra) == 0:
            print("No spectra were generated. Export cancelled.")
            return
        
        # Export all spectra to CSV
        with open(filename, 'w') as f:
            # Write header
            f.write(f"# All averaged spectra from HSI: {selhsi}\n")
            f.write(f"# Wavelength range: {self.wlstart} - {self.wlend} nm\n")
            f.write(f"# Generated: {len(generated_spectra)} spectral data types\n")
            f.write("#\n")
            
            # Get data type mapping for column headers
            type_names = []
            for spec_name in generated_spectra.keys():
                # Extract data type from name (format: HSI0_PLB_avg -> PLB)
                parts = spec_name.split('_')
                if len(parts) >= 2:
                    data_key = parts[-2]  # Get the data type key
                    # Find human-readable name
                    for name, key in self.speckeys.items():
                        if key == data_key:
                            type_names.append(name)
                            break
                    else:
                        type_names.append(data_key)
                else:
                    type_names.append(spec_name)
            
            # Write column headers
            headers = ['Wavelength (nm)'] + type_names
            f.write(','.join(headers) + '\n')
            
            # Get wavelength array from first spectrum
            first_spec = list(generated_spectra.values())[0]
            wavelengths = first_spec.WL
            
            # Write data rows
            for i, wl in enumerate(wavelengths):
                row = [str(wl)]
                for spec in generated_spectra.values():
                    if i < len(spec.Spec):
                        row.append(str(spec.Spec[i]))
                    else:
                        row.append('NaN')
                f.write(','.join(row) + '\n')
        
        print(f"Exported {len(generated_spectra)} averaged spectra to: {filename}")
        print("Export complete!")

    def correctSpectrum(self, specname):
        # correct the selected spectrum with the entered correction spectrum
        if self.correctionspecname == '':
            self.correctionspecname = tkfd.askopenfilename(filetypes=[('Correction spectrum', '*')])
        self.correctionWL, self.correctionSpec = deflib.loadexpspec(self.correctionspecname)

        # clear plt.figure
        plt.clf()
        # plot the correction spectrum
        plt.plot(self.correctionWL, self.correctionSpec)
        plt.plot(self.disspecs[specname].WL, self.disspecs[specname].Spec)
        plt.show()

        # correct the selected spectrum with the correction spectrum. Make sure to match the WL by interpolation
        corrected_spec = deflib.correct_spectrum(self.disspecs[specname].Spec, self.disspecs[specname].WL, self.correctionSpec, self.correctionWL)
        spname = '{}_corrected'.format(self.createdisspecname())
        # Create new Spectra object more efficiently - avoid deepcopy
        source_spec = self.disspecs[specname]
        self.disspecs[spname] = PMlib.Spectra(
            corrected_spec,
            source_spec.WL,
            source_spec.metadata.copy() if isinstance(source_spec.metadata, dict) else source_spec.metadata,
            source_spec.parenthsi
        )
        # Copy derivative data if present
        if hasattr(source_spec, 'Spec_d1') and source_spec.Spec_d1 is not None:
            self.disspecs[spname].Spec_d1 = source_spec.Spec_d1.copy() if hasattr(source_spec.Spec_d1, 'copy') else source_spec.Spec_d1
        if hasattr(source_spec, 'Spec_d2') and source_spec.Spec_d2 is not None:
            self.disspecs[spname].Spec_d2 = source_spec.Spec_d2.copy() if hasattr(source_spec.Spec_d2, 'copy') else source_spec.Spec_d2
            
        self._enforce_disspecs_memory_limit()
        
        # update the combobox with the new spectrum name
        # todo: fix this
        self.specselect['values'] = list(self.disspecs.keys())
        self.specselect.set(list(self.disspecs.keys())[-1])
        #self.disspecs['{}_corrected'.format(self.createdisspecname())].Spec = corrected_spec
    
    def loadcorrectionSpectrum(self):
        # load a correction spectrum
        with open(self.correctionspecname, 'r') as file:
            _ = file.readline()
            lines = file.readlines()
        WL = []
        spec = []
        for line in lines:
            line = line.strip()
            if line:
                parts = line.split('\n')[0].split('\t')
                WL.append(float(parts[0]))
                spec.append(float(parts[1]))
        self.correctionWL = WL
        self.correctionSpec = spec
        del lines
    
    def _enforce_disspecs_memory_limit(self):
        """Limit the number of stored spectra to prevent memory leaks."""
        MAX_SPECTRA_CACHE = 50
        while len(self.disspecs) > MAX_SPECTRA_CACHE:
            oldest_key = next(iter(self.disspecs))
            del self.disspecs[oldest_key]
            
    def averageHSItoSpecData(self):
        # average all HSI to a single SpectrumData
        self.hsiselected = self.hsiselect.get()
        # get wlstart and wlend from self.hsiselected
        meta = getattr(self.PMdict[self.hsiselected], 'metadata', {})
        if not isinstance(meta, dict):
            meta = {}
            
        if 'wlstart' in meta and 'wlend' in meta and 'aqpixstart' in meta and 'aqpixend' in meta:
            self.wlstart = meta['wlstart']
            self.wlend = meta['wlend']
            self.aqpixstart = meta['aqpixstart']
            self.aqpixend = meta['aqpixend']
        else:
            self.updatewl()
        
        # Get selected data source
        data_source_name = self.selectspecbox.get()
        # Fallback to PLB if something is wrong or key not found
        data_key = self.speckeys.get(data_source_name, 'PLB')
        
        # Check if we are trying to average something that isn't an array (like WL)
        if data_key == 'WL':
             print("Cannot average Wavelength axis. Using Spectrum (PL-BG) instead.")
             data_key = 'PLB'
             data_source_name = 'Spectrum (PL-BG)'
        
        # update wlstart and wlend entries on GUI
        metadata = {'wlstart': self.wlstart, 'wlend': self.wlend, 
                    'countthresh': self.countthreshv, 'aqpixstart': self.aqpixstart, 
                    'aqpixend': self.aqpixend, 'averaged_data_type': data_source_name}
        
        WL = self.WL[self.aqpixstart: self.aqpixend]
        # Initialize with zeros instead of copying WL - more efficient and correct
        accumulated_data = np.zeros_like(WL, dtype=float)
        speccount = 0

        # Optimized averaging: vectorized inner loop
        aqstart = int(self.aqpixstart)
        aqend = int(self.aqpixend)
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                if np.isnan(self.PMdict[self.hsiselected].PixMatrix[i][j]) == False:
                    # Check if attribute exists and is valid
                    pixel_obj = self.SpecDataMatrix[i][j]
                    if hasattr(pixel_obj, data_key):
                         data_array = getattr(pixel_obj, data_key)
                         
                         # Ensure data array is compatible and exists
                         if data_array is not None and len(data_array) >= aqend:
                            speccount += 1
                            # Vectorized addition
                            accumulated_data += data_array[aqstart:aqend]
        
        if speccount > 0:
            accumulated_data = np.divide(accumulated_data, speccount)
        
        new_spec = PMlib.Spectra(accumulated_data, WL, metadata, self.hsiselected)
        
        # Calculate derivatives using optimized function
        if self.derivative_polynomarray and len(self.derivative_polynomarray) >= 4:
            try:
                # Handle both Tkinter variables and direct values
                d1_val = self.derivative_polynomarray[0]
                calc_d1 = d1_val.get() if hasattr(d1_val, 'get') else d1_val
                
                d2_val = self.derivative_polynomarray[1]
                calc_d2 = d2_val.get() if hasattr(d2_val, 'get') else d2_val
                
                if calc_d1 or calc_d2:
                    p_order_val = self.derivative_polynomarray[2]
                    poly_order = int(p_order_val.get()) if hasattr(p_order_val, 'get') else int(p_order_val)
                    
                    n_points_val = self.derivative_polynomarray[3]
                    N_fitpoints = int(n_points_val.get()) if hasattr(n_points_val, 'get') else int(n_points_val)
                    
                    # Use optimized calc_derivative function instead of manual loop
                    derivative_config = [calc_d1, calc_d2, poly_order, N_fitpoints]
                    PMlib.calc_derivative(new_spec, derivative_config)
            except Exception as e:
                print(f"Error calculating derivatives for averaged spectrum: {e}")

        self.disspecs[self.createdisspecname()] = new_spec
        self._enforce_disspecs_memory_limit()

        # update the selectbox for spectral data
        self.specselect['values'] = list(self.disspecs.keys())
        self.specselect.set(list(self.disspecs.keys())[-1])
    
    def averageHSItoSpecDataMultiple(self, data_types=None):
        """
        Average HSI to multiple spectral data types simultaneously.
        
        Parameters:
        -----------
        data_types : list, optional
            List of data type keys to average. If None, defaults to all available types:
            ['PLB', 'Specdiff1', 'Specdiff2', 'Specdiff1_norm', 'Specdiff2_norm', 'Specdiff1_norm_intensity', 'Specdiff2_norm_intensity', 'Specdiff1_norm_counts', 'Specdiff2_norm_counts']
        
        Returns:
        --------
        dict : Dictionary of generated averaged spectra with their names as keys
        """
        # Default to all available data types if not specified
        if data_types is None:
            data_types = ['PLB', 'Specdiff1', 'Specdiff2', 'Specdiff1_norm', 'Specdiff2_norm', 'Specdiff1_norm_intensity', 'Specdiff2_norm_intensity', 'Specdiff1_norm_counts', 'Specdiff2_norm_counts']
        
        # Get current HSI selection
        self.hsiselected = self.hsiselect.get()
        meta = getattr(self.PMdict[self.hsiselected], 'metadata', {})
        if not isinstance(meta, dict):
            meta = {}
            
        if 'wlstart' in meta and 'wlend' in meta and 'aqpixstart' in meta and 'aqpixend' in meta:
            self.wlstart = meta['wlstart']
            self.wlend = meta['wlend']
            self.aqpixstart = meta['aqpixstart']
            self.aqpixend = meta['aqpixend']
        else:
            self.updatewl()
        
        # Prepare metadata
        base_metadata = {
            'wlstart': self.wlstart, 
            'wlend': self.wlend, 
            'countthresh': self.countthreshv, 
            'aqpixstart': self.aqpixstart, 
            'aqpixend': self.aqpixend
        }
        
        WL = self.WL[self.aqpixstart: self.aqpixend]
        aqstart = int(self.aqpixstart)
        aqend = int(self.aqpixend)
        
        # Check if derivatives need to be calculated for individual pixels
        needs_derivatives = any(dt in data_types for dt in ['Specdiff1', 'Specdiff2', 'Specdiff1_norm', 'Specdiff2_norm', 'Specdiff1_norm_intensity', 'Specdiff2_norm_intensity', 'Specdiff1_norm_counts', 'Specdiff2_norm_counts'])
        
        if needs_derivatives:
            # Ensure derivatives are calculated for all spectra in SpecDataMatrix
            print("Checking if derivatives are calculated for individual spectra...")
            derivatives_exist = False
            
            # Check if derivatives exist in SpecDataMatrix
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    pixel_obj = self.SpecDataMatrix[i][j]
                    if pixel_obj is not None:
                        if hasattr(pixel_obj, 'Specdiff1') and pixel_obj.Specdiff1 is not None:
                            derivatives_exist = True
                            break
                if derivatives_exist:
                    break
            
            if not derivatives_exist:
                print("Derivatives not found in SpecDataMatrix. Calculating derivatives...")
                self.calculate_derivatives()
            else:
                print("Derivatives already exist in SpecDataMatrix.")
        
        generated_spectra = {}
        
        for data_key in data_types:
            # Get the human-readable name for this data type
            data_source_name = None
            for name, key in self.speckeys.items():
                if key == data_key:
                    data_source_name = name
                    break
            
            if data_source_name is None:
                print(f"Warning: Unknown data type '{data_key}', skipping...")
                continue
            
            print(f"Averaging {data_source_name}...")
            
            # Initialize accumulator
            accumulated_data = np.zeros_like(WL, dtype=float)
            speccount = 0
            
            # Average over all valid pixels in SpecDataMatrix
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    if np.isnan(self.PMdict[self.hsiselected].PixMatrix[i][j]) == False:
                        pixel_obj = self.SpecDataMatrix[i][j]
                        if hasattr(pixel_obj, data_key):
                            data_array = getattr(pixel_obj, data_key)
                            
                            if data_array is not None and len(data_array) >= aqend:
                                speccount += 1
                                accumulated_data += data_array[aqstart:aqend]
            
            if speccount > 0:
                accumulated_data = np.divide(accumulated_data, speccount)
            else:
                print(f"Warning: No valid spectra found for {data_source_name}")
                continue
            
            # Create metadata for this specific data type
            metadata = base_metadata.copy()
            metadata['averaged_data_type'] = data_source_name
            
            # Create the spectrum object
            new_spec = PMlib.Spectra(accumulated_data, WL, metadata, self.hsiselected)
            
            # Calculate derivatives for the averaged spectrum itself
            if self.derivative_polynomarray and len(self.derivative_polynomarray) >= 4:
                try:
                    d1_val = self.derivative_polynomarray[0]
                    calc_d1 = d1_val.get() if hasattr(d1_val, 'get') else d1_val
                    
                    d2_val = self.derivative_polynomarray[1]
                    calc_d2 = d2_val.get() if hasattr(d2_val, 'get') else d2_val
                    
                    if calc_d1 or calc_d2:
                        p_order_val = self.derivative_polynomarray[2]
                        poly_order = int(p_order_val.get()) if hasattr(p_order_val, 'get') else int(p_order_val)
                        
                        n_points_val = self.derivative_polynomarray[3]
                        N_fitpoints = int(n_points_val.get()) if hasattr(n_points_val, 'get') else int(n_points_val)
                        
                        derivative_config = [calc_d1, calc_d2, poly_order, N_fitpoints]
                        PMlib.calc_derivative(new_spec, derivative_config)
                except Exception as e:
                    print(f"Error calculating derivatives for averaged spectrum: {e}")
            
            # Create a descriptive name for this averaged spectrum
            import re
            match = re.search(r'(.*)(\(roi\d+\))(.*)$', self.hsiselected)
            if match:
                spec_name = f"{match.group(1)}_{data_key}_avg{match.group(2)}{match.group(3)}"
            else:
                spec_name = f"{self.hsiselected}_{data_key}_avg"
            
            # Store in disspecs
            self.disspecs[spec_name] = new_spec
            self._enforce_disspecs_memory_limit()
            generated_spectra[spec_name] = new_spec
            
            print(f"  Created averaged spectrum: {spec_name}")
        
        # Update the GUI combobox
        if hasattr(self, 'specselect'):
            self.specselect['values'] = list(self.disspecs.keys())
            if len(self.disspecs) > 0:
                self.specselect.set(list(self.disspecs.keys())[-1])
        
        print(f"Generated {len(generated_spectra)} averaged spectra")
        return generated_spectra
    
    def createdisspecname(self): # create a new spectral data name
        import re
        HSIname = self.hsiselect.get()
        if HSIname == '':
            base_name = 'SpectrumData'
            roi_suffix = ''
        else: 
            match = re.search(r'(.*)(\(roi\d+\))(.*)$', HSIname)
            if match:
                base_name = match.group(1)
                roi_suffix = match.group(2) + match.group(3)
            else:
                base_name = HSIname
                roi_suffix = ''

        i = 0
        specname = f"{base_name}{i}{roi_suffix}"
        while specname in self.disspecs.keys():
            i += 1
            specname = f"{base_name}{i}{roi_suffix}"
            
        return specname 
    
    def saveSpectrum(self, specname):
        savename = tkfd.asksaveasfilename(defaultextension='.txt', filetypes=[('Text files', '*.txt')])
        if savename:
            self.disspecs[specname].save(savename)
    
    def delSpecData(self, specname): # self.disspecs is type dict
        # delete the selected spectral data
        try:
            del self.disspecs[specname]
            if len(self.disspecs) == 0:
                pass
            else:
                self.specselect['values'] = list(self.disspecs.keys())
                self.specselect.set(list(self.disspecs.keys())[-1])
        except Exception as e:
            print('Error deleting spectral data.', e)
    
    def plotSpectral(self, specname):
        # plot the selected spectral data
        try:
            selspec = self.specselect.get()
            selspecdataclass = self.disspecs[selspec] # Use selspec (the name) to get the object
            
            # Check if we want to plot derivatives based on selection in selectspecpixbox (or similar)
            # But wait, plotSpectral is triggered by button b7 in build_roi_frame.
            # It plots the spectrum selected in self.specselect.
            
            # We need to know WHAT to plot (PLB, d1, d2). 
            # Currently, it seems to default to 'Averaged HSI Spectrum' (PLB).
            # Let's check if we can use the selection from selectspecpixbox or similar if available,
            # OR just plot what is available.
            
            # For now, let's plot the main spectrum. 
            # If the user wants derivatives, they might need to select "first derivative" in a dropdown.
            # However, the user request says: "Plot Spectrum with Selected Data Set second derivative must display the 2nd derivative spectrum"
            
            # This implies we should look at self.selectspecpixbox (from build_button_frame) or similar?
            # No, build_button_frame is for Pixel Spectrum.
            # build_roi_frame has "Plot Spectrum" button.
            
            # Let's look at where "first derivative" / "second derivative" are selected.
            # They are in self.speckeys, used in self.selectspecpixbox (Pixel Spectrum) and self.selectspecbox (Fit Window).
            
            # If the user wants to see the derivative of the AVERAGED spectrum (disspecs),
            # we need to check if disspecs has those attributes.
            
            # Let's try to plot what is requested if possible, or default to Spec.
            
            # Actually, the user might be referring to PlotPixelSpectrum (for a single pixel) OR plotSpectral (for averaged).
            # The request says "Hyperspectra Notebook, where the user cliks on plot spectrum".
            # This usually refers to the averaged spectrum from ROI.
            
            # Let's modify plotSpectral to check if derivatives exist and plot them if they do?
            # Or better, let's make sure PlotPixelSpectrum works for derivatives (it uses Specdiff1/2).
            
            # But wait, the user said "Plot Spectrum with Selected Data Set second derivative".
            # This likely refers to the "Select Data Set" combobox in the "Pixel Spectrum" area (build_button_frame).
            # That combobox is self.selectspecpixbox.
            # The button is "Plot Spectrum" (b1 in build_button_frame), which calls self.PlotPixelSpectrum.
            
            # So let's look at PlotPixelSpectrum.
            
            # Get averaged data type from metadata if available
            caption = 'Averaged HSI Spectrum'
            if hasattr(selspecdataclass, 'metadata') and 'averaged_data_type' in selspecdataclass.metadata:
                caption = f"Averaged HSI ({selspecdataclass.metadata['averaged_data_type']})"
            
            self.PlotSpectrum(selspecdataclass.Spec, selspecdataclass.WL, caption)
            
            # If the averaged spectrum has derivatives, maybe we should plot them too?
            if hasattr(selspecdataclass, 'Spec_d1') and selspecdataclass.Spec_d1 is not None:
                 self.PlotSpectrum(selspecdataclass.Spec_d1, selspecdataclass.WL, '1st Derivative')
            if hasattr(selspecdataclass, 'Spec_d2') and selspecdataclass.Spec_d2 is not None:
                 self.PlotSpectrum(selspecdataclass.Spec_d2, selspecdataclass.WL, '2nd Derivative')

        except Exception as e:
            print('Error plotting spectral data.', e)
    
    def build_plotfitparamframe(self, parframe):
        frame = tk.Frame(parframe, border=5, relief="sunken")
        frame.grid(row=0, column=6, rowspan=6, sticky=tk.NW)
        tk.Label(frame, text="Plot Fit Parameters".format(self.DataSpecMin)).grid(row=0, column=0)
        # add selectbox for fit parameters
        self.selectfitparambox = ttk.Combobox(frame, values=['FWHM', 'Center X', 'Center Y'])
    
    def fit2dgausstopixmatrix(self):
        try: 
            self.maxiter = int(self.fitmaxiter.get())
        except:
            print('Maxiter must be type int. Using default 1000.')
            self.maxiter = 1000
        try:
            PMselect = self.getPixMatrixSelection(self.hsiselect.get())
            matl.fitgaussian2dtomatrix(self.PMdict[PMselect].PixMatrix, True, self.gdx, self.gdy, self.colormap.get(), maxfev=self.maxiter)
        except:
            print('2D Gaussian Fit failed.')
            
    def fit2drotgausstopixmatrix(self):
        try: 
            self.maxiter = int(self.fitmaxiter.get())
        except:
            print('Maxiter must be type int. Using default 1000.')
            self.maxiter = 1000
        try:
            PMselect = self.getPixMatrixSelection(self.hsiselect.get())
            matl.fitgaussiand2dtomatrixrot(self.PMdict[PMselect].PixMatrix, True, self.gdx, self.gdy, self.colormap.get(), maxfev=self.maxiter)
        except Exception as Error:
            print('2D rotational Gaussian Fit failed.', Error)
        
    # Matrix with Pixels to obtain spectrum
    def build_button_frame(self, placeframe, width=600, height=600):
        # create new subframe

        parframe = tk.Frame(placeframe)
        #parframe.pack(side=tk.TOP, anchor=tk.W, expand=False)#, fill=tk.Y)
        parframe.grid(row=0, column=0, sticky=tk.NW)

        # create input GUI, ButtonMatrix is created in buildButtonMatrix in seperate frame
        # Build full GUI even if no data is loaded yet
        
        # Determine if we have data loaded
        has_data = len(self.PMdict) > 0
        
        # Get dimensions for button matrix if data exists
        if has_data:
            firstPM = list(self.PMdict.keys())[0]
            n = len(firstPM)
            m = len(firstPM[0])
        else:
            n = 0
            m = 0
        
        plotframe = tk.Frame(parframe, border=5, relief="raised")
        plotframe.grid(row=0, column=0, sticky=tk.NW)
        
        # Show appropriate status message
        if has_data:
            self.loadinfolabel = tk.Label(plotframe, text='Pixel Loaded: {} x {}'.format(len(self.SpecDataMatrix[0]), len(self.SpecDataMatrix))).pack(side=tk.TOP, anchor=tk.W)
        else:
            self.loadinfolabel = tk.Label(plotframe, text='No data loaded - Load HSI data', ).pack(side=tk.TOP, anchor=tk.W)
        
        tk.Label(plotframe, text='selected Pixel: ').pack(side=tk.TOP, anchor=tk.W)
        xyframe = tk.Frame(plotframe)
        xyframe.pack(side=tk.TOP, anchor=tk.W)
        # inside xyframe
        tk.Label(xyframe, text='X =').grid(row=0, column=0)
        self.selectPixX = tk.Entry(xyframe)
        self.selectPixX.grid(row=0, column=1)
        self.selectPixX.insert(0, self.defentries['selected_pixel_x'])
        tk.Label(xyframe, text='Y =').grid(row=1, column=0)
        self.selectPixY = tk.Entry(xyframe)#, text="0")
        self.selectPixY.grid(row=1, column=1)
        self.selectPixY.insert(0, self.defentries['selected_pixel_y'])
        # inside xyframe
        self.newselx = 0
        self.newsely = 0

        tk.Label(plotframe, text="Select Data Set".format(self.DataSpecMax)).pack(side=tk.TOP, anchor=tk.W)
        self.selectspecpixbox = ttk.Combobox(plotframe, values=list(self.speckeys.keys()), width=40)
        if self.default_dataset in self.speckeys:
            self.selectspecpixbox.set(self.default_dataset)
        else:
            self.selectspecpixbox.set(list(self.speckeys.keys())[-1])
        self.selectspecpixbox.pack(side=tk.TOP, anchor=tk.W)

        # add a combobox to select wavelength / nm or energy / eV axis
        tk.Label(plotframe, text="Select WL Axis").pack(side=tk.TOP, anchor=tk.W)
        self.selectspecxbox = ttk.Combobox(plotframe, values=['Wavelength (nm)', 'Energy (eV)'])
        self.selectspecxbox.set(self.defentries['selected_WL_axis'])
        self.selectspecxbox.pack(side=tk.TOP, anchor=tk.W)

        # bind event to selectspecxbox when selection changes
        self.selectspecxbox.bind('<<ComboboxSelected>>', lambda event: self.WL_selection.set(self.selectspecxbox.get()))

        b1 = tk.Button(plotframe, text="Plot Spectrum", command=self.PlotPixelSpectrum)
        b1.pack(side=tk.TOP, anchor=tk.W)
        # Disable button if no data
        if not has_data:
            b1.config(state='disabled')

        if self.defentries['enable_buttonmatrix'] == True and has_data:
            b2 = tk.Button(parframe, text="Create Button Matrix", command= lambda: self.buildButtonMatrix(parframe, n, m)) # do not use the buttonmatrix !!!
            b2.pack(side=tk.BOTTOM, anchor=tk.W)

        # add new frame for Fit Parameters
        fitframe = tk.Frame(placeframe, border=5, relief="raised")
        fitframe.grid(row=0, column=1, sticky=tk.NW)

        # fit to spectrum
        tk.Label(fitframe, text="Select Fit function".format(self.DataSpecMax)).pack(side=tk.TOP, anchor=tk.W)
        self.selectwindowbox = ttk.Combobox(fitframe, values=list(self.windowfunctions))
        self.selectwindowbox.set(self.defentries['selected_fit_function'])
        self.selectwindowbox.pack(side=tk.TOP, anchor=tk.W)
        b3 = tk.Button(fitframe, text="Fit Window to Spectrum", command=lambda: self.fitwindowtospec('fitmaxX', newfit=True))
        b3.pack(side=tk.TOP, anchor=tk.W)
        if not has_data:
            b3.config(state='disabled')
            
        self.sepfitfunct = tk.BooleanVar()
        b4 = tk.Button(fitframe, text="plot existing fit and spectrum", command=lambda: self.runPlotFitSpectrum())
        b4.pack(side=tk.TOP, anchor=tk.W)
        if not has_data:
            b4.config(state='disabled')
            
        self.sepfitfunct.set(False)
        self.sepfitfunctsbut = tk.Checkbutton(fitframe, text="Seperate Fit Functions", variable=self.sepfitfunct)
        self.sepfitfunctsbut.pack(side=tk.TOP, anchor=tk.W)

        # create a select box for the parameters of each fit
        tk.Label(fitframe, text="HSI from Fit Parameter").pack(side=tk.TOP, anchor=tk.W)
        self.selectfitparambox = ttk.Combobox(fitframe, values=self.allfpnamesinone, width=27)
        self.selectfitparambox.pack(side=tk.TOP, anchor=tk.W)
        self.selectfitparambox.set(self.allfpnamesinone[0]) # set default value
        b5 = tk.Button(fitframe, text="Plot HSI from Fit Parameter", command= lambda: self.plotHSIfromfitparam())
        b5.pack(side=tk.TOP, anchor=tk.W)
        if not has_data:
            b5.config(state='disabled')
            
        # add a bool var to check and uncheck HSI_from_fitparam_useROI
        b6 = tk.Checkbutton(fitframe, text="Use ROI for parameter plot", variable=self.HSI_from_fitparam_useROI).pack(side=tk.TOP, anchor=tk.W)

        self.build_roi_frame(placeframe)
        self.build_plot_options_frame(placeframe)
        buttons = []
        self.SpecButtons = []
        return buttons
        
        # fram = Tkinter frame, n = len(self.PMdict[i].PixMatrix), m = len(self.PMdict[i].PixMatrix[0])
    def buildButtonMatrix(self, frame, n, m):

        # create new subframe
        butmatframe = tk.Frame(frame)
        butmatframe.pack(side=tk.TOP, anchor=tk.W, fill=tk.BOTH, expand=True)

        # create buttons
        buttons = []
        self.vvars = []
        for row in range(n):
            row_buttons = []
            row_vars = []
            for col in range(m):
                var = [col, row]
                button = tk.Button(butmatframe, bg="red", command=lambda v=var: self.button_click(v))
                button.grid(row=row+1, column=col+1, sticky=tk.NSEW)
                row_buttons.append(button)
                row_vars.append(var)
            buttons.append(row_buttons)
            self.vvars.append(row_vars)
        # Configure row and column weights to make buttons fill the butmatframe
        for row in range(1, n+1):
            butmatframe.rowconfigure(row, weight=1)
        for col in range(1, m+1):
            butmatframe.columnconfigure(col, weight=1)
        # Add row axis
        for row in range(n):
            label = tk.Label(butmatframe, text=str(round(row*self.gdx, 10)), relief=tk.RAISED)
            label.grid(row=row+1, column=0, sticky=tk.NSEW)
        # Add column axis
        for col in range(m):
            label = tk.Label(butmatframe, text=str(round(col*self.gdy, 10)), relief=tk.RAISED)
            label.grid(row=0, column=col+1, sticky=tk.NSEW)
        self.SpecButtons = buttons
        self.buttonframe_updateColor()
        #return buttons

    def button_click(self, var):
        self.newselx = var[0]
        self.newsely = var[1]
        self.updateselectionentries()
    
    def buttonframe_updateColor(self):
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                if type(self.SpecDataMatrix[i][j]) == SpectrumData:
                    if self.SpecDataMatrix[i][j].dataokay == True:
                        self.SpecButtons[i][j].config(bg="green")
                    else:
                        self.SpecButtons[i][j].config(bg="red")
                else:
                    self.SpecButtons[i][j].config(bg="red")
    
    def get_normalization_matrix(self):
        """
        Generate normalization matrix based on current GUI settings.
        
        Returns:
            2D numpy array: Normalization matrix or None if normalization is disabled
        """
        # Check if normalization is enabled
        if not self.normalize_HSI_var.get():
            return None
        
        # Get normalization method
        method_name = self.normalize_method_var.get()
        
        # Map GUI names to internal method names
        method_map = {
            'None': 'none',
            'Integrated Counts': 'integrated_counts',
            'Max Intensity': 'max_intensity',
            'Counts at Wavelength': 'counts_at_wavelength',
            'Normalize intern': 'normalize_intern'
        }
        
        method = method_map.get(method_name, 'none')
        
        # If method is 'none', return None (no normalization)
        if method == 'none':
            return None
        
        # Get selected data key from normalize_dataset_box (independent of plot dataset)
        # This allows normalizing e.g., derivatives by counts from the original signal
        selected_data_name = self.normalize_dataset_var.get()
        data_key = self.speckeys.get(selected_data_name, 'PLB')
        
        # Prepare parameters based on method
        params = {'data_key': data_key}
        
        if method == 'integrated_counts' or method == 'max_intensity':
            # Use current wavelength range
            params['wl_start'] = self.wlstart # type: ignore
            params['wl_end'] = self.wlend # type: ignore
        elif method == 'counts_at_wavelength':
            # Get wavelength from entry
            try:
                params['wavelength'] = float(self.normalize_wl_var.get()) # type: ignore
            except ValueError:
                print(f"Warning: Invalid normalization wavelength '{self.normalize_wl_var.get()}'. Using mean wavelength.")
                params['wavelength'] = np.mean(self.WL) # type: ignore
        
        # Create normalization object and generate matrix
        try:
            normalizer = hsi_normalization.HSINormalization(self.SpecDataMatrix, np.array(self.WL))
            norm_matrix = normalizer.generate_normalization_matrix(method, params)
            return norm_matrix
        except Exception as e:
            print(f"Error generating normalization matrix: {e}")
            import traceback
            traceback.print_exc()
            return None

    def updateselectionentries(self):
        self.selectPixX.delete(0, tk.END)
        self.selectPixX.insert(0, str(self.newselx))
        self.selectPixY.delete(0, tk.END)
        self.selectPixY.insert(0, str(self.newsely))

    # Max Counts Colormap
    def buildandPlotIntCmap(self, savetoimage='False', plot=True, datatype='Spectrum (PL-BG)', wlstart=None, wlend=None):
        hsi_name = self.hsiselect.get()
        if not hsi_name or hsi_name not in self.PMdict:
            print("No valid HSI selected to plot.")
            return
            
        self.readfontsize()
        self.updatecountthresh()
        # update spec min and max values # self.proc_spec_min
        # if wlstart and wlend are provided, write them to self.proc_spec_min and self.proc_spec_max as entries
        if wlstart is not None:
            self.proc_spec_min.delete(0, tk.END)
            self.proc_spec_min.insert(0, str(wlstart))
        if wlend is not None:
            self.proc_spec_max.delete(0, tk.END)
            self.proc_spec_max.insert(0, str(wlend))
        self.updatewl()
        # create a new colormap by using the selected HSI
        lastpm = copy.deepcopy(self.PMdict[hsi_name].PixMatrix)
        
        metadata = {}
        # set metadata for the new colormap
        if datatype in deflib.datatype2unit:
            metadata['unit'] = deflib.datatype2unit[datatype]
            metadata['quantity'] = 'signal'
            metadata['fitmodel'] = None
            metadata['fitparameter'] = None
        else: 
            metadata['unit'] = 'A. U.'
            metadata['quantity'] = 'signal'
            metadata['fitmodel'] = None
            metadata['fitparameter'] = None
        newpm = self.writetopixmatrix(lastpm, None, metadata=metadata, method='Integrated')
        self.getPLpixelIntervalMaxIndex(self.PMdict[newpm].PixMatrix, False)
        
        # Apply normalization if enabled
        norm_matrix = self.get_normalization_matrix()
        if norm_matrix is not None:
            self.PMdict[newpm].PixMatrix = hsi_normalization.HSINormalization.apply_normalization(
                self.PMdict[newpm].PixMatrix, norm_matrix
            )
        
        self.UpdateHSIselect()
        if plot:
            self.plotPixelMatrix(self.hsiselect.get(), savetoimage=savetoimage)
        
    # Spectral Maximum Colormap
    def buildandPlotSpecCmap(self):
        hsi_name = self.hsiselect.get()
        if not hsi_name or hsi_name not in self.PMdict:
            print("No valid HSI selected to plot.")
            return

        self.updatewl()
        self.updatecountthresh()
        self.readfontsize()
        lastpm = copy.deepcopy(self.PMdict[hsi_name].PixMatrix)
        # before create newpm: create meatadata
        metadata = {}
        if self.selectwindowbox.get() in matl.fitkeys:
            metadata['unit'] = 'eV'
            metadata['quantity'] = 'wavelength'
            metadata['fitmodel'] = self.selectwindowbox.get()
            metadata['fitparameter'] = None
        else:
            metadata['unit'] = 'A. U.'
            metadata['quantity'] = 'signal'
            metadata['fitmodel'] = None
            metadata['fitparameter'] = None
        newpm = self.writetopixmatrix(lastpm, metadata=metadata, method='spectral_shift')#str(self.selectspecpixbox.get()))
        if self.HSI_fit_useROI.get() == False:
            self.fittoMatrixfitparams(self.PMdict[newpm].PixMatrix, 'fitmaxX', mode='fullHSI', roi=None)
        else:
            # get the selected ROI mask
            if self.roisel.get() == '':
                print('Error', 'No ROI selected. Please select a ROI first.')
                return
            elif self.roisel.get() not in self.roihandler.roilist.keys():
                print('Error', 'Selected ROI not found. Please select a valid ROI.')
                return
            else: # if the ROI is valid, copy it onto roi
                roi = self.roihandler.roilist[self.roisel.get()][:]
            # use the ROI mask to fit the pixel matrix
            self.fittoMatrixfitparams(self.PMdict[newpm].PixMatrix, 'fitmaxX', mode='roi', roi=roi)
        self.getPLpixelSpecMax(self.PMdict[newpm].PixMatrix)
        # update matadata of newpm: 'unit', 'quantity', 'fitmodel', 'fitparams'
        # 1. set unit according to selected fit function mathlib.fitkeys[self.selectwindowbox.get()][6][0]
        if self.selectwindowbox.get() in matl.fitkeys:
            self.PMdict[newpm].metadata['unit'] = 'eV'
            self.PMdict[newpm].metadata['quantity'] = 'wavelength'
            self.PMdict[newpm].metadata['fitmodel'] = self.selectwindowbox.get()
            self.PMdict[newpm].metadata['fitparameter'] = None
        else:
            self.PMdict[newpm].metadata['unit'] = 'A. U.'
        
        # Apply normalization if enabled
        norm_matrix = self.get_normalization_matrix()
        if norm_matrix is not None:
            self.PMdict[newpm].PixMatrix = hsi_normalization.HSINormalization.apply_normalization(
                self.PMdict[newpm].PixMatrix, norm_matrix
            )
        
        self.UpdateHSIselect()
    
        try:
            self.plotPixelMatrix(self.hsiselect.get())
        except Exception as e:
            print('Plot failed. {}'.format(str(e)))
            plt.cla()
    
    def updateandPlotSpecCmap(self, variable):
        #self.updatePixelMatrix(variable)
        self.plotPixelMatrix(self.hsiselect.get())
        self.UpdateHSIselect()

    def fitwindowtospec(self, variable, newfit=False):
        self.updatewl()
        x, y, valid = self.validpixelinput()
        data = None
        if valid[0] == True and valid[1] == True:
            if type(self.SpecDataMatrix[y][x]) == SpectrumData:
                if self.SpecDataMatrix[y][x].dataokay == True:
                    self.selectwindowboxVari = self.selectwindowbox.get()
                    self.selectspecboxVari = self.selectspecbox.get()
                    if self.speckeys[self.selectspecboxVari] == 'WL': #Wavelength
                        data = self.SpecDataMatrix[y][x].WL
                        self.PlotSpectrum(data, self.SpecDataMatrix[y][x].WL_eV, 'Wavelength', xunit='eV')
                    elif self.speckeys[self.selectspecboxVari] == 'BG': #Background
                        data = self.SpecDataMatrix[y][x].BG
                        self.PlotSpectrum(data, self.SpecDataMatrix[y][x].WL_eV, 'Background Counts', xunit='eV')
                    elif self.speckeys[self.selectspecboxVari] == 'PL': # Counts
                        data = self.SpecDataMatrix[y][x].PL
                        self.PlotSpectrum(data, self.SpecDataMatrix[y][x].WL_eV, 'Spectrometer Counts', xunit='eV')
                    elif self.speckeys[self.selectspecboxVari] == 'PLB': #Spectrum
                        data = self.SpecDataMatrix[y][x].PLB[self.aqpixstart: self.aqpixend]
                    try: # fit function to spectrum  
                        if newfit == True or self.SpecDataMatrix[y][x].fitdata == [None]:
                            if self.SpecDataMatrix[y][x].fitdata == [None]:
                                print('No fit data found. Fitting new data.')
                            try: 
                                self.maxiter = int(self.fitmaxiter.get())
                            except:
                                print('Maxiter must be type int. Using default 1000.')
                                self.maxiter = 1000
                            self.SpecDataMatrix[y][x].fitdata = self.fitkeys[self.selectwindowboxVari][1](self.aqpixstart, self.aqpixend, self.SpecDataMatrix[y][x].WL_eV, self.SpecDataMatrix[y][x].PLB, self.maxiter)
                            x_window = self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart:self.aqpixend]
                            fit_xmin = float(np.min(x_window))
                            fit_xmax = float(np.max(x_window))
                            self.SpecDataMatrix[y][x].fitmaxX, self.SpecDataMatrix[y][x].fitmaxY = self.fitkeys[self.selectwindowboxVari][2](fit_xmin, fit_xmax, *self.SpecDataMatrix[y][x].fitdata[:-1])
                        
                            #plt.scatter(self.SpecDataMatrix[y][x].fitmaxX, self.SpecDataMatrix[y][x].fitmaxY, color='red')
                        self.PlotFitSpectrum(self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart: self.aqpixend], data, ['Spectrometer counts', self.fitkeys[self.selectwindowboxVari][3]], [self.SpecDataMatrix[y][x].fitdata[:-1]], [self.fitkeys[self.selectwindowboxVari][0]])
                        
                    except Exception as e:
                        print('Fit failed. {}'.format(str(e)))
        else:
            print(self.SpecDataMatrix[y][x])  

    def execute_fit_with_error_handling(self, spec_obj, fit_func, aqpixstart, aqpixend, maxiter, fitbackup):
        """
        Execute a fit function with proper error handling and return status.
        
        Returns:
            tuple: (fitdata, fit_status) where fit_status is:
                   0 = not fitted
                   1 = fitted successfully
                   2 = fit failed (convergence error or other exception)
        """
        try:
            fitdata = fit_func(
                aqpixstart,
                aqpixend,
                spec_obj.WL_eV,
                spec_obj.PLB,
                maxiter,
                fitbackup
            )
            # Check if fitdata is valid (not None and not [None])
            if fitdata is None or fitdata == [None]:
                return [None], 2  # Fit failed
            return fitdata, 1  # Fit succeeded
        except Exception as e:
            # Catch curve_fit convergence errors and other exceptions
            print(f'Fit failed for spectrum: {e}')
            return [None], 2  # Fit failed

    def fittoMatrixfitparams(self, PixMatrix, variable='fitmaxX', incmin=0.01, incmax=-0.01, nmin=20, nmax=20, mode='fullHSI', roi=None):
        # init empty self.fitbakup
        self.fitbackup = None

        if mode not in ['fullHSI', 'roi']:
            raise ValueError("Mode must be 'fullHSI' or 'roi'.")
        if mode == 'roi' and roi is None:
            # check if roi is 2D array
            if not isinstance(roi, np.ndarray) or roi.ndim != 2:
                raise ValueError("roi must be a 2D numpy array. Got: {}".format(type(roi)))        

        # fill matrix with data of the selected enry:
        self.updatecountthresh()
        # update wlstart and wlend
        self.selectwindowboxVari = self.selectwindowbox.get()
        self.selectspecboxVari = self.selectspecbox.get()
        
        # Precompute lookups and indices outside the loop
        # Check if maxiter is available, otherwise use default
        try:
            precomputed_maxiter = int(self.fitmaxiter.get())
        except:
            print('Maxiter must be int. Using default 1000.')
            precomputed_maxiter = 1000
            
        a_index = list(matl.fitkeys.keys()).index(self.selectwindowboxVari)
        num_fit_params = matl.fitkeys[list(matl.fitkeys.keys())[a_index]][4]
        
        idx_max_x = matl.addtofitparms.index('max_x') - len(matl.addtofitparms)
        idx_max_y = matl.addtofitparms.index('max_y') - len(matl.addtofitparms)
        idx_pixstart = matl.addtofitparms.index('pixstart') - len(matl.addtofitparms)
        idx_pixend = matl.addtofitparms.index('pixend') - len(matl.addtofitparms)
        idx_wlstart = matl.addtofitparms.index('wlstart') - len(matl.addtofitparms)
        idx_wlend = matl.addtofitparms.index('wlend') - len(matl.addtofitparms)
        idx_fwhm = matl.addtofitparms.index('fwhm') - len(matl.addtofitparms)
        idx_r_squared = matl.addtofitparms.index('r_squared') - len(matl.addtofitparms)
        idx_ss_res = matl.addtofitparms.index('ss_res') - len(matl.addtofitparms)
        idx_ss_tot = matl.addtofitparms.index('ss_tot') - len(matl.addtofitparms)
        idx_fit_status = matl.addtofitparms.index('fit_status') - len(matl.addtofitparms)
        
        fit_func = self.fitkeys[self.selectwindowboxVari][1]
        fit_calc_max = self.fitkeys[self.selectwindowboxVari][2]
        fit_calc_fwhm = matl.fitkeys[list(matl.fitkeys.keys())[a_index]][7]
        fit_eval = self.fitkeys[self.selectwindowboxVari][0]

        # Update GUI elements definitively once before threading
        self.updatewl()
        base_aqstart = self.aqpixstart
        base_aqend = self.aqpixend
        
        def process_pixel(i, j):
            if type(self.SpecDataMatrix[i][j]) != SpectrumData or not self.SpecDataMatrix[i][j].dataokay:
                return (i, j, False)
                
            self.SpecDataMatrix[i][j].dofit = True
            tries = 1
            worked = False
            adjmin = True
            
            # Start with the base values to prevent race conditions on self
            local_aqstart = base_aqstart
            local_aqend = base_aqend
            local_aqstart_ev = self.SpecDataMatrix[i][j].WL_eV[local_aqstart]
            local_aqend_ev = self.SpecDataMatrix[i][j].WL_eV[local_aqend]
            local_fitbackup = None  # Using None avoids thread contention and bad p0 jumping
            
            while tries < nmin+nmax and not worked and self.SpecDataMatrix[i][j].dofit:
                tries += 1
                skip_fit = False
                if mode == 'fullHSI':
                    if np.sum(self.SpecDataMatrix[i][j].PLB[int(local_aqstart):int(local_aqend)]) < self.countthreshv:
                        skip_fit = True
                elif mode == 'roi':
                    if np.isnan(roi[i][j]): # type: ignore
                        skip_fit = True
                else:
                    print('Mode not implemented yet.')
                    skip_fit = True
                    
                if skip_fit:
                    PixMatrix[i][j] = np.nan
                    self.SpecDataMatrix[i][j].dofit = False
                    break
                
                # Execute fit with local subset to ensure thread safety
                fitdata, fit_status = self.execute_fit_with_error_handling(
                    self.SpecDataMatrix[i][j],
                    fit_func,
                    int(local_aqstart),
                    int(local_aqend),
                    precomputed_maxiter,
                    local_fitbackup
                )
                self.SpecDataMatrix[i][j].fitdata = fitdata
                
                # Only process if fit succeeded
                if fit_status == 1 and self.SpecDataMatrix[i][j].fitdata != [None]:
                    x_window = self.SpecDataMatrix[i][j].WL_eV[int(local_aqstart):int(local_aqend)]
                    fit_xmin = float(np.min(x_window))
                    fit_xmax = float(np.max(x_window))
                    fitmaxX, fitmaxY = fit_calc_max(fit_xmin, fit_xmax, *self.SpecDataMatrix[i][j].fitdata[:-1])
                    self.SpecDataMatrix[i][j].fitmaxX = fitmaxX
                    self.SpecDataMatrix[i][j].fitmaxY = fitmaxY
                    r_squared, ss_res, ss_tot = matl.calc_r_squared(
                        self.SpecDataMatrix[i][j].PLB[int(local_aqstart):int(local_aqend)], 
                        fit_eval(self.SpecDataMatrix[i][j].WL_eV[int(local_aqstart):int(local_aqend)], *self.SpecDataMatrix[i][j].fitdata[:-1])
                    )
                    
                    try:
                        for k in range(num_fit_params):
                            self.SpecDataMatrix[i][j].fitparams[a_index][k] = self.SpecDataMatrix[i][j].fitdata[k]
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_max_x] = fitmaxX
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_max_y] = fitmaxY
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_pixstart] = local_aqstart
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_pixend] = local_aqend
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_wlstart] = self.wlstart
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_wlend] = self.wlend
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_fwhm] = fit_calc_fwhm(self.SpecDataMatrix[i][j].fitdata[:num_fit_params])
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_r_squared] = r_squared
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_ss_res] = ss_res
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_ss_tot] = ss_tot
                        self.SpecDataMatrix[i][j].fitparams[a_index][idx_fit_status] = fit_status
                        local_fitbackup = self.SpecDataMatrix[i][j].fitdata
                    except Exception as e:
                        print('Fit parameter update failed in function XYMap.fittoMatrixfitparams:', e)
                else:
                    self.SpecDataMatrix[i][j].fitparams[a_index][idx_fit_status] = fit_status
                    PixMatrix[i][j] = np.nan
                    
                if self.SpecDataMatrix[i][j].fitdata == [None] or self.SpecDataMatrix[i][j].fitdata is None:
                    PixMatrix[i][j] = np.nan
                else:
                    wl_min = float(np.min(self.SpecDataMatrix[i][j].WL_eV[int(local_aqstart):int(local_aqend)]))
                    wl_max = float(np.max(self.SpecDataMatrix[i][j].WL_eV[int(local_aqstart):int(local_aqend)]))

                    if wl_min <= getattr(self.SpecDataMatrix[i][j], 'fitmaxX', -1) <= wl_max:
                        worked = True
                        PixMatrix[i][j] = self.SpecDataMatrix[i][j].get_attribute(variable)
                
                if not worked:
                    try:
                        if adjmin:
                            local_aqstart += incmin
                            adjmin = False
                        else:
                            local_aqend += incmax
                            adjmin = True
                    except Exception:
                        worked = True
                        
            return (i, j, worked)

        # Build thread pool to evaluate all pixels concurrently
        workers = min(32, (os.cpu_count() or 1) + 4)
        futures = []
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    futures.append(executor.submit(process_pixel, i, j))
        
        # Await completion to ensure matrix is completely populated before proceeding
        for future in as_completed(futures):
            future.result()

    # function currently not in use
    def fittoSpecfitparams(self, PixMatrix, variable='fitmaxX', incmin=2, incmax=-2, nmin=20, nmax=20):
        # fill matrix with data of the selected enry
        self.updatewl()
        x, y, valid = self.validpixelinput()
        if valid[0] == True and valid[1] == True:
            if type(self.SpecDataMatrix[y][x]) == SpectrumData:
                if self.SpecDataMatrix[y][x].dataokay == True:
                    self.selectwindowboxVari = self.selectwindowbox.get()
                    self.selectspecboxVari = self.selectspecbox.get()
                    tries = 1
                    worked = False
                    adjmin = True
                    # if fit does not work, adjust the window size
                    while tries < nmin+nmax and worked == False:
                        tries += 1
                        try:
                            if self.speckeys[self.selectspecboxVari] == 'PLB': #Spectrum
                                if np.sum(self.SpecDataMatrix[y][x].PLB[self.aqpixstart:self.aqpixend]) < self.countthreshv:
                                    PixMatrix[y][x] = np.nan
                                else:
                                    try:
                                        self.maxiter = int(self.fitmaxiter.get())
                                    except:
                                        print('Maxiter must be int. Using default 1000.')
                                        self.maxiter = 1000
                                    self.SpecDataMatrix[y][x].fitdata = self.fitkeys[self.selectwindowboxVari][1](self.aqpixstart, self.aqpixend, self.SpecDataMatrix[y][x].WL, self.SpecDataMatrix[y][x].PLB, self.maxiter)
                                    x_window = self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart:self.aqpixend]
                                    fit_xmin = float(np.min(x_window))
                                    fit_xmax = float(np.max(x_window))
                                    self.SpecDataMatrix[y][x].fitmaxX, self.SpecDataMatrix[y][x].fitmaxY = self.fitkeys[self.selectwindowboxVari][2](fit_xmin, fit_xmax, *self.SpecDataMatrix[y][x].fitdata[:-1])#[1]
                                    r_squared, ss_res, ss_tot = matl.calc_r_squared(self.SpecDataMatrix[y][x].PLB[self.aqpixstart:self.aqpixend], self.fitkeys[self.selectwindowboxVari][0](self.SpecDataMatrix[y][x].WL[self.aqpixstart:self.aqpixend], *self.SpecDataMatrix[y][x].fitdata[:-1]))
                                    a =  list(matl.fitkeys.keys()).index(self.selectwindowbox.get())
                                    # put fitdata into fitparams
                                    try:
                                        for k in range(matl.fitkeys[list(matl.fitkeys.keys())[a]][4]):
                                            self.SpecDataMatrix[y][x].fitparams[a][k] = self.SpecDataMatrix[y][x].fitdata[k]
                                        # store fitmaxX and fitmaxY in fitparams[-1] and [-2]
                                        self.SpecDataMatrix[y][x].fitparams[a][-1] = self.SpecDataMatrix[y][x].fitmaxX
                                        self.SpecDataMatrix[y][x].fitparams[a][-2] = self.SpecDataMatrix[y][x].fitmaxY
                                        # store aqpixstart and aqpixend in fitparams[-3] and [-4]
                                        self.SpecDataMatrix[y][x].fitparams[a][-3] = self.aqpixstart
                                        self.SpecDataMatrix[y][x].fitparams[a][-4] = self.aqpixend
                                        # store r_squared, ss_res, ss_tot in fitparams[-5] and [-6] and [-7]
                                        self.SpecDataMatrix[y][x].fitparams[a][matl.addtofitparms.index('r_squared')-len(matl.addtofitparms)+1] = r_squared
                                        self.SpecDataMatrix[y][x].fitparams[a][matl.addtofitparms.index('ss_res')-len(matl.addtofitparms)+1] = ss_res
                                        self.SpecDataMatrix[y][x].fitparams[a][matl.addtofitparms.index('ss_tot')-len(matl.addtofitparms)+1] = ss_tot
                                    except Exception as e:
                                        print('Fit parameter update failed in new fitline in function {}. {}'.format(self.__class__.__name__, str(e)))
                        
                                if self.SpecDataMatrix[y][x].fitdata == [None]:
                                    PixMatrix[y][x] = np.nan 
                                else:
                                    # check if maximum is within the window and set to Intmatrix
                                    if self.SpecDataMatrix[y][x].fitmaxX >= float(np.min(self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart:self.aqpixend])) and self.SpecDataMatrix[y][x].fitmaxX <= float(np.max(self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart:self.aqpixend])):
                                        worked = True
                                        PixMatrix[y][x] = self.SpecDataMatrix[y][x].get_attribute(variable)
                                    else:
                                        pass
                                        #print(self.SpecDataMatrix[y][x].fitmaxX, self.SpecDataMatrix[y][x].fitmaxY)
                        except Exception as e:
                            try:
                                if adjmin == True:
                                    self.aqpixstart += incmin
                                    adjmin = False
                                else:
                                    self.aqpixend += incmax
                                    adjmin = True
                            except:
                                print("Fit Window ran out of Data. Fit to Matrix Failed at element {}, {}.\n{}".format(x, y, self.__class__.__name__, str(e)))
                                worked = True
                            print("Fit Window ran out of Data. Fit to Matrix Failed at element {}, {}.\n{}".format(x, y, str(e)))
                        self.updatewl()   

        # fill matrix with data of the selected enry:
        self.updatecountthresh()
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                if type(self.SpecDataMatrix[i][j]) == SpectrumData:
                    if self.SpecDataMatrix[i][j].dataokay == True:
                        self.selectwindowboxVari = self.selectwindowbox.get()
                        self.selectspecboxVari = self.selectspecbox.get()
                        tries = 1
                        worked = False
                        adjmin = True
                        # if fit does not work, adjust the window size
                        while tries < nmin+nmax and worked == False:
                            tries += 1
                            try:
                                if self.speckeys[self.selectspecboxVari] == 'PLB': #Spectrum
                                    if np.sum(self.SpecDataMatrix[i][j].PLB[self.aqpixstart:self.aqpixend]) < self.countthreshv:
                                        PixMatrix[i][j] = np.nan
                                    else:
                                        try:
                                            self.maxiter = int(self.fitmaxiter.get())
                                        except:
                                            print('Maxiter must be int. Using default 1000.')
                                            self.maxiter = 1000
                                        self.SpecDataMatrix[i][j].fitdata = self.fitkeys[self.selectwindowboxVari][1](self.aqpixstart, self.aqpixend, self.SpecDataMatrix[i][j].WL, self.SpecDataMatrix[i][j].PLB, self.maxiter)
                                        self.SpecDataMatrix[i][j].fitmaxX, self.SpecDataMatrix[i][j].fitmaxY = self.fitkeys[self.selectwindowboxVari][2](self.aqpixstart, self.aqpixend, *self.SpecDataMatrix[i][j].fitdata[:-1])#[1]
                                        r_squared, ss_res, ss_tot = matl.calc_r_squared(self.SpecDataMatrix[i][j].PLB[self.aqpixstart:self.aqpixend], self.fitkeys[self.selectwindowboxVari][0](self.SpecDataMatrix[i][j].WL[self.aqpixstart:self.aqpixend], *self.SpecDataMatrix[i][j].fitdata[:-1]))
                                        a =  list(matl.fitkeys.keys()).index(self.selectwindowbox.get())
                                        # put fitdata into fitparams
                                        try:
                                            for k in range(matl.fitkeys[list(matl.fitkeys.keys())[a]][4]):
                                                self.SpecDataMatrix[i][j].fitparams[a][k] = self.SpecDataMatrix[i][j].fitdata[k]
                                            # store fitmaxX and fitmaxY in fitparams[-1] and [-2]
                                            self.SpecDataMatrix[i][j].fitparams[a][-1] = self.SpecDataMatrix[i][j].fitmaxX
                                            self.SpecDataMatrix[i][j].fitparams[a][-2] = self.SpecDataMatrix[i][j].fitmaxY
                                            # store aqpixstart and aqpixend in fitparams[-3] and [-4]
                                            self.SpecDataMatrix[i][j].fitparams[a][-3] = self.aqpixstart
                                            self.SpecDataMatrix[i][j].fitparams[a][-4] = self.aqpixend
                                            # store r_squared in fitparams [-7], [-8], [-9]
                                            self.SpecDataMatrix[i][j].fitparams[a][matl.addtofitparms.index('r_squared')-len(matl.addtofitparms)+1] = r_squared
                                            self.SpecDataMatrix[i][j].fitparams[a][matl.addtofitparms.index('ss_res')-len(matl.addtofitparms)+1] = ss_res
                                            self.SpecDataMatrix[i][j].fitparams[a][matl.addtofitparms.index('ss_tot')-len(matl.addtofitparms)+1] = ss_tot
                                        except Exception as e:
                                            print('Fit parameter update failed in new fitline. {}'.format(str(e)))

                                if self.SpecDataMatrix[i][j].fitdata == [None]:
                                    PixMatrix[i][j] = np.nan 
                                else:
                                    # check if maximum is within the window and set to Intmatrix
                                    if self.SpecDataMatrix[i][j].fitmaxX >= float(np.min(self.SpecDataMatrix[i][j].WL_eV[self.aqpixstart:self.aqpixend])) and self.SpecDataMatrix[i][j].fitmaxX <= float(np.max(self.SpecDataMatrix[i][j].WL_eV[self.aqpixstart:self.aqpixend])):
                                        worked = True
                                        PixMatrix[i][j] = self.SpecDataMatrix[i][j].get_attribute(variable)
                                    else:
                                        pass
                                        #print(self.SpecDataMatrix[i][j].fitmaxX, self.SpecDataMatrix[i][j].fitmaxY)
                            except Exception as e:
                                try:
                                    if adjmin == True:
                                        self.aqpixstart += incmin
                                        adjmin = False
                                    else:
                                        self.aqpixend += incmax
                                        adjmin = True
                                except:
                                    print("Fit Window ran out of Data. Fit to Matrix Failed at element {}, {} in function {}.\n{}".format(i, j, 'XYMap.updatecountthresh', str(e)))
                                    worked = True
                                print("Fit to Matrix Failed at element {}, {} in function {}.\n{}".format(i, j, 'XYMap.updatecountthresh', str(e)))     
                            self.updatewl()
                                
    # function curently not in use
    def updatePixelMatrix(self, PixMatrix, variable, nonecase=np.nan):
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                try:
                    PixMatrix[i][j] = self.SpecDataMatrix[i][j].get_attribute(variable)
                except Exception as e:
                    PixMatrix[i][j] = nonecase
                    print("Update Pixel Matrix Failed at element {}, {}.\n{}".format(i, j, str(e)))

    def readfontsize(self):
        try:
            self.fontsize =abs(float(self.CMFont.get()))
        except Exception as e:
            print("Error", '{} Font Size must be Number.'.format(str(e)))

    def validpixelinput(self):
        x = self.selectPixX.get()
        y = self.selectPixY.get()
        valid = [False, False]
        try:
            x = int(x)
            if x < len(self.SpecDataMatrix[0]):
                valid[0] = True
            else:
                print("Error", "No Pixel on X-Position.")
        except Exception as e:
            print("Error", '{} Insert valid X-Position.'.format(str(e)))
        try:
            y = int(y)
            if y < len(self.SpecDataMatrix):
                valid[1] = True
            else:
                print("Error", "No Pixel on Y-Position.")
        except Exception as e:
            print("Error", '{} Insert valid Y-Position.'.format(str(e)))
        return(int(x), int(y), valid)

    # Function to handle click events of the image
    def on_click(self, event, image):
        if event.inaxes:
            self.newselx = round(event.xdata)
            self.newsely = round(event.ydata)
            self.updateselectionentries()
    
    def WL2selectedunit(self, i, j):
        if self.WL_selection.get() == 'Wavelength (nm)':
            self.WLunit = 'nm'
            return self.SpecDataMatrix[i][j].WL
        elif self.WL_selection.get() == 'Energy (eV)':
            self.WLunit = 'eV'
            #return deflib.wl_array_to_ev(self.SpecDataMatrix[i][j].WL)
            return self.SpecDataMatrix[i][j].WL_eV
        else: 
            self.WLunit = 'nm'
            print('No valid Wavelength Unit selected. Using nm as default.')
            return self.SpecDataMatrix[i][j].WL

    def PlotPixelSpectrum(self):
        x, y, valid = self.validpixelinput()        
        if valid[0] == True and valid[1] == True:
            if type(self.SpecDataMatrix[y][x]) == SpectrumData:
                if self.SpecDataMatrix[y][x].dataokay == True:
                    wl = self.WL2selectedunit(y, x)
                    self.selectdataboxVari = self.selectspecpixbox.get()
                    if self.speckeys[self.selectdataboxVari] == 'WL': #Wavelength
                        self.PlotSpectrum(wl, wl, 'Wavelength', xunit=self.WLunit)
                    elif self.speckeys[self.selectdataboxVari] == 'BG': #Background
                        self.PlotSpectrum(self.SpecDataMatrix[y][x].BG, wl, 'Background Counts', xunit=self.WLunit)
                    elif self.speckeys[self.selectdataboxVari] == 'PL': # Counts
                        self.PlotSpectrum(self.SpecDataMatrix[y][x].PL, wl, 'Spectrometer Counts', xunit=self.WLunit)
                    elif self.speckeys[self.selectdataboxVari] == 'PLB': #Spectrum
                        self.PlotSpectrum(self.SpecDataMatrix[y][x].PLB, wl, 'PL Spectrum', xunit=self.WLunit)
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff1': # First Derivative
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff1') and self.SpecDataMatrix[y][x].Specdiff1 is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff1, wl, '1st Derivative', xunit=self.WLunit)
                        else:
                            print("First derivative not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff2': # Second Derivative
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff2') and self.SpecDataMatrix[y][x].Specdiff2 is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff2, wl, '2nd Derivative', xunit=self.WLunit)
                        else:
                            print("Second derivative not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff1_norm': # First Derivative (normalized)
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff1_norm') and self.SpecDataMatrix[y][x].Specdiff1_norm is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff1_norm, wl, '1st Derivative (normalized)', xunit=self.WLunit, yunit='1/nm')
                        else:
                            print("Normalized first derivative not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff2_norm': # Second Derivative (normalized)
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff2_norm') and self.SpecDataMatrix[y][x].Specdiff2_norm is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff2_norm, wl, '2nd Derivative (normalized)', xunit=self.WLunit, yunit='1/nm²')
                        else:
                            print("Normalized second derivative not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff1_norm_intensity': # First Derivative (norm on intensity, then derive)
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff1_norm_intensity') and self.SpecDataMatrix[y][x].Specdiff1_norm_intensity is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff1_norm_intensity, wl, '1st Derivative (norm on intensity, then derive)', xunit=self.WLunit, yunit='1/nm')
                        else:
                            print("First derivative (norm on intensity) not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff2_norm_intensity': # Second Derivative (norm on intensity, then derive)
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff2_norm_intensity') and self.SpecDataMatrix[y][x].Specdiff2_norm_intensity is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff2_norm_intensity, wl, '2nd Derivative (norm on intensity, then derive)', xunit=self.WLunit, yunit='1/nm²')
                        else:
                            print("Second derivative (norm on intensity) not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff1_norm_counts': # First Derivative (norm on counts, then derive)
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff1_norm_counts') and self.SpecDataMatrix[y][x].Specdiff1_norm_counts is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff1_norm_counts, wl, '1st Derivative (norm on counts, then derive)', xunit=self.WLunit, yunit='1/nm')
                        else:
                            print("First derivative (norm on counts) not available for this pixel.")
                    elif self.speckeys[self.selectdataboxVari] == 'Specdiff2_norm_counts': # Second Derivative (norm on counts, then derive)
                        if hasattr(self.SpecDataMatrix[y][x], 'Specdiff2_norm_counts') and self.SpecDataMatrix[y][x].Specdiff2_norm_counts is not None:
                             self.PlotSpectrum(self.SpecDataMatrix[y][x].Specdiff2_norm_counts, wl, '2nd Derivative (norm on counts, then derive)', xunit=self.WLunit, yunit='1/nm²')
                        else:
                            print("Second derivative (norm on counts) not available for this pixel.")
                    else:
                        print('No valid Data set selected for the Plot.')
                    print('x-position: {}, y-position: {}'.format(x, y))

    def PlotSpectrum(self, x, y, label, xunit='nm', yunit='counts'):
        self.readfontsize()
        plt.figure(figsize=(10, 6))
        plt.plot(y, x, label=label, color='blue')       
        plt.xlabel(f'Wavelength / {xunit}', fontsize=self.fontsize)
        plt.ylabel(f'Intensity / {yunit}', fontsize=self.fontsize)
        plt.title('Spectrograph Data', fontsize=self.fontsize)
        plt.legend(fontsize=self.fontsize)
        plt.tick_params(axis='both', which='major', labelsize=self.fontsize)
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def runPlotFitSpectrum(self):
        # collect entries and run PlotFitSpectrum
        self.updatewl()
        data = None
        self.selectspecboxVari = self.selectspecbox.get()
        x, y, valid = self.validpixelinput()
        self.selectwindowboxVari = self.selectwindowbox.get()
        a_index = list(matl.fitkeys.keys()).index(self.selectwindowboxVari)
        # get the fitdate of the selected selectspecboxVari and the selected pixel (x, y) 
        """
        a_index = list(matl.fitkeys.keys()).index(self.selectwindowboxVari)
        self.SpecDataMatrix[i][j].fitparams[a_index][k] = self.SpecDataMatrix[i][j].fitdata[k]
        """
        #fitdata = self.SpecDataMatrix[y][x].fitparams[self.selectwindowboxVari]
        # wl is required in eV:
        self.WLunit = 'eV'
        if valid[0] == True and valid[1] == True:
            if self.speckeys[self.selectspecboxVari] == 'WL': #Wavelength
                data = self.SpecDataMatrix[y][x].WL_eV
                self.PlotSpectrum(data, self.SpecDataMatrix[y][x].WL_eV, 'Wavelength', xunit='eV')
            elif self.speckeys[self.selectspecboxVari] == 'BG': #Background
                data = self.SpecDataMatrix[y][x].BG
                self.PlotSpectrum(data, self.SpecDataMatrix[y][x].WL_eV, 'Background Counts', xunit='eV')
            elif self.speckeys[self.selectspecboxVari] == 'PL': # Counts
                data = self.SpecDataMatrix[y][x].PL
                self.PlotSpectrum(data, self.SpecDataMatrix[y][x].WL_eV, 'Spectrometer Counts', xunit='eV')
            elif self.speckeys[self.selectspecboxVari] == 'PLB': #Spectrum
                data = self.SpecDataMatrix[y][x].PLB[self.aqpixstart: self.aqpixend]
            ''' old version, not working with new fitdata structure
            self.PlotFitSpectrum(self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart: self.aqpixend], 
                                 data, 
                                 ['', self.fitkeys[self.selectwindowboxVari][3]], 
                                 [self.SpecDataMatrix[y][x].fitdata[:-1]], 
                                 [self.fitkeys[self.selectwindowboxVari][0]]
                                 )
            new version, working with new fitdata structure: use fitparameters of the selected windowfunction
            ''' 
            print('Plotting fit spectrum for pixel ({}, {}) with fitdata: {}'.format(x, y, self.SpecDataMatrix[y][x].fitparams[a_index][:-1]))
            self.PlotFitSpectrum(self.SpecDataMatrix[y][x].WL_eV[self.aqpixstart: self.aqpixend], 
                                 data, 
                                 ['', self.fitkeys[self.selectwindowboxVari][3]], 
                                 [self.SpecDataMatrix[y][x].fitparams[a_index][:-1]], 
                                 [self.fitkeys[self.selectwindowboxVari][0]]
                                 )
            
    def PlotFitSpectrum(self, x, y, label, fitdata, fitfunc):
        self.readfontsize()
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, label=label[0], color='blue')
        #plt.plot(x, fitfunc(*fitdata), label='Fitted function', color='red')
        #plt.plot(x, fitfunc(x, *fitdata), label='Fitted function', color='red')
        self.selectwindowboxVari = self.selectwindowbox.get()
        # update self.fontsize: 
        if self.sepfitfunct.get() == True:
            # plot double window function seperately
            if self.selectwindowboxVari == 'double gaussian':
                plt.plot(x, matl.gaussianwind(x, fitdata[0][0], fitdata[0][1], fitdata[0][2]), label='Gaussian 1', color='red')
                plt.plot(x, matl.gaussianwind(x, fitdata[0][3], fitdata[0][4], fitdata[0][5]), label='Gaussian 2', color='green')
            elif self.selectwindowboxVari == 'double lorentz':
                plt.plot(x, matl.lorentzwind(x, fitdata[0][0], fitdata[0][1], fitdata[0][2]), label='Lorentz 1', color='red')
                plt.plot(x, matl.lorentzwind(x, fitdata[0][3], fitdata[0][4], fitdata[0][5]), label='Lorentz 2', color='green')
            elif self.selectwindowboxVari == 'double voigt':
                # hight version: voigt_height(x, amp, cen, fwhm, eta)
                # width version: voigtwind(x, amp, cen, fwhm, eta)
                plt.plot(x, matl.voigt_height(x, fitdata[0][0], fitdata[0][1], fitdata[0][2], fitdata[0][3]), label='Voigt 1', color='red')
                plt.plot(x, matl.voigt_height(x, fitdata[0][4], fitdata[0][5], fitdata[0][6], fitdata[0][7]), label='Voigt 2', color='green')
            else:
                try:
                    if len(fitdata) == 1:
                        plt.plot(x, fitfunc[0](x, *fitdata[0][:self.fitkeys[self.selectwindowboxVari][4]]), label=label[1], color='red')
                    elif len(fitdata) > 1:
                        for i in range(len(fitdata)):
                            plt.plot(x, fitfunc[i](x, *fitdata[i][:self.fitkeys[self.selectwindowboxVari][4]]), label=label[i+1], color='red')  
                    else:
                        print('No fitdata available for plotting: ', fitdata)
                except:
                    print('Error in plotting fitdata (Seperated Fitkeys): ', fitdata)
            plt.legend(fontsize=self.fontsize)
        else:  
            try:
                if len(fitdata) == 1:
                    plt.plot(x, fitfunc[0](x, *fitdata[0][:self.fitkeys[self.selectwindowboxVari][4]]), label=label[1], color='red')
                elif len(fitdata) > 1:
                    for i in range(len(fitdata)):
                        plt.plot(x, fitfunc[i](x, *fitdata[i][:self.fitkeys[self.selectwindowboxVari][4]]), label=label[i+1], color='red')
                else:
                    print('No fitdata available for plotting: ', fitdata)
            except:
                print('Error in plotting fitdata: ', fitdata)
            plt.legend(fontsize=self.fontsize)
        plt.xlabel('Energy / eV', fontsize=self.fontsize)
        plt.ylabel('Intensity / counts', fontsize=self.fontsize)
        plt.title('Spectrograph Data', fontsize=self.fontsize)
        plt.tick_params(axis='both', which='major', labelsize=self.fontsize)
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def _get_hsi_title(self, HSIname):
        meta = getattr(self.PMdict[HSIname], 'metadata', {})
        if not isinstance(meta, dict):
            meta = {}
        pl_aqpixstart = meta.get('aqpixstart', self.aqpixstart)
        pl_aqpixend = meta.get('aqpixend', self.aqpixend)
        pl_wlstart = meta.get('wlstart', self.wlstart)
        pl_wlend = meta.get('wlend', self.wlend)
        return 'Integrated {0} pixels from {1} nm to {2} nm.'.format(int(pl_aqpixend-pl_aqpixstart+1), pl_wlstart, pl_wlend)

    def plotPixelMatrix(self, HSIname, leglabel='Spectrometer Counts', savetoimage='False'):
        fig, ax = plt.subplots()
        HSIimage = self.PMdict[HSIname].PixMatrix       
        # Display the data as an image with a colormap
        cax = ax.imshow(HSIimage, cmap=self.colormap.get()) # aspect='auto' for cubic image
        # Add a colorbar to the image
        cbar = fig.colorbar(cax, ax=ax)
        # Set the colorbar label
        leglabel = ''
        if 'quantity' in self.PMdict[HSIname].metadata:
            leglabel = self.PMdict[HSIname].metadata['quantity']
        if 'unit' in self.PMdict[HSIname].metadata:
            #leglabel += self.PMdict[HSIname].metadata['unit']
            leglabel = f'{leglabel} / {self.PMdict[HSIname].metadata["unit"]}'
        elif 'unit' not in self.PMdict[HSIname].metadata and leglabel == '':
            leglabel = 'Counts'
        cbar.set_label(leglabel, fontsize=self.fontsize)
        # Set the font size of the colorbar ticks
        cbar.ax.tick_params(labelsize=self.fontsize)
        # Set the plot title
        ax.set_title(self._get_hsi_title(HSIname), fontsize=self.fontsize)
        # Set the x and y axis labels
        ax.set_xlabel('Nanostage X Axis in \u03bcm', fontsize=self.fontsize)
        ax.set_ylabel('Nanostage Y Axis in \u03bcm', fontsize=self.fontsize)
        # set the axis ticks
        xticks = np.arange(0, self.PixAxX[-1]-self.PixAxX[0], self.gdx)
        yticks = np.arange(0, self.PixAxY[-1]-self.PixAxY[0], self.gdy)
        ax.set_xticks(np.arange(len(xticks)))
        ax.set_yticks(np.arange(len(yticks)))
        ax.set_xticklabels(np.round(xticks, 4))
        ax.set_yticklabels(np.round(yticks, 4))
        # automatic displayed ticks with autolocator
        ax.xaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_major_locator(AutoLocator())
        # Set the font size of the ticks on both axes
        ax.tick_params(axis='both', which='major', labelsize=self.fontsize)
        plt.tight_layout()
        # if savetoimage is not False, save the image to the given path
        if savetoimage != 'False':
            try:
                plt.savefig(savetoimage, dpi=600)
                # discard the plot after saving
                plt.close('all')
                print('Image saved to {}'.format(savetoimage))
            except Exception as e:
                print('Could not save image to {}. {}'.format(savetoimage, str(e)))

        else:
            # click event and update "selected Pixels X and Y"
            cid = fig.canvas.mpl_connect('button_press_event', lambda event: self.on_click(event, self.PMdict[self.getPixMatrixSelection(self.hsiselect.get())].PixMatrix))
            self.updateselectionentries()
            fig.canvas.mpl_connect('motion_notify_event', lambda event: deflib.fig_on_hoverevent(event, ax, fig, self.PMdict[self.getPixMatrixSelection(self.hsiselect.get())].PixMatrix, (self.PixAxX[0], self.PixAxX[-1]), (self.PixAxY[0], self.PixAxY[-1])))

            plt.show()

    def getPixMatrixSelection(self, PMname=None):
        if PMname == None:
            PMname = self.hsiselect.get()
        try:
            if PMname not in list(self.PMdict.keys()):
                print('No Pixel Matrix with the given name for PlotPixelMatrixSpectral.')
                PMname = list(self.PMdict.keys())[0]
        except Exception as e:
            print('No Pixel Matrix with the given name for PlotPixelMatrixSpectral. {}'.format(str(e)))
            PMname = list(self.PMdict.keys())[0]
        return PMname

    def plotPixelMatrixSpectral(self):
        PMname = self.getPixMatrixSelection(self.hsiselect.get())
        fig, ax = plt.subplots()
        # Display the data as an image with a colormap
        cax = ax.imshow(self.PMdict[PMname].PixMatrix, cmap=self.colormap.get())#'viridis')
        # Add a colorbar to the image
        cbar = fig.colorbar(cax, ax=ax)
        # Set the colorbar label
        cbar.set_label('Highest Intensity in nm', fontsize=self.fontsize)
        # Set the font size of the colorbar ticks
        cbar.ax.tick_params(labelsize=self.fontsize)
        # Get metadata for correct wavelength representation
        meta = getattr(self.PMdict[PMname], 'metadata', {})
        if not isinstance(meta, dict):
            meta = {}
        pl_wlstart = meta.get('wlstart', self.wlstart)
        pl_wlend = meta.get('wlend', self.wlend)
        # Set the plot title
        ax.set_title('Nanostage Intensity {} nm to {} nm.'.format(pl_wlstart, pl_wlend), fontsize=self.fontsize)
        # Set the x and y axis labels
        ax.set_xlabel('Nanostage X Axis in \u03bcm', fontsize=self.fontsize)
        ax.set_ylabel('Nanostage Y Axis in \u03bcm', fontsize=self.fontsize)
        # set the axis ticks
                # set the axis ticks
        xticks = np.arange(0, self.PixAxX[-1]-self.PixAxX[0], self.gdx)
        yticks = np.arange(0, self.PixAxY[-1]-self.PixAxY[0], self.gdy)
        ax.set_xticks(np.arange(len(xticks)))
        ax.set_yticks(np.arange(len(yticks)))
        ax.set_xticklabels(np.round(xticks, 4))
        ax.set_yticklabels(np.round(yticks, 4))

        # automatic displayed ticks with autolocator
        ax.xaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_major_locator(AutoLocator())
        # Set the font size of the ticks on both axes
        ax.tick_params(axis='both', which='major', labelsize=self.fontsize)
        # Create a legend
        plt.tight_layout()

        # click event
        cid = fig.canvas.mpl_connect('button_press_event', lambda event: deflib.on_click(event, self.PMdict[self.getPixMatrixSelection(self.hsiselect.get())].PixMatrix))
        self.updateselectionentries()

        plt.show()

    def getPLpixelIntervalMaxIndex(self, PixMatrix, makenan=False):
        # fill matrix with data of the selected enry:
        self.updatewl()
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                if type(self.SpecDataMatrix[i][j]) == SpectrumData:
                    try:
                        self.selectspecboxVari = self.selectspecbox.get()
                        if self.speckeys[self.selectspecboxVari] == 'WL': #Wavelength
                            PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].WL[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'BG': #Background
                            PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].BG[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'PL': # Counts
                            PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].PL[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'PLB': #Spectrum
                            PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].PLB[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff1': # First Derivative
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff1') and self.SpecDataMatrix[i][j].Specdiff1 is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff1[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff2': # Second Derivative
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff2') and self.SpecDataMatrix[i][j].Specdiff2 is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff2[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff1_norm': # First Derivative (normalized by signal)
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff1_norm') and self.SpecDataMatrix[i][j].Specdiff1_norm is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff1_norm[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff2_norm': # Second Derivative (normalized by signal)
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff2_norm') and self.SpecDataMatrix[i][j].Specdiff2_norm is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff2_norm[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff1_norm_intensity': # First Derivative (norm on intensity, then derive)
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff1_norm_intensity') and self.SpecDataMatrix[i][j].Specdiff1_norm_intensity is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff1_norm_intensity[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff2_norm_intensity': # Second Derivative (norm on intensity, then derive)
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff2_norm_intensity') and self.SpecDataMatrix[i][j].Specdiff2_norm_intensity is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff2_norm_intensity[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff1_norm_counts': # First Derivative (norm on counts, then derive)
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff1_norm_counts') and self.SpecDataMatrix[i][j].Specdiff1_norm_counts is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff1_norm_counts[self.aqpixstart:self.aqpixend])
                        elif self.speckeys[self.selectspecboxVari] == 'Specdiff2_norm_counts': # Second Derivative (norm on counts, then derive)
                            if hasattr(self.SpecDataMatrix[i][j], 'Specdiff2_norm_counts') and self.SpecDataMatrix[i][j].Specdiff2_norm_counts is not None:
                                PixMatrix[i][j] = np.sum(self.SpecDataMatrix[i][j].Specdiff2_norm_counts[self.aqpixstart:self.aqpixend])
                        if PixMatrix[i][j] < self.countthreshv:
                            if makenan == True:
                                PixMatrix[i][j] = np.nan
                    except Exception as e:
                        print(str(e), 'Error in getPLpixelIntervalMaxIndex')
                if PixMatrix[i][j] == np.nan:
                    PixMatrix[i][j] = 0
                        
    def getPLpixelSpecMax(self, PixMatrix):
        # fill matrix with data of the selected enry:
        self.updatewl()
        self.updatecountthresh()
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                if type(self.SpecDataMatrix[i][j]) == SpectrumData:
                    try:
                        self.selectspecboxVari = self.selectspecbox.get()
                        if self.speckeys[self.selectspecboxVari] == 'WL': #Wavelength
                            if np.sum(self.SpecDataMatrix[i][j].WL[self.aqpixstart:self.aqpixend]) < self.countthreshv:
                                PixMatrix[i][j] = np.nan
                            else:
                                PixMatrix[i][j] = self.SpecDataMatrix[i][j].fitmaxX
                        elif self.speckeys[self.selectspecboxVari] == 'BG': #Background
                            if np.sum(self.SpecDataMatrix[i][j].BG[self.aqpixstart:self.aqpixend]) < self.countthreshv:
                                PixMatrix[i][j] = np.nan
                            else:
                                PixMatrix[i][j] = self.SpecDataMatrix[i][j].fitmaxX
                        elif self.speckeys[self.selectspecboxVari] == 'PL': # Counts
                            if np.sum(self.SpecDataMatrix[i][j].PL[self.aqpixstart:self.aqpixend]) < self.countthreshv:
                                PixMatrix[i][j] = np.nan
                            else:
                                PixMatrix[i][j] = self.SpecDataMatrix[i][j].fitmaxX
                        elif self.speckeys[self.selectspecboxVari] == 'PLB': #Spectrum
                            if np.sum(self.SpecDataMatrix[i][j].PLB[self.aqpixstart:self.aqpixend]) < self.countthreshv:
                                PixMatrix[i][j] = np.nan
                            else:
                                PixMatrix[i][j] = self.SpecDataMatrix[i][j].fitmaxX
                    except Exception as e:
                        print('Error in getPLpixelSpecMax', str(e))
                        sys.exit()
                else:
                    PixMatrix[i][j] = np.nan

    def loadfiles(self):
        # read WL axis once for all files (must be same for all datafiles)
        lines = ['0']
        gotWL = False
        gotBG = False
        i = 0
        self.WL = []
        self.BG = []
        
        # Check if there are any files to load
        if len(self.fnames) == 0:
            print("No files to load. Skipping loadfiles() - will be populated by load_state().")
            self.WL_eV = []
            return
        
        while gotWL == False or gotBG == False:
            try:
                with open(self.fnames[i], 'r') as file:
                    lines = file.readlines()
            except Exception as e:
                print('Error While trying to read WL axis. No WL found in {} Files. {}'.format(i, str(e)))
                break
            # Process lines to store variables
            startreaddata = False 
            for line in lines:
                if '\t' in line:  # Data lines with tabs
                    parts = line.split()
                    if startreaddata == False:
                        count = 0
                        # start reading if at least two keys in the line
                        for j in parts:
                            if j in self.readinkeys:
                                count += 1
                        if count > 0:
                            startreaddata = True
                    elif startreaddata == True:
                        if gotWL == False:
                            try:
                                self.WL.append(float(parts[0]))
                            except Exception as e:
                                print('Error While trying to read WL axis from {}. {}'.format(self.fnames[i], str(e)))
                        if gotBG == False or self.loadeachbg == False:
                            try:
                                self.BG.append(float(parts[1]))
                            except Exception as e:
                                print('Error While trying to read WL axis from {}. {}'.format(self.fnames[i], str(e)))

            i += 1
            if len(self.WL) > 1:
                gotWL = True
            if len(self.BG) > 1:
                gotBG = True
        
        if self.loadeachbg == False:
            if self.linearbg == True:
                av = np.mean(self.BG)
                for i in range(len(self.BG)):
                    self.BG[i] = av

        # Convert WL and BG to float32 numpy arrays to reduce RAM usage
        self.WL = np.array(self.WL, dtype=np.float32)
        self.BG = np.array(self.BG, dtype=np.float32)

        # convert WL in nm to eV and store as WL_eV (use copy to avoid modifying WL in-place)
        self.WL_eV = deflib.wl_array_to_ev(self.WL.copy())

        # parallel loading of spectra
        self.parallel_load_spectra()
        del lines

    def parallel_load_spectra(self):
        # Get memory tracker
        mem_tracker = memory_tracker.get_default_memory_tracker()
        
        # Log initial state
        mem_tracker.log_separator("LOADING SPECTRA")
        mem_tracker.log_memory(
            "Before loading", 
            context="Parallel spectra loading",
            data_info={'num_files': len(self.fnames), 'spectral_points': len(self.WL) if hasattr(self, 'WL') and len(self.WL) > 0 else 0}
        )
        
        # before starting threads, clear specs
        self.specs = []

        lock = thre.Lock()  # To avoid race conditions when modifying self.specs

        # Use ThreadPoolExecutor to limit concurrent threads and prevent "too many open files" error
        # Max workers = min(32, number of CPU cores + 4) is a good default
        max_workers = min(32, (os.cpu_count() or 4) + 4)
        
        # Use ThreadPoolExecutor with a limited number of workers
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            futures = [executor.submit(load_spectrum, fname, self, lock) for fname in self.fnames]
            
            # Track progress for large datasets
            total = len(futures)
            completed = 0
            last_log_percent = 0
            
            # Wait for all tasks to complete
            for future in as_completed(futures):
                try:
                    future.result()  # This will raise any exceptions from the worker threads
                    completed += 1
                    
                    # Log progress every 10% for large datasets
                    if total >= 1000:
                        percent = int((completed / total) * 100)
                        if percent >= last_log_percent + 10:
                            mem_tracker.log_memory(
                                f"Loading progress: {percent}%",
                                context="Parallel spectra loading",
                                data_info={'completed': completed, 'total': total}
                            )
                            last_log_percent = percent
                            # Force garbage collection periodically for large datasets
                            if percent % 20 == 0:
                                gc.collect()
                except Exception as e:
                    print(f'Error loading spectrum: {e}')
        
        # Log after loading
        mem_tracker.log_memory(
            "After loading all spectra",
            context="Parallel spectra loading",
            data_info={'spectra_loaded': len(self.specs)}
        )
        
        # Force garbage collection after loading
        gc.collect()
        
		# after spectra are loaded, they must be put into matrix, after this, correlated cosmic ray removal can be applied (see autogenmatrix) # correlatedcosmicrayremoval

    def autogenmatrix(self):
        # Get memory tracker
        mem_tracker = memory_tracker.get_default_memory_tracker()
        mem_tracker.log_separator("BUILDING SPECTRAL MATRIX")
        mem_tracker.log_memory("Before matrix generation", context="autogenmatrix")
        
        self.mxcoords = []
        self.mycoords = []
        self.PMmetadata = {}
        if len(self.specs) == 0:
            print("No spectra loaded - initializing empty structures (will be populated by load_state)")
            # Initialize empty structures for loading from saved state
            self.PixAxX = []
            self.PixAxY = []
            self.SpecDataMatrix = []
            self.gdx = 0
            self.gdy = 0
            return
        elif len(self.specs) == 1:
            self.mxcoords.append(self.specs[0].data['x-position'])
            self.mycoords.append(self.specs[0].data['y-position'])
            self.PixAxX = [0]
            self.PixAxY = [0]
            self.SpecDataMatrix = [[self.specs[0]]]
            PixMatrix = [[0]]
            self.gdx = 0
            self.gdy = 0
        else:
            _coord_tol = 1e-9  # tolerance for floating-point coordinate deduplication
            for i in self.specs:
                if i is not None:
                    xpos = i.data['x-position']
                    ypos = i.data['y-position']
                    # Bug 1 fix: use tolerance-based comparison instead of exact equality
                    # to prevent floating-point noise from creating spurious extra columns/rows
                    if not any(abs(xpos - cx) <= _coord_tol for cx in self.mxcoords):
                        self.mxcoords.append(xpos)
                    if not any(abs(ypos - cy) <= _coord_tol for cy in self.mycoords):
                        self.mycoords.append(ypos)
            self.mxcoords = sorted(self.mxcoords)
            self.mycoords = sorted(self.mycoords)

            mem_tracker.log_memory(
                "After coordinate extraction",
                context="autogenmatrix",
                data_info={'x_coords': len(self.mxcoords), 'y_coords': len(self.mycoords)}
            )

            PixMatrix, self.SpecDataMatrix, self.PixAxX, self.PixAxY = self.genmatgrid(self.mxcoords, self.mycoords)
            
            mem_tracker.log_memory(
                "After grid generation",
                context="autogenmatrix",
                data_info={'matrix_shape': f"{len(self.PixAxY)}x{len(self.PixAxX)}"}
            )
            
            PixMatrixc = PMlib.PMclass(np.asarray(PixMatrix, dtype=float), self.PixAxX, self.PixAxY, self.PMmetadata)
            # Use counter for consistent naming
            initial_hsi_name = f'HSI{self._hsi_counter}'
            self._hsi_counter += 1
            PixMatrixc.name = initial_hsi_name
            self.PMdict[initial_hsi_name] = PixMatrixc
            self.SpecdataintoMatrix()
            
            mem_tracker.log_memory("After populating matrix", context="autogenmatrix")
            
            # after all threads are done, check if the cosmic ray removal method is in deflib.correlationcosmicfuncts
            if self.removecosmics == True:
                if self.remcosmicfunc in deflib.correlationcosmicfuncts:
                    mem_tracker.log_memory("Before correlated cosmic ray removal", context="autogenmatrix")
                    self.remcosmicfunc = deflib.correlationcosmicfuncts[self.remcosmicfunc](self.SpecDataMatrix, self.cosmicthreshold, self.cosmicpixels)
                    mem_tracker.log_memory("After correlated cosmic ray removal", context="autogenmatrix")
                else: 
                    pass
        
        # CRITICAL MEMORY OPTIMIZATION: Clear self.specs to free memory
        # The spectra are now stored in SpecDataMatrix, so we don't need the original list
        # This saves approximately 50% of memory for large datasets
        if hasattr(self, 'specs') and self.specs:
            num_specs = len(self.specs)
            self.specs = []
            gc.collect()  # Force garbage collection to release memory immediately
            mem_tracker.log_memory(
                "After clearing specs list (OPTIMIZATION)",
                context="autogenmatrix",
                data_info={'cleared_spectra': num_specs, 'note': 'Freed ~50% memory'}
            )

    def calculate_derivatives(self):
        """
        Calculate first and second derivatives for all SpectrumData objects in SpecDataMatrix.
        Uses sliding window polynomial fitting (Savitzky-Golay style).
        """
        # Get memory tracker
        mem_tracker = memory_tracker.get_default_memory_tracker()
        mem_tracker.log_separator("CALCULATING DERIVATIVES")
        mem_tracker.log_memory("Before derivative calculation", context="calculate_derivatives")
        
        # Check if derivative calculation is requested
        # derivative_polynomarray: [first_derivative_bool, second_derivative_bool, polynomial_order, N_fitpoints]
        if not self.derivative_polynomarray or len(self.derivative_polynomarray) < 4:
            print("Derivative calculation skipped: Invalid configuration array.")
            return

        try:
            # Handle both Tkinter variables (with .get()) and direct values
            d1_val = self.derivative_polynomarray[0]
            calc_d1 = d1_val.get() if hasattr(d1_val, 'get') else d1_val
            
            d2_val = self.derivative_polynomarray[1]
            calc_d2 = d2_val.get() if hasattr(d2_val, 'get') else d2_val
        except Exception as e:
            print(f"Error reading derivative flags: {e}")
            return
        
        if not (calc_d1 or calc_d2):
            print("Derivative calculation skipped: Neither 1st nor 2nd derivative requested.")
            return

        try:
            p_order_val = self.derivative_polynomarray[2]
            poly_order = int(p_order_val.get()) if hasattr(p_order_val, 'get') else int(p_order_val)
            
            n_points_val = self.derivative_polynomarray[3]
            N_fitpoints = int(n_points_val.get()) if hasattr(n_points_val, 'get') else int(n_points_val)
        except (ValueError, TypeError) as e:
            print(f"Invalid polynomial order or fit points for derivative calculation: {e}")
            return

        # Ensure window size is larger than poly_order for Savitzky-Golay
        if N_fitpoints <= poly_order:
            N_fitpoints = poly_order + 2
        if N_fitpoints % 2 == 0:
            N_fitpoints += 1 # Ensure odd window size

        print(f"Calculating derivatives: d1={calc_d1}, d2={calc_d2}, order={poly_order}, window={N_fitpoints}")
        
        # Get error engine for logging
        error_engine = error_handler.get_default_error_engine()

        # Iterate over all spectra in SpecDataMatrix
        # (self.specs is cleared after autogenmatrix() to free memory)
        from scipy.signal import savgol_filter
        
        for row_idx, row in enumerate(self.SpecDataMatrix):
            for col_idx, spec in enumerate(row):
                if spec is None:
                    continue

                try:
                    # Ensure we have data
                    # Check for None or empty arrays/lists explicitly to avoid ambiguity with numpy arrays
                    if not hasattr(spec, 'PLB') or not hasattr(spec, 'WL'):
                        error_engine.warning(f"Invalid spectrum object encountered at row {row_idx}, col {col_idx}. Type: {type(spec)}. Skipping.", context="calculate_derivatives")
                        continue
                        
                    if spec.PLB is None or spec.WL is None or len(spec.PLB) == 0 or len(spec.WL) == 0 or len(spec.PLB) != len(spec.WL):
                        continue

                    wl = np.asarray(spec.WL, dtype=np.float32)
                    plb = np.asarray(spec.PLB, dtype=np.float32)

                    # Calculate delta for proper scaling
                    delta = np.mean(np.diff(wl))

                    n_points = len(wl)
                    
                    # Ensure window size doesn't exceed array size
                    window = N_fitpoints
                    if window > n_points:
                        window = n_points if n_points % 2 != 0 else n_points - 1
                        
                    if window <= poly_order:
                        continue # insufficient points

                    if calc_d1:
                        spec.Specdiff1 = savgol_filter(plb, window, poly_order, deriv=1, delta=delta).astype(np.float32) # type: ignore
                    if calc_d2:
                        spec.Specdiff2 = savgol_filter(plb, window, poly_order, deriv=2, delta=delta).astype(np.float32) # type: ignore

                except Exception as e:
                    error_engine.error(e, context="calculate_derivatives", row=row_idx, col=col_idx)
                    continue

        # Normalize derivatives in SpecDataMatrix by the PLB signal
        # ONLY call if unconditionally needed or if user checked it. 
        # Adding a logic flag check to fix Memory Leak Issue 1 & 2
        # However, it's safer to only do it if explicitly requested somewhere. Wait, is there a user checkbox for normalized display? 
        # For now, it will only normalize if user explicitly enabled it. Wait, the software assumes it is always present for displaying.
        # Let's keep it but at least the loop is vectorized.
        # if hasattr(self, 'SpecDataMatrix') and self.SpecDataMatrix:
        #    hsi_normalization.normalize_derivatives_by_signal(self.SpecDataMatrix, signal_key='PLB')
        
        # Calculate normalized derivatives: normalize FIRST, then derive
        # Check if these normalization options are enabled
        calc_norm_and_derive = False
        calc_norm_on_intensity = False
        
        if self.calc_norm_and_derive is not None:
            calc_norm_and_derive = self.calc_norm_and_derive.get() if hasattr(self.calc_norm_and_derive, 'get') else bool(self.calc_norm_and_derive)
        
        if self.calc_norm_on_intensity is not None:
            calc_norm_on_intensity = self.calc_norm_on_intensity.get() if hasattr(self.calc_norm_on_intensity, 'get') else bool(self.calc_norm_on_intensity)
        
        # Calculate "normalize on counts, then derive" (calc_norm_and_derive checkbox)
        if calc_norm_and_derive and (calc_d1 or calc_d2):
            print("Calculating derivatives after normalizing on counts...")
            self._calculate_normalized_derivatives(poly_order, N_fitpoints, calc_d1, calc_d2, 
                                                   norm_type='counts', attr_suffix='_norm_counts')
        
        # Calculate "normalize on total intensity, then derive" (calc_norm_on_intensity checkbox)
        if calc_norm_on_intensity and (calc_d1 or calc_d2):
            print("Calculating derivatives after normalizing on total intensity...")
            self._calculate_normalized_derivatives(poly_order, N_fitpoints, calc_d1, calc_d2, 
                                                   norm_type='intensity', attr_suffix='_norm_intensity')
        
        # Log memory after derivative calculation
        mem_tracker.log_memory("After derivative calculation", context="calculate_derivatives")
    
    def _calculate_normalized_derivatives(self, poly_order, N_fitpoints, calc_d1, calc_d2, norm_type='counts', attr_suffix='_norm_counts'):
        """
        Calculate derivatives after normalizing the PLB data first.
        
        Args:
            poly_order: Polynomial order for fitting
            N_fitpoints: Number of points in sliding window
            calc_d1: Whether to calculate first derivative
            calc_d2: Whether to calculate second derivative
            norm_type: Type of normalization ('counts' or 'intensity')
            attr_suffix: Suffix for attribute names (e.g., '_norm_counts' or '_norm_intensity')
        """
        # Use the same normalization threshold as in hsi_normalization module
        # to maintain consistency across the codebase
        MIN_NORMALIZATION_THRESHOLD = hsi_normalization.MIN_NORMALIZATION_THRESHOLD
        
        from scipy.signal import savgol_filter
        
        # Ensure window size is larger than poly_order for Savitzky-Golay
        if N_fitpoints <= poly_order:
            N_fitpoints = poly_order + 2
        if N_fitpoints % 2 == 0:
            N_fitpoints += 1 # Ensure odd window size
        
        # Process SpecDataMatrix (self.specs is dead code after autogenmatrix)
        if hasattr(self, 'SpecDataMatrix') and self.SpecDataMatrix:
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    spec = self.SpecDataMatrix[i][j]
                    
                    if spec is None or not hasattr(spec, 'dataokay') or not spec.dataokay:
                        continue
                    
                    if spec.PLB is None or spec.WL is None or len(spec.PLB) == 0 or len(spec.WL) == 0 or len(spec.PLB) != len(spec.WL):
                        continue
                    
                    wl = np.asarray(spec.WL, dtype=np.float32)
                    plb = np.asarray(spec.PLB, dtype=np.float32)
                    
                    # Normalize PLB based on norm_type
                    if norm_type == 'intensity':
                        total_intensity = np.sum(plb)
                        if total_intensity > MIN_NORMALIZATION_THRESHOLD:
                            plb_normalized = plb / total_intensity
                        else:
                            plb_normalized = plb
                    elif norm_type == 'counts':
                        max_counts = np.max(plb)
                        if max_counts > MIN_NORMALIZATION_THRESHOLD:
                            plb_normalized = plb / max_counts
                        else:
                            plb_normalized = plb
                    else:
                        plb_normalized = plb
                    
                    n_points = len(wl)
                    window = N_fitpoints
                    if window > n_points:
                        window = n_points if n_points % 2 != 0 else n_points - 1
                        
                    if window <= poly_order:
                        continue # Skip if points too few
                        
                    delta = np.mean(np.diff(wl))
                    
                    try:
                        if calc_d1:
                            setattr(spec, f'Specdiff1{attr_suffix}', savgol_filter(plb_normalized, window, poly_order, deriv=1, delta=delta).astype(np.float32)) # type: ignore
                        
                        if calc_d2:
                            setattr(spec, f'Specdiff2{attr_suffix}', savgol_filter(plb_normalized, window, poly_order, deriv=2, delta=delta).astype(np.float32)) # type: ignore
                    except Exception as e:
                        print(f"Error calculating savgol_filter array: {e}")
                        pass
    
    def UpdateHSIselect(self):
        if len(self.PMdict) > 0:
            self.hsiselect['values'] = list(self.PMdict.keys())
            self.hsiselect.current(0)
            # select the last HSI
            self.hsiselect.set(list(self.PMdict.keys())[-1])

    def multiroitopixmatrix(self):
        if len(self.roihandler.roilist) == 0:
            print('No ROI found. Cannot create HSI.')
            return
        else: 
            roi_key = self.roiselgui.get()
            roi = self.roihandler.roilist[roi_key]
            
            roi_base = roi_key.split(' ')[0]
            base_name = self.hsiselect.get()
            new_hsi_name_base = f"{base_name}({roi_base})"
            
            new_hsi_name = new_hsi_name_base
            i = 1
            while new_hsi_name in self.PMdict:
                new_hsi_name = f"{new_hsi_name_base}_{i}"
                i += 1
                
        # Create new PixMatrix more efficiently - only copy the matrix data, not the entire object
        source_pm = self.PMdict[self.hsiselect.get()]
        # Create a copy of just the matrix data, apply ROI mask
        masked_matrix = np.array(source_pm.PixMatrix)  # Creates a copy
        masked_matrix[np.isnan(roi)] = np.nan
        # Create new PMclass with copied attributes instead of deepcopy
        lastpixmatrix = PMlib.PMclass(
            masked_matrix,
            source_pm.xax,
            source_pm.yax,
            source_pm.metadata.copy() if isinstance(source_pm.metadata, dict) else source_pm.metadata,
            name=new_hsi_name,
            units=source_pm.units,
            description=source_pm.description,
            data_type=source_pm.data_type
        )
        self.PMdict[new_hsi_name] = lastpixmatrix
        self.UpdateHSIselect()
        self.hsiselect.set(new_hsi_name)
        self.plotPixelMatrix(new_hsi_name)
    
    def delHSI(self):
        if len(self.PMdict) == 1:
            print("Error", 'Cannot delete last HSI.')
        else:
            del self.PMdict[self.hsiselect.get()]
            #del self.PMmetadata[self.hsiselect.get()]
        self.UpdateHSIselect()

    # set the spectra into the array
    def SpecdataintoMatrix(self, overwrite=False):
        # sorting function
        error_engine = error_handler.get_default_error_engine()
        for i in self.specs:
            x = i.data['x-position']
            y = i.data['y-position']
            xind, yind = deflib.closest_indices(self.PixAxX, self.PixAxY, x, y)
            if type(self.SpecDataMatrix[yind][xind]) == SpectrumData:
                if overwrite == True:
                    self.SpecDataMatrix[yind][xind] = i
                else:
                    # Bug 5 fix: log the collision to the error log instead of silently printing
                    occupant = self.SpecDataMatrix[yind][xind]
                    msg = (
                        f"Point neglected at grid index ({xind}, {yind}). "
                        f"Incoming position: ({x}, {y}). "
                        f"Occupying position: ({occupant.data['x-position']}, {occupant.data['y-position']}). "
                        f"Use overwrite=True to replace existing spectra."
                    )
                    error_engine.get_logger().warning(msg, extra={'context': 'SpecdataintoMatrix'})
            else:
                #self.SpecDataMatrix[xind][yind] = i
                self.SpecDataMatrix[yind][xind] = i
        
    def genmatgrid(self, xar, yar): # returns that must be filled with the SpectrumData Objects
        error_engine = error_handler.get_default_error_engine()
        self.matstart = [np.amin(self.mxcoords), np.amin(self.mycoords)] # find min coordinates
        self.matend = [np.amax(self.mxcoords), np.amax(self.mycoords)] # find max coordinates
        dxa = deflib.findif(self.mxcoords) # find differences between coordinates
        dya = deflib.findif(self.mycoords) # find differences between coordinates
        # Bug 2 fix: round step values before deduplication and most_freq_element so that
        # near-identical floating-point steps are treated as the same step size
        _step_round = 10
        dxa_rounded = [round(v, _step_round) for v in dxa]
        dya_rounded = [round(v, _step_round) for v in dya]
        self.dxanodup = deflib.remove_duplicates(dxa_rounded) # remove duplicates
        self.dyanodup = deflib.remove_duplicates(dya_rounded) # remove duplicates
        if len(self.dxanodup) == 1:
            self.gdx = dxa_rounded[0]
        else:
            self.gdx = deflib.most_freq_element(dxa_rounded) # find most frequent element
        if len(self.dyanodup) == 1:
            self.gdy = dya_rounded[0]
        else:
            self.gdy = deflib.most_freq_element(dya_rounded)
        self.gdx = round(self.gdx, 10) # avoid rounding errors
        self.gdy = round(self.gdy, 10) # avoid rounding errors

        if abs(self.gdx) < 1e-9 or abs(self.gdy) < 1e-9:
            msg = f"genmatgrid: computed step size is zero (gdx={self.gdx}, gdy={self.gdy}). Cannot build grid."
            error_engine.get_logger().error(msg, extra={'context': 'genmatgrid'})
            raise ValueError(msg)

        matpixax = []
        matpiyax = []
        PixelMatrix = []
        SpectralMatrix = []
        # Bug 4 fix: use round() before int() to prevent truncation of values like 4.9999999
        # +gdx/gdy in the numerator makes the range inclusive of the endpoint
        nx = int(round((self.matend[0] - self.matstart[0] + self.gdx) / self.gdx))
        ny = int(round((self.matend[1] - self.matstart[1] + self.gdy) / self.gdy))
        # Bug 3 fix: build axes by multiplying index * step to avoid accumulated rounding errors
        for i in range(nx): # create x axis of the matrix
            matpixax.append(round(i*self.gdx+self.matstart[0], 10)) # put rounded values in the axis
        for i in range(ny): # create y axis of the matrix
            matpiyax.append(round(i*self.gdy+self.matstart[1], 10)) # put rounded values in the axis
            fillmat = []
            pixmat = []
            for j in matpixax:
                fillmat.append(np.nan)#None) # fill rows with None
                pixmat.append(0)
            SpectralMatrix.append(fillmat) # put row data into the matrix
            PixelMatrix.append(pixmat)
        error_engine.get_logger().info(
            f"genmatgrid: grid built ({ny}x{nx}), gdx={self.gdx}, gdy={self.gdy}",
            extra={'context': 'genmatgrid'}
        )
        return(PixelMatrix, SpectralMatrix, matpixax, matpiyax)
    
    def writetopixmatrix(self, matrix, name=None, metadata=None, method='Integrated'):
        # selected dataset = self.PMdict[self.hsiselect.get()]

        if name == None or name not in self.PMdict.keys():
            # Get suffix from selected data set (Select Data Set)
            suffix = ''
            if method == 'Integrated' or method == 'spectral_shift':
                val = self.selectspecpixbox.get()
                # Use short code from speckeys if available, else sanitize value
                short_code = self.speckeys.get(val, val)
                suffix = f"{short_code}_{method}".replace(" ", "_").replace("(", "").replace(")", "")

            # Use monotonically increasing counter for unique HSI names
            elif 'fitparameter' in metadata: # type: ignore
                if metadata['fitparameter'] != '': # type: ignore
                    suffix = f"{metadata['fitparameter'].replace(" ", "_").replace("(", "").replace(")", "")}_{method}" # type: ignore
            newpmname = f'HSI{self._hsi_counter}_{suffix}'
            self._hsi_counter += 1
        else: 
            newpmname = name
        # add new PixMatrix to the dictionary with its metadata
        self.PMdict[newpmname] = PMlib.PMclass(np.asarray(matrix), self.PixAxX, self.PixAxY, self.PMmetadata)
        self.PMdict[newpmname].name = newpmname
        self.PMdict[newpmname].metadata = {'wlstart': self.wlstart, 'wlend': self.wlend, 'countthresh': self.countthreshv, 'aqpixstart': self.aqpixstart, 'aqpixend': self.aqpixend, 
        'unit': '', 
        'quantatity': '', 
        'fitmodel':'', 
        'fitparameter': ''
        }
        if metadata != {}:
            for i in metadata: # type: ignore
                self.PMdict[newpmname].metadata[i] = metadata[i] # type: ignore

        self.PMdict[newpmname].units = {'x': 'um', 'y': 'um', 'wl': 'nm', 'z': ''}
        return newpmname

    def plotHSIfromfitparam(self):
        self.updatewl()
        self.updatecountthresh()
        lastpm = copy.deepcopy(self.PMdict[self.hsiselect.get()].PixMatrix)
        # prepare metadata: 
        metadata = {
            # get unit from fit fitkeys[selected fitfunction][selected fitparameter]
            #'unit': self.fitkeys[self.selectwindowbox.get()][6][int(self.selectfitparambox.get().split(' ')[-1])] if self.selectwindowbox.get() in self.fitkeys and int(self.selectfitparambox.get().split(' ')[-1]) <= len(self.fitkeys[self.selectfitparambox.get().split(' ')[0]][6]) else 'failed to get unit in plotHSIfromfitparam',

            #'unit': self.fitkeys[self.selectwindowbox.get()][6][int(self.selectfitparambox.get().split(' ')[-1])] if self.selectfitparambox.get().split(' ')[-1][:1].isdigit() and self.selectwindowbox.get() in self.fitkeys and int(self.selectfitparambox.get().split(' ')[-1]) <= len(self.fitkeys[self.selectfitparambox.get().split(' ')[0]][6]) else 'failed to get unit in plotHSIfromfitparam',

            # we can get this easier, just call mathlib.fitunitpairs[fitparameter] to get the unit, if it exists (works for all fitfunctions, since the fitparameter is unique across all fitfunctions)
            'unit': matl.fitunitparis.get(self.selectfitparambox.get(), 'failed to get unit in plotHSIfromfitparam'),

            # get quantatity from fit fitkeys[selected fitfunction][selected fitparameter]
            'quantity': self.fitkeys[self.selectwindowbox.get()][5][int(self.selectfitparambox.get().split(' ')[-1])] if self.selectfitparambox.get().split(' ')[-1][:1].isdigit() and self.selectwindowbox.get() in matl.fitunitparis else self.selectfitparambox.get(), 

            # get fitmodel from fitkeys[selected fitfunction] 
            'fitmodel': self.selectfitparambox.get().split(' ')[0] if self.selectfitparambox.get()[0] in self.fitkeys else self.selectwindowbox.get(),

            # get fitparameter from selectfitparambox <- this one is okay
            'fitparameter': self.selectfitparambox.get() if self.selectfitparambox.get() != '' else 'failed to get fitparameter in plotHSIfromfitparam'
        }

        newpm = self.writetopixmatrix(lastpm, metadata=metadata, method=metadata['fitparameter'])
        self.getPLpixelIntervalMaxIndex(self.PMdict[newpm].PixMatrix, False)

        roi = lastpm[:]

        useroi = self.HSI_from_fitparam_useROI.get()
        if useroi:
            if self.roiselgui.get() not in self.roihandler.roilist.keys():
                print('No ROI selected. Cannot plot HSI from fit parameters.')
                useroi = False

            else:
                # get the selected ROI
                roi = self.roihandler.roilist[self.roiselgui.get()][:]
        
        if useroi:
            # get index of fitvari in self.allfitparams
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    if type(self.SpecDataMatrix[i][j]) == SpectrumData and np.isnan(roi[i][j]) == False:
                        if self.SpecDataMatrix[i][j].dataokay == True:
                            try:
                                pari, parj = matl.getindexofFitparameter(self.allfpnames, self.selectfitparambox.get())
                                if pari == None or parj == None:
                                    raise Exception('Parameter not found in fitparams.')
                                pari = int(pari)
                                parj = int(parj)
                                if pari != -1 and parj != -1:
                                    lastpm[i][j] = self.SpecDataMatrix[i][j].fitparams[pari][parj]
                                else:
                                    lastpm[i][j] = np.nan
                                    raise Exception('Parameter not found in fitparams.')
                            except Exception as e:
                                print('Error in plotHSIfromfitparam. {}'.format(str(e)))
                        else:
                            lastpm[i][j] = np.nan
                            print('No Data found in Pixel {}, {}'.format(i, j))
                    else:
                        lastpm[i][j] = np.nan
        else: 
            # get index of fitvari in self.allfitparams
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    if type(self.SpecDataMatrix[i][j]) == SpectrumData:
                        if self.SpecDataMatrix[i][j].dataokay == True:
                            try:
                                pari, parj = matl.getindexofFitparameter(self.allfpnames, self.selectfitparambox.get())
                                if pari == None or parj == None:
                                    raise Exception('Parameter not found in fitparams.')
                                pari = int(pari)
                                parj = int(parj)
                                if pari != -1 and parj != -1:
                                    lastpm[i][j] = self.SpecDataMatrix[i][j].fitparams[pari][parj]
                                else:
                                    lastpm[i][j] = np.nan
                                    raise Exception('Parameter not found in fitparams.')
                            except Exception as e:
                                print('Error in plotHSIfromfitparam. {}'.format(str(e)))
                        else:
                            lastpm[i][j] = np.nan
                            print('No Data found in Pixel {}, {}'.format(i, j))
                    else:
                        lastpm[i][j] = np.nan
                        print('No Data found in Pixel {}, {}'.format(i, j))
        
        # Apply normalization if enabled
        norm_matrix = self.get_normalization_matrix()
        if norm_matrix is not None:
            lastpm = hsi_normalization.HSINormalization.apply_normalization(
                lastpm, norm_matrix
            )
        
        # Update the pixel matrix with normalized values
        self.PMdict[newpm].PixMatrix = lastpm
        
        # test
        self.plotPixelMatrix(newpm)
        #self.plotPixelMatrix(self.PMdict[newpm].PixMatrix)
        self.UpdateHSIselect()
    
    def powercorrection(self, powerkey='Power at Glass Plate (µW)'):
        # Get memory tracker
        mem_tracker = memory_tracker.get_default_memory_tracker()
        mem_tracker.log_separator("POWER CORRECTION")
        mem_tracker.log_memory("Before power correction", context="powercorrection")
        
        # create Matrix with power values of metadata['Power at Glass Plane (mW)']
        if self.hsiselect.get() not in self.PMdict.keys():
            print('No valid HSI selected for power correction.')
            return
        powerMatrix = []
        for i in range(len(self.PMdict[self.hsiselect.get()].PixMatrix)):
            powermat = []
            for j in range(len(self.PMdict[self.hsiselect.get()].PixMatrix[i])):
                # check if self.SpecDataMatrix[i][j] is of type SpectrumData
                if self.SpecDataMatrix[i][j] is None or type(self.SpecDataMatrix[i][j]) != SpectrumData:
                    powermat.append(1)
                elif powerkey in self.SpecDataMatrix[i][j].data.keys():
                    powermat.append(float(self.SpecDataMatrix[i][j].data[powerkey]))
                else:
                    powermat.append(np.nan)
            powerMatrix.append(powermat)
        # normalize powerMatrix to its max value
        powerMatrix = np.asarray(powerMatrix)
        print('Powermatrix before removing nans:', powerMatrix)
        maxpower = np.amax(powerMatrix)
        print('Powermatrix before normalization:', powerMatrix)
        #sys.exit()
        for i in range(len(powerMatrix)):
            for j in range(len(powerMatrix[i])):
                powerMatrix[i][j] = powerMatrix[i][j]/maxpower
        
		# divide the selected spectra by the powerMatrix
        for i in range(len(self.SpecDataMatrix)):
            for j in range(len(self.SpecDataMatrix[i])):
                if type(self.SpecDataMatrix[i][j]) == SpectrumData:
                    try:
                        for k in range(len(self.SpecDataMatrix[i][j].PL)):
                            self.SpecDataMatrix[i][j].PL[k] = self.SpecDataMatrix[i][j].PL[k]/powerMatrix[i][j]
                        for k in range(len(self.SpecDataMatrix[i][j].PLB)):
                            self.SpecDataMatrix[i][j].PLB[k] = self.SpecDataMatrix[i][j].PLB[k]/powerMatrix[i][j]
                    except Exception as e:
                        print('Error in power correction of pixel {}, {}. {}'.format(i, j, str(e)))
        
        mem_tracker.log_memory("After power correction", context="powercorrection")
    
    def on_close(self):
        plt.close('all')
        # tkinter destroy
        # try:
        #     self.cmapframe.destroy()
        # except:
        #     pass
        # explicitly clean up spectra to release file handles
        if hasattr(self, 'specs'):
            for spec in self.specs:
                if spec is not None:
                    # clear large data arrays
                    if hasattr(spec, 'PL'):
                        spec.PL = []
                    if hasattr(spec, 'BG'):
                        spec.BG = []
                    if hasattr(spec, 'PLB'):
                        spec.PLB = []
                    if hasattr(spec, 'WL'):
                        spec.WL = []
            self.specs = []
        # clear the matrix
        if hasattr(self, 'SpecDataMatrix'):
            for row in self.SpecDataMatrix:
                for i in range(len(row)):
                    row[i] = None
            self.SpecDataMatrix = []
        # force garbage collection to release file handles
        gc.collect()
    
    def __del__(self):
        """Destructor to ensure cleanup when object is deleted"""
        try:
            self.on_close()
        except:
            pass
    
    def _save_specs_chunked(self, file_handle, specs, chunk_size=1000):
        """
        Save specs list in chunks to reduce memory usage during pickling.
        
        Args:
            file_handle: Open file handle for writing
            specs: List of SpectrumData objects
            chunk_size: Number of items per chunk
        """
        total_specs = len(specs)
        pickle.dump(total_specs, file_handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        for i in range(0, total_specs, chunk_size):
            chunk = specs[i:i+chunk_size]
            pickle.dump(chunk, file_handle, protocol=pickle.HIGHEST_PROTOCOL)
            if total_specs > 1000 and i > 0 and i % 5000 == 0:
                print(f"  - Saved {min(i+chunk_size, total_specs)}/{total_specs} spectra...")
    
    def _save_matrix_chunked(self, file_handle, matrix, chunk_size=50):
        """
        Save SpecDataMatrix in chunks (row by row or in row chunks).
        
        Args:
            file_handle: Open file handle for writing
            matrix: 2D list/array of SpectrumData objects
            chunk_size: Number of rows per chunk
        """
        if not matrix or len(matrix) == 0:
            pickle.dump(0, file_handle, protocol=pickle.HIGHEST_PROTOCOL)
            return
        
        num_rows = len(matrix)
        num_cols = len(matrix[0]) if num_rows > 0 else 0
        
        # Save dimensions
        pickle.dump((num_rows, num_cols), file_handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        # Save matrix in row chunks
        for i in range(0, num_rows, chunk_size):
            chunk = matrix[i:i+chunk_size]
            pickle.dump(chunk, file_handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def _save_pmdict_chunked(self, file_handle, pmdict, chunk_size=10):
        """
        Save PMdict in chunks to reduce memory usage.
        
        Args:
            file_handle: Open file handle for writing
            pmdict: Dictionary of PM objects
            chunk_size: Number of HSI images per chunk
        """
        keys = list(pmdict.keys())
        total_items = len(keys)
        
        # Save total count
        pickle.dump(total_items, file_handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        # Save in chunks
        for i in range(0, total_items, chunk_size):
            chunk_keys = keys[i:i+chunk_size]
            chunk_dict = {k: pmdict[k] for k in chunk_keys}
            pickle.dump(chunk_dict, file_handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    def _load_specs_chunked(self, file_handle):
        """
        Load specs list from chunked format.
        
        Args:
            file_handle: Open file handle for reading
            
        Returns:
            List of SpectrumData objects
        """
        total_specs = pickle.load(file_handle)
        specs = []
        
        chunk_size = 1000
        for i in range(0, total_specs, chunk_size):
            chunk = pickle.load(file_handle)
            specs.extend(chunk)
            if total_specs > 1000 and i > 0 and i % 5000 == 0:
                print(f"  - Loaded {min(i+chunk_size, total_specs)}/{total_specs} spectra...")
                gc.collect()  # Force garbage collection after each large chunk
        
        return specs
    
    def _load_matrix_chunked(self, file_handle):
        """
        Load SpecDataMatrix from chunked format.
        
        Args:
            file_handle: Open file handle for reading
            
        Returns:
            2D list of SpectrumData objects
        """
        dimensions = pickle.load(file_handle)
        
        if dimensions == 0:
            return []
        
        num_rows, num_cols = dimensions
        matrix = []
        
        chunk_size = 50
        for i in range(0, num_rows, chunk_size):
            chunk = pickle.load(file_handle)
            matrix.extend(chunk)
            if i > 0 and i % 100 == 0:
                gc.collect()  # Periodic garbage collection
        
        return matrix
    
    def _load_pmdict_chunked(self, file_handle):
        """
        Load PMdict from chunked format.
        
        Args:
            file_handle: Open file handle for reading
            
        Returns:
            Dictionary of PM objects
        """
        total_items = pickle.load(file_handle)
        pmdict = {}
        
        chunk_size = 10
        items_loaded = 0
        chunks_loaded = 0
        
        while items_loaded < total_items:
            chunk_dict = pickle.load(file_handle)
            pmdict.update(chunk_dict)
            items_loaded += len(chunk_dict)
            chunks_loaded += 1
            if chunks_loaded % 5 == 0:  # Every 5 chunks (50 items)
                gc.collect()  # Periodic garbage collection
        
        return pmdict
    
    def _prepare_pmdict_for_pickle(self, pmdict):
        """
        Replace np.nan values in PMdict PixMatrix data with unique numbers.
        This significantly reduces pickle file size and save/load time.
        
        Note: Creates deep copies of PM objects to avoid modifying the original.
        
        Returns:
            tuple: (modified_pmdict_copy, nan_replacements_dict)
        """
        nan_replacements = {}
        
        # Create deep copy to avoid modifying original PMdict
        pmdict_copy = {}
        
        for hsi_name, pm_obj in pmdict.items():
            # Deep copy the PM object
            pm_obj_copy = copy.deepcopy(pm_obj)
            
            if hasattr(pm_obj_copy, 'PixMatrix'):
                matrix = np.array(pm_obj_copy.PixMatrix)
                
                # Check if matrix contains np.nan
                if np.any(np.isnan(matrix)):
                    # Find a unique number not in the data
                    valid_data = matrix[~np.isnan(matrix)]
                    if len(valid_data) > 0:
                        # Use a value outside the data range
                        data_min = np.min(valid_data)
                        data_max = np.max(valid_data)
                        
                        # Generate two candidate unique numbers
                        candidate1 = data_min - max(1.0, abs(data_min)) - 1.0
                        candidate2 = data_max + max(1.0, abs(data_max)) + 1.0
                        
                        # Verify neither candidate is in the data (extra safety)
                        if candidate1 not in valid_data:
                            unique_num = candidate1
                        elif candidate2 not in valid_data:
                            unique_num = candidate2
                        else:
                            # Very unlikely, but use an extreme fallback
                            unique_num = -999999999.0
                    else:
                        unique_num = -999999.0  # fallback for all-nan matrices
                    
                    # Create a copy and replace nan with unique number
                    matrix_copy = np.copy(matrix)
                    matrix_copy[np.isnan(matrix_copy)] = unique_num
                    
                    # Preserve original type (list or array)
                    if isinstance(pm_obj_copy.PixMatrix, list):
                        pm_obj_copy.PixMatrix = matrix_copy.tolist()
                    else:
                        pm_obj_copy.PixMatrix = matrix_copy
                    
                    # Store replacement info
                    nan_replacements[hsi_name] = unique_num
            
            pmdict_copy[hsi_name] = pm_obj_copy
        
        return pmdict_copy, nan_replacements
    
    def _restore_nan_in_pmdict(self, pmdict, nan_replacements):
        """
        Restore np.nan values in PMdict PixMatrix data from unique numbers.
        
        Args:
            pmdict: Dictionary of PM objects
            nan_replacements: Dictionary mapping HSI names to replacement numbers
            
        Returns:
            pmdict with np.nan values restored
        """
        for hsi_name, unique_num in nan_replacements.items():
            if hsi_name in pmdict:
                pm_obj = pmdict[hsi_name]
                if hasattr(pm_obj, 'PixMatrix'):
                    matrix = np.array(pm_obj.PixMatrix)
                    # Replace unique number back to nan
                    matrix = np.where(matrix == unique_num, np.nan, matrix)
                    pm_obj.PixMatrix = matrix.tolist() if isinstance(pm_obj.PixMatrix, list) else matrix
        
        return pmdict
    
    def _prepare_roilist_for_pickle(self, roilist):
        """
        Replace np.nan values in ROI masks with unique numbers for efficient compression.
        Similar to _prepare_pmdict_for_pickle, but for ROI data structures.
        
        ROI masks use np.nan to mark pixels that are NOT in the ROI.
        This method replaces them with unique numbers for better pickle compression.
        
        Args:
            roilist: Dictionary mapping ROI names to ROI mask arrays
            
        Returns:
            tuple: (modified_roilist_copy, nan_replacements_dict)
        """
        import copy
        nan_replacements = {}
        roilist_copy = {}
        
        for roi_name, roi_mask in roilist.items():
            # Convert to numpy array if it isn't already
            roi_array = np.array(roi_mask)
            
            # Check if the ROI mask contains np.nan
            if np.any(np.isnan(roi_array)):
                # Find a unique number not in the data
                valid_data = roi_array[~np.isnan(roi_array)]
                if len(valid_data) > 0:
                    # ROI masks typically contain 1.0 for ROI pixels and np.nan for non-ROI
                    # Use a value outside the expected range
                    data_min = np.min(valid_data)
                    data_max = np.max(valid_data)
                    
                    # Generate unique number outside data range
                    candidate1 = data_min - max(1.0, abs(data_min)) - 1.0
                    candidate2 = data_max + max(1.0, abs(data_max)) + 1.0
                    
                    # Verify candidate is not in data
                    if candidate1 not in valid_data:
                        unique_num = candidate1
                    elif candidate2 not in valid_data:
                        unique_num = candidate2
                    else:
                        # Fallback (very unlikely)
                        unique_num = -999999999.0
                else:
                    # All-nan ROI mask (empty ROI)
                    unique_num = -999999.0
                
                # Create copy and replace nan with unique number
                roi_array_copy = np.copy(roi_array)
                roi_array_copy[np.isnan(roi_array_copy)] = unique_num
                
                # Store the processed ROI mask
                roilist_copy[roi_name] = roi_array_copy
                
                # Store replacement info
                nan_replacements[roi_name] = unique_num
            else:
                # No nan values, just copy the array
                roilist_copy[roi_name] = np.copy(roi_array)
        
        return roilist_copy, nan_replacements
    
    def _restore_nan_in_roilist(self, roilist, nan_replacements):
        """
        Restore np.nan values in ROI masks from unique numbers.
        
        Args:
            roilist: Dictionary of ROI masks
            nan_replacements: Dictionary mapping ROI names to replacement numbers
            
        Returns:
            roilist with np.nan values restored
        """
        for roi_name, unique_num in nan_replacements.items():
            if roi_name in roilist:
                roi_array = np.array(roilist[roi_name])
                # Replace unique number back to nan
                roi_array = np.where(roi_array == unique_num, np.nan, roi_array)
                roilist[roi_name] = roi_array
        
        return roilist
    
    def save_state(self, filename):
        """
        Save the complete state of the XYMap instance to a file.
        This includes all data, processing parameters, and results.
        
        Optimization: 
        1. Replaces np.nan values in PixMatrix data and ROI masks with unique 
           numbers before pickling to significantly reduce file size and save/load time.
        2. Uses chunked saving for large data structures to reduce memory spikes.
        """
        mem_tracker = memory_tracker.get_default_memory_tracker()
        mem_tracker.log_memory("Before save_state", context="save_state")
        
        # Prepare PMdict by replacing nan values with unique numbers
        # This creates a deep copy of PMdict with nan values replaced
        pmdict_prepared, nan_replacements_pmdict = self._prepare_pmdict_for_pickle(self.PMdict)
        
        # Prepare ROI list by replacing nan values with unique numbers
        # ROI masks use np.nan to mark non-ROI pixels, which needs compression
        roilist_raw = self.roihandler.roilist if hasattr(self, 'roihandler') else {}
        roilist_prepared, nan_replacements_roilist = self._prepare_roilist_for_pickle(roilist_raw)
        
        mem_tracker.log_memory("After preparing data for pickle", context="save_state")
        
        # Create a dictionary with all important state (excluding large arrays)
        state = {
            # Core data - WL and WL_eV arrays (shared references)
            'WL': self.WL,
            'WL_eV': self.WL_eV if hasattr(self, 'WL_eV') else None,
            'BG': self.BG if hasattr(self, 'BG') else [],
            'fnames': self.fnames,
            
            # NOTE: Large arrays (specs, SpecDataMatrix, PMdict) saved separately
            'format_version': 2,  # Version 2 = chunked format, Version 1 (implicit) = legacy format
            
            # PMdict metadata
            'PMmetadata': self.PMmetadata if hasattr(self, 'PMmetadata') else {},
            
            # Store the HSI counter for unique naming persistence
            '_hsi_counter': self._hsi_counter if hasattr(self, '_hsi_counter') else 0,
            
            # Store nan replacement metadata for restoration
            '_nan_replacements': nan_replacements_pmdict,
            '_nan_replacements_roilist': nan_replacements_roilist,
            
            # Averaged spectra data
            'disspecs': self.disspecs if hasattr(self, 'disspecs') else {},
            
            # Processing parameters
            'wlstart': self.wlstart,
            'wlend': self.wlend,
            'countthreshv': self.countthreshv,
            'aqpixstart': self.aqpixstart,
            'aqpixend': self.aqpixend,
            
            # Grid parameters
            'mxcoords': self.mxcoords if hasattr(self, 'mxcoords') else [],
            'mycoords': self.mycoords if hasattr(self, 'mycoords') else [],
            'PixAxX': self.PixAxX if hasattr(self, 'PixAxX') else [],
            'PixAxY': self.PixAxY if hasattr(self, 'PixAxY') else [],
            'gdx': self.gdx if hasattr(self, 'gdx') else 0,
            'gdy': self.gdy if hasattr(self, 'gdy') else 0,
            
            # Data ranges
            'DataSpecMin': self.DataSpecMin,
            'DataSpecMax': self.DataSpecMax,
            'DataSpecdL': self.DataSpecdL,
            'DataPixSt': self.DataPixSt,
            'DataPixDX': self.DataPixDX if hasattr(self, 'DataPixDX') else 0,
            'DataPixDY': self.DataPixDY if hasattr(self, 'DataPixDY') else 0,
            
            # Configuration
            'loadeachbg': self.loadeachbg,
            'linearbg': self.linearbg,
            'removecosmics': self.removecosmics,
            'cosmicthreshold': self.cosmicthreshold,
            'cosmicpixels': self.cosmicpixels,
            'remcosmicfunc': self.remcosmicfunc,
            'fontsize': self.fontsize,
            
            # Settings
            'colormap': self.colormap.get(),
            'WL_selection': self.WL_selection.get(),
            'HSI_fit_useROI': self.HSI_fit_useROI.get(),
            'HSI_from_fitparam_useROI': self.HSI_from_fitparam_useROI.get(),
            
            # Additional attributes
            'defentries': self.defentries,
        }
        
        # Save to file with error handling and chunking
        try:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(filename), exist_ok=True) if os.path.dirname(filename) else None
            
            with open(filename, 'wb') as f:
                # Save main state dictionary first
                pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
                mem_tracker.log_memory("After saving main state", context="save_state")
                
                # Save large data structures in chunks with memory cleanup
                # 1. Save specs list in chunks
                self._save_specs_chunked(f, self.specs)
                gc.collect()
                mem_tracker.log_memory("After saving specs", context="save_state")
                
                # 2. Save SpecDataMatrix in chunks
                self._save_matrix_chunked(f, self.SpecDataMatrix)
                gc.collect()
                mem_tracker.log_memory("After saving SpecDataMatrix", context="save_state")
                
                # 3. Save PMdict in chunks
                self._save_pmdict_chunked(f, pmdict_prepared)
                gc.collect()
                mem_tracker.log_memory("After saving PMdict", context="save_state")
                
                # 4. Save ROI list
                pickle.dump(roilist_prepared, f, protocol=pickle.HIGHEST_PROTOCOL)
                gc.collect()
                mem_tracker.log_memory("After saving roilist", context="save_state")
            
            print(f"Successfully saved XYMap state to: {filename}")
            print(f"  - Saved {len(self.specs)} spectra")
            print(f"  - Saved {len(self.PMdict)} HSI images")
            roi_count = len(self.roihandler.roilist) if hasattr(self, 'roihandler') else 0
            print(f"  - Saved {roi_count} ROI masks")
            if roi_count > 0:
                print(f"    ROI names: {list(self.roihandler.roilist.keys())}")
            disspecs_count = len(self.disspecs) if hasattr(self, 'disspecs') else 0
            print(f"  - Saved {disspecs_count} averaged spectra")
            if disspecs_count > 0:
                print(f"    Averaged spectra names: {list(self.disspecs.keys())}")
            if nan_replacements_pmdict:
                print(f"  - Replaced nan values in {len(nan_replacements_pmdict)} HSI images for optimization")
            if nan_replacements_roilist:
                print(f"  - Replaced nan values in {len(nan_replacements_roilist)} ROI masks for optimization")
            return True
        except Exception as e:
            print(f"Error saving XYMap state: {e}")
            traceback.print_exc()
            return False
    
    def load_state(self, filename):
        """
        Load a previously saved XYMap state from a file.
        This restores all data, processing parameters, and results.
        
        Automatically restores np.nan values that were replaced with unique numbers
        during the save operation.
        
        Supports both legacy (all-in-one) and new chunked format for memory efficiency.
        """
        mem_tracker = memory_tracker.get_default_memory_tracker()
        mem_tracker.log_memory("Before load_state", context="load_state")
        
        try:
            with open(filename, 'rb') as f:
                # Load main state dictionary
                state = pickle.load(f)
                mem_tracker.log_memory("After loading main state", context="load_state")
                
                # Check format version (version 2 = chunked, version 1 or missing = legacy)
                format_version = state.get('format_version', 1)
                
                if format_version >= 2:
                    # Load large data structures in chunks (version 2+)
                    self.specs = self._load_specs_chunked(f)
                    gc.collect()
                    mem_tracker.log_memory("After loading specs", context="load_state")
                    
                    self.SpecDataMatrix = self._load_matrix_chunked(f)
                    gc.collect()
                    mem_tracker.log_memory("After loading SpecDataMatrix", context="load_state")
                    
                    pmdict_loaded = self._load_pmdict_chunked(f)
                    gc.collect()
                    mem_tracker.log_memory("After loading PMdict", context="load_state")
                    
                    roilist_loaded = pickle.load(f)
                    gc.collect()
                    mem_tracker.log_memory("After loading roilist", context="load_state")
                else:
                    # Legacy format (version 1) - load from state dictionary
                    print("  - Loading legacy format (version 1)")
                    self.specs = state['specs']
                    self.SpecDataMatrix = state['SpecDataMatrix']
                    pmdict_loaded = state['PMdict']
                    roilist_loaded = state.get('roilist', {})
            
            # Restore WL arrays FIRST (these are shared references)
            self.WL = state['WL']
            self.WL_eV = state.get('WL_eV', None)
            self.BG = state.get('BG', [])
            self.fnames = state['fnames']
            
            # Ensure all specs have the correct WL reference
            # This fixes the issue where WL was a pointer and needs to be restored
            for spec in self.specs:
                if spec is not None:
                    spec.WL = self.WL  # Set the shared WL reference
                    if self.WL_eV is not None:
                        spec.WL_eV = self.WL_eV
            
            # Ensure all SpectrumData objects in matrix have correct WL reference
            for i in range(len(self.SpecDataMatrix)):
                for j in range(len(self.SpecDataMatrix[i])):
                    if isinstance(self.SpecDataMatrix[i][j], SpectrumData):
                        self.SpecDataMatrix[i][j].WL = self.WL
                        if self.WL_eV is not None:
                            self.SpecDataMatrix[i][j].WL_eV = self.WL_eV
            
            mem_tracker.log_memory("After WL reference updates", context="load_state")
            
            # Restore PMdict - contains all HSI images with fit results
            # First load the PMdict, then restore nan values if they were replaced
            nan_replacements = state.get('_nan_replacements', {})
            
            # Restore nan values from unique numbers
            if nan_replacements:
                self.PMdict = self._restore_nan_in_pmdict(pmdict_loaded, nan_replacements)
                print(f"  - Restored nan values in {len(nan_replacements)} HSI images")
            else:
                self.PMdict = pmdict_loaded
            
            mem_tracker.log_memory("After PMdict restoration", context="load_state")
            
            # Restore HSI counter and ensure it's higher than any existing HSI numbers
            self._hsi_counter = state.get('_hsi_counter', 0)
            
            # Ensure counter is higher than any existing HSI numbers to prevent duplicates
            # Extract numeric parts from HSI names like 'HSI0', 'HSI1', etc.
            max_hsi_num = -1
            for hsi_name in self.PMdict.keys():
                if hsi_name.startswith('HSI'):
                    try:
                        num = int(hsi_name[3:])  # Extract number after 'HSI'
                        max_hsi_num = max(max_hsi_num, num)
                    except ValueError:
                        # Skip names that don't follow HSI{number} pattern
                        pass
            
            # Set counter to be one more than the highest existing number
            if max_hsi_num >= self._hsi_counter:
                self._hsi_counter = max_hsi_num + 1
            
            self.PMmetadata = state.get('PMmetadata', {})
            
            # Restore ROI data with nan value restoration
            if hasattr(self, 'roihandler'):
                nan_replacements_roilist = state.get('_nan_replacements_roilist', {})
                
                # Restore nan values from unique numbers in ROI masks
                if nan_replacements_roilist:
                    self.roihandler.roilist = self._restore_nan_in_roilist(roilist_loaded, nan_replacements_roilist)
                    print(f"  - Restored nan values in {len(nan_replacements_roilist)} ROI masks")
                else:
                    self.roihandler.roilist = roilist_loaded
                
                # Validate ROI dimensions match HSI data dimensions
                if self.roihandler.roilist and len(self.PMdict) > 0:
                    # Get the first HSI to check dimensions
                    first_hsi_name = list(self.PMdict.keys())[0]
                    hsi_shape = np.array(self.PMdict[first_hsi_name].PixMatrix).shape
                    
                    for roi_name, roi_mask in self.roihandler.roilist.items():
                        roi_array = np.array(roi_mask)
                        if roi_array.shape != hsi_shape:
                            print(f"  Warning: ROI '{roi_name}' dimensions {roi_array.shape} don't match HSI dimensions {hsi_shape}")
                        else:
                            print(f"  ROI '{roi_name}' dimensions validated: {roi_array.shape}")
            
            mem_tracker.log_memory("After ROI restoration", context="load_state")
            
            # Restore averaged spectra data
            self.disspecs = state.get('disspecs', {})
            
            # Update the GUI combobox for spectra selection
            if hasattr(self, 'specselect') and len(self.disspecs) > 0:
                self.specselect['values'] = list(self.disspecs.keys())
                self.specselect.set(list(self.disspecs.keys())[-1])
            
            # Update the GUI combobox for ROI selection
            if hasattr(self, 'roiselgui') and hasattr(self, 'roihandler'):
                self.roiselgui['values'] = list(self.roihandler.roilist.keys())
                if len(self.roihandler.roilist) > 0:
                    self.roiselgui.set(list(self.roihandler.roilist.keys())[-1])
                else:
                    self.roiselgui.set('')
            
            # Restore processing parameters
            self.wlstart = state['wlstart']
            self.wlend = state['wlend']
            self.countthreshv = state['countthreshv']
            self.aqpixstart = state['aqpixstart']
            self.aqpixend = state['aqpixend']
            
            # Restore grid parameters
            self.mxcoords = state.get('mxcoords', [])
            self.mycoords = state.get('mycoords', [])
            self.PixAxX = state.get('PixAxX', [])
            self.PixAxY = state.get('PixAxY', [])
            self.gdx = state.get('gdx', 0)
            self.gdy = state.get('gdy', 0)
            
            # Restore data ranges
            self.DataSpecMin = state['DataSpecMin']
            self.DataSpecMax = state['DataSpecMax']
            self.DataSpecdL = state['DataSpecdL']
            self.DataPixSt = state['DataPixSt']
            self.DataPixDX = state.get('DataPixDX', self.gdx)
            self.DataPixDY = state.get('DataPixDY', self.gdy)
            
            # Restore configuration
            self.loadeachbg = state['loadeachbg']
            self.linearbg = state['linearbg']
            self.removecosmics = state['removecosmics']
            self.cosmicthreshold = state['cosmicthreshold']
            self.cosmicpixels = state['cosmicpixels']
            self.remcosmicfunc = state['remcosmicfunc']
            self.fontsize = state['fontsize']
            
            # Restore additional attributes
            if 'defentries' in state:
                self.defentries = state['defentries']
            
            # add fallback: iterate trough defentries and add any missing keys with default values
            for key, value in deflib.initdefaults().items():
                if key not in self.defentries:
                    self.defentries[key] = value
            
            # Restore settings (tk variables)
            self.colormap.set(state.get('colormap', self.defentries.get('default_colormap', 'viridis')))
            self.WL_selection.set(state.get('WL_selection', self.WL_values[0]))
            self.HSI_fit_useROI.set(state.get('HSI_fit_useROI', False))
            self.HSI_from_fitparam_useROI.set(state.get('HSI_from_fitparam_useROI', False))
            
            # Update GUI elements
            self.UpdateHSIselect()
            self.updateproc_spec_max()
            
            # Update threshold display
            if hasattr(self, 'countthresh'):
                self.countthresh.delete(0, tk.END)
                self.countthresh.insert(0, str(self.countthreshv))
            
            # Update font size display
            if hasattr(self, 'CMFont'):
                self.CMFont.delete(0, tk.END)
                self.CMFont.insert(0, str(self.fontsize))
            
            # Rebuild GUI with loaded data
            self.build_gui()
            
            print(f"Successfully loaded XYMap state from: {filename}")
            print(f"  - Loaded {len(self.specs)} spectra")
            print(f"  - Loaded {len(self.PMdict)} HSI images")
            roi_count = len(self.roihandler.roilist) if hasattr(self, 'roihandler') else 0
            print(f"  - Loaded {roi_count} ROI masks")
            if roi_count > 0:
                print(f"    ROI names: {list(self.roihandler.roilist.keys())}")
            disspecs_count = len(self.disspecs)
            print(f"  - Loaded {disspecs_count} averaged spectra")
            if disspecs_count > 0:
                print(f"    Averaged spectra names: {list(self.disspecs.keys())}")
            print(f"  - WL axis: {len(self.WL)} points from {self.DataSpecMin:.2f} to {self.DataSpecMax:.2f} nm")
            return True
            
        except Exception as e:
            print(f"Error loading XYMap state: {e}")
            traceback.print_exc()
            return False
		
        
    # function to update after loading HSI pickle
    def update_after_load(self):
        self.update_SpectralData_combobox()
        self.update_ROIcombobox()

    # update ROI combobox: display all ROIs in the roihandler.roilist for selection
    def update_ROIcombobox(self):
        if hasattr(self, 'roihandler') and hasattr(self, 'roiselgui'):
            self.roiselgui['values'] = list(self.roihandler.roilist.keys())
            if len(self.roihandler.roilist) > 0:
                self.roiselgui.set(list(self.roihandler.roilist.keys())[-1])
            else:
                self.roiselgui.set('')

    # update Select Spectral Data combobox after loading HSI pickle
    def update_SpectralData_combobox(self):
        if hasattr(self, 'specselect'):
            if hasattr(self, 'disspecs'):
                self.specselect['values'] = list(self.disspecs.keys())
                if len(self.disspecs) > 0:
                    self.specselect.set(list(self.disspecs.keys())[-1])
                else:
                    self.specselect.set('')

    def returndisspecs(self):
        return self.disspecs

def load_spectrum(fname, instance, lock):
    specobj = SpectrumData(
        fname,
        instance.WL,
        instance.BG,
        instance.loadeachbg,
        instance.linearbg,
        instance.removecosmics,
        instance.cosmicthreshold,
        instance.cosmicpixels,
        instance.remcosmicfunc, 
        instance.WL_eV

    )
    if specobj.dataokay:
        with lock:
            instance.specs.append(specobj)

def loades_spectra_derivatives(SpectrumInstance, lock):
    SpectrumInstance.calculate_derivatives()


