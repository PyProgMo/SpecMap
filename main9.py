import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from matplotlib.figure import Figure
#from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import os, sys, pickle, re
from PIL import Image, ImageTk
import lib9 as lib # type: ignore
import numpy as np
import deflib1 as deflib
import claralib1 as claralib
import export2 as xplib
import newtonspeclib1 as newtonlib
import threading as thr
import matplotlib.pyplot as plt
import HSI_debugger as DBG
import TCSPClib as tcspclib
import shutil, gc
import error_handler  # Centralized error handling and logging
import datetime as datet
import plotspecs
import cube2image
try:
    import marwin_specplot.marwin_specplot1 as marwin_specplotlib
    importmarwin = True
    print("Marwins specplotter imported successfully.")
except:
    importmarwin = False
    print("Marwins specplotter not imported, check if marwin_specplot.py is in the same directory as main.py")

class FileProcessorApp:
    def __init__(self, root, defaults):
        self.defaults = defaults
        self.root = root
        self.root.title("File Processor")
        self.multiple_BG = tk.IntVar()
        self.linearBG = tk.IntVar()
        self.removecosmicsBool = tk.IntVar()
        # Track threads with stop events for proper cleanup
        self.managed_threads = []
        self.derivative_polynomarray = [None, None, None, None]  # Placeholder for derivative settings
        self._is_closing = False
        self.createmenue()
        self.windownotebook(deflib.Notebooks)
        # init XYMap GUI components
        
        self.createbuttons(self.nodeframes['Load Data'])
        
        # Add scrollbars to Hyperspectra notebook tab
        self.hyper_canvas = tk.Canvas(self.nodeframes['Hyperspectra'])
        self.hyper_canvas.grid(row=0, column=0, sticky='nsew')
        
        # Configure grid weights for Hyperspectra frame
        self.nodeframes['Hyperspectra'].grid_rowconfigure(0, weight=1)
        self.nodeframes['Hyperspectra'].grid_columnconfigure(0, weight=1)
        
        # Add vertical scrollbar
        self.hyper_vsb = tk.Scrollbar(self.nodeframes['Hyperspectra'], orient="vertical", command=self.hyper_canvas.yview)
        self.hyper_vsb.grid(row=0, column=1, sticky='ns')
        
        # Add horizontal scrollbar
        self.hyper_hsb = tk.Scrollbar(self.nodeframes['Hyperspectra'], orient="horizontal", command=self.hyper_canvas.xview)
        self.hyper_hsb.grid(row=1, column=0, sticky='ew')
        
        # Configure canvas scrolling
        self.hyper_canvas.configure(yscrollcommand=self.hyper_vsb.set, xscrollcommand=self.hyper_hsb.set)
        
        # Create content frame inside canvas for Hyperspectra
        self.hyper_content_frame = tk.Frame(self.hyper_canvas)
        self.hyper_canvas_window = self.hyper_canvas.create_window((0, 0), window=self.hyper_content_frame, anchor='nw')
        
        # Bind configure event to update scroll region
        self.hyper_content_frame.bind('<Configure>', lambda e: self.hyper_canvas.configure(scrollregion=self.hyper_canvas.bbox('all')))
        
        # Bind canvas resize to adjust frame width (minimum 1200 pixels for Hyperspectra)
        self.hyper_canvas.bind('<Configure>', self._on_hyper_canvas_configure)
        
        # Enable mousewheel scrolling
        self.hyper_canvas.bind('<Enter>', self._bind_hyper_mousewheel)
        self.hyper_canvas.bind('<Leave>', self._unbind_hyper_mousewheel)
        
        # Create frames for Nanomap object inside hyper_content_frame
        self.cmapframe = tk.Frame(self.hyper_content_frame, width=100, height=50, borderwidth=5, relief="raised")
        self.cmapframe.pack(fill=tk.BOTH)
        self.specframe = tk.Frame(self.hyper_content_frame, borderwidth=5, relief="sunken")
        self.specframe.pack(fill=tk.BOTH, expand=True)
        
        # Create Nanomap object with GUI (no data loaded yet)
        self.Nanomap = lib.XYMap(
            [], self.cmapframe, self.specframe, 
            False, False, False, 
            20, 3, list(deflib.cosmicfuncts.keys())[0], 
            self.defaults, self.derivative_polynomarray,
            self.calc_norm_and_deriveBool, self.calc_norm_on_intensityBool
        )
        #self.Nanomap.build_gui()  # Build GUI components now
        
        # construct Exportframe to create a fancy Image from the HSI
        self.Exporter = xplib.Exportframe(self.nodeframes['HSI Plot'], self.Nanomap)

        # build the Plot Spectra notebook
        self.Specplotter = plotspecs.Specplottergui(
            Specdata={}, 
            Specdiffsets={}, 
            guiroot=self.nodeframes['Plot Spectra'],
            disspecs=self.Nanomap.disspecs
            )
        
        # build the cube2image notebook
        self.Cube2Imager = cube2image.Cube2Image(Nanomap=self.Nanomap, guiroot=self.nodeframes['Cube2Image'])

    def createmenue(self):
        # Create the menu bar
        menu_bar = tk.Menu(root)
        self.root.config(menu=menu_bar)
        deflib.create_menu(self.root, menu_bar)
    
    def windownotebook(self, notebookentries):
        # Create a notebook widget
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        # Create the frames for the notebook
        self.nodeframes = {}
        for i in range(len(notebookentries)):
            # add a new frame to the notebook and store it by its name in the frames dictionary
            frame = ttk.Frame(self.notebook)
            self.notebook.add(frame, text=notebookentries[i])
            self.nodeframes[notebookentries[i]] = frame
        
        if 'Marwins Spekplotter' in notebookentries and importmarwin:
            # check if self.marwinspecplot is already created, if not create it
            if hasattr(self, 'marwin_specplot'):
                self.marwin_specplot.clear()
            self.marwin_specplot = marwin_specplotlib.HSIExplorer(self.nodeframes['Marwins Spekplotter'])

    def createbuttons(self, Notebook):
        # Create a canvas with scrollbars for the Load Data tab
        self.load_canvas = tk.Canvas(Notebook)
        self.load_canvas.grid(row=0, column=0, sticky='nsew')
        
        # Configure grid weights to allow expansion
        Notebook.grid_rowconfigure(0, weight=1)
        Notebook.grid_columnconfigure(0, weight=1)
        
        # Add vertical scrollbar
        self.load_vsb = tk.Scrollbar(Notebook, orient="vertical", command=self.load_canvas.yview)
        self.load_vsb.grid(row=0, column=1, sticky='ns')
        
        # Add horizontal scrollbar
        self.load_hsb = tk.Scrollbar(Notebook, orient="horizontal", command=self.load_canvas.xview)
        self.load_hsb.grid(row=1, column=0, sticky='ew')
        
        # Configure canvas scrolling
        self.load_canvas.configure(yscrollcommand=self.load_vsb.set, xscrollcommand=self.load_hsb.set)
        
        # Create a frame inside the canvas to hold all content
        self.load_content_frame = tk.Frame(self.load_canvas)
        self.canvas_window = self.load_canvas.create_window((0, 0), window=self.load_content_frame, anchor='nw')
        
        # Bind configure event to update scroll region
        self.load_content_frame.bind('<Configure>', lambda e: self.load_canvas.configure(scrollregion=self.load_canvas.bbox('all')))
        
        # Bind canvas resize to adjust frame width
        self.load_canvas.bind('<Configure>', self._on_canvas_configure)
        
        # Enable mousewheel scrolling
        self.load_canvas.bind('<Enter>', self._bind_mousewheel)
        self.load_canvas.bind('<Leave>', self._unbind_mousewheel)
        
        # Specmap Load Frame (now inside load_content_frame)
        self.open_frame = tk.Frame(self.load_content_frame, width=100, height=100, borderwidth=5, relief="ridge")
        self.open_frame.grid(row=0, column=0, padx=5, pady=5)
        # Folder selection
        self.SpecMapLoad_label = tk.Label(self.open_frame, text="Select folder with spectra for Hyperspectral processing", width=100)
        self.SpecMapLoad_label.pack()
        self.folder_label = tk.Label(self.open_frame, text="Select Folder")
        self.folder_label.pack()
        
        self.folder_entry = tk.Entry(self.open_frame)
        self.folder_entry.pack(fill=tk.X)
        self.folder_entry.insert(0, defaults['data_file'])
        self.folder_button = tk.Button(self.open_frame, text="Browse", command=lambda: deflib.browse_folder(self.folder_entry))
        self.folder_button.pack()
        # Filename input
        self.filename_label = tk.Label(self.open_frame, text="Enter Filename")
        self.filename_label.pack()
        self.filename_entry = tk.Entry(self.open_frame)
        self.filename_entry.pack(fill=tk.X)
        self.filename_entry.insert(0, defaults['filename'])
        # File format
        self.fileformat_label = tk.Label(self.open_frame, text="Enter File format")
        self.fileformat_label.pack()
        self.fileformat_entry = tk.Entry(self.open_frame)
        self.fileformat_entry.pack(fill=tk.X)
        self.fileformat_entry.insert(0, defaults['file_extension'])
        # Process button
        self.process_button = tk.Button(self.open_frame, text="Load HSI data", command=self.init_spec_loadfiles)
        self.process_button.pack()
        # space between frames
        tk.Frame(self.open_frame, height=10).pack()

        # Load frame for background selection and cosmic removal place subframes inside
        self.loadframe = tk.Frame(self.open_frame)
        self.loadframe.pack()
        # BG selection
        # add extra frame for background selection onto loadframe
        self.bgframe = tk.Frame(self.loadframe, borderwidth=2, relief="sunken")
        self.bgframe.grid(row=0, column=0)
        self.multiple_BG.set(defaults['multiple_Background'])
        self.chkmultiple = tk.Checkbutton(self.bgframe, text="Multiple Backgrounds", variable=self.multiple_BG)
        self.chkmultiple.grid(row=0, column=0)
		# linear background
        self.linearBG.set(defaults['linear_Background'])
        self.linearBGcheck = tk.Checkbutton(self.bgframe, text="Linear Background", variable=self.linearBG)
        self.linearBGcheck.grid(row=1, column=0)
        self.linearBG.set(defaults['linear_Background'])
        # perform power correction 
        self.powercorrectionBool = tk.IntVar()
        self.powercorrectionBool.set(defaults['power_correction'])
        self.powercorrectioncheck = tk.Checkbutton(self.bgframe, text="Power Correction", variable=self.powercorrectionBool)
        self.powercorrectioncheck.grid(row=2, column=0)

        # below the frame add another frame for making multiple HSIs
        # spacings
        self.multiple_HSIs_inp_frame = tk.Frame(self.loadframe)
        self.multiple_HSIs_inp_frame.grid(row=4, column=0, columnspan=2)
        tk.Label(self.multiple_HSIs_inp_frame, text="Make Multiple HSIs from folders in selected dir").grid(row=0, column=0, columnspan=2)
        self.make_multiple_HSIsbool = tk.IntVar()
        self.make_multiple_HSIsbool.set(defaults['process_multiple_HSIs_bool'])
        self.process_multiple_HSIs = tk.Checkbutton(self.multiple_HSIs_inp_frame, text="Process Multiple HSIs", variable=self.make_multiple_HSIsbool).grid(row=1, column=0)
        # Add a checkbox: save HSI objects via pickle for later use
        self.save_multiple_HSIs_bool = tk.IntVar()
        self.save_multiple_HSIs_bool.set(defaults['save_multiple_HSIs_bool'])
        self.save_multiple_HSIs_check = tk.Checkbutton(self.multiple_HSIs_inp_frame, text="Save each HSI object", variable=self.save_multiple_HSIs_bool).grid(row=1, column=1)
        # main directory for multiple HSIs entry
        tk.Label(self.multiple_HSIs_inp_frame, text="File Main Directory:").grid(row=2, column=0)
        self.multiple_HSIs_dir_entry = tk.Entry(self.multiple_HSIs_inp_frame, width=80)
        self.multiple_HSIs_dir_entry.grid(row=2, column=1)
        self.multiple_HSIs_dir_entry.insert(0, defaults['hsifilesorter_maindir'])
        self.Browse_multiple_HSIs_dir_button = tk.Button(self.multiple_HSIs_inp_frame, text="Browse", command=lambda: deflib.browse_folder(self.multiple_HSIs_dir_entry))
        self.Browse_multiple_HSIs_dir_button.grid(row=2, column=2)
        # save directory to save the HSI output images
        tk.Label(self.multiple_HSIs_inp_frame, text="Save Images Directory:").grid(row=3, column=0)
        self.multiple_HSIs_save_dir_entry = tk.Entry(self.multiple_HSIs_inp_frame, width=80)
        self.multiple_HSIs_save_dir_entry.grid(row=3, column=1)
        self.multiple_HSIs_save_dir_entry.insert(0, defaults['hsifilesorter_savedir'])
        self.Browse_multiple_HSIs_save_dir_button = tk.Button(self.multiple_HSIs_inp_frame, text="Browse", command=lambda: deflib.browse_folder(self.multiple_HSIs_save_dir_entry))
        self.Browse_multiple_HSIs_save_dir_button.grid(row=3, column=2)
        # save directory for multiple HSI objects if save multiple HSIs checkbox is selected
        tk.Label(self.multiple_HSIs_inp_frame, text="Save HSI objects Directory:").grid(row=4, column=0)
        self.multiple_HSIs_save_pkl_dir_entry = tk.Entry(self.multiple_HSIs_inp_frame, width=80)
        self.multiple_HSIs_save_pkl_dir_entry.grid(row=4, column=1)
        self.Browse_hsiobject_save_dir_button = tk.Button(self.multiple_HSIs_inp_frame, text="Browse", command=lambda: deflib.browse_folder(self.multiple_HSIs_save_pkl_dir_entry)).grid(row=4, column=2)

        # Cosmic removal
        # add extra frame for cosmic removal onto loadframe
        self.cosmicframe = tk.Frame(self.loadframe, borderwidth=2, relief="sunken")
        self.cosmicframe.grid(row=0, column=1)

        self.removecosmicsBool.set(defaults['remove_cosmics'])
        self.removecosmics = tk.Checkbutton(self.cosmicframe, text="Remove Cosmics", variable=self.removecosmicsBool)
        self.removecosmics.grid(row=0, column=0)
        tk.Label(self.cosmicframe, text="Cosmic Removal function:").grid(row=1, column=0)
        self.cosmicremoval = ttk.Combobox(self.cosmicframe, values=list(deflib.cosmicfuncts.keys()), width=20)
        self.cosmicremoval.grid(row=1, column=1)
        self.cosmicremoval.set(list(deflib.cosmicfuncts.keys())[0])
        tk.Label(self.cosmicframe, text="Cosmic Threshold:").grid(row=2, column=0)
        self.cosmicthresholdentry = tk.Entry(self.cosmicframe, width=10)
        self.cosmicthresholdentry.grid(row=3, column=0)
        self.cosmicthresholdentry.insert(0, defaults['cosmic_threshold'])
        tk.Label(self.cosmicframe, text="Cosmic Width:").grid(row=2, column=1)
        self.cosmicwidthentry = tk.Entry(self.cosmicframe, width=10)
        self.cosmicwidthentry.grid(row=3, column=1)
        self.cosmicwidthentry.insert(0, defaults['cosmic_width'])
        tk.Label(self.cosmicframe, text="Laser Spotsize (nm):").grid(row=2, column=2)
        self.laserspotsizeentry = tk.Entry(self.cosmicframe, width=10)
        self.laserspotsizeentry.grid(row=3, column=2)
        self.laserspotsizeentry.insert(0, defaults['laser_spotsize_nm'])
        # input checkbox: calculate derivatives: first and second
        self.calculate_firstderivativeBool = tk.IntVar()
        self.calculate_firstderivativeBool.set(defaults['calculate_first_derivative'])
        self.calculate_firstderivative_check = tk.Checkbutton(self.cosmicframe, text="calc 1th derivative", variable=self.calculate_firstderivativeBool)
        self.calculate_firstderivative_check.grid(row=4, column=0)
        self.calculate_secondderivativeBool = tk.IntVar()
        self.calculate_secondderivativeBool.set(defaults['calculate_second_derivative'])
        self.calculate_secondderivative_check = tk.Checkbutton(self.cosmicframe, text="calc 2nd derivative", variable=self.calculate_secondderivativeBool)
        self.calculate_secondderivative_check.grid(row=4, column=1)
        # for derivative: select polynomial order
        tk.Label(self.cosmicframe, text="Polynom Order Derivative fit:").grid(row=5, column=0)
        defaultentries = [1, 2, 3, 4, 5]
        selection = defaults['derivative_polynomial_order'] if defaults['derivative_polynomial_order'] in defaultentries else 2
        self.derivative_polyorder_entry = ttk.Combobox(self.cosmicframe, values=[1, 2, 3, 4, 5], width=10) # type: ignore
        self.derivative_polyorder_entry.grid(row=5, column=1)
        self.derivative_polyorder_entry.insert(0, str(selection))
        self.derivative_fitpoints_label = tk.Label(self.cosmicframe, text="use n points:")
        self.derivative_fitpoints_label.grid(row=5, column=2)
        self.derivative_fitpoints_entry = tk.Entry(self.cosmicframe, width=10)
        self.derivative_fitpoints_entry.grid(row=5, column=3)
        if 'derivative_Nfitpoints' in defaults:
            self.derivative_fitpoints_entry.insert(0, str(defaults['derivative_Nfitpoints']))
        else:
            print("No default for derivative N fit points found, using 5")
            self.derivative_fitpoints_entry.insert(0, str(5))
        # construct polynomarray: [ first_derivative_bool, second_derivative_bool, polynomial_order ]
        self.derivative_polynomarray = [ self.calculate_firstderivativeBool, self.calculate_secondderivativeBool, self.derivative_polyorder_entry, self.derivative_fitpoints_entry]
        
        # Calculate first and second derivatives but normalize before derivative is calculated.
        self.calc_norm_and_deriveBool = tk.IntVar()
        self.calc_norm_and_deriveBool.set(defaults['calc_norm_and_derive'])
        self.calc_norm_and_derive_check = tk.Checkbutton(self.cosmicframe, text="norm counts, then derive", variable=self.calc_norm_and_deriveBool)
        self.calc_norm_and_derive_check.grid(row=4, column=2)

        # Normalize on total intensity
        self.calc_norm_on_intensityBool = tk.IntVar()
        self.calc_norm_on_intensityBool.set(defaults['calc_norm_on_intensity'])
        self.calc_norm_on_intensity_check = tk.Checkbutton(self.cosmicframe, text="norm on I, then derive", variable=self.calc_norm_on_intensityBool)
        self.calc_norm_on_intensity_check.grid(row=4, column=3)

        # if Marwins specplotter importet
        if importmarwin:
            self.marvs = marwin_specplotlib.HSIExplorer(self.nodeframes['Marwins Spekplotter'])
            # buld frame to select dir and open Marwins specplotter
            self.mwframe = tk.Frame(self.load_content_frame, width=60, height=100, borderwidth=5, relief="ridge")
            self.mwframe.grid(row=1, column=1, padx=5, pady=5)

            # Folder selection
            self.marwin_label = tk.Label(self.mwframe, text="Select folder with spectra for Marwins Specplotter", width=100)
            self.marwin_label.pack()
            self.marwin_label_folder_label = tk.Label(self.mwframe, text="Select Folder")
            self.marwin_label_folder_label.pack()
            self.marwin_folder_entry = tk.Entry(self.mwframe)
            self.marwin_folder_entry.pack(fill=tk.X)

            # add Browse and load buttons
            self.marwin_folder_button = tk.Button(self.mwframe, text="Browse", command=lambda: deflib.browse_folder(self.marwin_folder_entry))
            self.marwin_folder_button.pack(side=tk.LEFT, padx=5)
            self.marwin_load_button = tk.Button(self.mwframe, text="Load Marwins Plotter", command=lambda: self.marvs.load_dataset(self.marwin_folder_entry.get()))
            self.marwin_load_button.pack(side=tk.LEFT, padx=5)

        # Clara load frame (now inside load_content_frame)
        self.claraloadframe = tk.Frame(self.load_content_frame, width=60, height=100, borderwidth=5, relief="ridge")
        self.claraloadframe.grid(row=1, column=0, padx=5, pady=5)
        # Folder selection
        self.clara_label = tk.Label(self.claraloadframe, text="Select folder with spectra for Clara processing", width=100)
        self.clara_label.pack()
        self.clara_file_label = tk.Label(self.claraloadframe, text="Select File")
        self.clara_file_label.pack()

        # Clara process frame
        self.cl_file_entrystr = tk.StringVar()
        self.cl_file_entry = tk.Entry(self.claraloadframe, textvariable=self.cl_file_entrystr, width=113)
        self.cl_file_entry.pack()
        self.cl_file_entry.insert(0, defaults['clara_image'])
        self.cl_process_button = tk.Button(self.claraloadframe, text="Load Clara data", command=self.cl_loadfiles)
        self.cl_process_button.pack(side=tk.RIGHT, anchor='center')
        spacer = tk.Frame(self.claraloadframe, width=10)
        spacer.pack(side=tk.RIGHT, anchor='center')
        # add combobox to select the clara scaling factors (20x, 50x, 100x)
        self.cl_scaling_label = tk.Label(self.claraloadframe, text="Scaling Factor").pack(side=tk.LEFT, anchor='center')
        self.cl_scaling = ttk.Combobox(self.claraloadframe, values=list(claralib.cl_scalings.keys()), width=10)
        self.cl_scaling.pack(side=tk.LEFT, anchor='center')
        self.cl_scaling.set(list(claralib.cl_scalings.keys())[0])
        # Process button
        self.cl_folder_button = tk.Button(self.claraloadframe, text="Browse", command= lambda: deflib.select_file(self.cl_file_entrystr))
        self.cl_folder_button.pack(side=tk.RIGHT, anchor='center')

        # frame to save the current hsi object (now inside load_content_frame)
        self.saveframe = tk.Frame(self.load_content_frame, width=100, height=100, borderwidth=5, relief="ridge")
        self.saveframe.grid(row=2, column=0, padx=5, pady=5)
        
        # save the current hsi object
        self.save_label = tk.Label(self.saveframe, text="Save the current Data object")
        self.save_label.grid(row=0, column=0, columnspan=3, pady=(5,2))
        
        # Save path entry and browse button
        save_entry_frame = tk.Frame(self.saveframe)
        save_entry_frame.grid(row=1, column=0, columnspan=3, padx=5, pady=2)
        self.savehsipath = tk.StringVar()
        self.save_entry = tk.Entry(save_entry_frame, textvariable=self.savehsipath, width=100)
        self.save_entry.pack(side=tk.LEFT, padx=(0,5))
        self.save_entry.insert(0, defaults['save_hsi'])
        tk.Button(save_entry_frame, text="Browse", 
                 command=lambda: self.browse_save_path()).pack(side=tk.LEFT, padx=(0,5))
        self.save_button = tk.Button(save_entry_frame, text="Save", 
                                     command=lambda: self.saveNanomap(self.savehsipath.get()),
                                     width=5)
        self.save_button.pack(side=tk.LEFT)
        
        # Separator
        tk.Frame(self.saveframe, height=2, bg='gray').grid(row=2, column=0, columnspan=3, sticky='ew', padx=5, pady=10)
        
        # load saved spectra to current object
        self.load_label = tk.Label(self.saveframe, text="Load a saved Data object")
        self.load_label.grid(row=3, column=0, columnspan=3, pady=(5,2))
        
        # Load path entry and browse button
        load_entry_frame = tk.Frame(self.saveframe)
        load_entry_frame.grid(row=4, column=0, columnspan=3, padx=5, pady=2)
        self.loadhsipath = tk.StringVar()
        self.load_entry = tk.Entry(load_entry_frame, textvariable=self.loadhsipath, width=100)
        self.load_entry.pack(side=tk.LEFT, padx=(0,5))
        self.load_entry.insert(0, defaults['load_hsi_saved'])
        tk.Button(load_entry_frame, text="Browse", 
                 command=lambda: self.browse_load_path()).pack(side=tk.LEFT, padx=(0,5))
        self.load_button = tk.Button(load_entry_frame, text="Load", 
                                     command=lambda: self.loadhsisaved(self.loadhsipath.get()),
                                     width=5)
        self.load_button.pack(side=tk.LEFT)

        # frame to load Newton spectrum (now inside load_content_frame)
        self.newtonframe = tk.Frame(self.load_content_frame, width=100, height=100, borderwidth=5, relief="ridge")
        self.newtonframe.grid(row=3, column=0, padx=5, pady=5)
        # Newton spectrum load
        self.newton_label = tk.Label(self.newtonframe, text="Select Newton spectrum file", width=85)
        self.newton_label.grid(row=0, column=1)
        self.newton_file_label = tk.Label(self.newtonframe, text="Select File")
        self.newton_file_label.grid(row=1, column=0)
        self.newton_file_entrystr = tk.StringVar()
        self.newton_file_entry = tk.Entry(self.newtonframe, textvariable=self.newton_file_entrystr, width=97)
        self.newton_file_entry.grid(row=1, column=1)
        self.newton_file_entry.insert(0, defaults['newton_spectrum'])
        self.newton_process_button = tk.Button(self.newtonframe, text="Load Newton data", command=self.newtonloadfiles)
        self.newton_process_button.grid(row=2, column=1)
        self.newton_folder_button = tk.Button(self.newtonframe, text="Browse", command= lambda: deflib.select_file(self.newton_file_entrystr))
        self.newton_folder_button.grid(row=1, column=3)

        # frame to load TCSPC data (now inside load_content_frame)
        self.tcspcframe = tk.Frame(self.load_content_frame, width=100, height=100, borderwidth=5, relief="ridge")
        self.tcspcframe.grid(row=4, column=0, padx=5, pady=5)
        # TCSPC load
        self.tcspc_label = tk.Label(self.tcspcframe, text="Select TCSPC directory", width=86)
        self.tcspc_label.grid(row=0, column=1)
        self.tcspc_file_label = tk.Label(self.tcspcframe, text="TCSPC dir")
        self.tcspc_file_label.grid(row=1, column=0)
        self.tcspc_maindir_entrystr = tk.StringVar()
        self.tcspc_maindir_entry = tk.Entry(self.tcspcframe, textvariable=self.tcspc_maindir_entrystr, width=90)
        self.tcspc_maindir_entry.grid(row=1, column=1)
        self.tcspc_maindir_entry.insert(0, defaults['tcspc_maindir'])
        self.tcspc_subdir_label = tk.Label(self.tcspcframe, text="save dir")
        self.tcspc_subdir_label.grid(row=2, column=0)
        self.tcspc_subdir_entrystr = tk.StringVar()
        self.tcspc_subdir_entry = tk.Entry(self.tcspcframe, textvariable=self.tcspc_subdir_entrystr, width=90)
        self.tcspc_subdir_entry.grid(row=2, column=1)
        self.tcspc_process_button = tk.Button(self.tcspcframe, text="Load TCSPC data", command=self.tcspcloadfiles)
        self.tcspc_process_button.grid(row=3, column=1)
        self.tcspc_folder_button = tk.Button(self.tcspcframe, text="Browse", command= lambda: deflib.select_folder(self.tcspc_maindir_entrystr))
        self.tcspc_folder_button.grid(row=1, column=3)

        self.TCSPC_Processor = tcspclib.TCSPCprocessor(self.nodeframes['TCSPC'], self.tcspc_maindir_entry, self.tcspc_subdir_entry)
        self.TCSPC_Processor.build_frame()

        # load HSI on start ater constructing the GUI
        if defaults['loadonstart'] == True:
            self.spec_loadfiles()
        
        # build specfilesorter frame on the specfilesorter notebook
        self.specfilesorterframe = specfilesorter(
            self.nodeframes['HSI File Sorter'], 
            defaults['hsifilesorter_maindir'], 
            defaults['hsifilesorter_filename'], 
            defaults['hsifilesorter_fileend'], 
            defaults['hsifilesorter_savedir'], 
            defaults['hsifilesorter_processdir']
            )
        
    def _on_canvas_configure(self, event):
        """Update canvas window width when canvas is resized"""
        # Ensure minimum width of 500 pixels for proper horizontal scrolling
        canvas_width = max(event.width, 500)
        self.load_canvas.itemconfig(self.canvas_window, width=canvas_width)
    
    def _on_hyper_canvas_configure(self, event):
        """Update canvas window width for Hyperspectra canvas when resized"""
        # Ensure minimum width of 1200 pixels for Hyperspectra frame
        canvas_width = max(event.width, 1200)
        self.hyper_canvas.itemconfig(self.hyper_canvas_window, width=canvas_width)
    
    def _on_mousewheel(self, event):
        """Handle mousewheel scrolling"""
        if event.delta:
            self.load_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        else:
            if event.num == 4:
                self.load_canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                self.load_canvas.yview_scroll(1, "units")
    
    def _bind_mousewheel(self, event):
        """Bind mousewheel events when mouse enters canvas"""
        self.load_canvas.bind_all("<MouseWheel>", self._on_mousewheel)
        self.load_canvas.bind_all("<Button-4>", self._on_mousewheel)
        self.load_canvas.bind_all("<Button-5>", self._on_mousewheel)
    
    def _unbind_mousewheel(self, event):
        """Unbind mousewheel events when mouse leaves canvas"""
        self.load_canvas.unbind_all("<MouseWheel>")
        self.load_canvas.unbind_all("<Button-4>")
        self.load_canvas.unbind_all("<Button-5>")
    
    def _on_hyper_mousewheel(self, event):
        """Handle mousewheel scrolling for Hyperspectra canvas"""
        if event.delta:
            self.hyper_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        else:
            if event.num == 4:
                self.hyper_canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                self.hyper_canvas.yview_scroll(1, "units")
    
    def _bind_hyper_mousewheel(self, event):
        """Bind mousewheel events when mouse enters Hyperspectra canvas"""
        self.hyper_canvas.bind_all("<MouseWheel>", self._on_hyper_mousewheel)
        self.hyper_canvas.bind_all("<Button-4>", self._on_hyper_mousewheel)
        self.hyper_canvas.bind_all("<Button-5>", self._on_hyper_mousewheel)
    
    def _unbind_hyper_mousewheel(self, event):
        """Unbind mousewheel events when mouse leaves Hyperspectra canvas"""
        self.hyper_canvas.unbind_all("<MouseWheel>")
        self.hyper_canvas.unbind_all("<Button-4>")
        self.hyper_canvas.unbind_all("<Button-5>")
    
    def browse_save_path(self):
        """Open file dialog to select save path"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".pkl",
            filetypes=[("Pickle files", "*.pkl"), ("All files", "*.*")],
            title="Select save location for Data object"
        )
        if filename:
            self.savehsipath.set(filename)
    
    def browse_load_path(self):
        """Open file dialog to select load path"""
        filename = filedialog.askopenfilename(
            defaultextension=".pkl",
            filetypes=[("Pickle files", "*.pkl"), ("All files", "*.*")],
            title="Select Data object to load"
        )
        if filename:
            self.loadhsipath.set(filename)
    
    def tcspcloadfiles(self):
        file = self.tcspc_maindir_entry.get()
        # check, if file is valid
        if not file:
            print('Please select a TCSPC main directory')
            return
        else: 
            # check if file exists
            if not os.path.exists(file):
                print('TCSPC main directory does not exist:', file)
                return
        
        self.TCSPC_Processor.filepath = file
        self.TCSPC_Processor.load_tcspc()

    def newtonloadfiles(self):
        file = self.newton_file_entry.get()        
        self.newtonclass = newtonlib.newtonspecopener(self.nodeframes['Newton Spectrum'], file)

    def init_spec_loadfiles(self):
        if self.make_multiple_HSIsbool.get() == True:
            # try to get File main directory
            filemaindir = self.multiple_HSIs_dir_entry.get()
            savedir = self.multiple_HSIs_save_dir_entry.get()
            pklsavedir = self.multiple_HSIs_save_pkl_dir_entry.get()

            if not filemaindir:
                print("Error while loading HSI data, please select a main directory for multiple HSIs")
                self.spec_loadfiles()
            if not savedir:
                print("Error while loading HSI data, please select a save directory for multiple HSIs")
                self.spec_loadfiles()
            else:
                foldersinmaindir = [f for f in os.listdir(filemaindir) if os.path.isdir(os.path.join(filemaindir, f))]
                for folder in foldersinmaindir:
                    # debug: print all open files
                    fullfolderpath =  os.path.join(filemaindir, folder)
                    self.folder_entry.delete(0, tk.END)
                    self.folder_entry.insert(0, fullfolderpath)
                    filename = str(folder+"_HSI.png")
                    imagesavefolder = os.path.join(savedir, filename)
                    try:
                        # before loading: delete self.Nanomap if it exists
                        if hasattr(self, 'Nanomap'):
                            self.Nanomap.on_close()
                            del self.Nanomap
                            # force garbage collection to release file handles immediately
                            gc.collect()
                        self.spec_loadfiles()
                        # create intensity colormap and save spectra
                        self.Nanomap.buildandPlotIntCmap(savetoimage=imagesavefolder)

                        # save the HSI object if the checkbox is selected
                        if self.save_multiple_HSIs_bool.get() == True:
                            hsi_object_save_path = os.path.join(pklsavedir, folder+"_HSI_object.pkl")
                            self.saveNanomap(hsi_object_save_path)
                    except:
                        print("Error processing folder:", fullfolderpath)
                        continue
                
        else:
            self.spec_loadfiles()
        
    # testcomment
    def spec_loadfiles(self):
        # Close pyplot-managed figures individually so a destroyed Tk toolbar
        # cannot abort the load path with a TclError.
        for fig_num in list(plt.get_fignums()):
            try:
                plt.close(fig_num)
            except tk.TclError:
                continue
            except Exception:
                continue
        
        # Stop managed threads safely
        self._stop_managed_threads()
        
        # Fallback: try to stop any remaining threads with stop_event
        for thread in thr.enumerate():
            if thread.name != 'MainThread' and thread.is_alive():
                # Use getattr with default instead of hasattr to avoid exception suppression
                stop_event = getattr(thread, 'stop_event', None)
                if stop_event is not None:
                    try:
                        stop_event.set()
                    except Exception as e:
                        print(f"Warning: Could not stop thread {thread.name}: {e}")
		
        # kill any existing Nanomap object
        try:
            if hasattr(self, 'Nanomap'):
                nanomap = getattr(self, 'Nanomap', None)
                if nanomap is not None:
                    try:
                        nanomap.on_close()
                    except Exception as e:
                        print(f"Error calling on_close on Nanomap: {e}")
                del self.Nanomap
        except Exception as e:
            print(f"Error closing/deleting Nanomap: {e}")

        # Destroy cmapframe if it exists
        if hasattr(self, 'cmapframe'):
            try:
                if self.cmapframe.winfo_exists():
                    self.cmapframe.destroy()
            except tk.TclError:
                pass
            except Exception as e:
                print(f"Error destroying cmapframe: {e}")

        # Destroy specframe if it exists
        if hasattr(self, 'specframe'):
            try:
                if self.specframe.winfo_exists():
                    self.specframe.destroy()
            except tk.TclError:
                pass
            except Exception as e:
                print(f"Error destroying specframe: {e}")

        try:
            if hasattr(self, 'Exporter'):
                del self.Exporter
        except:
             pass
        
        # force garbage collection to release file handles immediately
        try:
            gc.collect()
        except:
            pass
        
        # Re-initialize defaults and read input fields
        self.defaults = deflib.initdefaults()
        folder = self.folder_entry.get()
        filename = self.filename_entry.get()
        fileend = self.fileformat_entry.get()
        
        if not folder:
            print("Error while loading HSI data, please select a folder")
            return
        if not filename:
            print("Error while loading HSI data, please select a file")
            return
            
        files_processed = []
        for dirpath, _, filenames in os.walk(folder):
            for i in self.filter_substring(filenames, filename):
                if fileend in i:
                    file_path = os.path.join(dirpath, i)
                    files_processed.append(file_path) # can be accessed to process the files

        # frames for the colormap and spectral buttons
        if files_processed:
            # Double check frames are cleared
            try:
                if hasattr(self, 'cmapframe') and self.cmapframe.winfo_exists():
                    self.cmapframe.destroy()
                if hasattr(self, 'specframe') and self.specframe.winfo_exists():
                    self.specframe.destroy()
            except:
                pass
                
            # Recreate frames inside hyper_content_frame
            self.cmapframe = tk.Frame(self.hyper_content_frame, width=100, height=50, borderwidth=5, relief="raised")
            self.cmapframe.pack(fill=tk.BOTH)
            self.specframe = tk.Frame(self.hyper_content_frame, borderwidth=5, relief="sunken")
            self.specframe.pack(fill=tk.BOTH, expand=True)
            
            # get the cosmic threshold and width
            try:
                self.cosmicthreshold = int(self.cosmicthresholdentry.get())
            except:
                self.cosmicthreshold = 20
                self.cosmicthresholdentry.insert(0, '20')
            try:
                self.cosmicwidth = int(self.cosmicwidthentry.get())
            except:
                self.cosmicwidth = 3
                self.cosmicwidthentry.insert(0, '3')
            self.Nanomap = lib.XYMap(
                files_processed, self.cmapframe, self.specframe, 
                bool(self.multiple_BG.get()), bool(self.linearBG.get()), bool(self.removecosmicsBool.get()), 
                self.cosmicthreshold, self.cosmicwidth, self.cosmicremoval.get(), 
                self.defaults, self.derivative_polynomarray,
                self.calc_norm_and_deriveBool, self.calc_norm_on_intensityBool
                )
            if self.powercorrectionBool.get() == 1:
                self.Nanomap.powercorrection()
            self.Exporter = xplib.Exportframe(self.nodeframes['HSI Plot'], self.Nanomap)
            self.Specplotter = plotspecs.Specplottergui(
                Specdata={}, 
                Specdiffsets={}, 
                guiroot=self.nodeframes['Plot Spectra'],
                disspecs=self.Nanomap.disspecs
                )

            # self.Cube2Imager: update 
            if hasattr(self, 'Cube2Imager'):
                # update CubeImager's reference to the new Nanomap instance
                if hasattr(self.Cube2Imager, 'gui') and self.Cube2Imager.gui is not None:
                    self.Cube2Imager.gui.Nanomap = self.Nanomap
                self.Cube2Imager.update_bounds()

        else:
            print("No files found with the specified name.")

    def filter_substring(self, a, b):
        return [element for element in a if b in element]
    
    def _stop_managed_threads(self):
        """
        Safely stop all managed threads by setting their stop_event.
        Removes dead threads from the list.
        """
        threads_to_remove = []
        
        for thread in self.managed_threads:
            if thread.is_alive():
                stop_event = getattr(thread, 'stop_event', None)
                if stop_event is not None:
                    try:
                        stop_event.set()
                        # Optionally wait for thread to finish (with timeout)
                        thread.join(timeout=0.5)
                    except Exception as e:
                        print(f"Warning: Error stopping thread {thread.name}: {e}")
            else:
                threads_to_remove.append(thread)
        
        # Clean up dead threads from the list
        for thread in threads_to_remove:
            self.managed_threads.remove(thread)
    
    def cl_loadfiles(self):
        file = self.cl_file_entry.get()
        if not file:
            print("Error", "Please select a file")
            return
        try:
            # clara 100x scaling factor 56.8 nm /pixel
            # adjust dx and dy for different scaling factors
            dx = 0.0568*2#0.0568 
            dy = 0.0568*2#0.0568
            scaling = self.cl_scaling.get()
            if scaling in claralib.cl_scalings.keys():
                dx = claralib.cl_scalings[scaling]
                dy = claralib.cl_scalings[scaling]
            self.claraimage = claralib.imageprocessor(self.nodeframes['Clara Image'], file, deflib.loadclaraimage, None, dx, dy)# 0.568, 0.568) 100x scaling
        except Exception as error:
            print("Error", "Could not load Clara image. {}".format(error))
    
    def on_closing(self):
        if self._is_closing:
            return
        self._is_closing = True

        # Tear down Cube2Image first so embedded matplotlib/tk canvas does not outlive root.
        if hasattr(self, 'Cube2Imager') and self.Cube2Imager is not None:
            try:
                self.Cube2Imager.destroy()
            except Exception as e:
                print(f"Cube2Image cleanup failed: {e}")

        # Close Nanomap resources if available.
        if hasattr(self, 'Nanomap') and self.Nanomap is not None:
            try:
                if hasattr(self.Nanomap, 'on_close'):
                    self.Nanomap.on_close()
            except Exception as e:
                print(f"Nanomap cleanup failed: {e}")

        # Cleanup tracked worker threads/events.
        try:
            self._stop_managed_threads()
        except Exception as e:
            print(f"Thread cleanup failed: {e}")

    def saveNanomap(self, filename):
        """
        Save the complete XYMap state to a file.
        This saves all data, processing parameters, HSI images, and ROI masks.
        """
        # check if the filename is empty
        if not filename:
            print("Error: Please enter a filename")
            return
        
        # check if the Nanomap object exists
        if not hasattr(self, 'Nanomap'):
            print("Error: No data loaded to save")
            return
        
        # Ensure .pkl extension
        if not filename.endswith('.pkl'):
            filename += '.pkl'
        
        # Handle file already exists
        if os.path.exists(filename):
            filename = deflib.increment_filename(filename)
            print(f"Saving to: {filename}")
        
        # Call the XYMap save_state method
        success = self.Nanomap.save_state(filename)
        
        if success:
            print("Success", f"Data saved successfully to:\n{filename}")
            # Update the entry field with the actual saved path
            self.savehsipath.set(filename)
        else:
            print("Error", "Failed to save data. Check console for details.")
    
    def loadhsisaved(self, filename):
        """
        Load a previously saved XYMap state from a file.
        This restores all data, processing parameters, HSI images, and ROI masks.
        """
        # check if the filename is empty
        if not filename:
            print('Error: Please enter a filename to load')
            return
        
        # check if filename is a valid path
        if not os.path.exists(filename):
            print(f'Error: File not found: {filename}')
            return
        
        # Check if the Nanomap object exists; if not, we need to create one
        if not hasattr(self, 'Nanomap'):
            print("Creating new XYMap instance to load data into...")
            # Create minimal frames for the XYMap if they don't exist
            try:
                self.cmapframe = tk.Frame(self.hyper_content_frame, width=100, height=50, borderwidth=5, relief="raised")
                self.cmapframe.pack(fill=tk.BOTH)
                self.specframe = tk.Frame(self.hyper_content_frame, borderwidth=5, relief="sunken")
                self.specframe.pack(fill=tk.BOTH, expand=True)
            except:
                pass
            
            # Create a temporary XYMap instance with minimal parameters
            # Skip GUI building - it will be built after load_state() populates the data
            self.Nanomap = lib.XYMap(
                [], self.cmapframe, self.specframe, 
                False, False, False, 
                20, 3, list(deflib.cosmicfuncts.keys())[0], 
                self.defaults, self.derivative_polynomarray,
                self.calc_norm_and_deriveBool, self.calc_norm_on_intensityBool,
                skip_gui_build=True  # Don't build GUI yet - will build after loading data
            )
            
            # Create Exporter
            self.Exporter = xplib.Exportframe(self.nodeframes['HSI Plot'], self.Nanomap)
            # build the Plot Spectra notebook
            self.Specplotter = plotspecs.Specplottergui(
                Specdata={}, 
                Specdiffsets={}, 
                guiroot=self.nodeframes['Plot Spectra'],
                disspecs=self.Nanomap.disspecs
                )
        
            # self.Cube2Imager: updte 
            if hasattr(self, 'Cube2Imager'):
                # update CubeImager's reference to the new Nanomap instance
                if hasattr(self.Cube2Imager, 'gui') and self.Cube2Imager.gui is not None:
                    self.Cube2Imager.gui.Nanomap = self.Nanomap
                self.Cube2Imager.update_bounds()

        # Load the state
        success = self.Nanomap.load_state(filename)

        if success:
            print("Success", f"Data loaded successfully from:\n{filename}")
            # Update the entry field
            self.loadhsipath.set(filename)
            self.Nanomap.update_after_load()
            
            # Create Exporter
            self.Exporter = xplib.Exportframe(self.nodeframes['HSI Plot'], self.Nanomap)
            # build the Plot Spectra notebook
            self.Specplotter = plotspecs.Specplottergui(
                Specdata={}, 
                Specdiffsets={}, 
                guiroot=self.nodeframes['Plot Spectra'],
                disspecs=self.Nanomap.disspecs
                )
            
            # self.Cube2Imager: updte 
            if hasattr(self, 'Cube2Imager'):
                # update CubeImager's reference to the new Nanomap instance
                if hasattr(self.Cube2Imager, 'gui') and self.Cube2Imager.gui is not None:
                    self.Cube2Imager.gui.Nanomap = self.Nanomap
                self.Cube2Imager.update_bounds()      

        else:
            print("Error", "Failed to load data. Check console for details.")

# sort files and multiple HSI processings
# what is the taks of specfilesorter: 
# 1. Step:
# select a dir with multiple measurement days, on each day HSIs were measured, but some of them went for multiple days. Therefore find the "deepest" folder instance for each folder name. 
# then create a list of all days on which HSIs were measured. For each day, on which on the following day a measurement appears, check if they must be merged. Merge them into a single folder. If no measurement followed, just copy the folder to the savedir (no merge needed)
# 2. Step: 
# for each HSI in savedir, process the HSI with the FileProcessorApp class, create a HSI image, then save it to the savedir as a image. by the HSI name. 

class specfilesorter:
    def __init__(self, tkroot, maindir, filename, fileend, savedir, processdir):
        self.tkroot = tkroot
        self.maindir = maindir
        self.filename = filename
        self.fileend = fileend
        self.savedir = savedir
        self.maindir_entrystr = tk.StringVar()
        self.maindir_entrystr.set(maindir)
        self.filename_entrystr = tk.StringVar()
        self.filename_entrystr.set(filename)
        self.fileend_entrystr = tk.StringVar()
        self.fileend_entrystr.set(fileend)
        self.savedir_entrystr = tk.StringVar()
        self.savedir_entrystr.set(savedir)
        self.processdir_entrystr = tk.StringVar()
        self.processdir_entrystr.set(processdir)
        
        self.files_processed = []
        # UI state
        self.scan_results = []
        self.selected_items = []
        self.merge_var = tk.IntVar(value=0)
        self.buildGUI()

    def buildGUI(self):
        # Main container
        self.sortframe = ttk.LabelFrame(self.tkroot, text="HSI File Sorter", padding=(8, 8))
        self.sortframe.pack(fill='both', expand=True)

        # Left: inputs and controls
        left = tk.Frame(self.sortframe)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

        # Main dir
        tk.Label(left, text="Main directory:").grid(row=0, column=0, sticky='w')
        self.maindir_entry = tk.Entry(left, textvariable=self.maindir_entrystr, width=63)
        self.maindir_entry.grid(row=1, column=0, columnspan=4, sticky='w')
        tk.Button(left, text="Browse...", command=lambda: deflib.browse_folder(self.maindir_entry)).grid(row=1, column=2, padx=4)

        # Filename pattern
        tk.Label(left, text="Filename contains:").grid(row=2, column=0, sticky='w', pady=(8, 0))
        self.filename_entry = tk.Entry(left, textvariable=self.filename_entrystr, width=20)
        self.filename_entry.grid(row=3, column=0, sticky='w')

        # Fileend pattern
        tk.Label(left, text="File extension contains:").grid(row=2, column=1, sticky='w', pady=(8, 0))
        self.fileend_entry = tk.Entry(left, textvariable=self.fileend_entrystr, width=12)
        self.fileend_entry.grid(row=3, column=1, sticky='w')

        # Save dir
        tk.Label(left, text="Save directory:").grid(row=4, column=0, sticky='w', pady=(8, 0))
        self.savedir_entry = tk.Entry(left, textvariable=self.savedir_entrystr, width=63)
        self.savedir_entry.grid(row=5, column=0, columnspan=4, sticky='w')
        tk.Button(left, text="Browse...", command=lambda: deflib.browse_folder(self.savedir_entry)).grid(row=5, column=2, padx=4)

        # Process dir (where merged/processed folders are written)
        tk.Label(left, text="Process directory:").grid(row=6, column=0, sticky='w', pady=(8, 0))
        self.processdir_entry = tk.Entry(left, textvariable=self.processdir_entrystr, width=63)
        self.processdir_entry.grid(row=7, column=0, columnspan=4, sticky='w')
        tk.Button(left, text="Browse...", command=lambda: deflib.browse_folder(self.processdir_entry)).grid(row=7, column=2, padx=4)

        # Options
        tk.Checkbutton(left, text="Merge consecutive days", variable=self.merge_var).grid(row=8, column=0, columnspan=2, sticky='w', pady=(8, 0))

        # Control buttons
        self.btnframe = tk.Frame(left)
        self.btnframe.grid(row=9, column=0, columnspan=3, pady=(12, 0), sticky='w')
        self.scan_button = tk.Button(self.btnframe, text="Scan", width=12, command=self.scan_maindir)
        self.scan_button.pack(side=tk.LEFT)
        self.preview_button = tk.Button(self.btnframe, text="Preview Selected", width=14, command=self.preview_selected)
        self.preview_button.pack(side=tk.LEFT, padx=6)
        self.process_button_sf = tk.Button(self.btnframe, text="Process Selected", width=14, command=self.sort_and_process)
        self.process_button_sf.pack(side=tk.LEFT)
        self.clear_button = tk.Button(self.btnframe, text="Clear", width=8, command=self.clear_list)
        self.clear_button.pack(side=tk.LEFT, padx=6)
        # cancel button for background copy (disabled until a job runs)
        self.cancel_button = tk.Button(self.btnframe, text="Cancel", width=8, command=self.cancel_copy, state='disabled')
        self.cancel_button.pack(side=tk.LEFT, padx=6)

        # Right: results treeview
        right = tk.Frame(self.sortframe)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=6, pady=6)

        cols = ("#1", "#2", "#3")
        self.tree = ttk.Treeview(right, columns=cols, show='headings', selectmode='extended')
        self.tree.heading('#1', text='Folder')
        self.tree.heading('#2', text='Files')
        self.tree.heading('#3', text='Path')
        self.tree.column('#1', width=180)
        self.tree.column('#2', width=60, anchor='center')
        self.tree.column('#3', width=380)

        vsb = ttk.Scrollbar(right, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.LEFT, fill=tk.Y)

        self.tree.bind('<<TreeviewSelect>>', self.on_tree_select)

        # Progress bar
        # below buttonframe, ad a progress barframe
        self.progressfame = tk.Frame(left)
        self.progressfame.grid(row=10, column=0, columnspan=3, pady=(12, 0), sticky='we')
        lablel = tk.Label(self.progressfame, text="Progress:").pack(side=tk.LEFT)
        self.progress = ttk.Progressbar(self.progressfame, orient='horizontal', mode='determinate')
        self.progress.pack(side=tk.LEFT, padx=8, pady=(0, 8), fill=tk.X, expand=True)
        #self.progress.pack(fill=tk.X, padx=8, pady=(0, 8))
        #self.progress.pack(side=tk.LEFT, padx=8, pady=(12, 0), fill=tk.X, expand=True)
    
    def sort_and_process(self):
        # Prepare processing tasks and launch background worker thread to copy files safely
        # Operate on the currently scanned folders (self.scan_results)
        # 1) For each unique folder name, create a folder with that name inside the Process directory
        # 2) If duplicates (same folder name) exist among scanned results, only the first instance is used
        # 3) Copy files matching filename and file extension filters into the created process folder
        # 4) leave a `pass` placeholder for further per-folder processing

        if not self.scan_results:
            print("No data", "No scanned folders available. Run Scan first.")
            return

        savedir = self.savedir_entry.get()
        if not savedir:
            print("No process directory", "Please select a Process directory before processing.")
            return
        # ensure processdir exists
        try:
            os.makedirs(savedir, exist_ok=True)
        except Exception as e:
            print("Error", f"Could not create/access Process directory: {e}")
            return

        # build a mapping of unique folder name -> first seen path
        unique = {}
        for name, count, path in self.scan_results:
            if name not in unique:
                unique[name] = path

    # collect files to copy per unique folder and compute total file count
        filename_filter = self.filename_entry.get()
        fileend_filter = self.fileend_entry.get()
        tasks = []  # list of (folder_name, src_path, [files])
        total_files = 0
        for name, src in unique.items():
            try:
                entries = os.listdir(src)
            except Exception as e:
                print('Could not list directory', src, e)
                continue
            matched = []
            for f in entries:
                full = os.path.join(src, f)
                if not os.path.isfile(full):
                    continue
                if (not filename_filter or filename_filter in f) and (not fileend_filter or fileend_filter in f):
                    matched.append(f)
            if matched:
                tasks.append((name, src, matched))
                total_files += len(matched)

        if not tasks:
            print("No files", "No files matching the filters were found in scanned folders.")
            return

        # perform copying with progress
        copied = 0
        self.progress['maximum'] = total_files
        for name, src, files in tasks:
            target_dir = os.path.join(savedir, name)
            os.makedirs(target_dir, exist_ok=True)
            # copy each matching file
            for f in files:
                srcfile = os.path.join(src, f)
                dstfile = os.path.join(target_dir, f)
                try:
                    shutil.copy2(srcfile, dstfile)
                except Exception as e:
                    print('Failed to copy', srcfile, '->', dstfile, e)
                copied += 1
                # update progress
                try:
                    self.progress['value'] = copied
                    self.sortframe.update_idletasks()
                except Exception:
                    if not tasks:
                        print("No files", "No files matching the filters were found in scanned folders.")
                        return

                    # store task info for worker
                    self._copy_tasks = tasks
                    self._savedir = savedir
                    self._total_files = total_files
                    self._copied_count = 0
                    self._stop_event = thr.Event()
                    self._lock = thr.Lock()

                    # disable controls that should not be used while copying
                    self.scan_button.config(state='disabled')
                    self.preview_button.config(state='disabled')
                    self.process_button_sf.config(state='disabled')
                    self.clear_button.config(state='disabled')
                    self.cancel_button.config(state='normal')

                    # prepare progress
                    self.progress['maximum'] = total_files
                    self.progress['value'] = 0

                    # launch background thread (daemon so it won't block exit)
                    self._copy_thread = thr.Thread(target=self._copy_worker, daemon=True)
                    self._copy_thread.start()
                    return

        # enforce merge consecutive days if the option is selected (synchronous copy path)
        if self.merge_var.get() == 1:
            self._merge_consecutive_days(savedir)

    def _merge_consecutive_days(self, savedir):
        """Merge folders that represent the same measurement split across consecutive days.
        
        Folders are identified as belonging to the same measurement if they share the same
        name except for the 8-digit date portion (YYYYMMDD).  All files from a later-day
        folder are moved into the earliest-day folder, and the now-empty later-day folder
        is removed.
        """
        folderlist = [f for f in os.listdir(savedir) if os.path.isdir(os.path.join(savedir, f))]

        # Group folders by (prefix, suffix) key to identify measurements split across days
        groups = {}
        for folder in folderlist:
            # Regex to find 8-digit date (starting with 19 or 20) in the folder name
            # Pattern matches: Prefix + YYYYMMDD + Suffix
            match = re.search(r'^(.*?)((?:19|20)\d{6})(.*)$', folder)
            if match:
                prefix, date_str, suffix = match.groups()
                try:
                    date_obj = datet.datetime.strptime(date_str, "%Y%m%d").date()
                    # The key is the identifier (Name + Info) without the date
                    key = (prefix, suffix)
                    if key not in groups:
                        groups[key] = []
                    groups[key].append({'date': date_obj, 'folder': folder, 'path': os.path.join(savedir, folder)})
                except ValueError:
                    continue  # Skip if date parsing fails

        # Process each group to find and merge consecutive sequences
        for key, items in groups.items():
            # Sort folders by date so we can find sequences naturally
            items.sort(key=lambda x: x['date'])

            if len(items) < 2:
                continue

            # 'base_idx' points to the folder that starts the current chain (earliest day)
            base_idx = 0

            # Iterate starting from the second item to check consecutiveness
            for i in range(1, len(items)):
                curr_item = items[i]
                prev_item = items[i - 1]

                # Check strictly for consecutive days (current day == previous day + 1)
                if curr_item['date'] - prev_item['date'] == datet.timedelta(days=1):
                    # It is consecutive; merge 'curr' into the 'base' folder
                    base_folder = items[base_idx]['path']
                    src_folder = curr_item['path']

                    try:
                        # Move all contents from the later day to the first day
                        for item_name in os.listdir(src_folder):
                            src_path = os.path.join(src_folder, item_name)
                            dst_path = os.path.join(base_folder, item_name)

                            # Handle filename collisions: append date to avoid overwriting
                            if os.path.exists(dst_path):
                                f_name, f_ext = os.path.splitext(item_name)
                                # e.g., image.png -> image_20250102.png
                                new_name = f"{f_name}_{curr_item['date'].strftime('%Y%m%d')}{f_ext}"
                                dst_path = os.path.join(base_folder, new_name)

                            shutil.move(src_path, dst_path)

                        # Remove the now-empty source directory (the later day)
                        os.rmdir(src_folder)
                        print(f"Merged consecutive day: {src_folder} -> {base_folder}")

                    except Exception as e:
                        print(f"Error merging {src_folder} into {base_folder}: {e}")

                else:
                    # Sequence broken (gap in days); current becomes the new base
                    base_idx = i

    def _copy_worker(self):
        """Background worker that copies files and updates progress via the main thread."""
        try:
            for name, src, files in self._copy_tasks:
                if self._stop_event.is_set():
                    break
                target_dir = os.path.join(self._savedir, name)
                os.makedirs(target_dir, exist_ok=True)
                for f in files:
                    if self._stop_event.is_set():
                        break
                    srcfile = os.path.join(src, f)
                    dstfile = os.path.join(target_dir, f)
                    try:
                        shutil.copy2(srcfile, dstfile)
                    except Exception as e:
                        print('Failed to copy', srcfile, '->', dstfile, e)
                    # increment counter in thread-safe manner
                    with self._lock:
                        self._copied_count += 1
                        copied = self._copied_count
                    # schedule a GUI update on the main thread
                    try:
                        self.sortframe.after(0, self._update_progress, copied)
                    except Exception:
                        pass
                # placeholder for further processing per folder
                # pass

            # after the loop, enforce merge consecutive days if the option is selected
            if self.merge_var.get() == 1:
                self._merge_consecutive_days(self._savedir)

        finally:
            # schedule finalizer on main thread
            self.sortframe.after(0, self._copy_finished)

    def _update_progress(self, value):
        try:
            self.progress['value'] = value
        except Exception:
            pass

    def cancel_copy(self):
        # signal worker to stop
        if hasattr(self, '_stop_event') and self._stop_event:
            self._stop_event.set()
            # disable cancel to avoid repeated clicks
            try:
                self.cancel_button.config(state='disabled')
            except Exception:
                pass

    def _copy_finished(self):
        # re-enable controls
        try:
            self.scan_button.config(state='normal')
            self.preview_button.config(state='normal')
            self.process_button_sf.config(state='normal')
            self.clear_button.config(state='normal')
            self.cancel_button.config(state='disabled')
        except Exception:
            pass
        # inform user
        copied = getattr(self, '_copied_count', 0)
        total = getattr(self, '_total_files', 0)
        if getattr(self, '_stop_event', None) and self._stop_event.is_set():
            print("Cancelled", f"Copy cancelled after {copied} of {total} files.")
        else:
            print("Done", f"Copied {copied} files into {getattr(self, '_savedir', '')}")
        # reset progress
        try:
            self.progress['value'] = 0
        except Exception:
            pass
        self.scan_results.sort()
        self.populate_list()
    
    def scan_maindir(self):
        folder = self.maindir_entry.get()
        filename_filter = self.filename_entry.get()
        fileend_filter = self.fileend_entry.get()

        if not folder:
            print("Error", "Please select a main directory to scan.")
            return
        if not os.path.exists(folder):
            print("Error", "The selected main directory does not exist.")
            return

        self.scan_results = []
        for dirpath, dirnames, filenames in os.walk(folder):
            matched_files = []
            for f in filenames:
                if (not filename_filter or filename_filter in f) and (not fileend_filter or fileend_filter in f):
                    matched_files.append(f)
            if matched_files:
                folder_name = os.path.basename(dirpath)
                file_count = len(matched_files)
                self.scan_results.append((folder_name, file_count, dirpath))

        self.scan_results.sort()
        self.populate_list()

    def populate_list(self):
        self.tree.delete(*self.tree.get_children())
        for name, count, path in self.scan_results:
            self.tree.insert('', 'end', values=(name, count, path))

    def on_tree_select(self, event):
        sel = self.tree.selection()
        self.selected_items = [self.tree.set(i, '#3') for i in sel]

    def preview_selected(self):
        # Open the folder(s) in file explorer where possible
        if not self.selected_items:
            print("No selection", "Please select a folder to preview.")
            return
        for p in self.selected_items:
            try:
                if sys.platform == 'win32':
                    os.startfile(p)
                else:
                    # cross-platform fallback
                    import subprocess
                    subprocess.run(['xdg-open', p], check=False)
            except Exception as e:
                print('Could not open folder', p, e)

    def clear_list(self):
        self.tree.delete(*self.tree.get_children())
        self.scan_results = []
        self.selected_items = []

    def make_HSI_from_folder(self, folder):
        # set self.folder_entry, self.filename_entry and self.fileformat_entry accordingly
        pass

        #self.folder_entry.delete(0, tk.END)
        #self.folder_entry.insert(0, folder)


def pressclose(root, app):
    # Get all running threads
    for thread in thr.enumerate():
        # Skip the main thread
        if thread is not thr.main_thread():
            # Set daemon to True so thread terminates when main program exits
            thread.daemon = True
    
    # Run app cleanup first so widgets/canvases are destroyed while Tk is still valid.
    app.on_closing()
    # Then terminate tkinter loop and root window.
    try:
        root.quit()
    except Exception:
        pass
    root.destroy()

if __name__ == "__main__":
    # Initialize centralized error handling and logging
    # This creates a rotating log file at logs/specmap.log (10MB max, 5 backups)
    # Provides thread-safe logging and user-friendly error messages via tkinter messagebox
    error_engine = error_handler.ErrorEngine(
        log_file='logs/specmap.log',
        max_bytes=10*1024*1024,  # 10MB
        backup_count=5
    )
    
    # Set as default error engine for application-wide use
    # This allows lib9.py and other modules to use the same error handler
    error_handler._default_error_engine = error_engine
    
    # Get logger for application-wide use
    logger = error_engine.get_logger()
    logger.info("SpecMap application starting", extra={'context': 'Application Startup'})
    
    # init debugger
    debugger = DBG.main_Debugger()

    # Create the main window

    # default values for GUI
    defaults = deflib.initdefaults()
    print('Starting...')

    root = tk.Tk()
    root.geometry('{}x{}'.format(int(defaults['windowsize_X']), int(defaults['windowsize_Y'])))
    frame = tk.Frame(root)
    
    # Store error_engine in app for use throughout application
    app = FileProcessorApp(root, defaults)
    app.error_engine = error_engine  # Make error engine accessible to app

    # Set the protocol for closing the window
    root.protocol("WM_DELETE_WINDOW", lambda: pressclose(root, app))

    logger.info("GUI initialized successfully", extra={'context': 'Application Startup'})
    
    # Run the application
    root.mainloop()
    
    logger.info("SpecMap application shutting down", extra={'context': 'Application Shutdown'})