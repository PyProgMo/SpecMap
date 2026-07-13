import tkinter as tk
from tkinter import ttk, filedialog
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button
import numpy as np
from math import floor, ceil
import re
import os

matplotlib.use('TkAgg')

class PlotHSI:
    """
    Interactive HSI Plotter with Tkinter GUI.
    
    Features:
    - Entry GUI for all plot parameters
    - Separate plot window with drag-and-drop elements
    - Update capability to refresh plot with new parameters
    - Transparent structure for easy element management
    
    Usage:
        root = tk.Tk()
        plotter = PlotHSI(root)
        root.mainloop()
    """
    
    def __init__(self, root, get_hsi_func=None):
        """
        Initialize the PlotHSI GUI.
        
        Args:
            root: Tkinter root window
            get_hsi_func: Optional function to load HSI data (signature: func(filepath) -> (metadata, data))
        """
        self.root = root
        #self.root.title("HSI Plotter")
        self.get_hsi_func = get_hsi_func
        
        # Data storage
        self.data = None
        self.metadata = None
        self.current_file = None
        
        # Plot window and elements
        self.plot_window = None
        self.fig = None
        self.ax = None
        self.canvas = None
        
        # Plot element storage (transparent structure for easy access)
        self.elements = {
            'image': None,          # imshow object
            'colorbar': None,       # colorbar object
            'colorbar_ax': None,    # colorbar axes
            'scalebar': None,       # Rectangle patch
            'scalebar_text': None,  # Text artist
            'image_label': None,    # Text artist
            'save_button': None,    # Button widget
            'save_button_ax': None  # Button axes
        }
        
        # Tk variables for all parameters
        self._init_variables()
        
        # Create GUI
        self._create_gui()
    
    def _init_variables(self):
        """Initialize all Tkinter variables with default values."""
        # Colormap settings
        self.cmap_var = tk.StringVar(value='twilight')
        self.vmin_var = tk.StringVar(value='')
        self.vmax_var = tk.StringVar(value='')
        
        # Scalebar settings
        self.scalebarlength_var = tk.IntVar(value=20)
        self.scalebarwidth_var = tk.IntVar(value=2)
        self.scalebarpos_x_var = tk.IntVar(value=50)
        self.scalebarpos_y_var = tk.IntVar(value=50)
        self.dx_var = tk.DoubleVar(value=1.0)
        self.unit_var = tk.StringVar(value='$\\mu m$')
        
        # Figure settings
        self.figsize_w_var = tk.IntVar(value=7)
        self.figsize_h_var = tk.IntVar(value=6)
        self.title_var = tk.StringVar(value='')
        self.xlabel_var = tk.StringVar(value='')
        self.ylabel_var = tk.StringVar(value='')
        self.title_color_var = tk.StringVar(value='black')
        
        # Colorbar settings
        self.show_colorbar_var = tk.BooleanVar(value=True)
        self.cbar_unit_var = tk.StringVar(value='')
        self.cmapcaption_var = tk.StringVar(value='')
        self.cbar_divisor_var = tk.DoubleVar(value=1.0)
        self.cbar_caption_color_var = tk.StringVar(value='white')
        self.cbar_border_color_var = tk.StringVar(value='black')
        self.cbar_upper_num_color_var = tk.StringVar(value='black')
        self.cbar_lower_num_color_var = tk.StringVar(value='black')
        
        # Scalebar color settings
        self.scalebar_face_color_var = tk.StringVar(value='white')
        self.scalebar_edge_color_var = tk.StringVar(value='black')
        self.scalebar_text_color_var = tk.StringVar(value='black')
        
        # Font size settings (separate for each element)
        self.title_fontsize_var = tk.IntVar(value=12)
        self.axis_label_fontsize_var = tk.IntVar(value=10)
        self.colorbar_fontsize_var = tk.IntVar(value=27)
        self.colorbar_caption_fontsize_var = tk.IntVar(value=27)
        self.scalebar_fontsize_var = tk.IntVar(value=27)
        self.image_label_fontsize_var = tk.IntVar(value=27)
        
        # Label settings
        self.image_label_var = tk.StringVar(value='(525$\\pm$25)nm')
        self.image_label_color_var = tk.StringVar(value='black')
        
        # Visibility toggles for easy element management
        self.show_scalebar_var = tk.BooleanVar(value=True)
        self.show_image_label_var = tk.BooleanVar(value=True)
        self.enable_drag_var = tk.BooleanVar(value=True)
    
    def _create_gui(self):
        """Create the main GUI with scrollable parameter entries in 2 columns."""
        # Create canvas for scrolling (both vertical and horizontal)
        canvas = tk.Canvas(self.root)
        v_scrollbar = ttk.Scrollbar(self.root, orient="vertical", command=canvas.yview)
        h_scrollbar = ttk.Scrollbar(self.root, orient="horizontal", command=canvas.xview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        # Layout with grid
        canvas.grid(row=0, column=0, sticky=tk.N+tk.S+tk.E+tk.W)
        v_scrollbar.grid(row=0, column=1, sticky=tk.N+tk.S)
        h_scrollbar.grid(row=1, column=0, sticky=tk.E+tk.W)
        
        # Configure root grid weights
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        
        # Main container with 2 columns
        main_container = ttk.Frame(scrollable_frame, padding="10")
        main_container.grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)
        
        # Left column frame
        left_frame = ttk.Frame(main_container, padding="5")
        left_frame.grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S, padx=5)
        
        # Right column frame
        right_frame = ttk.Frame(main_container, padding="5")
        right_frame.grid(row=0, column=1, sticky=tk.W+tk.E+tk.N+tk.S, padx=5)
        
        # LEFT COLUMN
        row = 0
        
        # === File Loading Section ===
        ttk.Label(left_frame, text="=== File Loading ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(0, 5), sticky=tk.W)
        row += 1
        
        ttk.Button(left_frame, text="Load HSI Data", command=self._load_data).grid(
            row=row, column=0, columnspan=2, pady=5, sticky=tk.W+tk.E)
        row += 1
        
        self.file_label = ttk.Label(left_frame, text="No file loaded", foreground="gray")
        self.file_label.grid(row=row, column=0, columnspan=2, pady=(0, 10), sticky=tk.W)
        row += 1
        
        # === Colormap Section ===
        ttk.Label(left_frame, text="=== Colormap Settings ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(5, 5), sticky=tk.W)
        row += 1
        
        self._add_entry(left_frame, row, "Colormap:", self.cmap_var)
        row += 1
        self._add_entry(left_frame, row, "V Min:", self.vmin_var)
        row += 1
        self._add_entry(left_frame, row, "V Max:", self.vmax_var)
        row += 1
        
        # === HSI Image Settings ===
        ttk.Label(left_frame, text="=== HSI Image Settings ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1

        hsi_btn_frame = ttk.Frame(left_frame)
        hsi_btn_frame.grid(row=row, column=0, columnspan=2, pady=5, sticky=tk.W+tk.E)
        ttk.Button(hsi_btn_frame, text="Mirror Horizontally", command=self._mirror_horizontal).pack(side=tk.LEFT, padx=2)
        ttk.Button(hsi_btn_frame, text="Mirror Vertically", command=self._mirror_vertical).pack(side=tk.LEFT, padx=2)
        ttk.Button(hsi_btn_frame, text="Rotate 90°", command=self._rotate_90).pack(side=tk.LEFT, padx=2)
        row += 1

        # === Scalebar Section ===
        ttk.Label(left_frame, text="=== Scalebar Settings ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1
        
        ttk.Checkbutton(left_frame, text="Show Scalebar", variable=self.show_scalebar_var).grid(
            row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1
        
        self._add_entry(left_frame, row, "Length (px):", self.scalebarlength_var)
        row += 1
        self._add_entry(left_frame, row, "Width:", self.scalebarwidth_var)
        row += 1
        self._add_entry(left_frame, row, "Position X:", self.scalebarpos_x_var)
        row += 1
        self._add_entry(left_frame, row, "Position Y:", self.scalebarpos_y_var)
        row += 1
        self._add_entry(left_frame, row, "dx (scale):", self.dx_var)
        row += 1
        self._add_entry(left_frame, row, "Unit:", self.unit_var)
        row += 1
        
        # === Figure Settings ===
        ttk.Label(left_frame, text="=== Figure Settings ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1
        
        self._add_entry(left_frame, row, "Width:", self.figsize_w_var)
        row += 1
        self._add_entry(left_frame, row, "Height:", self.figsize_h_var)
        row += 1
        self._add_entry(left_frame, row, "Title:", self.title_var)
        row += 1
        self._add_entry(left_frame, row, "Title Color:", self.title_color_var)
        row += 1
        self._add_entry(left_frame, row, "X Label:", self.xlabel_var)
        row += 1
        self._add_entry(left_frame, row, "Y Label:", self.ylabel_var)
        row += 1
        self._add_entry(left_frame, row, "Title Font Size:", self.title_fontsize_var)
        row += 1
        self._add_entry(left_frame, row, "Axis Label Font Size:", self.axis_label_fontsize_var)
        row += 1
        
        # RIGHT COLUMN
        row = 0
        
        # === Colorbar Section ===
        ttk.Label(right_frame, text="=== Colorbar Settings ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(0, 5), sticky=tk.W)
        row += 1
        
        ttk.Checkbutton(right_frame, text="Show Colorbar", variable=self.show_colorbar_var).grid(
            row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1
        
        self._add_entry(right_frame, row, "Unit:", self.cbar_unit_var)
        row += 1
        self._add_entry(right_frame, row, "Caption:", self.cmapcaption_var)
        row += 1
        self._add_entry(right_frame, row, "Divisor:", self.cbar_divisor_var)
        row += 1
        self._add_entry(right_frame, row, "Numbers Font Size:", self.colorbar_fontsize_var)
        row += 1
        self._add_entry(right_frame, row, "Caption Font Size:", self.colorbar_caption_fontsize_var)
        row += 1
        self._add_entry(right_frame, row, "Caption Color:", self.cbar_caption_color_var)
        row += 1
        self._add_entry(right_frame, row, "Border Color:", self.cbar_border_color_var)
        row += 1
        self._add_entry(right_frame, row, "Upper Number Color:", self.cbar_upper_num_color_var)
        row += 1
        self._add_entry(right_frame, row, "Lower Number Color:", self.cbar_lower_num_color_var)
        row += 1
        
        # === Scalebar Color Settings ===
        ttk.Label(right_frame, text="=== Scalebar Colors ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1
        
        self._add_entry(right_frame, row, "Bar Fill Color:", self.scalebar_face_color_var)
        row += 1
        self._add_entry(right_frame, row, "Bar Edge Color:", self.scalebar_edge_color_var)
        row += 1
        self._add_entry(right_frame, row, "Text Color:", self.scalebar_text_color_var)
        row += 1
        self._add_entry(right_frame, row, "Font Size:", self.scalebar_fontsize_var)
        row += 1
        
        # === Labels Section ===
        ttk.Label(right_frame, text="=== Label Settings ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1
        
        ttk.Checkbutton(right_frame, text="Show Image Label", variable=self.show_image_label_var).grid(
            row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1
        
        self._add_entry(right_frame, row, "Image Label:", self.image_label_var)
        row += 1
        self._add_entry(right_frame, row, "Label Color:", self.image_label_color_var)
        row += 1
        self._add_entry(right_frame, row, "Font Size:", self.image_label_fontsize_var)
        row += 1
        
        # === Configuration Management ===
        ttk.Label(right_frame, text="=== Configuration ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1
        
        config_frame = ttk.Frame(right_frame)
        config_frame.grid(row=row, column=0, columnspan=2, pady=5, sticky=tk.W+tk.E)
        ttk.Button(config_frame, text="Save Config", command=self._save_config).pack(side=tk.LEFT, padx=2)
        ttk.Button(config_frame, text="Load Config", command=self._load_config).pack(side=tk.LEFT, padx=2)
        row += 1
        
        # === Interaction Settings ===
        ttk.Label(right_frame, text="=== Interaction ===", font=('Arial', 10, 'bold')).grid(
            row=row, column=0, columnspan=2, pady=(10, 5), sticky=tk.W)
        row += 1
        
        ttk.Checkbutton(right_frame, text="Enable Drag & Drop", variable=self.enable_drag_var).grid(
            row=row, column=0, columnspan=2, sticky=tk.W)
        row += 1
        
        # === Action Buttons (span both columns) ===
        button_frame = ttk.Frame(main_container)
        button_frame.grid(row=1, column=0, columnspan=2, pady=15)
        
        ttk.Button(button_frame, text="Plot HSI", command=self._plot_hsi, 
                  style='Accent.TButton').pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Update Plot", command=self._update_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Close Plot", command=self._close_plot).pack(side=tk.LEFT, padx=5)
        
        # Configure column weights for both frames
        left_frame.columnconfigure(1, weight=1)
        right_frame.columnconfigure(1, weight=1)
    
    def _add_entry(self, parent, row, label_text, variable):
        """Helper to add a labeled entry field."""
        ttk.Label(parent, text=label_text).grid(row=row, column=0, sticky=tk.W, pady=2)
        ttk.Entry(parent, textvariable=variable, width=25).grid(
            row=row, column=1, sticky=(tk.W, tk.E), pady=2)
    
    def _load_data(self):
        """Load HSI data from file."""
        filename = filedialog.askopenfilename(
            title="Select HSI file",
            filetypes=[("All files", "*.*")]
        )
        if not filename:
            return
        
        try:
            if self.get_hsi_func:
                self.metadata, self.data = self.get_hsi_func(filename)
            else:
                # use 20x20 dummy data if no function provided
                self.metadata = {'wlstart': '400', 'wlend': '1000'}
                self.data = np.random.rand(100, 200)
            
            self.current_file = filename
            self.file_label.config(
                text=f"Loaded: {os.path.basename(filename)} | Shape: {self.data.shape}",
                foreground="green"
            )
            print(f"Successfully loaded: {filename}")
            print(f"Data shape: {self.data.shape}")
            print(f"Metadata: {self.metadata}")
            
        except Exception as e:
            self.file_label.config(text=f"Error: {str(e)}", foreground="red")
            print(f"Error loading data: {e}")
    
    def _mirror_horizontal(self):
        """Mirror the loaded HSI data horizontally (left-right) and refresh the plot."""
        if self.data is None:
            print("No data loaded!")
            return
        self.data = np.fliplr(self.data)
        if self.plot_window is not None:
            self._update_plot()

    def _mirror_vertical(self):
        """Mirror the loaded HSI data vertically (up-down) and refresh the plot."""
        if self.data is None:
            print("No data loaded!")
            return
        self.data = np.flipud(self.data)
        if self.plot_window is not None:
            self._update_plot()

    def _rotate_90(self):
        """Rotate the loaded HSI data by 90° counter-clockwise and refresh the plot."""
        if self.data is None:
            print("No data loaded!")
            return
        self.data = np.rot90(self.data)
        if self.plot_window is not None:
            self._update_plot()

    def _plot_hsi(self):
        """Create initial plot in separate window."""
        if self.data is None:
            print("No data loaded! Please load HSI data first.")
            return
        
        # Close existing plot window
        if self.plot_window is not None:
            try:
                self.plot_window.destroy()
            except:
                pass
        
        # Create new window
        self.plot_window = tk.Toplevel(self.root)
        self.plot_window.title("HSI Plot - Interactive")
        
        # Create the plot
        self._create_plot()
    
    def _update_plot(self):
        """Update existing plot with new parameters."""
        if self.plot_window is None or self.data is None:
            print("No active plot to update! Create a plot first.")
            return
        
        # Clear old plot
        if self.fig is not None:
            plt.close(self.fig)
        
        # Recreate plot
        self._create_plot()
    
    def _close_plot(self):
        """Close the plot window."""
        if self.plot_window is not None:
            try:
                if self.fig is not None:
                    plt.close(self.fig)
                self.plot_window.destroy()
                self.plot_window = None
            except:
                pass
    
    def _create_plot(self):
        """Create the matplotlib plot with all interactive features."""
        # Get parameters
        vmin = None if self.vmin_var.get() == '' else float(self.vmin_var.get())
        vmax = None if self.vmax_var.get() == '' else float(self.vmax_var.get())
        figsize = (self.figsize_w_var.get(), self.figsize_h_var.get())
        
        # Create figure
        self.fig, self.ax = plt.subplots(figsize=figsize)
        
        # Calculate extent for wavelength data
        extent, wl0, wl1 = self._calculate_extent()
        
        # Create image
        data_array = np.array(self.data, dtype=float)
        self.elements['image'] = self.ax.imshow(
            data_array, origin='lower', aspect='auto',
            cmap=self.cmap_var.get(), vmin=vmin, vmax=vmax, extent=extent
        )
        
        # Set labels with font sizes
        self.ax.set_title(self.title_var.get(), color=self.title_color_var.get(), 
                         fontsize=self.title_fontsize_var.get())
        self.ax.set_xlabel(self.xlabel_var.get(), fontsize=self.axis_label_fontsize_var.get())
        self.ax.set_ylabel(self.ylabel_var.get(), fontsize=self.axis_label_fontsize_var.get())
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        
        # Add elements based on visibility settings
        if self.show_colorbar_var.get():
            self._add_colorbar()
        
        if self.show_scalebar_var.get():
            self._add_scalebar(extent, wl0, wl1)
        
        if self.show_image_label_var.get():
            self._add_image_label()
        
        # Add save button
        self._add_save_button()
        
        # Setup drag handlers
        if self.enable_drag_var.get():
            self._setup_drag_handlers()
        
        plt.tight_layout()
        
        # Embed in Tkinter
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_window)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def _calculate_extent(self):
        """Calculate extent for wavelength data."""
        extent = None
        wl0 = wl1 = None
        try:
            if self.metadata and 'wlstart' in self.metadata and 'wlend' in self.metadata:
                wl0 = round(float(self.metadata['wlstart']), 1)
                wl1 = round(float(self.metadata['wlend']), 1)
                extent = (wl0, wl1, 0, self.data.shape[0])
        except Exception:
            pass
        return extent, wl0, wl1
    
    def _add_colorbar(self):
        """Add inset colorbar to the plot."""
        scale_factor = self.colorbar_fontsize_var.get() / 14.0
        ax_pos = self.ax.get_position()
        
        # Calculate colorbar dimensions
        cax_w = ax_pos.width * 0.08 * scale_factor
        cax_h = ax_pos.height * 0.28 * (1.0 + 0.15 * (scale_factor - 1.0))
        cax_x = ax_pos.x0 + ax_pos.width * (1.0 - 0.08 * scale_factor - 0.10)
        cax_y = ax_pos.y0 + ax_pos.height * 0.06
        
        if self.fig is not None:
            self.elements['colorbar_ax'] = self.fig.add_axes((cax_x, cax_y, cax_w, cax_h))  # type: ignore
        cax = self.elements['colorbar_ax']
        
        # Style
        cax.patch.set_facecolor('white')
        cax.patch.set_edgecolor('white')
        cax.patch.set_linewidth(2.6 * scale_factor)
        
        self.elements['colorbar'] = self.fig.colorbar(
            self.elements['image'], cax=cax, orientation='vertical'
        )
        
        # Apply border color from settings
        border_color = self.cbar_border_color_var.get()
        for spine in cax.spines.values():
            spine.set_edgecolor(border_color)
            spine.set_linewidth(0.8 * scale_factor)
        
        # Add labels
        self._add_colorbar_labels(cax, scale_factor)
    
    def _add_colorbar_labels(self, cax, scale_factor):
        """Add labels to colorbar."""
        try:
            clim = self.elements['image'].get_clim()
            vminr = floor(clim[0])
            vmaxr = ceil(clim[1])
            if vminr == vmaxr:
                vminr = int(round(clim[0]))
                vmaxr = int(round(clim[1]))
            
            self.elements['colorbar'].set_ticks([])
            cax.yaxis.set_ticks([])
            
            # Apply divisor
            divisor = self.cbar_divisor_var.get()
            if divisor == 0:
                divisor = 1.0  # Avoid division by zero
            
            vminr_display = vminr / divisor
            vmaxr_display = vmaxr / divisor
            
            # Format labels
            if self.cbar_unit_var.get():
                top_label = f"{int(vmaxr_display)} {self.cbar_unit_var.get()}"
                bot_label = f"{int(vminr_display)} {self.cbar_unit_var.get()}"
            else:
                top_label = f"{int(vmaxr_display)}"
                bot_label = f"{int(vminr_display)}"
            
            # Get font sizes
            numbers_fontsize = self.colorbar_fontsize_var.get()
            caption_fontsize = self.colorbar_caption_fontsize_var.get()
            
            # Get colors from settings
            upper_color = self.cbar_upper_num_color_var.get()
            lower_color = self.cbar_lower_num_color_var.get()
            caption_color = self.cbar_caption_color_var.get()
            border_color = self.cbar_border_color_var.get()
            
            cax.text(0.5, 0.94, top_label, transform=cax.transAxes,
                    ha='center', va='top', fontsize=numbers_fontsize, color=upper_color, weight='bold')
            cax.text(0.5, 0.06, bot_label, transform=cax.transAxes,
                    ha='center', va='bottom', fontsize=numbers_fontsize, color=lower_color, weight='bold')
            cax.text(0.5, 1.05, self.cmapcaption_var.get(), transform=cax.transAxes,
                    ha='center', va='bottom', fontsize=caption_fontsize, color=caption_color, weight='bold')
            
            self.elements['colorbar'].outline.set_edgecolor(border_color)
            self.elements['colorbar'].outline.set_linewidth(0.8 * scale_factor)
        except Exception as e:
            print(f"Warning: Could not add colorbar labels: {e}")
    
    def _add_scalebar(self, extent, wl0, wl1):
        """Add scalebar to the plot."""
        fontsize = self.scalebar_fontsize_var.get()
        scale_factor = fontsize / 14.0
        sbx = self.scalebarpos_x_var.get()
        sby = self.scalebarpos_y_var.get()
        scalebarlength = self.scalebarlength_var.get()
        scaled_width = self.scalebarwidth_var.get() * scale_factor
        
        # Get colors from settings
        face_color = self.scalebar_face_color_var.get()
        edge_color = self.scalebar_edge_color_var.get()
        text_color = self.scalebar_text_color_var.get()
        
        # Get image dimensions in pixels (rows, cols)
        nrows, ncols = self.data.shape
        
        # The scalebar position and length are in pixel coordinates
        # These correspond to the actual image array indices
        # Since imshow displays with origin='lower', y=0 is at bottom, y=nrows at top
        # x ranges from 0 to ncols (or extent if specified)
        
        # If extent is specified (wavelength mode), convert pixel positions to extent coordinates
        if extent is not None and wl0 is not None and wl1 is not None:
            # extent = (wl0, wl1, 0, nrows)
            # Convert pixel x-position to wavelength coordinate
            x_scale = (wl1 - wl0) / ncols
            sbx_display = wl0 + sbx * x_scale
            sb_width_display = scalebarlength * x_scale
            label_x = sbx_display + (sb_width_display / 2)
            # Y coordinate stays in pixel space (0 to nrows)
            sby_display = sby
        else:
            # No extent - use direct pixel coordinates
            sbx_display = sbx
            sby_display = sby
            sb_width_display = scalebarlength
            label_x = sbx + (scalebarlength / 2)
        
        # Create scalebar rectangle
        if self.ax is not None:
            self.elements['scalebar'] = Rectangle(
                (sbx_display, sby_display), sb_width_display, scaled_width,
                facecolor=face_color, edgecolor=edge_color,
                linewidth=1.1 * scale_factor, zorder=5,
                transform=self.ax.transData  # type: ignore
            )
            self.ax.add_patch(self.elements['scalebar'])
        
        # Add label
        phys_len = scalebarlength * self.dx_var.get()
        formatted_len = f"{phys_len:.2f}".rstrip('0').rstrip('.')
        label_text = f"{formatted_len} {self.unit_var.get()}"
        
        if self.ax is not None:
            self.elements['scalebar_text'] = self.ax.text(
                label_x, sby_display + scaled_width + 1 * scale_factor, label_text,
                color=text_color, ha='center', va='bottom',
                fontsize=fontsize, weight='bold', zorder=6
            )
    
    def _add_image_label(self):
        """Add draggable image label."""
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        xpos = xlim[1] - 0.02 * (xlim[1] - xlim[0])
        ypos = ylim[1] - 0.02 * (ylim[1] - ylim[0])
        
        label_color = self.image_label_color_var.get()
        fontsize = self.image_label_fontsize_var.get()
        
        self.elements['image_label'] = self.ax.text(
            xpos, ypos, self.image_label_var.get(),
            color=label_color, ha='right', va='top',
            fontsize=fontsize, weight='bold', zorder=6
        )
    
    def _add_save_button(self):
        """Add save button to figure."""
        try:
            self.elements['save_button_ax'] = self.fig.add_axes([0.01, 0.01, 0.13, 0.04])
            self.elements['save_button'] = Button(self.elements['save_button_ax'], 'Save 600 dpi')
            self.elements['save_button'].on_clicked(self._on_save)
        except Exception as e:
            print(f"Warning: Could not create save button: {e}")
    
    def _on_save(self, event=None):
        """Save the figure at 600 DPI."""
        save_outline = None
        bax = self.elements['save_button_ax']
        cax = self.elements['colorbar_ax']
        
        try:
            # Hide save button
            if bax is not None:
                bax.set_visible(False)
            
            # Temporarily modify colorbar for saving
            if cax is not None:
                orig_face = cax.patch.get_facecolor()
                orig_edge = cax.patch.get_edgecolor()
                orig_linewidth = cax.patch.get_linewidth()
                orig_spines = {k: v.get_visible() for k, v in cax.spines.items()}
                
                cax.patch.set_facecolor('none')
                cax.patch.set_edgecolor('none')
                cax.patch.set_linewidth(0)
                for spine in cax.spines.values():
                    spine.set_visible(False)
                
                pos = cax.get_position()
                save_outline = Rectangle(
                    (pos.x0, pos.y0), pos.width, pos.height,
                    transform=self.fig.transFigure, fill=False,
                    edgecolor='white', linewidth=max(1.5, orig_linewidth),
                    zorder=10, clip_on=False
                )
                self.fig.add_artist(save_outline)
            
            self.fig.canvas.draw()
        except Exception as e:
            print(f"Warning during save preparation: {e}")
        
        # Ask for filename
        fname = filedialog.asksaveasfilename(
            defaultextension='.png',
            filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg')]
        )
        
        if fname:
            self.fig.savefig(fname, dpi=600, bbox_inches='tight', pad_inches=0)
            print(f'Saved: {fname}')
        
        # Restore visibility
        try:
            if bax is not None:
                bax.set_visible(True)
            if cax is not None:
                cax.patch.set_facecolor(orig_face)
                cax.patch.set_edgecolor(orig_edge)
                cax.patch.set_linewidth(orig_linewidth)
                for k, v in cax.spines.items():
                    v.set_visible(orig_spines[k])
                if save_outline is not None:
                    save_outline.remove()
            self.fig.canvas.draw()
        except Exception as e:
            print(f"Warning during restore: {e}")
    
    def _setup_drag_handlers(self):
        """Setup interactive drag-and-drop handlers."""
        self.drag_state = {'target': None, 'offset': (0, 0)}
        
        self.fig.canvas.mpl_connect('button_press_event', self._on_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)
    
    def _on_press(self, event):
        """Handle mouse press for dragging."""
        # Check scalebar
        if self.elements['scalebar'] is not None and event.inaxes == self.ax:
            contains, _ = self.elements['scalebar'].contains(event)
            if contains and event.xdata and event.ydata:
                self.drag_state['target'] = 'scalebar'
                self.drag_state['offset'] = (
                    event.xdata - self.elements['scalebar'].get_x(),
                    event.ydata - self.elements['scalebar'].get_y()
                )
                return
        
        # Check image label
        if self.elements['image_label'] is not None and event.inaxes == self.ax:
            try:
                contains, _ = self.elements['image_label'].contains(event)
                if contains and event.xdata and event.ydata:
                    self.drag_state['target'] = 'image_label'
                    pos = self.elements['image_label'].get_position()
                    self.drag_state['offset'] = (event.xdata - pos[0], event.ydata - pos[1])
                    return
            except:
                pass
        
        # Check colorbar
        if self.elements['colorbar_ax'] is not None:
            fx, fy = self.fig.transFigure.inverted().transform((event.x, event.y))
            cpos = self.elements['colorbar_ax'].get_position()
            if 0 <= fx <= 1 and 0 <= fy <= 1 and cpos.contains(fx, fy):
                self.drag_state['target'] = 'colorbar'
                self.drag_state['offset'] = (fx - cpos.x0, fy - cpos.y0)
    
    def _on_motion(self, event):
        """Handle mouse motion for dragging."""
        if self.drag_state['target'] is None:
            return
        
        if self.drag_state['target'] == 'scalebar':
            if event.inaxes != self.ax or not event.xdata or not event.ydata:
                return
            offx, offy = self.drag_state['offset']
            self.elements['scalebar'].set_x(event.xdata - offx)
            self.elements['scalebar'].set_y(event.ydata - offy)
            if self.elements['scalebar_text']:
                self.elements['scalebar_text'].set_x(
                    event.xdata - offx + self.elements['scalebar'].get_width() / 2
                )
                self.elements['scalebar_text'].set_y(
                    event.ydata - offy + self.elements['scalebar'].get_height() + 1
                )
            self.fig.canvas.draw_idle()
        
        elif self.drag_state['target'] == 'image_label':
            if event.inaxes != self.ax or not event.xdata or not event.ydata:
                return
            offx, offy = self.drag_state['offset']
            self.elements['image_label'].set_x(event.xdata - offx)
            self.elements['image_label'].set_y(event.ydata - offy)
            self.fig.canvas.draw_idle()
        
        elif self.drag_state['target'] == 'colorbar':
            fx, fy = self.fig.transFigure.inverted().transform((event.x, event.y))
            offx, offy = self.drag_state['offset']
            cax = self.elements['colorbar_ax']
            new_x = max(0.0, min(fx - offx, 1.0 - cax.get_position().width))
            new_y = max(0.0, min(fy - offy, 1.0 - cax.get_position().height))
            cax.set_position([new_x, new_y, cax.get_position().width, cax.get_position().height])
            self.fig.canvas.draw_idle()
    
    def _on_release(self, event):
        """Handle mouse release."""
        self.drag_state['target'] = None
    
    def show_element(self, element_name, visible=True):
        """
        Show or hide a plot element.
        
        Args:
            element_name: One of 'image', 'colorbar', 'scalebar', 'image_label', 'save_button'
            visible: True to show, False to hide
        """
        element = self.elements.get(element_name)
        if element is None:
            print(f"Element '{element_name}' not found or not created yet.")
            return
        
        try:
            if element_name == 'colorbar':
                if self.elements['colorbar_ax']:
                    self.elements['colorbar_ax'].set_visible(visible)
            elif element_name == 'scalebar':
                element.set_visible(visible)
                if self.elements['scalebar_text']:
                    self.elements['scalebar_text'].set_visible(visible)
            elif element_name == 'save_button':
                if self.elements['save_button_ax']:
                    self.elements['save_button_ax'].set_visible(visible)
            else:
                element.set_visible(visible)
            
            if self.canvas:
                self.canvas.draw()
        except Exception as e:
            print(f"Error toggling element visibility: {e}")
    
    def _save_config(self):
        """Save current configuration to a JSON file."""
        import json
        
        filename = filedialog.asksaveasfilename(
            defaultextension='.json',
            filetypes=[('JSON files', '*.json'), ('All files', '*.*')],
            title='Save Configuration'
        )
        
        if not filename:
            return
        
        try:
            config = {
                # Colormap settings
                'cmap': self.cmap_var.get(),
                'vmin': self.vmin_var.get(),
                'vmax': self.vmax_var.get(),
                
                # Scalebar settings
                'scalebar_length': self.scalebarlength_var.get(),
                'scalebar_width': self.scalebarwidth_var.get(),
                'scalebar_pos_x': self.scalebarpos_x_var.get(),
                'scalebar_pos_y': self.scalebarpos_y_var.get(),
                'dx': self.dx_var.get(),
                'unit': self.unit_var.get(),
                
                # Figure settings
                'figsize_w': self.figsize_w_var.get(),
                'figsize_h': self.figsize_h_var.get(),
                'title': self.title_var.get(),
                'xlabel': self.xlabel_var.get(),
                'ylabel': self.ylabel_var.get(),
                'title_color': self.title_color_var.get(),
                
                # Font sizes
                'title_fontsize': self.title_fontsize_var.get(),
                'axis_label_fontsize': self.axis_label_fontsize_var.get(),
                'colorbar_fontsize': self.colorbar_fontsize_var.get(),
                'colorbar_caption_fontsize': self.colorbar_caption_fontsize_var.get(),
                'scalebar_fontsize': self.scalebar_fontsize_var.get(),
                'image_label_fontsize': self.image_label_fontsize_var.get(),
                
                # Colorbar settings
                'show_colorbar': self.show_colorbar_var.get(),
                'cbar_unit': self.cbar_unit_var.get(),
                'cbar_caption': self.cmapcaption_var.get(),
                'cbar_divisor': self.cbar_divisor_var.get(),
                'cbar_caption_color': self.cbar_caption_color_var.get(),
                'cbar_border_color': self.cbar_border_color_var.get(),
                'cbar_upper_num_color': self.cbar_upper_num_color_var.get(),
                'cbar_lower_num_color': self.cbar_lower_num_color_var.get(),
                
                # Scalebar colors
                'scalebar_face_color': self.scalebar_face_color_var.get(),
                'scalebar_edge_color': self.scalebar_edge_color_var.get(),
                'scalebar_text_color': self.scalebar_text_color_var.get(),
                
                # Label settings
                'image_label': self.image_label_var.get(),
                'image_label_color': self.image_label_color_var.get(),
                
                # Visibility toggles
                'show_scalebar': self.show_scalebar_var.get(),
                'show_image_label': self.show_image_label_var.get(),
                'enable_drag': self.enable_drag_var.get()
            }
            
            with open(filename, 'w') as f:
                json.dump(config, f, indent=4)
            
            print(f"Configuration saved to: {filename}")
            
        except Exception as e:
            print(f"Error saving configuration: {e}")
    
    def _load_config(self):
        """Load configuration from a JSON file."""
        import json
        
        filename = filedialog.askopenfilename(
            filetypes=[('JSON files', '*.json'), ('All files', '*.*')],
            title='Load Configuration'
        )
        
        if not filename:
            return
        
        try:
            with open(filename, 'r') as f:
                config = json.load(f)
            
            # Colormap settings
            if 'cmap' in config:
                self.cmap_var.set(config['cmap'])
            if 'vmin' in config:
                self.vmin_var.set(config['vmin'])
            if 'vmax' in config:
                self.vmax_var.set(config['vmax'])
            
            # Scalebar settings
            if 'scalebar_length' in config:
                self.scalebarlength_var.set(config['scalebar_length'])
            if 'scalebar_width' in config:
                self.scalebarwidth_var.set(config['scalebar_width'])
            if 'scalebar_pos_x' in config:
                self.scalebarpos_x_var.set(config['scalebar_pos_x'])
            if 'scalebar_pos_y' in config:
                self.scalebarpos_y_var.set(config['scalebar_pos_y'])
            if 'dx' in config:
                self.dx_var.set(config['dx'])
            if 'unit' in config:
                self.unit_var.set(config['unit'])
            
            # Figure settings
            if 'figsize_w' in config:
                self.figsize_w_var.set(config['figsize_w'])
            if 'figsize_h' in config:
                self.figsize_h_var.set(config['figsize_h'])
            if 'title' in config:
                self.title_var.set(config['title'])
            if 'xlabel' in config:
                self.xlabel_var.set(config['xlabel'])
            if 'ylabel' in config:
                self.ylabel_var.set(config['ylabel'])
            if 'title_color' in config:
                self.title_color_var.set(config['title_color'])
            
            # Font sizes
            if 'title_fontsize' in config:
                self.title_fontsize_var.set(config['title_fontsize'])
            if 'axis_label_fontsize' in config:
                self.axis_label_fontsize_var.set(config['axis_label_fontsize'])
            if 'colorbar_fontsize' in config:
                self.colorbar_fontsize_var.set(config['colorbar_fontsize'])
            if 'colorbar_caption_fontsize' in config:
                self.colorbar_caption_fontsize_var.set(config['colorbar_caption_fontsize'])
            if 'scalebar_fontsize' in config:
                self.scalebar_fontsize_var.set(config['scalebar_fontsize'])
            if 'image_label_fontsize' in config:
                self.image_label_fontsize_var.set(config['image_label_fontsize'])
            
            # Colorbar settings
            if 'show_colorbar' in config:
                self.show_colorbar_var.set(config['show_colorbar'])
            if 'cbar_unit' in config:
                self.cbar_unit_var.set(config['cbar_unit'])
            if 'cbar_caption' in config:
                self.cmapcaption_var.set(config['cbar_caption'])
            if 'cbar_divisor' in config:
                self.cbar_divisor_var.set(config['cbar_divisor'])
            if 'cbar_caption_color' in config:
                self.cbar_caption_color_var.set(config['cbar_caption_color'])
            if 'cbar_border_color' in config:
                self.cbar_border_color_var.set(config['cbar_border_color'])
            if 'cbar_upper_num_color' in config:
                self.cbar_upper_num_color_var.set(config['cbar_upper_num_color'])
            if 'cbar_lower_num_color' in config:
                self.cbar_lower_num_color_var.set(config['cbar_lower_num_color'])
            
            # Scalebar colors
            if 'scalebar_face_color' in config:
                self.scalebar_face_color_var.set(config['scalebar_face_color'])
            if 'scalebar_edge_color' in config:
                self.scalebar_edge_color_var.set(config['scalebar_edge_color'])
            if 'scalebar_text_color' in config:
                self.scalebar_text_color_var.set(config['scalebar_text_color'])
            
            # Label settings
            if 'image_label' in config:
                self.image_label_var.set(config['image_label'])
            if 'image_label_color' in config:
                self.image_label_color_var.set(config['image_label_color'])
            
            # Visibility toggles
            if 'show_scalebar' in config:
                self.show_scalebar_var.set(config['show_scalebar'])
            if 'show_image_label' in config:
                self.show_image_label_var.set(config['show_image_label'])
            if 'enable_drag' in config:
                self.enable_drag_var.set(config['enable_drag'])
            
            print(f"Configuration loaded from: {filename}")
            
        except Exception as e:
            print(f"Error loading configuration: {e}")


def plot_HSI(data, metadata=None, cmap='twilight', vmin=None, vmax=None,
             scalebarlength=20, scalebarwidth=2,
             scalebarpos=(50, 50), figsize=(7,6), title='',
             xlabel=None, ylabel=None, show_colorbar=True, cbar_unit='', scalebarfontsize=14,
             enable_drag=True, cmapcaption='', dx=1, image_label='', unit=''):
    """
    Plot a 2D HSI-like array as an image with an inset colorbar and a publication-style scalebar.

    Interactive features when enable_drag=True:
      - Drag-and-drop the scalebar (click+drag on the white bar) to reposition it.
      - Drag-and-drop the inset colorbar (click+drag inside the colorbar box) to reposition it.
      - A "Save 600 dpi" button that does not overlap the image and is hidden during saving.

    Returns:
        fig, ax: matplotlib Figure and Axes objects
    """

    data = np.array(data, dtype=float)
    if xlabel is None:
        xlabel = ''
    if ylabel is None:
        ylabel = ''

    fig, ax = plt.subplots(figsize=figsize)

    # If metadata contains wlstart/wlend, set extent so x-axis shows wavelengths
    extent = None
    wl0 = wl1 = None
    try:
        if metadata is not None and 'wlstart' in metadata and 'wlend' in metadata:
            wl0 = round(float(metadata['wlstart']), 1)
            wl1 = round(float(metadata['wlend']), 1)
            ncols = data.shape[1]
            extent = (wl0, wl1, 0, data.shape[0])
    except Exception:
        extent = None

    im = ax.imshow(data, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks([])
    ax.set_yticks([])

    # Inset colorbar inside the image (lower-right)
    cax = None
    cb = None
    if show_colorbar:
        # Scale colorbar size based on fontsize (base fontsize is 14)
        scale_factor = scalebarfontsize / 14.0
        ax_pos = ax.get_position()  # in figure coordinates
        cax_w = ax_pos.width * 0.08 * scale_factor  # reduced from 0.12 to 0.08 for narrower colorbar
        # Height scales only slightly with fontsize (reduced scaling factor)
        cax_h = ax_pos.height * 0.28 * (1.0 + 0.15 * (scale_factor - 1.0))
        cax_x = ax_pos.x0 + ax_pos.width * (1.0 - 0.08 * scale_factor - 0.10)  # adjusted for new width
        cax_y = ax_pos.y0 + ax_pos.height * 0.06
        cax = fig.add_axes([cax_x, cax_y, cax_w, cax_h])
        # keep a visible border/background for on-screen editing; we hide/change these when saving
        cax.patch.set_facecolor('white') # set to white
        cax.patch.set_edgecolor('white') 
        cax.patch.set_linewidth(2.6 * scale_factor)
        cb = fig.colorbar(im, cax=cax, orientation='vertical')
        for spine in cax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(0.8 * scale_factor)
        try:
            clim = im.get_clim()
            vminr = floor(clim[0])
            vmaxr = ceil(clim[1])
            if vminr == vmaxr:
                vminr = int(round(clim[0]))
                vmaxr = int(round(clim[1]))

            # Remove default ticks/labels and draw the two labels inside the colorbar axes
            try:
                cb.set_ticks([])
            except Exception:
                pass
            try:
                cax.yaxis.set_ticks([])
            except Exception:
                pass

            # Format colorbar label text
            if cbar_unit:
                top_label = f"{int(vmaxr)} {cbar_unit}" 
                bot_label = f"{int(vminr)} {cbar_unit}" 
            else:
                top_label = f"{int(vmaxr/1000)}" # scale down by 1000 for display
                bot_label = f"{int(vminr/1000)}" # scale down by 1000 for display

            # Choose a readable fontsize (fallback if scalebarfontsize not provided)
            try:
                labfs = max(8, int(scalebarfontsize * 0.8))
            except Exception:
                labfs = 10

            # Place labels inside the colorbar using axis-relative coords
            # top (near top of the cax) and bottom (near bottom)
            cax.text(0.5, 0.94, top_label, transform=cax.transAxes,
                     ha='center', va='top', fontsize=scalebarfontsize, color='black', weight='bold') # top text in black
            cax.text(0.5, 0.06, bot_label, transform=cax.transAxes,
                     ha='center', va='bottom', fontsize=scalebarfontsize, color='black', weight='bold') # bottom text in white
            # add a label above the colormap that says "Intensity" or similar if desired
            cax.text(0.5, 1.05, cmapcaption, transform=cax.transAxes,
                     ha='center', va='bottom', fontsize=scalebarfontsize, color='white', weight='bold')

            try:
                cb.outline.set_edgecolor('black')
                cb.outline.set_linewidth(0.8 * scale_factor)
            except Exception:
                pass
        except Exception:
            pass

    # Draw scalebar (only width scales with fontsize, length is bound to pixel size)
    scale_factor = scalebarfontsize / 14.0
    sbx, sby = scalebarpos
    # scalebarlength is NOT scaled - it represents actual pixels and is tied to physical size via dx
    scaled_width = scalebarwidth * scale_factor
    if extent is not None and (wl0 is not None and wl1 is not None):
        try:
            ncols = data.shape[1]
            x_per_pixel = (wl1 - wl0) / max(1, ncols)
            sb_width_display = scalebarlength * x_per_pixel  # use unscaled scalebarlength
            sb_height_display = scaled_width
            rect = Rectangle((sbx, sby), sb_width_display, sb_height_display,
                             facecolor='white', edgecolor='black', linewidth=1.1 * scale_factor, zorder=5)
        except Exception:
            rect = Rectangle((sbx, sby), scalebarlength, scaled_width,  # use unscaled scalebarlength
                             facecolor='white', edgecolor='black', linewidth=1.1 * scale_factor, zorder=5)
    else:
        rect = Rectangle((sbx, sby), scalebarlength, scaled_width,  # use unscaled scalebarlength
                         facecolor='white', edgecolor='black', linewidth=1.1 * scale_factor, zorder=5)
    ax.add_patch(rect)

    # label text
    try:
        phys_len = scalebarlength * float(dx)
        # Format with 2 decimal places, then remove trailing zeros
        formatted_len = f"{phys_len:.2f}"
        # Remove up to 2 trailing zeros after decimal point
        if '.' in formatted_len:
            formatted_len = formatted_len.rstrip('0').rstrip('.')
        label_text = f"{formatted_len} {unit}"
    except Exception:
        label_text = f"{scalebarlength} px"

    if extent is not None and (wl0 is not None and wl1 is not None):
        label_x = sbx + (sb_width_display / 2)
    else:
        label_x = sbx + (scalebarlength / 2)  # use unscaled scalebarlength

    text_artist = ax.text(label_x, sby + scaled_width + 1 * scale_factor, label_text,
            color='black', ha='center', va='bottom', fontsize=scalebarfontsize, weight='bold', zorder=6,
            #bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', linewidth=0.6))
    )

    # Image label (draggable) shown near top-right of the axes and saved with the figure.
    try:
        img_label_text = image_label  # uses module-level variable if present
    except Exception:
        img_label_text = ''
    try:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xpos = xlim[1] - 0.02 * (xlim[1] - xlim[0])
        ypos = ylim[1] - 0.02 * (ylim[1] - ylim[0])
    except Exception:
        xpos, ypos = 0.98, 0.98
    img_label_artist = ax.text(xpos, ypos, img_label_text,
                               color='black', ha='right', va='top',
                               fontsize=scalebarfontsize, weight='bold', zorder=6)

    # Storage for save button (will be created after tight_layout)
    button_data = {'bax': None, 'button': None}

    def on_save(event=None):
        # Hide the save button so it is not included in the saved image
        save_outline = None
        bax = button_data['bax']
        try:
            if bax is not None:
                bax.set_visible(False)
            # For saving: remove the filled white patch that can create a filled block/stripe
            # but add a temporary outline rectangle (no fill) in figure coords so the saved
            # image retains a white border around the colorbar without adding a filled white area.
            if cax is not None:
                orig_face = cax.patch.get_facecolor()
                orig_edge = cax.patch.get_edgecolor()
                orig_linewidth = cax.patch.get_linewidth()
                orig_spines = {k: v.get_visible() for k, v in cax.spines.items()}

                # hide the patch and spines to avoid a filled rectangle being saved over the image
                cax.patch.set_facecolor('none')
                cax.patch.set_edgecolor('none')
                cax.patch.set_linewidth(0)
                for spine in cax.spines.values():
                    spine.set_visible(False)

                # add a temporary outline rectangle in figure coordinates (no fill)
                try:
                    pos = cax.get_position()
                    save_outline = Rectangle((pos.x0, pos.y0), pos.width, pos.height,
                                             transform=fig.transFigure, fill=False,
                                             edgecolor='white', linewidth=max(1.5, orig_linewidth), zorder=10, clip_on=False)
                    fig.add_artist(save_outline)
                except Exception:
                    save_outline = None
            # redraw
            fig.canvas.draw()
        except Exception:
            pass

        # ask filename and save
        root = tk.Tk()
        root.withdraw()
        fname = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[('PNG', '*.png')])
        root.destroy()
        if fname:
            # Save without extra white border
            fig.savefig(fname, dpi=600, bbox_inches='tight', pad_inches=0)
            print(f'Saved: {fname}')

        # restore visibility and colorbar styling
        try:
            if bax is not None:
                bax.set_visible(True)
            if cax is not None:
                # restore patch and spines
                cax.patch.set_facecolor(orig_face)
                cax.patch.set_edgecolor(orig_edge)
                cax.patch.set_linewidth(orig_linewidth)
                for k, v in cax.spines.items():
                    v.set_visible(orig_spines[k])
                # remove temporary outline if present
                try:
                    if save_outline is not None:
                        save_outline.remove()
                except Exception:
                    pass
            fig.canvas.draw()
        except Exception:
            pass

    # Interactive drag handlers
    if enable_drag:
        drag_state = {'target': None, 'offset': (0, 0)}

        def on_press(event):
            if event.inaxes == ax:
                contains_rect, _ = rect.contains(event)
                if contains_rect:
                    drag_state['target'] = 'rect'
                    if event.xdata is not None and event.ydata is not None:
                        drag_state['offset'] = (event.xdata - rect.get_x(), event.ydata - rect.get_y())
                    return
                # check if user clicked the image label
                try:
                    contains_label, _ = img_label_artist.contains(event)
                except Exception:
                    contains_label = False
                if contains_label:
                    drag_state['target'] = 'img_label'
                    if event.xdata is not None and event.ydata is not None:
                        drag_state['offset'] = (event.xdata - img_label_artist.get_position()[0],
                                                event.ydata - img_label_artist.get_position()[1])
                    return
            if cax is not None:
                fx, fy = fig.transFigure.inverted().transform((event.x, event.y))
                cpos = cax.get_position()
                if 0 <= fx <= 1 and 0 <= fy <= 1 and cpos.contains(fx, fy):
                    drag_state['target'] = 'cax'
                    drag_state['offset'] = (fx - cpos.x0, fy - cpos.y0)
                    return

        def on_motion(event):
            if drag_state['target'] is None:
                return
            if drag_state['target'] == 'rect':
                if event.inaxes != ax:
                    return
                if event.xdata is None or event.ydata is None:
                    return
                offx, offy = drag_state['offset']
                new_x = event.xdata - offx
                new_y = event.ydata - offy
                rect.set_x(new_x)
                rect.set_y(new_y)
                # also move label
                text_artist.set_x(new_x + rect.get_width() / 2)
                text_artist.set_y(new_y + rect.get_height() + 1)
                fig.canvas.draw_idle()
            elif drag_state['target'] == 'img_label':
                if event.inaxes != ax:
                    return
                if event.xdata is None or event.ydata is None:
                    return
                offx, offy = drag_state['offset']
                new_x = event.xdata - offx
                new_y = event.ydata - offy
                try:
                    img_label_artist.set_x(new_x)
                    img_label_artist.set_y(new_y)
                except Exception:
                    pass
                fig.canvas.draw_idle()
            elif drag_state['target'] == 'cax':
                fx, fy = fig.transFigure.inverted().transform((event.x, event.y))
                offx, offy = drag_state['offset']
                new_x = fx - offx
                new_y = fy - offy
                new_x = max(0.0, min(new_x, 1.0 - cax.get_position().width))
                new_y = max(0.0, min(new_y, 1.0 - cax.get_position().height))
                cax.set_position([new_x, new_y, cax.get_position().width, cax.get_position().height])
                fig.canvas.draw_idle()

        def on_release(event):
            drag_state['target'] = None

        fig.canvas.mpl_connect('button_press_event', on_press)
        fig.canvas.mpl_connect('motion_notify_event', on_motion)
        fig.canvas.mpl_connect('button_release_event', on_release)

    # Apply tight_layout first, then create the save button so it doesn't get repositioned
    plt.tight_layout()
    
    # Create the Save button after tight_layout in the bottom-left corner
    # The button will be hidden during saving so it doesn't appear in the saved image
    try:
        # Position button in bottom-left corner of figure (inside the figure area)
        button_data['bax'] = fig.add_axes([0.01, 0.01, 0.13, 0.04])
        button_data['button'] = Button(button_data['bax'], 'Save 600 dpi')
        button_data['button'].on_clicked(on_save)
    except Exception as e:
        print(f"Warning: Could not create save button: {e}")
    
    plt.show()
    return fig, ax


# Helper function to load HSI data
def get_HSI(file_path):
    """
    Reads a spectral file with metadata header and semicolon-separated data.
    
    Returns:
        tuple: (metadata (dict), data (np.ndarray))
    """
    metadata = {}
    data_rows = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if ':' in line and ';' not in line:
            k, v = line.split(':', 1)
            k = k.strip()
            v = v.strip()
            try:
                if re.match(r'^-?\d+$', v):
                    metadata[k] = int(v)
                else:
                    metadata[k] = float(v)
            except Exception:
                metadata[k] = v
        elif ';' in line:
            parts = [p.strip() for p in line.split(';') if p.strip()!='']
            row = []
            for p in parts:
                try:
                    row.append(float(p))
                except Exception:
                    row.append(np.nan)
            data_rows.append(row)

    if len(data_rows) == 0:
        data_arr = np.empty((0, 0), dtype=float)
    else:
        max_len = max(len(r) for r in data_rows)
        padded = [r + [np.nan] * (max_len - len(r)) for r in data_rows]
        data_arr = np.array(padded, dtype=float)

    return metadata, data_arr

def hsi_asarray(array2d):
    """Convert a 2D array to a NumPy array of floats."""
    return np.array(array2d, dtype=float)

def main():
    """Main function to run the PlotHSI GUI."""
    root = tk.Tk()
    app = PlotHSI(root, get_hsi_func=get_HSI)
    root.mainloop()


if __name__ == "__main__":
    main()
    #print("os.getcwd():", os.getcwd())
