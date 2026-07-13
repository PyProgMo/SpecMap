# Visualization of the Select Sepctral Data as spectra
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import ttk
try:
    from . import deflib1 as deflib
except ImportError:
    import deflib1 as deflib

import traceback


class Specplottergui:
    # Fallback defaults in case file reading fails
    FALLBACK_DEFAULTS = {
        'x_min': 'Auto',
        'x_max': 'Auto',
        'y_min': 'Auto',
        'y_max': 'Auto',
        'x_tick': '50',
        'y_tick': '100',
        'label_font_size': '24',
        'tick_font_size': '22',
        'legend_font_size': '22',
        'font_family': 'Arial',
        'x_label': 'Wavelength (nm)',
        'y_label': 'norm. Intensity',
        'title': 'Spectra Plot',
        'show_grid': 'True',
        'normalize_point': 'False',
        'normalize_max': 'False',
        'use_nm_to_ev': 'False',
        'save_path': 'C:\\Users\\volib\\Desktop',
        'file_name': 'specs1'
    }
    
    def __init__(self, Specdata={}, Specdiffsets={}, guiroot=None, disspecs=None, defaultsfile='defaults_specplotter.txt'):
        self.Specdata = Specdata
        self.Specdiffsets = Specdiffsets
        self.guiroot = guiroot
        self.defaultsfile = defaultsfile
        self.defaults = self._load_defaults()
        
        # clear guiroot before building to avoid duplicates when reloading data
        if self.guiroot is not None:
            for widget in self.guiroot.winfo_children():
                widget.destroy()
        self.disspecs = disspecs

        if self.guiroot is not None:
            self.buildgui()
    
    def _load_defaults(self):
        """
        Load defaults from file, with fallback to FALLBACK_DEFAULTS.
        """
        defaults = self.FALLBACK_DEFAULTS.copy()
        try:
            with open(self.defaultsfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip comments and empty lines
                    if not line or line.startswith('#'):
                        continue
                    # Parse key = value
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        # Convert boolean strings to actual values
                        if value.lower() == 'true':
                            value = 'True'
                        elif value.lower() == 'false':
                            value = 'False'
                        defaults[key] = value
        except FileNotFoundError:
            print(f"Defaults file '{self.defaultsfile}' not found. Using fallback defaults.")
        except Exception as e:
            print(f"Error loading defaults from '{self.defaultsfile}': {e}. Using fallback defaults.")
        
        return defaults

    def save_defaults(self):
        """
        Save current plot options to the defaults file.
        """
        try:
            with open(self.defaultsfile, 'w') as f:
                f.write('# Default settings for Spectra Plotter GUI\n')
                f.write('# Format: key = value\n\n')
                
                # Axes Limits
                f.write('# Axes Limits\n')
                f.write(f"x_min = {self.opt_vars['x_min'].get()}\n")
                f.write(f"x_max = {self.opt_vars['x_max'].get()}\n")
                f.write(f"y_min = {self.opt_vars['y_min'].get()}\n")
                f.write(f"y_max = {self.opt_vars['y_max'].get()}\n\n")
                
                # Tick spacing
                f.write('# Tick spacing\n')
                f.write(f"x_tick = {self.opt_vars['x_tick'].get()}\n")
                f.write(f"y_tick = {self.opt_vars['y_tick'].get()}\n\n")
                
                # Font sizes
                f.write('# Font sizes\n')
                f.write(f"label_font_size = {self.opt_vars['label_font_size'].get()}\n")
                f.write(f"tick_font_size = {self.opt_vars['tick_font_size'].get()}\n")
                f.write(f"legend_font_size = {self.opt_vars['legend_font_size'].get()}\n")
                f.write(f"font_family = {self.opt_vars['font_family'].get()}\n\n")
                
                # Axis labels
                f.write('# Axis labels\n')
                f.write(f"x_label = {self.opt_vars['x_label'].get()}\n")
                f.write(f"y_label = {self.opt_vars['y_label'].get()}\n")
                f.write(f"title = {self.opt_vars['title'].get()}\n\n")
                
                # Display options
                f.write('# Display options\n')
                f.write(f"show_grid = {self.opt_vars['show_grid'].get()}\n")
                f.write(f"normalize_point = {self.opt_vars['normalize_point'].get()}\n")
                f.write(f"normalize_max = {self.opt_vars['normalize_max'].get()}\n")
                f.write(f"use_nm_to_ev = {self.opt_vars['use_nm_to_ev'].get()}\n\n")
                
                # File settings
                f.write('# File settings\n')
                f.write(f"save_path = {self.opt_vars['save_path'].get()}\n")
                f.write(f"file_name = {self.opt_vars['file_name'].get()}\n")
            
            print(f"Defaults saved to '{self.defaultsfile}'")
        except Exception as e:
            print(f"Error saving defaults to '{self.defaultsfile}': {e}")

    def buildgui(self):
        # build 2 tabs, one to select what spectral data to plot, and one to select the plot type
        self.tab_control = ttk.Notebook(self.guiroot)
        self.tab1 = ttk.Frame(self.tab_control)
        self.tab2 = ttk.Frame(self.tab_control)
        
        self.tab_control.add(self.tab1, text='Data Selection')
        self.tab_control.add(self.tab2, text='Plot Options')
        self.tab_control.pack(expand=1, fill="both")

        self.build_data_selection_tab()
        self.build_plot_options_tab()

    def build_data_selection_tab(self):
        # Matrix to select spectra: HSI names (rows) x Datatypes (columns)
        # Using short names to avoid excessive column width
        self.specdatatypes = [
            'PL-BG', 'D1', 'D2', 'D1-N', 'D2-N', 'D1-NI', 'D2-NI', 'D1-NC', 'D2-NC'
        ]
        
        self.dataset_names = list(self.disspecs.keys()) if self.disspecs else []

        if not self.dataset_names:
            self.dataset_names = ['Dataset_1', 'Dataset_2']

        self.plot_vars = {}
        self.col_headers = {}  # Store column header references
        self.row_headers = {}  # Store row header references

        # Button frame for Select All / Deselect All / Transfer
        btn_frame = ttk.Frame(self.tab1)
        btn_frame.pack(fill="x", padx=10, pady=10)
        ttk.Button(btn_frame, text="Select All", command=self.select_all_specs).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Deselect All", command=self.deselect_all_specs).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Transfer to Plot Options", command=self.transfer_to_plot_options).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Generate Plot", command=self.generate_plot).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Refresh Data", command=self.refresh_data).pack(side="left", padx=5)

        # Create a separate frame for the grid-based matrix (cannot mix pack and grid on same frame)
        self.matrix_frame = ttk.Frame(self.tab1)
        self.matrix_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Create headers
        ttk.Label(self.matrix_frame, text="Datasets", font=('Arial', 10, 'bold')).grid(row=0, column=0, padx=5, pady=5)
        for col, dtype in enumerate(self.specdatatypes):
            col_label = ttk.Label(self.matrix_frame, text=dtype.upper(), font=('Arial', 10, 'bold'), foreground='blue')
            col_label.grid(row=0, column=col+1, padx=5, pady=5)
            col_label.bind('<Button-1>', lambda e, c=col: self.toggle_column(c))
            col_label.config(cursor='hand2')
            self.col_headers[col] = col_label
        
        self.checkbox_widgets = []
        self.build_plot_checkbox()
    
    def build_plot_checkbox(self):
        # Destroy existing widgets before building new ones
        self.destroy_plot_checkboxes()
        
        for row, ds_name in enumerate(self.dataset_names):
            label = ttk.Label(self.matrix_frame, text=ds_name, foreground='blue')
            label.grid(row=row+1, column=0, sticky='w', padx=5, pady=2)
            label.bind('<Button-1>', lambda e, r=row: self.toggle_row(r))
            label.config(cursor='hand2')
            self.row_headers[row] = label
            self.checkbox_widgets.append(label)
            
            self.plot_vars[ds_name] = {}
            for col, dtype in enumerate(self.specdatatypes):
                var = tk.BooleanVar()
                chk = ttk.Checkbutton(self.matrix_frame, variable=var)
                chk.grid(row=row+1, column=col+1, padx=5, pady=2)
                self.plot_vars[ds_name][dtype] = var
                self.checkbox_widgets.append(chk)

    def destroy_plot_checkboxes(self):
        """Destroy all checkbox and label widgets."""
        for widget in self.checkbox_widgets:
            widget.destroy()
        self.checkbox_widgets = []

        
    def refresh_data(self):
        """Refresh dataset names from self.disspecs and rebuild the checkbox matrix."""
        if not self.disspecs:
            print("No spectral data available (disspecs is empty)")
            # even if empty, we should clear the old view
            self.dataset_names = []
            self.plot_vars = {}
            self.build_plot_checkbox()
            return
        
        # Update dataset names and reset plot variables
        self.dataset_names = list(self.disspecs.keys())
        self.plot_vars = {}
        
        # Rebuild the checkboxes
        self.build_plot_checkbox()
        
    def select_all_specs(self):
        """Set all checkboxes to True"""
        for ds in self.plot_vars:
            for dtype in self.plot_vars[ds]:
                self.plot_vars[ds][dtype].set(True)

    def deselect_all_specs(self):
        """Set all checkboxes to False"""
        for ds in self.plot_vars:
            for dtype in self.plot_vars[ds]:
                self.plot_vars[ds][dtype].set(False)

    def toggle_column(self, col_idx):
        """Toggle selection of all checkboxes in a column"""
        # Determine the current state of the column
        # If any checkbox is unchecked, check all; otherwise uncheck all
        column_values = []
        for ds_name in self.dataset_names:
            dtype = self.specdatatypes[col_idx]
            if dtype in self.plot_vars[ds_name]:
                column_values.append(self.plot_vars[ds_name][dtype].get())
        
        # Determine new state: if all are True, set to False; otherwise set to True
        new_state = not all(column_values) if column_values else True
        
        # Apply new state to all checkboxes in this column
        for ds_name in self.dataset_names:
            dtype = self.specdatatypes[col_idx]
            if dtype in self.plot_vars[ds_name]:
                self.plot_vars[ds_name][dtype].set(new_state)

    def toggle_row(self, row_idx):
        """Toggle selection of all checkboxes in a row"""
        if row_idx >= len(self.dataset_names):
            return
        
        ds_name = self.dataset_names[row_idx]
        
        # Determine the current state of the row
        row_values = [self.plot_vars[ds_name][dtype].get() for dtype in self.specdatatypes if dtype in self.plot_vars[ds_name]]
        
        # Determine new state: if all are True, set to False; otherwise set to True
        new_state = not all(row_values) if row_values else True
        
        # Apply new state to all checkboxes in this row
        for dtype in self.specdatatypes:
            if dtype in self.plot_vars[ds_name]:
                self.plot_vars[ds_name][dtype].set(new_state)

    def transfer_to_plot_options(self):
        """
        Transfer selected specs to Plot Options tab for configuration.
        """
        if not self.disspecs:
            print("Error: No spectral data available")
            return

        # Collect selected spectra
        self.selected_plots = []
        for ds_name, types_dict in self.plot_vars.items():
            for dtype, is_selected in types_dict.items():
                if is_selected.get():
                    if ds_name in self.disspecs:
                        spec_obj = self.disspecs[ds_name]
                        data = self._get_spectra_array(spec_obj, dtype)
                        if data is not None:
                            self.selected_plots.append({
                                'dataset': ds_name,
                                'datatype': dtype,
                                'wl': spec_obj.WL,
                                'data': data,
                                'spec_obj': spec_obj
                            })
        
        if not self.selected_plots:
            print("Error: No spectra selected")
            return

        # Dynamically build plot options in the second tab
        self._build_dynamic_plot_options()
        
        # Switch to the Plot Options tab
        self.tab_control.select(self.tab2)
        
        print(f"Transferred {len(self.selected_plots)} spectra to Plot Options.")

    def generate_plot(self):
        """
        Generates and displays the plot based on selections and configurations.
        """
        if not hasattr(self, 'selected_plots') or not self.selected_plots:
            print("Error: No plots selected to generate. Please select data and transfer first.")
            return

        # Create figure and plot
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Define color palette (fallback)
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta', 'yellow', 'teal', 'navy', 'maroon', 'olive', 'lime', 'coral', 'indigo', 'gold', 'silver']
        
        # Get normalization and conversion options
        normalize_by_point = self.opt_vars['normalize_point'].get()
        normalize_by_max = self.opt_vars['normalize_max'].get()
        use_nm_to_ev = self.opt_vars['use_nm_to_ev'].get()
        
        # Process each spectrum with normalization
        processed_plots = []
        for plot_data in self.selected_plots:
            wl = np.array(plot_data['wl'])
            data = np.array(plot_data['data'])
            
            # Apply normalization by point
            if normalize_by_point:
                # Find max point (or could be configurable)
                max_idx = np.argmax(data)
                norm_value = data[max_idx]
                if norm_value != 0:
                    data = data / norm_value
            
            # Apply normalization by max
            elif normalize_by_max:
                max_val = np.nanmax(np.abs(data))
                if max_val != 0:
                    data = data / max_val
            
            # Convert wavelength axis if needed
            if use_nm_to_ev:
                # Convert nm to eV using E = 1240 / lambda
                wl = 1240.0 / wl
                # Reverse the order since eV increases as nm decreases
                wl = wl[::-1]
                data = data[::-1]
            
            processed_plots.append({
                'wl': wl,
                'data': data,
                'dataset': plot_data['dataset'],
                'datatype': plot_data['datatype']
            })
        
        # Plot each processed spectrum with configured options
        
        # Sort plots by their configured order
        if hasattr(self, 'series_options'):
            for idx, plot_data in enumerate(processed_plots):
                if idx in self.series_options:
                    try:
                        order_val = int(self.series_options[idx]['order'].get())
                    except:
                        order_val = idx + 1
                    plot_data['order'] = order_val
                else:
                    plot_data['order'] = idx + 1
                    
            processed_plots.sort(key=lambda x: x['order'])
            
            # Create a lookup mapping from sorted position back to original idx
            # so we can fetch the correct styles
            # Wait, sorting processed_plots changes their indices but they correspond to original indices?
            # Let's keep the original index inside the dict
            for idx, plot_data in enumerate(self.selected_plots):
                # add original index
                pass
                
        # Re-populate processed plots with original indices
        processed_plots = []
        for orig_idx, plot_data in enumerate(self.selected_plots):
            wl = np.array(plot_data['wl'])
            data = np.array(plot_data['data'])
            
            # Apply normalization by point
            if normalize_by_point:
                # Find max point (or could be configurable)
                max_idx = np.argmax(data)
                norm_value = data[max_idx]
                if norm_value != 0:
                    data = data / norm_value
            
            # Apply normalization by max
            elif normalize_by_max:
                max_val = np.nanmax(np.abs(data))
                if max_val != 0:
                    data = data / max_val
            
            # Convert wavelength axis if needed
            if use_nm_to_ev:
                # Convert nm to eV using E = 1240 / lambda
                wl = 1240.0 / wl
                # Reverse the order since eV increases as nm decreases
                wl = wl[::-1]
                data = data[::-1]
            
            order_val = orig_idx + 1
            if hasattr(self, 'series_options') and orig_idx in self.series_options:
                try:
                    order_val = int(self.series_options[orig_idx]['order'].get())
                except:
                    pass

            processed_plots.append({
                'orig_idx': orig_idx,
                'order': order_val,
                'wl': wl,
                'data': data,
                'dataset': plot_data['dataset'],
                'datatype': plot_data['datatype']
            })
            
        processed_plots.sort(key=lambda x: x['order'])

        for plot_data in processed_plots:
            orig_idx = plot_data['orig_idx']
            # Get custom options for this series
            if hasattr(self, 'series_options') and orig_idx in self.series_options:
                options = self.series_options[orig_idx]
                try:
                    line_thickness = float(options['line_thickness'].get())
                except:
                    line_thickness = 1.5
                
                try:
                    color = options['line_color'].get()
                except:
                    color = colors[idx % len(colors)]
                
                try:
                    label = options['legend_label'].get()
                except:
                    label = f"{plot_data['dataset']} - {plot_data['datatype']}"
                
                line_style_str = options['line_style'].get().lower()
                linestyle_map = {'solid': '-', 'dashed': '--', 'dotted': ':'}
                linestyle = linestyle_map.get(line_style_str, '-')
                
                is_scatter = options['scatter_graph'].get()
                
                try:
                    point_size = float(options['point_size'].get())
                except:
                    point_size = 15
            else:
                # Fallback to defaults if series_options not available
                line_thickness = 1.5
                color = colors[idx % len(colors)]
                label = f"{plot_data['dataset']} - {plot_data['datatype']}"
                linestyle = '-'
                is_scatter = False
                point_size = 15
            
            # Plot either as scatter or line
            if is_scatter:
                ax.scatter(plot_data['wl'], plot_data['data'], color=color, s=point_size, label=label, alpha=0.7)
            else:
                ax.plot(plot_data['wl'], plot_data['data'], color=color, linewidth=line_thickness, 
                       linestyle=linestyle, label=label)
        
        # Apply plot options
        try:
            # X and Y limits
            x_min = self.opt_vars['x_min'].get()
            x_max = self.opt_vars['x_max'].get()
            y_min = self.opt_vars['y_min'].get()
            y_max = self.opt_vars['y_max'].get()
            
            if x_min != 'Auto':
                x_min = float(x_min)
            else:
                x_min = None
            
            if x_max != 'Auto':
                x_max = float(x_max)
            else:
                x_max = None
            
            if y_min != 'Auto':
                y_min = float(y_min)
            else:
                y_min = None
            
            if y_max != 'Auto':
                y_max = float(y_max)
            else:
                y_max = None
            
            if x_min is not None or x_max is not None:
                ax.set_xlim(x_min, x_max)
            
            if y_min is not None or y_max is not None:
                ax.set_ylim(y_min, y_max)
            
            # Ticks
            try:
                x_tick = float(self.opt_vars['x_tick'].get())
                y_tick = float(self.opt_vars['y_tick'].get())
                ax.xaxis.set_major_locator(plt.MultipleLocator(x_tick))
                ax.yaxis.set_major_locator(plt.MultipleLocator(y_tick))
            except:
                pass
            
            # Labels and title
            x_label = self.opt_vars['x_label'].get()
            y_label = self.opt_vars['y_label'].get()
            title = self.opt_vars['title'].get()
            
            # Update x_label if converting to eV
            if use_nm_to_ev and 'nm' in x_label.lower():
                x_label = x_label.replace('nm', 'eV').replace('(nm)', '(eV)').replace('Nm', 'eV')
            elif use_nm_to_ev and 'eV' not in x_label.lower():
                x_label = x_label + ' (eV)'
            
            # Update y_label if normalizing
            if normalize_by_point or normalize_by_max:
                if 'norm' not in y_label.lower():
                    y_label = 'Normalized ' + y_label
            
            # Font sizes
            label_font_size = int(self.opt_vars['label_font_size'].get())
            tick_font_size = int(self.opt_vars['tick_font_size'].get())
            legend_font_size = int(self.opt_vars['legend_font_size'].get())
            
            ax.set_xlabel(x_label, fontsize=label_font_size)
            ax.set_ylabel(y_label, fontsize=label_font_size)
            ax.set_title(title, fontsize=label_font_size)
            
            ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
            
            # Grid
            if self.opt_vars['show_grid'].get():
                ax.grid(True, alpha=0.3)
            
            # Legend
            ax.legend(fontsize=legend_font_size)
            
            plt.tight_layout()
            plt.show()
                        
        except Exception as e:
            print(f"Error creating plot: {e}")
            traceback.print_exc()
    
    
    
    def _get_spectra_array(self, spec_obj, datatype):
        """
        Extract the appropriate spectral array from a Spectra object based on datatype.
        
        Datatypes:
        - PL-BG: Main spectrum
        - D1: First derivative
        - D2: Second derivative
        - D1-N: First derivative normalized by max
        - D2-N: Second derivative normalized by max
        - D1-NI: First derivative normalized inverted
        - D2-NI: Second derivative normalized inverted
        - D1-NC: First derivative normalized clipped
        - D2-NC: Second derivative normalized clipped
        """
        try:
            if datatype == 'PL-BG':
                return np.array(spec_obj.Spec)
            
            elif datatype == 'D1':
                if spec_obj.Spec_d1 is not None:
                    return np.array(spec_obj.Spec_d1)
                else:
                    print(f"Warning: D1 not available for {spec_obj.parenthsi}")
                    return None
            
            elif datatype == 'D2':
                if spec_obj.Spec_d2 is not None:
                    return np.array(spec_obj.Spec_d2)
                else:
                    print(f"Warning: D2 not available for {spec_obj.parenthsi}")
                    return None
            
            elif datatype == 'D1-N':
                if spec_obj.Spec_d1 is not None:
                    arr = np.array(spec_obj.Spec_d1)
                    max_val = np.nanmax(np.abs(arr))
                    return arr / max_val if max_val != 0 else arr
                else:
                    return None
            
            elif datatype == 'D2-N':
                if spec_obj.Spec_d2 is not None:
                    arr = np.array(spec_obj.Spec_d2)
                    max_val = np.nanmax(np.abs(arr))
                    return arr / max_val if max_val != 0 else arr
                else:
                    return None
            
            elif datatype == 'D1-NI':
                if spec_obj.Spec_d1 is not None:
                    arr = np.array(spec_obj.Spec_d1)
                    max_val = np.nanmax(np.abs(arr))
                    normalized = arr / max_val if max_val != 0 else arr
                    return -normalized
                else:
                    return None
            
            elif datatype == 'D2-NI':
                if spec_obj.Spec_d2 is not None:
                    arr = np.array(spec_obj.Spec_d2)
                    max_val = np.nanmax(np.abs(arr))
                    normalized = arr / max_val if max_val != 0 else arr
                    return -normalized
                else:
                    return None
            
            elif datatype == 'D1-NC':
                if spec_obj.Spec_d1 is not None:
                    arr = np.array(spec_obj.Spec_d1)
                    max_val = np.nanmax(np.abs(arr))
                    normalized = arr / max_val if max_val != 0 else arr
                    return np.clip(normalized, -1, 1)
                else:
                    return None
            
            elif datatype == 'D2-NC':
                if spec_obj.Spec_d2 is not None:
                    arr = np.array(spec_obj.Spec_d2)
                    max_val = np.nanmax(np.abs(arr))
                    normalized = arr / max_val if max_val != 0 else arr
                    return np.clip(normalized, -1, 1)
                else:
                    return None
            
            else:
                print(f"Unknown datatype: {datatype}")
                return None
        
        except Exception as e:
            print(f"Error extracting {datatype} from {spec_obj.parenthsi}: {e}")
            return None

    def build_plot_options_tab(self):
        # Variables for plot options, using loaded defaults
        self.opt_vars = {
            'x_min': tk.StringVar(value=self.defaults.get('x_min', 'Auto')),
            'x_max': tk.StringVar(value=self.defaults.get('x_max', 'Auto')),
            'y_min': tk.StringVar(value=self.defaults.get('y_min', 'Auto')),
            'y_max': tk.StringVar(value=self.defaults.get('y_max', 'Auto')),
            'x_tick': tk.StringVar(value=self.defaults.get('x_tick', '50')),
            'y_tick': tk.StringVar(value=self.defaults.get('y_tick', '100')),
            'label_font_size': tk.StringVar(value=self.defaults.get('label_font_size', '24')),
            'tick_font_size': tk.StringVar(value=self.defaults.get('tick_font_size', '22')),
            'legend_font_size': tk.StringVar(value=self.defaults.get('legend_font_size', '22')),
            'font_family': tk.StringVar(value=self.defaults.get('font_family', 'Arial')),
            'x_label': tk.StringVar(value=self.defaults.get('x_label', 'Wavelength (nm)')),
            'y_label': tk.StringVar(value=self.defaults.get('y_label', 'norm. Intensity')),
            'title': tk.StringVar(value=self.defaults.get('title', 'Spectra Plot')),
            'show_grid': tk.BooleanVar(value=self.defaults.get('show_grid', 'True').lower() == 'true'),
            'normalize_point': tk.BooleanVar(value=self.defaults.get('normalize_point', 'False').lower() == 'true'),
            'normalize_max': tk.BooleanVar(value=self.defaults.get('normalize_max', 'False').lower() == 'true'),
            'use_nm_to_ev': tk.BooleanVar(value=self.defaults.get('use_nm_to_ev', 'False').lower() == 'true'),
            'save_path': tk.StringVar(value=self.defaults.get('save_path', 'C:\\Users\\volib\\Desktop')),
            'file_name': tk.StringVar(value=self.defaults.get('file_name', 'specs1'))
        }

        # Create canvas with scrollbar for the options
        canvas = tk.Canvas(self.tab2)
        scrollbar = ttk.Scrollbar(self.tab2, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas)
        
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        scrollbar.pack(side="right", fill="y")

        # Axes Limits frame
        limits_frame = ttk.LabelFrame(self.scrollable_frame, text="Axes Limits")
        limits_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(limits_frame, text="X Min:").grid(row=0, column=0, padx=3, pady=3, sticky='e')
        ttk.Entry(limits_frame, textvariable=self.opt_vars['x_min'], width=10).grid(row=0, column=1, padx=3, pady=3)
        ttk.Button(limits_frame, text="Auto", width=6).grid(row=0, column=2, padx=3, pady=3)
        
        ttk.Label(limits_frame, text="X Max:").grid(row=0, column=3, padx=3, pady=3, sticky='e')
        ttk.Entry(limits_frame, textvariable=self.opt_vars['x_max'], width=10).grid(row=0, column=4, padx=3, pady=3)
        ttk.Button(limits_frame, text="Auto", width=6).grid(row=0, column=5, padx=3, pady=3)
        
        ttk.Label(limits_frame, text="Y Min:").grid(row=0, column=6, padx=3, pady=3, sticky='e')
        ttk.Entry(limits_frame, textvariable=self.opt_vars['y_min'], width=10).grid(row=0, column=7, padx=3, pady=3)
        ttk.Button(limits_frame, text="Auto", width=6).grid(row=0, column=8, padx=3, pady=3)
        
        ttk.Label(limits_frame, text="Y Max:").grid(row=0, column=9, padx=3, pady=3, sticky='e')
        ttk.Entry(limits_frame, textvariable=self.opt_vars['y_max'], width=10).grid(row=0, column=10, padx=3, pady=3)
        ttk.Button(limits_frame, text="Auto", width=6).grid(row=0, column=11, padx=3, pady=3)

        # Tick and Font frame
        format_frame = ttk.LabelFrame(self.scrollable_frame, text="Formatting")
        format_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(format_frame, text="X Tick:").grid(row=0, column=0, padx=3, pady=3, sticky='e')
        ttk.Entry(format_frame, textvariable=self.opt_vars['x_tick'], width=10).grid(row=0, column=1, padx=3, pady=3)
        
        ttk.Label(format_frame, text="Y Tick:").grid(row=0, column=2, padx=3, pady=3, sticky='e')
        ttk.Entry(format_frame, textvariable=self.opt_vars['y_tick'], width=10).grid(row=0, column=3, padx=3, pady=3)
        
        ttk.Label(format_frame, text="Label Font Size:").grid(row=0, column=4, padx=3, pady=3, sticky='e')
        ttk.Entry(format_frame, textvariable=self.opt_vars['label_font_size'], width=8).grid(row=0, column=5, padx=3, pady=3)
        
        ttk.Label(format_frame, text="Tick Font Size:").grid(row=0, column=6, padx=3, pady=3, sticky='e')
        ttk.Entry(format_frame, textvariable=self.opt_vars['tick_font_size'], width=8).grid(row=0, column=7, padx=3, pady=3)
        
        ttk.Label(format_frame, text="Legend Font Size:").grid(row=0, column=8, padx=3, pady=3, sticky='e')
        ttk.Entry(format_frame, textvariable=self.opt_vars['legend_font_size'], width=8).grid(row=0, column=9, padx=3, pady=3)
        
        ttk.Label(format_frame, text="Font:").grid(row=0, column=10, padx=3, pady=3, sticky='e')
        ttk.Combobox(format_frame, textvariable=self.opt_vars['font_family'], values=['Arial', 'Times', 'Courier'], width=10).grid(row=0, column=11, padx=3, pady=3)

        # Labels frame
        labels_frame = ttk.LabelFrame(self.scrollable_frame, text="Axes Labels")
        labels_frame.pack(fill="x", padx=5, pady=5)
        
        ttk.Label(labels_frame, text="X Label:").grid(row=0, column=0, padx=3, pady=3, sticky='e')
        ttk.Entry(labels_frame, textvariable=self.opt_vars['x_label'], width=20).grid(row=0, column=1, padx=3, pady=3)
        
        ttk.Label(labels_frame, text="Y Label:").grid(row=0, column=2, padx=3, pady=3, sticky='e')
        ttk.Entry(labels_frame, textvariable=self.opt_vars['y_label'], width=20).grid(row=0, column=3, padx=3, pady=3)

        # Options frame
        central_frame = ttk.Frame(self.scrollable_frame)
        central_frame.pack(fill="x", padx=5, pady=10)
        
        # Checkboxes
        chk_frame = ttk.Frame(central_frame)
        chk_frame.pack(pady=2)
        ttk.Checkbutton(chk_frame, text="Show Grid", variable=self.opt_vars['show_grid']).pack(side="left", padx=10)
        ttk.Checkbutton(chk_frame, text="Normalize by Point", variable=self.opt_vars['normalize_point']).pack(side="left", padx=10)
        ttk.Checkbutton(chk_frame, text="Normalize by Max", variable=self.opt_vars['normalize_max']).pack(side="left", padx=10)
        ttk.Checkbutton(chk_frame, text="Use nm to eV", variable=self.opt_vars['use_nm_to_ev']).pack(side="left", padx=10)

        # Save Path
        path_frame = ttk.Frame(central_frame)
        path_frame.pack(pady=2)
        ttk.Label(path_frame, text="Save Path: ").pack(side="left")
        ttk.Label(path_frame, textvariable=self.opt_vars['save_path']).pack(side="left")

        # Change Path button
        btn_path_frame = ttk.Frame(central_frame)
        btn_path_frame.pack(pady=2)
        ttk.Button(btn_path_frame, text="Change Path", width=12).pack()

        # Plot / File Name / Save Plot
        plot_btn_frame = ttk.Frame(central_frame)
        plot_btn_frame.pack(pady=5)
        ttk.Button(plot_btn_frame, text="Plot", command=self.generate_plot).pack(side="left", padx=5)
        ttk.Label(plot_btn_frame, text="File Name:").pack(side="left", padx=2)
        ttk.Entry(plot_btn_frame, textvariable=self.opt_vars['file_name'], width=20).pack(side="left", padx=2)
        ttk.Button(plot_btn_frame, text="Save Plot", command=self.save_plot).pack(side="left", padx=5)
        ttk.Button(plot_btn_frame, text="Save Defaults", command=self.save_defaults).pack(side="left", padx=5)

        # Container for the dynamic series rows
        self.series_container = ttk.Frame(self.scrollable_frame)
        self.series_container.pack(fill="x", padx=5, pady=10)
        
        # Initially, the dynamic plot options are not shown
        # self._add_mock_series_rows()

    def _build_dynamic_plot_options(self):
        # Clear existing dynamic rows
        for widget in self.series_container.winfo_children():
            widget.destroy()

        # Define color palette
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta', 'yellow', 'teal', 'navy', 'maroon', 'olive', 'lime', 'coral', 'indigo', 'gold', 'silver']

        # Initialize series options storage
        self.series_options = {}
        
        for idx, plot_data in enumerate(self.selected_plots):
            title = f"Plot "
            default_color = colors[idx % len(colors)]
            default_label = f"{plot_data['dataset']}-{plot_data['datatype']}"
            self.series_options[idx] = self._create_series_row(self.series_container, title, idx + 1, f": Wavelength (nm) vs {plot_data['dataset']} ({plot_data['datatype']})", default_color, default_label)

    def _create_series_row(self, parent, title_prefix, plot_num, title_suffix, default_color, default_label):
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=10)
        
        title_frame = ttk.Frame(frame)
        title_frame.grid(row=0, column=0, columnspan=16, sticky='w', pady=(0, 5))
        
        ttk.Label(title_frame, text=title_prefix).pack(side='left')
        
        e_order = ttk.Entry(title_frame, width=3)
        e_order.pack(side='left')
        e_order.insert(0, str(plot_num))
        
        ttk.Label(title_frame, text=title_suffix).pack(side='left')
        
        ttk.Label(frame, text="Line Thickness:").grid(row=1, column=0, padx=2)
        e_thick = ttk.Entry(frame, width=5)
        e_thick.grid(row=1, column=1, padx=2)
        e_thick.insert(0, "1.5")
        
        ttk.Label(frame, text="Line Color:").grid(row=1, column=2, padx=2)
        e_col = ttk.Entry(frame, width=10)
        e_col.grid(row=1, column=3, padx=2)
        e_col.insert(0, default_color)
        
        ttk.Label(frame, text="Legend Label:").grid(row=1, column=4, padx=2)
        e_leg = ttk.Entry(frame, width=20)
        e_leg.grid(row=1, column=5, padx=2)
        e_leg.insert(0, default_label)
        
        ttk.Label(frame, text="Line Style:").grid(row=1, column=6, padx=2)
        cb_style = ttk.Combobox(frame, values=["Solid", "Dashed", "Dotted"], width=8)
        cb_style.grid(row=1, column=7, padx=2)
        cb_style.set("Solid")
        
        ttk.Label(frame, text="Normalize Point:").grid(row=1, column=8, padx=2)
        e_norm = ttk.Entry(frame, width=5)
        e_norm.grid(row=1, column=9, padx=2)
        e_norm.insert(0, "1.0")
        
        var_scatter = tk.BooleanVar(value=False)
        chk_scatter = ttk.Checkbutton(frame, text="Scatter Graph", variable=var_scatter)
        chk_scatter.grid(row=1, column=10, padx=5)
        
        ttk.Label(frame, text="Point Size:").grid(row=1, column=11, padx=2)
        e_pt = ttk.Entry(frame, width=5)
        e_pt.grid(row=1, column=12, padx=2)
        e_pt.insert(0, "15")
        
        # Return a dictionary of options for this series
        return {
            'order': e_order,
            'line_thickness': e_thick,
            'line_color': e_col,
            'legend_label': e_leg,
            'line_style': cb_style,
            'normalize_point': e_norm,
            'scatter_graph': var_scatter,
            'point_size': e_pt
        }

    def save_plot(self):
        """Save plot (placeholder)"""
        pass

    def trigger_plot(self):
        print("Gathering data to plot based on:")
        for ds, types in self.plot_vars.items():
            for dtype, var in types.items():
                if var.get():
                    print(f" - Dataset: {ds}, Type: {dtype}")
        
        print("\nPlot Options:")
        for k, v in self.opt_vars.items():
            print(f" - {k}: {v.get()}")
    
    class Plotclass:
        def __init__(self, speclist):
            """
            Plotclass to handle advanced plotting of spectral data.
            
            Args:
                speclist: List of spectral data dictionaries to plot
            """
            self.speclist = speclist

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Spectra Plotter GUI")
    root.geometry("1200x700")
    
    # Test with dummy data
    app = Specplottergui(guiroot=root)
    root.mainloop()


