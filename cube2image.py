import tkinter as tk
from tkinter import ttk
from collections import OrderedDict
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class Cube2ImageGUI:
    # Number of datatype cubes to keep cached simultaneously.
    # Each cached cube costs rows*cols*n_wl*8 bytes (float64) x2 (cube + cumsum),
    # so keep this small for large hyperspectral maps.
    CUBE_CACHE_SIZE = 2

    def __init__(self, root, Nanomap=None, wlstart=0.0, wlend=1000.0, zoomlen=700.0):
        try: 
            self.wlstart = float(wlstart)
        except:
            self.wlstart = 0.0
        try:
            self.wlend = float(wlend)
        except:
            self.wlend = 1000.0
        try:
            self.zoomlen = float(zoomlen)
        except:
            self.zoomlen = 700.0
        self.wlcenter = tk.DoubleVar(value=(self.wlend + self.wlstart) / 2)
        self.wlwidth = tk.DoubleVar(value=10.0)
        self.root = root
        self.Nanomap = Nanomap
        self.colormaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds']
        # overwrite colormaps with matplotlib's registered colormaps, but keep the default list order
        self.colormaps = [cmap for cmap in self.colormaps if cmap in plt.colormaps()]
        self.colormap = 'viridis'
        self.default_colormap = 'viridis'

        # --- performance cache: dt_label -> (wl, cumsum) ---
        # cumsum has shape (rows, cols, n_wl+1), cumsum[:,:,k] = sum of band[0:k]
        # along the wavelength axis, so any window sum is a single slice/subtract.
        self._cube_cache = OrderedDict()
        self._cached_smat_id = None  # identifies which SpecDataMatrix the cache belongs to
        
        main_frame = ttk.Frame(self.root, padding='10')
        main_frame.grid(row=0, column=0, sticky='nsew')
        
        # Datatype combobox
        ttk.Label(main_frame, text='Select spectral datatype:').grid(row=0, column=0, sticky=tk.W)
        self.datatype_var = tk.StringVar()
        self.datatype_cb = ttk.Combobox(main_frame, textvariable=self.datatype_var)
        self.datatype_cb.grid(row=0, column=1, sticky='we')
        
        # Center and Width sliders
        ttk.Label(main_frame, text='Integration Central WL (nm):').grid(row=1, column=0, sticky=tk.W)
        self.center_slider = ttk.Scale(main_frame, from_=self.wlstart, to=self.wlend, variable=self.wlcenter, command=self.update_plot, length=self.zoomlen)
        self.center_slider.grid(row=1, column=1, sticky='we')
        self.center_val_label = ttk.Label(main_frame, text=f'{self.wlcenter.get():.1f}', width=8)
        self.center_val_label.grid(row=1, column=2, sticky=tk.W)
        
        ttk.Label(main_frame, text='Wavelength Width (corresponding unit):').grid(row=2, column=0, sticky=tk.W)
        self.width_slider = ttk.Scale(main_frame, from_=0.0, to=200, variable=self.wlwidth, command=self.update_plot, length=self.zoomlen)
        self.width_slider.grid(row=2, column=1, sticky='we')
        self.width_val_label = ttk.Label(main_frame, text='10.0', width=8)
        self.width_val_label.grid(row=2, column=2, sticky=tk.W)

        # Manual wavelength bounds
        manual_frame = ttk.Frame(main_frame)
        manual_frame.grid(row=3, column=0, columnspan=3, sticky='we', pady=(4, 0))
        ttk.Label(manual_frame, text='WL Center:').grid(row=0, column=0, sticky=tk.W)
        self.manual_wlcenter_var = tk.StringVar(value=f'{self.wlcenter.get():.2f}')
        self.manual_wlcenter_entry = ttk.Entry(manual_frame, textvariable=self.manual_wlcenter_var, width=10)
        self.manual_wlcenter_entry.grid(row=0, column=1, sticky=tk.W, padx=(4, 10))
        ttk.Label(manual_frame, text='WL width:').grid(row=0, column=2, sticky=tk.W)
        self.manual_wl_width_var = tk.StringVar(value='10.00')
        self.manual_wlwidth_entry = ttk.Entry(manual_frame, textvariable=self.manual_wl_width_var, width=10)
        self.manual_wlwidth_entry.grid(row=0, column=3, sticky=tk.W, padx=(4, 10))
        self.set_wl_button = ttk.Button(manual_frame, text='Set WL', command=self.set_manual_wavelengths)
        self.set_wl_button.grid(row=0, column=4, sticky=tk.W)

        # add colormap selection next to Set WL button
        ttk.Label(manual_frame, text='Colormap:').grid(row=0, column=5, sticky=tk.W, padx=(10, 0))
        self.colormap_var = tk.StringVar(value=self.colormap)
        # Use a Combobox for colormap selection
        self.colormap_cb = ttk.Combobox(manual_frame, textvariable=self.colormap_var, values=self.colormaps, width=12)
        self.colormap_cb.grid(row=0, column=6, sticky=tk.W, padx=(4, 0))
        # on colormap change, update self.colormap and redraw the plot
        self.colormap_cb.bind('<<ComboboxSelected>>', lambda e: self.change_colormap())
        
        # Create frame for Plot and Create HSI buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=4, column=0, columnspan=3, pady=10)
        ttk.Button(button_frame, text='Plot', command=self.update_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text='Create HSI', command=self.createHSI).pack(side=tk.LEFT, padx=5)
        
        # Matplotlib canvas
        self.fig = Figure(figsize=(5, 5))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=main_frame)
        self.canvas.get_tk_widget().grid(row=5, column=0, columnspan=3, sticky='nsew')

        # figure for the line profiles
        self.hl_profile_fig = Figure(figsize=(5, 2))
        self.vl_profile_fig = Figure(figsize=(2, 5))
        self.hl_profile_ax = self.hl_profile_fig.add_subplot(111)
        self.vl_profile_ax = self.vl_profile_fig.add_subplot(111)
        
        self.datatype_map = {
            'Wavelength axis': 'WL', 'Background (BG)': 'BG', 'Counts (PL)': 'PL', 'Spectrum (PL-BG)': 'PLB', 
            'first derivative': 'Specdiff1', 'second derivative': 'Specdiff2', 
            'first derivative (normalized)': 'Specdiff1_norm', 'second derivative (normalized)': 'Specdiff2_norm', 
            'first derivative (norm on intensity, then derive)': 'Specdiff1_norm_intensity', 
            'second derivative (norm on intensity, then derive)': 'Specdiff2_norm_intensity', 
            'first derivative (norm on counts, then derive)': 'Specdiff1_norm_counts', 
            'second derivative (norm on counts, then derive)': 'Specdiff2_norm_counts'
        }
        
        self.datatype_cb['values'] = list(self.datatype_map.keys())
        self.datatype_cb.current(3) # Default to 'Spectrum (PL-BG)'
        self.datatype_cb.bind('<<ComboboxSelected>>', self.update_plot)

        self.image_artist = None
        self.colorbar = None

        self.update_plot()  # Initial plot
    
    def change_colormap(self):
        selected_cmap = self.colormap_var.get()
        if selected_cmap in self.colormaps:
            self.colormap = selected_cmap
            self.update_plot()
        else:
            self.colormap = self.default_colormap
            self.update_plot()

    def _clear_plot(self):
        if self.image_artist is not None:
            try:
                self.image_artist.remove()
            except Exception:
                pass
            self.image_artist = None
        if self.colorbar is not None:
            try:
                self.colorbar.remove()
            except Exception:
                pass
            self.colorbar = None
        self.ax.clear()

    def _draw_image(self, img, title=None):
        self._clear_plot()
        self.image_artist = self.ax.imshow(img, cmap=self.colormap)
        if title:
            self.ax.set_title(title)
        self.canvas.draw_idle()

    def _centered_bounds(self, centerwl, width):
        width = max(0.0, min(float(width), self.wlend - self.wlstart))
        half_width = width / 2.0
        min_center = self.wlstart + half_width
        max_center = self.wlend - half_width
        if min_center <= max_center:
            centerwl = min(max(centerwl, min_center), max_center)
        else:
            centerwl = (self.wlstart + self.wlend) / 2.0
        start = centerwl - half_width
        end = centerwl + half_width
        return centerwl, width, start, end

    # Optimized Version for faster performance: cache the computed cubes for each datatype, and only recompute if the underlying SpecDataMatrix changes (detected by id()). Optimization procedure provided by claude. 
    def invalidate_cube_cache(self):
        """Call this whenever Nanomap.SpecDataMatrix contents change
        (e.g. after createHSI/buildandPlotIntCmap regenerates spectra),
        or when a new Nanomap is attached. Cheap - just drops cached arrays."""
        self._cube_cache.clear()
        self._cached_smat_id = None

    def _get_cube(self, dt):
        """Return (wl, cumsum) for datatype `dt`, building & caching it once.

        cumsum has shape (rows, cols, n_wl + 1) with cumsum[..., 0] = 0, so that
        the integral over wl indices [i0, i1) for every pixel simultaneously is
        just cumsum[..., i1] - cumsum[..., i0] -- one vectorized numpy op,
        independent of band width or number of pixels.
        """
        smat = self.Nanomap.SpecDataMatrix
        smat_id = id(smat)
        if smat_id != self._cached_smat_id:
            # underlying data object changed (new scan / new Nanomap) - drop stale cache
            self._cube_cache.clear()
            self._cached_smat_id = smat_id

        cached = self._cube_cache.get(dt)
        if cached is not None:
            # mark as most-recently-used
            self._cube_cache.move_to_end(dt)
            return cached

        wl = getattr(smat[0][0], 'WL', None)
        if wl is None:
            return None
        wl = np.asarray(wl, dtype=float)

        rows = len(smat)
        cols = len(smat[0])
        n_wl = wl.shape[0]

        cube = np.zeros((rows, cols, n_wl), dtype=float)
        for i in range(rows):
            row = smat[i]
            for j in range(cols):
                data = getattr(row[j], dt, None)
                if data is not None:
                    cube[i, j, :] = data

        cumsum = np.empty((rows, cols, n_wl + 1), dtype=float)
        cumsum[:, :, 0] = 0.0
        np.cumsum(cube, axis=2, out=cumsum[:, :, 1:])

        result = (wl, cumsum)
        self._cube_cache[dt] = result
        if len(self._cube_cache) > self.CUBE_CACHE_SIZE:
            self._cube_cache.popitem(last=False)  # evict least-recently-used
        return result
    
    def update_plot(self, *args):
        centerwl = float(self.wlcenter.get())
        width = float(self.wlwidth.get())
        centerwl, width, start, end = self._centered_bounds(centerwl, width)
        
        self.center_val_label.config(text=f"{centerwl:.1f}")
        self.width_val_label.config(text=f"{width:.1f}")
        self.wlcenter.set(centerwl)
        self.wlwidth.set(width)

        if not self.Nanomap:
            self._clear_plot()
            self.canvas.draw_idle()
            return
        dt_label = self.datatype_var.get()
        dt = self.datatype_map.get(dt_label)
        if not dt:
            self._clear_plot()
            self.canvas.draw_idle()
            return
        
        try:
            cached = self._get_cube(dt)
            if cached is not None:
                wl, cumsum = cached
                # wl assumed monotonically increasing (as in the original np.where mask)
                i0 = int(np.searchsorted(wl, start, side='left'))
                i1 = int(np.searchsorted(wl, end, side='right'))
                img = cumsum[:, :, i1] - cumsum[:, :, i0]

                self._draw_image(img, title=f'{dt_label}: {start:.1f} - {end:.1f} nm')
                return
        except Exception as e:
            self._clear_plot()
            self.ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center')
            self.canvas.draw_idle()
            return

        self._clear_plot()
        self.ax.text(0.5, 0.5, 'No wavelength data available', ha='center', va='center')
        self.canvas.draw_idle()


    def set_manual_wavelengths(self):
        try:
            wlcenter = round(float(self.manual_wlcenter_var.get()), 2)
            wlwidth = round(float(self.manual_wl_width_var.get()), 2)
        except Exception:
            self._clear_plot()
            self.ax.text(0.5, 0.5, 'Invalid wavelength bounds', ha='center', va='center')
            self.canvas.draw_idle()
            return

        wlcenter, wlwidth, _, _ = self._centered_bounds(wlcenter, wlwidth)

        self.manual_wlcenter_var.set(f'{wlcenter:.2f}')
        self.manual_wl_width_var.set(f'{wlwidth:.2f}')
        self.center_slider.set(wlcenter)
        self.width_slider.set(wlwidth)
        self.update_plot()
    
    def createHSI(self):

        centerwl = float(self.center_slider.get())
        width = float(self.width_slider.get())

        if not self.Nanomap: return
        dt_label = self.datatype_var.get()
        dt = self.datatype_map.get(dt_label)
        if not dt: return

        centerwl, width, wlstart, wlend = self._centered_bounds(centerwl, width)
        self.center_slider.set(centerwl)
        self.width_slider.set(width)
        self.manual_wlcenter_var.set(f'{centerwl:.2f}')
        self.manual_wl_width_var.set(f'{width:.2f}')

        if hasattr(self.Nanomap, 'selectspecbox'):
            try:
                self.Nanomap.selectspecbox.set(dt_label)
            except Exception:
                pass

        if hasattr(self.Nanomap, 'buildandPlotIntCmap'):
            try:
                self.Nanomap.buildandPlotIntCmap(savetoimage='False', plot=False, datatype=dt, wlstart=round(wlstart, 2), wlend=round(wlend, 2))
                # buildandPlotIntCmap may have mutated SpecDataMatrix entries in place
                # (same object id), so an id()-based cache check wouldn't catch it -
                # invalidate explicitly to avoid serving stale cached spectra.
                self.invalidate_cube_cache()
                self.update_plot()
            except Exception as e:
                self._clear_plot()
                self.ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center')
                self.canvas.draw_idle()
                return
    
    def update_bounds(self, wlstart, wlend):
        self.wlstart = wlstart
        self.wlend = wlend
        self.center_slider.config(from_=self.wlstart, to=self.wlend, length=self.zoomlen)
        self.center_slider.set((self.wlend + self.wlstart) / 2)
        self.manual_wlcenter_var.set(f'{(self.wlend + self.wlstart) / 2:.2f}')
        self.update_plot()
    
    def destroy(self):
        # clean up the GUI resources
        try:
            self.center_slider.destroy()
            self.width_slider.destroy()
            self.datatype_cb.destroy()
            self.center_val_label.destroy()
            self.width_val_label.destroy()
            self.manual_wlcenter_entry.destroy()
            self.manual_wlwidth_entry.destroy()
            self.set_wl_button.destroy()
        except Exception:
            pass
        # clean up the matplotlib resources
        try:
            self.fig.clear()
            self.canvas.get_tk_widget().destroy()
            plt.close(self.fig)
        except Exception:
            pass
        # drop cached cubes and references so tkinter/gc can close cleanly
        self._cube_cache.clear()
        self.Nanomap = None
        self.root = None

class lineprofiles:
    def __init__(self, xdata, ydata, label=None, mplwidget=None):
        self.xdata = xdata
        self.ydata = ydata
        self.label = label
        self.mplwidget = mplwidget

class Cube2Image:
    def __init__(self, Nanomap=None, guiroot=None):
        # befor passing guiroot: add scrollbar to guiroot if it doesn't have one
        frame = self.add_scrollbar(guiroot)
        self.gui = Cube2ImageGUI(frame, Nanomap)
    
    def add_scrollbar(self, root):
        # Check if the root already has a scrollbar
        if isinstance(root, tk.Tk) or isinstance(root, tk.Toplevel):
            # Create a canvas and a vertical scrollbar for scrolling
            canvas = tk.Canvas(root)
            scrollbar = ttk.Scrollbar(root, orient="vertical", command=canvas.yview)
            scrollable_frame = ttk.Frame(canvas)

            scrollable_frame.bind(
                "<Configure>",
                lambda e: canvas.configure(
                    scrollregion=canvas.bbox("all")
                )
            )

            canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)
            # add a horizontal scrollbar
            h_scrollbar = ttk.Scrollbar(root, orient="horizontal", command=canvas.xview)
            canvas.configure(xscrollcommand=h_scrollbar.set)
            h_scrollbar.pack(side="bottom", fill="x")

            # Pack the canvas and scrollbar
            canvas.pack(side="left", fill="both", expand=True)
            scrollbar.pack(side="right", fill="y")
            h_scrollbar.pack(side="bottom", fill="x")

            return scrollable_frame
        else:
            return root  # If it's not a Tk or Toplevel, just return it as is

    def destroy(self):
        if hasattr(self, 'gui') and self.gui is not None:
            self.gui.destroy()
            self.gui = None

    # Backward-compatible typo alias
    def destory(self):
        self.destroy()
    def update_bounds(self):
        if not hasattr(self, 'gui') or self.gui is None:
            print("Cube2Image GUI is not available.")
            return
        nanomap = getattr(self.gui, 'Nanomap', None)
        if nanomap is not None and hasattr(nanomap, 'DataSpecMax') and hasattr(nanomap, 'DataSpecMin'):
            self.gui.update_bounds(nanomap.DataSpecMin, nanomap.DataSpecMax)
            print(f"Updated Cube2Image bounds to wlstart={nanomap.DataSpecMin}, wlend={nanomap.DataSpecMax}")
        else:
            print("Nanomap does not have wlstart and wlend attributes.")

def testgui():
    root = tk.Tk()
    root.title('Cube2Image Test')
    #root.protocol("WM_DELETE_WINDOW", root.destroy)
    
    # Create a dummy Nanomap with necessary attributes for testing
    class DummySpec:
        def __init__(self, spec=None):
            self.WL = np.linspace(400, 700, 100)
            if spec is not None:
                self.PLB = np.random.rand(100)
            elif spec == 'gaussian':
                # gaussian peak for testing
                self.PLB = np.exp(-0.5 * ((self.WL - 550) / 20) ** 2)
            elif spec == 'lorentzian':
                # lorentzian peak for testing
                self.PLB = 1 / (1 + ((self.WL - 600) / 10) ** 2)
            elif spec == 'sine':
                # sine wave for testing
                self.PLB = 1
    
    class DummyNanomap:
        def __init__(self):
            self.speckeys = ['Spectrum1', 'Spectrum2', 'Spectrum3', 'Spectrum4']
            self.specs = [DummySpec(spec='gaussian'), DummySpec(spec='lorentzian'), DummySpec(spec='sine'), DummySpec()]
            # test: print the spectra
            self.SpecDataMatrix = [[spec for spec in self.specs] for _ in range(4)]
            # for testing: print the spectra of the first pixel
    
    nanomap = DummyNanomap()
    
    cube2image_gui = Cube2Image(Nanomap=nanomap, guiroot=root)
    cube2image_gui.update_bounds()  # Update bounds based on DummyNanomap
    
    # bind the close event to ensure proper cleanup
    def on_closing():
        print("Closing Cube2Image GUI...")
        cube2image_gui.destroy()
        root.quit()
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_closing)

    root.mainloop()

if __name__ == '__main__':
    testgui()