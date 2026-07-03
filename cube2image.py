import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class Cube2ImageGUI:
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
        self.root = root
        self.Nanomap = Nanomap
        self.colormaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds']
        # overwrite colormaps with matplotlib's registered colormaps, but keep the default list order
        self.colormaps = [cmap for cmap in self.colormaps if cmap in plt.colormaps()]
        self.colormap = 'viridis'
        self.default_colormap = 'viridis'
        
        main_frame = ttk.Frame(self.root, padding='10')
        main_frame.grid(row=0, column=0, sticky='nsew')
        
        # Datatype combobox
        ttk.Label(main_frame, text='Select spectral datatype:').grid(row=0, column=0, sticky=tk.W)
        self.datatype_var = tk.StringVar()
        self.datatype_cb = ttk.Combobox(main_frame, textvariable=self.datatype_var)
        self.datatype_cb.grid(row=0, column=1, sticky='we')
        
        # Start and Width sliders
        ttk.Label(main_frame, text='Integration Central WL (nm):').grid(row=1, column=0, sticky=tk.W)
        self.start_slider = ttk.Scale(main_frame, from_=self.wlstart, to=self.wlend, value=(self.wlend-self.wlstart)/2, command=self.update_plot, length=self.zoomlen)
        self.start_slider.grid(row=1, column=1, sticky='we')
        self.start_val_label = ttk.Label(main_frame, text='500.0', width=8)
        self.start_val_label.grid(row=1, column=2, sticky=tk.W)
        
        ttk.Label(main_frame, text='Wavelength Width (corresponding unit):').grid(row=2, column=0, sticky=tk.W)
        self.width_slider = ttk.Scale(main_frame, from_=0.0, to=200, value=10, command=self.update_plot, length=self.zoomlen)
        self.width_slider.grid(row=2, column=1, sticky='we')
        self.width_val_label = ttk.Label(main_frame, text='10.0', width=8)
        self.width_val_label.grid(row=2, column=2, sticky=tk.W)

        # Manual wavelength bounds
        manual_frame = ttk.Frame(main_frame)
        manual_frame.grid(row=3, column=0, columnspan=3, sticky='we', pady=(4, 0))
        ttk.Label(manual_frame, text='WL Start:').grid(row=0, column=0, sticky=tk.W)
        self.manual_wlstart_var = tk.StringVar(value=f'{self.wlstart:.2f}')
        self.manual_wlstart_entry = ttk.Entry(manual_frame, textvariable=self.manual_wlstart_var, width=10)
        self.manual_wlstart_entry.grid(row=0, column=1, sticky=tk.W, padx=(4, 10))
        ttk.Label(manual_frame, text='WL width:').grid(row=0, column=2, sticky=tk.W)
        self.manual_wl_width_var = tk.StringVar(value='10.00')
        self.manual_wlend_entry = ttk.Entry(manual_frame, textvariable=self.manual_wl_width_var, width=10)
        self.manual_wlend_entry.grid(row=0, column=3, sticky=tk.W, padx=(4, 10))
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
    
    def update_plot(self, *args):
        centerwl = float(self.start_slider.get())
        start = centerwl - float(self.width_slider.get()) / 2
        width = float(self.width_slider.get())
        
        self.start_val_label.config(text=f"{start:.1f}")
        self.width_val_label.config(text=f"{width:.1f}")

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
            wl = getattr(self.Nanomap.SpecDataMatrix[0][0], 'WL', None)
            if wl is not None:
                idx = np.where((wl >= start) & (wl <= start + width))[0]
                img = np.zeros((len(self.Nanomap.SpecDataMatrix), len(self.Nanomap.SpecDataMatrix[0])))
                for i in range(img.shape[0]):
                    for j in range(img.shape[1]):
                        data = getattr(self.Nanomap.SpecDataMatrix[i][j], dt, np.zeros_like(wl))
                        img[i, j] = np.sum(data[idx])
                
                self._draw_image(img, title=f'{dt_label}: {start:.1f} - {start + width:.1f} nm')
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
            wlstart = round(float(self.manual_wlstart_var.get()), 2)
            wlwidth = round(float(self.manual_wl_width_var.get()), 2)
        except Exception:
            self._clear_plot()
            self.ax.text(0.5, 0.5, 'Invalid wavelength bounds', ha='center', va='center')
            self.canvas.draw_idle()
            return

        wlstart = max(self.wlstart, min(self.wlend, wlstart))
        max_width = max(0.0, self.wlend - wlstart)
        wlwidth = max(0.0, min(max_width, wlwidth))

        self.manual_wlstart_var.set(f'{wlstart:.2f}')
        self.manual_wl_width_var.set(f'{wlwidth:.2f}')
        self.start_slider.set(wlstart)
        self.width_slider.set(wlwidth)
        self.update_plot()
    
    def createHSI(self):

        start = float(self.start_slider.get())
        width = float(self.width_slider.get())

        if not self.Nanomap: return
        dt_label = self.datatype_var.get()
        dt = self.datatype_map.get(dt_label)
        if not dt: return

        wlstart = start
        wlend = start + width

        if hasattr(self.Nanomap, 'selectspecbox'):
            try:
                self.Nanomap.selectspecbox.set(dt_label)
            except Exception:
                pass
        
        if hasattr(self.Nanomap, 'proc_spec_min') and hasattr(self.Nanomap, 'proc_spec_max'):
            try:
                self.Nanomap.proc_spec_min.delete(0, tk.END)
                self.Nanomap.proc_spec_min.insert(0, str(round(wlstart, 2)))
                self.Nanomap.proc_spec_max.delete(0, tk.END)
                self.Nanomap.proc_spec_max.insert(0, str(round(wlend, 2)))
            except Exception:
                pass

        if hasattr(self.Nanomap, 'buildandPlotIntCmap'):
            try:
                self.Nanomap.buildandPlotIntCmap(savetoimage='False', plot=False, datatype=dt, wlstart=wlstart, wlend=wlend)
                self.update_plot()
            except Exception as e:
                self._clear_plot()
                self.ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center')
                self.canvas.draw_idle()
                return
    
    def update_bounds(self, wlstart, wlend):
        self.wlstart = wlstart
        self.wlend = wlend
        self.start_slider.config(from_=self.wlstart, to=self.wlend, length=self.zoomlen)
        self.start_slider.set((self.wlend - self.wlstart) / 2)
        self.update_plot()
    
    def destroy(self):
        # clean up the GUI resources
        try:
            self.start_slider.destroy()
            self.width_slider.destroy()
            self.datatype_cb.destroy()
            self.start_val_label.destroy()
            self.width_val_label.destroy()
            self.manual_wlstart_entry.destroy()
            self.manual_wlend_entry.destroy()
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
        # break references so tkinter can close cleanly
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
        def __init__(self):
            self.WL = np.linspace(400, 700, 100)
            self.Spectrum1 = np.random.rand(100)
            self.Spectrum2 = np.random.rand(100)
    
    class DummyNanomap:
        def __init__(self):
            self.speckeys = ['Spectrum1', 'Spectrum2']
            self.SpecDataMatrix = [[DummySpec() for _ in range(5)] for _ in range(5)]
    
    nanomap = DummyNanomap()
    
    cube2image_gui = Cube2Image(Nanomap=nanomap, guiroot=root)
    
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
