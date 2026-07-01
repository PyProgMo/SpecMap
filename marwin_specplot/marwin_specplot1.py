"""
HSI Explorer for PLEM/HSI TXT exports
Start in Anaconda Prompt:
    python HSI_Explorer.py
Controls:
    Left click on map  = add pixel/spectrum
    Right click on map = remove last selected pixel
    c                  = clear all selected pixels
    s                  = save selected pixel table as CSV
    p                  = save current figure as PNG
Features:
    - Toggle normalization on/off in the GUI
    - Peak position and FWHM are printed, shown in the status bar/list, and exported to CSV
Edit the SETTINGS block below if needed.
"""
import glob
import os
import re
import sys
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.signal import savgol_filter
# ==========================================================
# SETTINGS
# ==========================================================
#'F:\Promotion\2026_MoSe2_BNNTs\MoSe2\HSI20260603_HSI_MoSe2_BNNT6T_M3_500nm'
DEFAULT_FOLDER = 'F:\\Promotion\\2026_MoSe2_BNNTs\\MoSe2\\HSI20260603_HSI_MoSe2_BNNT6T_M3_500nm'
#r"C:\Users\s416107\Desktop\Daten HSI\HSI20260603_HSI_MoSe2_BNNT6T_M3_500nm"
PEAK_MIN = 760
PEAK_MAX = 810
MAP_INTEGRATION_MIN = 760
MAP_INTEGRATION_MAX = 810
PIXEL_X_MAX = 50
PIXEL_Y_MAX = 56
SMOOTH_WINDOW = 11
SMOOTH_POLYORDER = 3
NORMALIZE_SPECTRA_DEFAULT = True
SPECTRUM_XLIM = (725, 850)
SPECTRUM_YLIM_NORMALIZED = (0, 1.05)
SPECTRUM_YLIM_RAW = None  # None = autoscale
COLORS = [
    "red", "blue", "green", "orange", "magenta", "cyan",
    "black", "purple", "brown", "lime", "navy", "gold"
]
# ==========================================================
# DATA LOADING
# ==========================================================
def find_data_start(text: str) -> int:
    """Return skiprows index for the numeric data after the WL/BG/PL header."""
    match = re.search(r"^\s*WL\s+BG\s+PL", text, flags=re.MULTILINE)
    if match is None:
        raise ValueError("No data header line 'WL BG PL' found.")
    return text[:match.end()].count("\n") + 1
def read_one_file(path: str):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    x_match = re.search(r"x-position:\s*([0-9.+\-Ee]+)", text)
    y_match = re.search(r"y-position:\s*([0-9.+\-Ee]+)", text)
    if x_match is None or y_match is None:
        raise ValueError(f"No x-position/y-position found in {path}")
    x = float(x_match.group(1))
    y = float(y_match.group(1))
    start = find_data_start(text)
    data = np.loadtxt(path, skiprows=start)
    wl = data[:, 0]
    bg = data[:, 1]
    pl = data[:, 2]
    pl_corr = pl - bg
    return x, y, wl, pl_corr
def load_dataset(folder: str):
    files = sorted(glob.glob(os.path.join(folder, "*.txt")))
    if not files:
        raise FileNotFoundError(f"No .txt files found in folder:\n{folder}")
    x_all = []
    y_all = []
    spectra = []
    wl_ref = None
    for i, file in enumerate(files):
        try:
            x, y, wl, pl_corr = read_one_file(file)
        except Exception as exc:
            print(f"Skipped {os.path.basename(file)}: {exc}")
            continue
        if wl_ref is None:
            wl_ref = wl
        elif len(wl) != len(wl_ref) or not np.allclose(wl, wl_ref, atol=1e-6):
            print(f"Skipped {os.path.basename(file)}: wavelength axis differs")
            continue
        if SMOOTH_WINDOW and SMOOTH_WINDOW >= 5:
            pl_corr = savgol_filter(pl_corr, SMOOTH_WINDOW, SMOOTH_POLYORDER)
        x_all.append(x)
        y_all.append(y)
        spectra.append(pl_corr)
    if wl_ref is None or len(spectra) == 0:
        raise RuntimeError("No spectra could be loaded.")
    x_all = np.array(x_all)
    y_all = np.array(y_all)
    spectra = np.array(spectra)
    map_region = (wl_ref >= MAP_INTEGRATION_MIN) & (wl_ref <= MAP_INTEGRATION_MAX)
    int_all = np.trapezoid(spectra[:, map_region], wl_ref[map_region], axis=1)
    peak_region = (wl_ref >= PEAK_MIN) & (wl_ref <= PEAK_MAX)
    peak_indices = np.argmax(spectra[:, peak_region], axis=1)
    peak_wl_axis = wl_ref[peak_region]
    peak_all = peak_wl_axis[peak_indices]
    fwhm_all = []
    for s in spectra:
        sig = s[peak_region]
        wl_peak = wl_ref[peak_region]
        half = np.max(sig) / 2
        above = sig >= half
        if np.sum(above) > 2:
            fwhm_all.append(wl_peak[above][-1] - wl_peak[above][0])
        else:
            fwhm_all.append(np.nan)
    fwhm_all = np.array(fwhm_all)
    # Pixel coordinates resembling the original HSI viewer, horizontally mirrored/oriented.
    px_all = (x_all - np.min(x_all)) / (np.max(x_all) - np.min(x_all)) * PIXEL_X_MAX
    py_all = (np.max(y_all) - y_all) / (np.max(y_all) - np.min(y_all)) * PIXEL_Y_MAX
    return {
        "folder": folder,
        "files": files,
        "wl": wl_ref,
        "x": x_all,
        "y": y_all,
        "px": px_all,
        "py": py_all,
        "spectra": spectra,
        "intensity": int_all,
        "peak": peak_all,
        "fwhm": fwhm_all,
    }
# ==========================================================
# APP
# ==========================================================
class HSIExplorer:
    def __init__(self, root):
        self.root = root
    
    def load_dataset(self, folder):
        # check if folder exists
        if not os.path.exists(folder):
            print(f"Folder does not exist: {folder}")
            return
        try:
            dataset = load_dataset(folder)
        except Exception as exc:
            raise
        print(f"Loaded spectra: {len(dataset['spectra'])}")
        print(f"Peak mean: {np.nanmean(dataset['peak']):.3f} nm")
        print(f"Peak std: {np.nanstd(dataset['peak']):.3f} nm")
        print(f"FWHM mean: {np.nanmean(dataset['fwhm']):.3f} nm")
        print(f"FWHM std: {np.nanstd(dataset['fwhm']):.3f} nm")
        self.buildgui(dataset)

    def clear(self):
        # Clear the GUI and data for a new dataset
        self.data = None
        self.selected = []
        self.markers = []
        self.texts = []
        self.fig = None
        self.ax_map = None
    
    def buildgui(self, dataset):
        self.data = dataset
        self.selected = []
        self.markers = []
        self.texts = []
        self.fig, (self.ax_map, self.ax_spec) = plt.subplots(1, 2, figsize=(14, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        # Control panel
        control_frame = tk.Frame(self.root)
        control_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.normalize_var = tk.BooleanVar(value=NORMALIZE_SPECTRA_DEFAULT)
        tk.Checkbutton(
            control_frame,
            text="Normalize spectra",
            variable=self.normalize_var,
            command=self.redraw_spectra
        ).pack(side=tk.LEFT, padx=6)
        tk.Button(control_frame, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=4)
        tk.Button(control_frame, text="Save CSV", command=self.save_csv).pack(side=tk.LEFT, padx=4)
        tk.Button(control_frame, text="Save PNG", command=self.save_png).pack(side=tk.LEFT, padx=4)
        self.status = tk.StringVar()
        self.status.set("Left click: add spectrum | Right click: remove last | c: clear | s: save CSV | p: save PNG")
        tk.Label(self.root, textvariable=self.status, anchor="w").pack(side=tk.BOTTOM, fill=tk.X)
        self.info_box = tk.Text(self.root, height=7, wrap="none")
        self.info_box.pack(side=tk.BOTTOM, fill=tk.X)
        self.info_box.insert("end", "Selected pixels will appear here.\n")
        self.info_box.configure(state="disabled")
        self.draw_initial()
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
    def draw_initial(self):
        d = self.data
        self.ax_map.clear()
        self.ax_spec.clear()
        sc = self.ax_map.scatter(d["px"], d["py"], c=d["intensity"], s=55, cmap="inferno")
        self.fig.colorbar(sc, ax=self.ax_map, label=f"Integrated PL intensity {MAP_INTEGRATION_MIN}–{MAP_INTEGRATION_MAX} nm")
        self.ax_map.set_title("HSI map: click pixels")
        self.ax_map.set_xlabel("Pixel X")
        self.ax_map.set_ylabel("Pixel Y")
        self.ax_map.invert_yaxis()
        self.ax_map.set_aspect("equal")
        self.setup_spectrum_axis()
        self.fig.tight_layout()
        self.canvas.draw_idle()
    def setup_spectrum_axis(self):
        normalized = self.normalize_var.get() if hasattr(self, "normalize_var") else NORMALIZE_SPECTRA_DEFAULT
        self.ax_spec.set_title("Selected spectra")
        self.ax_spec.set_xlabel("Wavelength (nm)")
        self.ax_spec.set_ylabel("Normalized intensity" if normalized else "PL - BG")
        self.ax_spec.set_xlim(*SPECTRUM_XLIM)
        if normalized:
            self.ax_spec.set_ylim(*SPECTRUM_YLIM_NORMALIZED)
        elif SPECTRUM_YLIM_RAW is not None:
            self.ax_spec.set_ylim(*SPECTRUM_YLIM_RAW)
        self.ax_spec.grid(True)
    def nearest_index(self, px, py):
        d = self.data
        dist = np.sqrt((d["px"] - px) ** 2 + (d["py"] - py) ** 2)
        return int(np.argmin(dist))
    def spectrum_for_plot(self, idx):
        d = self.data
        wl = d["wl"]
        s = d["spectra"][idx].copy()
        peak_region = (wl > PEAK_MIN) & (wl < PEAK_MAX)
        if self.normalize_var.get():
            max_val = np.max(s[peak_region])
            if max_val != 0:
                s = s / max_val
        return s
    def update_info_box(self):
        d = self.data
        self.info_box.configure(state="normal")
        self.info_box.delete("1.0", "end")
        if not self.selected:
            self.info_box.insert("end", "Selected pixels will appear here.\n")
        else:
            self.info_box.insert("end", "Nr | PixelX/Y | StageX/Y | Peak (nm) | FWHM (nm) | Integral\n")
            self.info_box.insert("end", "-" * 82 + "\n")
            for n, idx in enumerate(self.selected, start=1):
                self.info_box.insert(
                    "end",
                    f"{n:2d} | {d['px'][idx]:5.1f}/{d['py'][idx]:5.1f} | "
                    f"{d['x'][idx]:7.2f}/{d['y'][idx]:7.2f} | "
                    f"{d['peak'][idx]:8.2f} | {d['fwhm'][idx]:8.2f} | {d['intensity'][idx]:10.1f}\n"
                )
        self.info_box.configure(state="disabled")
    def redraw_spectra(self):
        d = self.data
        wl = d["wl"]
        self.ax_spec.clear()
        self.setup_spectrum_axis()
        for n, idx in enumerate(self.selected):
            color = COLORS[n % len(COLORS)]
            s = self.spectrum_for_plot(idx)
            peak = d["peak"][idx]
            self.ax_spec.plot(
                wl,
                s,
                color=color,
                linewidth=2.5,
                label=f"{n+1}: px {d['px'][idx]:.1f}/{d['py'][idx]:.1f}, peak {peak:.2f} nm, FWHM {d['fwhm'][idx]:.2f} nm"
            )
        if self.selected:
            self.ax_spec.legend(fontsize=8)
        if not self.normalize_var.get():
            self.ax_spec.relim()
            self.ax_spec.autoscale_view(scalex=False, scaley=True)
        self.update_info_box()
        self.canvas.draw_idle()
    def add_point(self, idx):
        d = self.data
        n = len(self.selected) + 1
        color = COLORS[(n - 1) % len(COLORS)]
        self.selected.append(idx)
        marker = self.ax_map.scatter(d["px"][idx], d["py"][idx], s=240, marker="x", color=color, linewidth=3)
        text = self.ax_map.text(d["px"][idx] + 0.5, d["py"][idx] + 0.5, str(n), color=color, fontsize=12, weight="bold")
        self.markers.append(marker)
        self.texts.append(text)
        msg = (
            f"Point {n}: Pixel {d['px'][idx]:.1f}/{d['py'][idx]:.1f} | "
            f"Stage {d['x'][idx]:.2f}/{d['y'][idx]:.2f} | "
            f"Peak {d['peak'][idx]:.2f} nm | FWHM {d['fwhm'][idx]:.2f} nm | "
            f"Integral {d['intensity'][idx]:.1f}"
        )
        print(msg)
        self.status.set(msg)
        self.redraw_spectra()
    def remove_last(self):
        if not self.selected:
            return
        self.selected.pop()
        if self.markers:
            self.markers.pop().remove()
        if self.texts:
            self.texts.pop().remove()
        self.status.set("Removed last point")
        self.redraw_spectra()
    def clear_all(self):
        self.selected.clear()
        while self.markers:
            self.markers.pop().remove()
        while self.texts:
            self.texts.pop().remove()
        self.status.set("Selection cleared")
        self.redraw_spectra()
    def save_csv(self):
        if not self.selected:
            print("No selected pixels to save.")
            return
        out_path = filedialog.asksaveasfilename(
            title="Save selected pixels",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile="selected_hsi_pixels.csv"
        )
        if not out_path:
            return
        d = self.data
        rows = []
        for n, idx in enumerate(self.selected, start=1):
            rows.append([n, d["px"][idx], d["py"][idx], d["x"][idx], d["y"][idx], d["peak"][idx], d["fwhm"][idx], d["intensity"][idx]])
        rows = np.array(rows)
        header = "Nr,PixelX,PixelY,StageX,StageY,Peak_nm,FWHM_nm,Integral"
        np.savetxt(out_path, rows, delimiter=",", header=header, comments="")
        self.status.set(f"Saved CSV: {out_path}")
        print(f"Saved CSV: {out_path}")
    def save_png(self):
        out_path = filedialog.asksaveasfilename(
            title="Save figure",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            initialfile="hsi_explorer_selection.png"
        )
        if not out_path:
            return
        self.fig.savefig(out_path, dpi=300, bbox_inches="tight")
        self.status.set(f"Saved PNG: {out_path}")
        print(f"Saved PNG: {out_path}")
    def on_click(self, event):
        if event.inaxes != self.ax_map:
            return
        if self.toolbar.mode:
            return
        if event.button == 1:
            idx = self.nearest_index(event.xdata, event.ydata)
            self.add_point(idx)
        elif event.button == 3:
            self.remove_last()
    def on_key(self, event):
        if event.key == "c":
            self.clear_all()
        elif event.key == "s":
            self.save_csv()
        elif event.key == "p":
            self.save_png()
            
def main():
    folder = DEFAULT_FOLDER
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    try:
        dataset = load_dataset(folder)
    except Exception as exc:
        print("HSI Explorer", str(exc))
        return
    print(f"Loaded spectra: {len(dataset['spectra'])}")
    print(f"Peak mean: {np.nanmean(dataset['peak']):.3f} nm")
    print(f"Peak std: {np.nanstd(dataset['peak']):.3f} nm")
    print(f"FWHM mean: {np.nanmean(dataset['fwhm']):.3f} nm")
    print(f"FWHM std: {np.nanstd(dataset['fwhm']):.3f} nm")
    root = tk.Tk()
    app = HSIExplorer(root)
    app.load_dataset(folder)
    root.mainloop()
if __name__ == "__main__":
    main()