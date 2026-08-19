import os, gc
import numpy as np
import deflib1 as deflib
import memory_tracker as memory_tracker
import threading as thre
from concurrent.futures import ThreadPoolExecutor, as_completed

loadingmethods = ['PLM Spectra', 'HDF5', 'ENVI', 'OME-TIFF', 'NetCDF', 'Zarr']
# if u have different data loading methods, feel free to add them to the list here and to the dict at the bottom. Dict name is: loadingmethodstofunctions
# I know one could just do this only with the dict (with the dict.keys()), but here it is listed one on the start. If extended do not forget to add the new method to the dict at the bottom of this file.

# start of the ''PLM Spectra' loading method, this one reads X*Y spectra multi-threaded

def loadPLMspecs(self):
    """
    Load PLM spectra from files and populate the XYMap object.
    """
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
    parallel_load_spectra(self)
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

def load_spectrum(fname, instance, lock):
    specobj = deflib.SpectrumData(
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

# end of the ''PLM Spectra' loading method
# start of the 'HDF5' loading method
def loadHDF5(self):
    """
    Load HDF5 data from files and populate the XYMap object.
    """
    # Implement HDF5 loading logic here
    pass
# end of the 'HDF5' loading method
# start of the 'ENVI' loading method
def loadENVI(self):
    """
    Load ENVI data from files and populate the XYMap object.
    """
    # Implement ENVI loading logic here
    pass
# end of the 'ENVI' loading method
# start of the 'OME-TIFF' loading method
def loadOMETIFF(self):
    """
    Load OME-TIFF data from files and populate the XYMap object.
    """
    # Implement OME-TIFF loading logic here
    pass
# end of the 'OME-TIFF' loading method
# start of the 'NetCDF' loading method
def loadNetCDF(self):
    """
    Load NetCDF data from files and populate the XYMap object.
    """
    # Implement NetCDF loading logic here
    pass
# end of the 'NetCDF' loading method
# start of the 'Zarr' loading method
def loadZarr(self):
    """
    Load Zarr data from files and populate the XYMap object.
    """
    # Implement Zarr loading logic here
    pass
# end of the 'Zarr' loading method
# if u have different data loading methods, feel free to add them to the list here and to the dict at the bottom. Dict name is: loadingmethodstofunctions

# last but not least, the dict that maps loading methods to their corresponding functions

loadingmethodstofunctions = {
    'PLM Spectra': loadPLMspecs, 
    'HDF5': ,
    'ENVI': , 
    'OME-TIFF': , 
    'NetCDF': ,
    'Zarr'
}
