import os, gc
from typing import cast
import numpy as np
import memory_tracker as memory_tracker
import threading as thre
from concurrent.futures import ThreadPoolExecutor, as_completed

loadingmethods = ['PLM Spectra', 'HDF5', 'ENVI', 'OME-TIFF', 'NetCDF', 'Zarr']
import deflib1 as deflib
import mathlib3 as matl
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

def threaded_3D_array2SpectrumData(
    self, array3D, wl_array=None, x_axis=None, y_axis=None, metadata=None
):
    """
    Convert a hyperspectral cube to SpectrumData objects using multithreading.

    The default cube layout is (Y, X, wavelength). When x_axis and y_axis
    are supplied, the cube layout is canonical (X, Y, wavelength).
    """
    array3D = np.asarray(array3D)
    if array3D.ndim != 3:
        raise ValueError("Input array must be 3D.")

    if x_axis is not None or y_axis is not None:
        if x_axis is None or y_axis is None:
            raise ValueError("x_axis and y_axis must be supplied together.")
        cube = array3D
        x_axis = np.asarray(x_axis, dtype=np.float32)
        y_axis = np.asarray(y_axis, dtype=np.float32)
        if cube.shape[:2] != (len(x_axis), len(y_axis)):
            raise ValueError("X and Y axis lengths must match the cube shape.")
        x_size, y_size, band_count = cube.shape
    else:
        spectral_axis = getattr(self, 'spectral_axis', -1)
        if spectral_axis not in (-1, 0, 2):
            raise ValueError("spectral_axis must be 0 or 2 for a 3D cube.")
        cube = np.moveaxis(array3D, spectral_axis, -1)
        cube = np.transpose(cube, (1, 0, 2))
        x_size, y_size, band_count = cube.shape
        x_axis = np.arange(x_size, dtype=np.float32)
        y_axis = np.arange(y_size, dtype=np.float32)

    wavelength = np.asarray(
        getattr(self, 'WL', []) if wl_array is None else wl_array,
        dtype=np.float32,
    )
    if wavelength.size == 0:
        wavelength = np.arange(band_count, dtype=np.float32)
    if wavelength.ndim != 1 or wavelength.size != band_count:
        raise ValueError(
            f"WL must contain one value per spectral band ({band_count}); "
            f"got shape {wavelength.shape}."
        )
    self.WL = wavelength
    self.WL_eV = deflib.wl_array_to_ev(wavelength.copy())
    if metadata is not None:
        self.HDF5metadata = metadata
    background = np.asarray(getattr(self, 'BG', np.zeros(band_count)), dtype=np.float32)
    if background.size == 0:
        background = np.zeros(band_count, dtype=np.float32)
    if background.ndim != 1 or background.size != band_count:
        raise ValueError("BG must contain one value per spectral band.")
    self.BG = background

    futures = {}
    with ThreadPoolExecutor(max_workers=min(32, (os.cpu_count() or 4) + 4)) as executor:
        for x in range(x_size):
            for y in range(y_size):
                future = executor.submit(
                    _spectrum_data_from_array,
                    cube[x, y, :], wavelength, background,
                    x_axis[x], y_axis[y], self
                )
                futures[future] = (x, y)

        ordered_specs = [None] * (x_size * y_size)
        for future in as_completed(futures):
            x, y = futures[future]
            specobj = future.result()
            if specobj.dataokay:
                ordered_specs[x * y_size + y] = specobj

    self.specs = [specobj for specobj in ordered_specs if specobj is not None]


def _spectrum_data_from_array(values, wavelength, background, x, y, instance):
    """Build one SpectrumData object without routing array data through a file parser."""
    values = np.asarray(values, dtype=np.float32)
    specobj = deflib.SpectrumData.__new__(deflib.SpectrumData)
    specobj.removecosmicsmethod = getattr(instance, 'remcosmicfunc', 'median')
    specobj.loadeachbg = getattr(instance, 'loadeachbg', False)
    specobj.linearbg = getattr(instance, 'linearbg', False)
    specobj.removecosmics = getattr(instance, 'removecosmics', False)
    specobj.cosmicthreshold = getattr(instance, 'cosmicthreshold', 20)
    specobj.cosmicpixels = getattr(instance, 'cosmicpixels', 3)
    specobj.WL = wavelength
    specobj.WL_eV = getattr(instance, 'WL_eV', None)
    specobj.filename = f"array[{y},{x}]"
    specobj.default_dataset = 'Spectrum (PL-BG)'
    specobj.dataokay = True
    specobj.data = {
        'x-position': float(x),
        'y-position': float(y),
        'Delta Wavelength (nm)': float(np.mean(np.diff(wavelength))) if len(wavelength) > 1 else 0.0,
    }
    specobj.roistore = {}
    specobj.PL = values + background
    specobj.BG = background
    specobj.PLB = values
    specobj.Specdiff1 = None
    specobj.Specdiff2 = None
    specobj.dofit = False
    specobj.fwhm = np.nan
    specobj.fitmaxX = np.nan
    specobj.fitmaxY = np.nan
    specobj.fitdata = [None]
    specobj.fitparams = matl.buildfitparas()
    specobj.fitparamunits = matl.buildfitparas()
    return specobj

# end of the ''PLM Spectra' loading method
# start of the 'HDF5' loading method
def loadHDF5(self):
    import h5py
    """
    Load HDF5 data from files and populate the XYMap object.
    """
    def _decode(value):
        if isinstance(value, bytes):
            return value.decode()
        if isinstance(value, np.ndarray):
            return [_decode(item) for item in value.tolist()]
        return value

    filenames = self.fnames
    if len(filenames) != 1:
        raise ValueError("HDF5 loading requires exactly one input file.")
    filepath = filenames[0]
    requested_path = getattr(self, 'hdf5_dataset_path', None)

    with h5py.File(filepath, "r") as f:
        dataset_path = requested_path
        if dataset_path is None and 'raw/data' in f:
            dataset_path = 'raw/data'
        if dataset_path is not None:
            if dataset_path not in f:
                raise KeyError(f"'{dataset_path}' not found in {filepath}")
            dset = f[dataset_path]
            if not isinstance(dset, h5py.Dataset):
                raise ValueError(f"'{dataset_path}' is not a dataset")
        else:
            candidates = []
            f.visititems(
                lambda name, obj: candidates.append(name)
                if isinstance(obj, h5py.Dataset) and obj.ndim == 3 else None
            )
            if not candidates:
                raise ValueError(f"No 3D datasets found in {filepath}")
            dataset_path = candidates[0]
            dset = f[dataset_path]
            if len(candidates) > 1:
                print(f"Warning: multiple 3D datasets found: {candidates}. Using '{dataset_path}'.")

        if not isinstance(dset, h5py.Dataset):
            raise ValueError(f"'{dataset_path}' is not a dataset")
        if dset.ndim != 3:
            raise ValueError(f"'{dataset_path}' is not 3D (shape={dset.shape})")
        array3D = dset[()]
        attrs = {key: _decode(value) for key, value in dset.attrs.items()}

        axis_order = tuple(attrs.get('axis_order', ('Y', 'X', 'Lambda')))
        if set(axis_order) != {'X', 'Y', 'Lambda'}:
            raise ValueError(f"Invalid axis_order in '{dataset_path}': {axis_order}")
        permutation = [axis_order.index(label) for label in ('X', 'Y', 'Lambda')]
        cube = np.transpose(array3D, axes=permutation)

        x_path = attrs.get('x_axis_path')
        y_path = attrs.get('y_axis_path')
        wavelength_path = attrs.get('wavelength_path')
        if not all((x_path, y_path, wavelength_path)):
            x_axis = np.arange(cube.shape[0], dtype=np.float32)
            y_axis = np.arange(cube.shape[1], dtype=np.float32)
            wavelength = np.asarray(attrs.get('wavelength', np.arange(cube.shape[2])), dtype=np.float32)
            metadata = attrs
        else:
            x_obj = f[x_path]
            y_obj = f[y_path]
            wavelength_obj = f[wavelength_path]
            if not isinstance(x_obj, h5py.Dataset):
                raise ValueError(f"'{x_path}' is not a dataset.")
            if not isinstance(y_obj, h5py.Dataset):
                raise ValueError(f"'{y_path}' is not a dataset.")
            if not isinstance(wavelength_obj, h5py.Dataset):
                raise ValueError("HDF5 axis paths must reference datasets.")
            x_dset = cast(h5py.Dataset, x_obj)
            y_dset = cast(h5py.Dataset, y_obj)
            wavelength_dset = cast(h5py.Dataset, wavelength_obj)
            x_axis = np.asarray(x_dset[()], dtype=np.float32)
            y_axis = np.asarray(y_dset[()], dtype=np.float32)
            wavelength = np.asarray(wavelength_dset[()], dtype=np.float32)
            metadata = dict(attrs)
            metadata.update({
                'x_units': _decode(x_dset.attrs.get('units')),
                'y_units': _decode(y_dset.attrs.get('units')),
                'wavelength_units': _decode(wavelength_dset.attrs.get('units')),
            })
        metadata['dataset_path'] = dataset_path
        metadata['axis_order'] = axis_order

    threaded_3D_array2SpectrumData(
        self, cube, wl_array=wavelength, x_axis=x_axis, y_axis=y_axis,
        metadata=metadata
    )

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
    'HDF5': loadHDF5,
    'ENVI': loadENVI, 
    'OME-TIFF': loadOMETIFF, 
    'NetCDF': loadNetCDF,
    'Zarr': loadZarr
}
