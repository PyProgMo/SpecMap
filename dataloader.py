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
    Load ENVI (.hdr + raw binary) hyperspectral data and populate the XYMap object.

    Requires the 'spectral' package (pip install spectral).
    self.fnames should contain the .hdr header file (or any of the paired
    .hdr/.img/.dat/.raw files — the matching .hdr will be located).
    """
    try:
        import spectral
    except ImportError as e:
        raise ImportError(
            "ENVI loading requires the 'spectral' package. Install it with: pip install spectral"
        ) from e

    filenames = self.fnames
    if len(filenames) == 0:
        raise ValueError("ENVI loading requires at least one input file.")

    # Resolve the header file path (spectral needs the .hdr, not the raw data file)
    hdr_path = None
    for fname in filenames:
        if fname.lower().endswith('.hdr'):
            hdr_path = fname
            break
    if hdr_path is None:
        candidate = os.path.splitext(filenames[0])[0] + '.hdr'
        if os.path.exists(candidate):
            hdr_path = candidate
        else:
            raise ValueError("Could not find an ENVI .hdr header file among the provided files.")

    img = spectral.open_image(hdr_path)
    cube_yxl = np.asarray(img.load(), dtype=np.float32)  # (rows=Y, cols=X, bands=Lambda)
    if cube_yxl.ndim != 3:
        raise ValueError(f"ENVI cube is not 3D (shape={cube_yxl.shape}).")

    cube = np.transpose(cube_yxl, (1, 0, 2))  # -> (X, Y, Lambda)
    x_size, y_size, band_count = cube.shape

    # Wavelength axis, if present in the header metadata
    bands = getattr(img, 'bands', None)
    if bands is not None and getattr(bands, 'centers', None):
        wavelength = np.array(bands.centers, dtype=np.float32)
    else:
        wavelength = np.arange(band_count, dtype=np.float32)

    # Pixel size from ENVI 'map info' field, if present, else fall back to index spacing
    img_metadata = getattr(img, 'metadata', {}) or {}
    map_info = img_metadata.get('map info')
    pixel_size_x = pixel_size_y = 1.0
    if map_info is not None and len(map_info) >= 7:
        try:
            pixel_size_x = float(map_info[5])
            pixel_size_y = float(map_info[6])
        except (ValueError, IndexError):
            pass
    x_axis = np.arange(x_size, dtype=np.float32) * pixel_size_x
    y_axis = np.arange(y_size, dtype=np.float32) * pixel_size_y

    metadata = dict(img_metadata)
    metadata['source_format'] = 'ENVI'
    metadata['hdr_path'] = hdr_path

    threaded_3D_array2SpectrumData(
        self, cube, wl_array=wavelength, x_axis=x_axis, y_axis=y_axis,
        metadata=metadata
    )
# end of the 'ENVI' loading method
# start of the 'OME-TIFF' loading method
def loadOMETIFF(self):
    """
    Load an OME-TIFF hyperspectral/multi-channel stack and populate the XYMap object.

    Requires the 'tifffile' package (pip install tifffile). Channel axis 'C'
    is treated as the spectral (Lambda) axis; singleton T/Z axes are squeezed out.
    """
    try:
        import tifffile
    except ImportError as e:
        raise ImportError(
            "OME-TIFF loading requires the 'tifffile' package. Install it with: pip install tifffile"
        ) from e

    filenames = self.fnames
    if len(filenames) != 1:
        raise ValueError("OME-TIFF loading requires exactly one input file.")
    filepath = filenames[0]

    with tifffile.TiffFile(filepath) as tif:
        array = np.asarray(tif.asarray())
        axes = tif.series[0].axes if tif.series else None  # e.g. 'TCZYX', 'CYX', 'ZYX'
        ome_xml = tif.ome_metadata

    # Reduce to 3D by squeezing out singleton T/Z axes
    remaining_axes = axes
    if axes is not None:
        squeeze_axes = tuple(
            i for i, ax in enumerate(axes) if ax in ('T', 'Z') and array.shape[i] == 1
        )
        if squeeze_axes:
            array = np.squeeze(array, axis=squeeze_axes)
        remaining_axes = ''.join(ax for i, ax in enumerate(axes) if i not in squeeze_axes)

    if array.ndim != 3:
        raise ValueError(
            f"OME-TIFF data is not reducible to 3D after squeezing singleton axes "
            f"(shape={array.shape}, axes={remaining_axes})."
        )

    # Determine axis order: prefer explicit axes string, else assume (C, Y, X)
    if remaining_axes and set(remaining_axes) == {'C', 'Y', 'X'}:
        order = [remaining_axes.index(ax) for ax in ('C', 'Y', 'X')]
        cyx = np.transpose(array, order)
    else:
        cyx = array  # assume already (C, Y, X)

    band_count, y_size, x_size = cyx.shape
    cube = np.transpose(cyx, (2, 1, 0))  # -> (X, Y, Lambda/C)

    # Try to pull channel wavelengths and pixel size from the OME-XML metadata
    wavelength = np.arange(band_count, dtype=np.float32)
    pixel_size_x = pixel_size_y = 1.0
    metadata = {'source_format': 'OME-TIFF'}
    if ome_xml:
        try:
            ome_dict = tifffile.xml2dict(ome_xml)
            image_node = ome_dict.get('OME', {}).get('Image', {})
            if isinstance(image_node, list):
                image_node = image_node[0]
            pixels = image_node.get('Pixels', {})
            channels = pixels.get('Channel', [])
            if isinstance(channels, dict):
                channels = [channels]
            centers = [
                float(ch['EmissionWavelength']) for ch in channels
                if isinstance(ch, dict) and 'EmissionWavelength' in ch
            ]
            if len(centers) == band_count:
                wavelength = np.array(centers, dtype=np.float32)
            pixel_size_x = float(pixels.get('PhysicalSizeX', 1.0))
            pixel_size_y = float(pixels.get('PhysicalSizeY', 1.0))
            metadata['ome_pixels'] = pixels
        except Exception as e:
            print(f"Warning: could not fully parse OME-XML metadata: {e}")

    x_axis = np.arange(x_size, dtype=np.float32) * pixel_size_x
    y_axis = np.arange(y_size, dtype=np.float32) * pixel_size_y

    threaded_3D_array2SpectrumData(
        self, cube, wl_array=wavelength, x_axis=x_axis, y_axis=y_axis,
        metadata=metadata
    )
# end of the 'OME-TIFF' loading method
# start of the 'NetCDF' loading method
def loadNetCDF(self):
    """
    Load a NetCDF hyperspectral cube and populate the XYMap object.

    Requires 'xarray' (and a NetCDF backend such as 'netCDF4').
    Install with: pip install xarray netCDF4

    Looks for a 3D data variable (set self.netcdf_variable to pick one
    explicitly). Dimensions are matched to X/Y/Lambda by name where
    possible, falling back to positional (Y, X, Lambda) order otherwise.
    """
    try:
        import xarray as xr
    except ImportError as e:
        raise ImportError(
            "NetCDF loading requires 'xarray' (and 'netCDF4'). Install with: pip install xarray netCDF4"
        ) from e

    filenames = self.fnames
    if len(filenames) != 1:
        raise ValueError("NetCDF loading requires exactly one input file.")
    filepath = filenames[0]

    requested_var = getattr(self, 'netcdf_variable', None)

    with xr.open_dataset(filepath) as ds:
        if requested_var is not None:
            if requested_var not in ds.data_vars:
                raise KeyError(f"'{requested_var}' not found in {filepath}")
            var_name = requested_var
        else:
            candidates = [name for name, da in ds.data_vars.items() if da.ndim == 3]
            if not candidates:
                raise ValueError(f"No 3D data variables found in {filepath}")
            var_name = candidates[0]
            if len(candidates) > 1:
                print(f"Warning: multiple 3D variables found: {candidates}. Using '{var_name}'.")

        da = ds[var_name]
        dims = list(da.dims)

        def _find_dim(*keywords):
            for dim in dims:
                if any(kw in dim.lower() for kw in keywords):
                    return dim
            return None

        x_dim = _find_dim('x', 'col')
        y_dim = _find_dim('y', 'row')
        lambda_dim = _find_dim('wavelength', 'lambda', 'band', 'wl', 'spectral')

        if not (x_dim and y_dim and lambda_dim) or len({x_dim, y_dim, lambda_dim}) != 3:
            if len(dims) != 3:
                raise ValueError(f"Could not resolve axis roles from dims {dims}.")
            # fall back to positional assumption: (Y, X, Lambda)
            y_dim, x_dim, lambda_dim = dims

        array = da.transpose(x_dim, y_dim, lambda_dim).values
        cube = np.asarray(array, dtype=np.float32)

        def _axis_values(dim):
            if dim in ds.coords:
                return np.asarray(ds.coords[dim].values, dtype=np.float32)
            return np.arange(da.sizes[dim], dtype=np.float32)

        x_axis = _axis_values(x_dim)
        y_axis = _axis_values(y_dim)
        wavelength = _axis_values(lambda_dim)

        metadata = dict(da.attrs)
        metadata['source_format'] = 'NetCDF'
        metadata['variable'] = var_name
        metadata['dims'] = {'x': x_dim, 'y': y_dim, 'lambda': lambda_dim}

    threaded_3D_array2SpectrumData(
        self, cube, wl_array=wavelength, x_axis=x_axis, y_axis=y_axis,
        metadata=metadata
    )
# end of the 'NetCDF' loading method
# start of the 'Zarr' loading method
def loadZarr(self):
    """
    Load a Zarr store and populate the XYMap object.

    Requires the 'zarr' package (pip install zarr). Follows the same
    axis_order / x_axis_path / y_axis_path / wavelength_path attribute
    convention documented in h5_format_spec.md for the project's HDF5
    format — a Zarr store written with that layout (array attrs instead
    of h5py dataset attrs) will load with no extra configuration.
    """
    try:
        import zarr
    except ImportError as e:
        raise ImportError(
            "Zarr loading requires the 'zarr' package. Install it with: pip install zarr"
        ) from e

    def _decode(value):
        if isinstance(value, bytes):
            return value.decode()
        if isinstance(value, (list, tuple)):
            return [_decode(item) for item in value]
        return value

    filenames = self.fnames
    if len(filenames) != 1:
        raise ValueError("Zarr loading requires exactly one input path (a .zarr store/directory).")
    filepath = filenames[0]
    requested_path = getattr(self, 'zarr_dataset_path', None)

    root = zarr.open(filepath, mode='r')

    dataset_path = requested_path
    if dataset_path is None and 'raw/data' in root:
        dataset_path = 'raw/data'

    if dataset_path is not None:
        if dataset_path not in root:
            raise KeyError(f"'{dataset_path}' not found in {filepath}")
        arr = root[dataset_path]
    else:
        candidates = []

        def _visit(name, obj):
            if hasattr(obj, 'ndim') and obj.ndim == 3:
                candidates.append(name)

        root.visititems(_visit)
        if not candidates:
            raise ValueError(f"No 3D arrays found in {filepath}")
        dataset_path = candidates[0]
        arr = root[dataset_path]
        if len(candidates) > 1:
            print(f"Warning: multiple 3D arrays found: {candidates}. Using '{dataset_path}'.")

    if arr.ndim != 3:
        raise ValueError(f"'{dataset_path}' is not 3D (shape={arr.shape})")

    array3D = np.asarray(arr[:])
    attrs = {key: _decode(value) for key, value in dict(arr.attrs).items()}

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
        x_arr = root[x_path]
        y_arr = root[y_path]
        wl_arr = root[wavelength_path]
        x_axis = np.asarray(x_arr[:], dtype=np.float32)
        y_axis = np.asarray(y_arr[:], dtype=np.float32)
        wavelength = np.asarray(wl_arr[:], dtype=np.float32)
        metadata = dict(attrs)
        metadata.update({
            'x_units': _decode(dict(x_arr.attrs).get('units')),
            'y_units': _decode(dict(y_arr.attrs).get('units')),
            'wavelength_units': _decode(dict(wl_arr.attrs).get('units')),
        })
    metadata['dataset_path'] = dataset_path
    metadata['axis_order'] = axis_order
    metadata['source_format'] = 'Zarr'

    threaded_3D_array2SpectrumData(
        self, cube, wl_array=wavelength, x_axis=x_axis, y_axis=y_axis,
        metadata=metadata
    )
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