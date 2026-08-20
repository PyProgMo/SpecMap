# Dataloader Format Documentation

This document describes every loading method supported by `dataloader.py` (`loadingmethodstofunctions`), what input(s) each one expects, which Python packages are required, and how the data is mapped into the project's canonical `(X, Y, Lambda)` layout before being converted into `SpectrumData` objects.

All loaders (except `PLM Spectra`) ultimately call `threaded_3D_array2SpectrumData(self, cube, wl_array, x_axis, y_axis, metadata)`, so regardless of source format, the result is the same: one `SpectrumData` object per `(x, y)` pixel, each holding the full spectrum along `Lambda`.

---

## Overview

| Method key | Function | Required package(s) | Input |
|---|---|---|---|
| `'PLM Spectra'` | `loadPLMspecs` | — | Multiple text files (one spectrum each) |
| `'HDF5'` | `loadHDF5` | `h5py` | One `.h5` file |
| `'ENVI'` | `loadENVI` | `spectral` | `.hdr` + raw binary pair |
| `'OME-TIFF'` | `loadOMETIFF` | `tifffile` | One `.ome.tiff` / `.ome.tif` file |
| `'NetCDF'` | `loadNetCDF` | `xarray`, `netCDF4` | One `.nc` file |
| `'Zarr'` | `loadZarr` | `zarr` | One `.zarr` store (directory) |

Set `self.fnames` to the input path(s) before calling the loader, exactly as with the existing `PLM Spectra` / `HDF5` methods.

---

## PLM Spectra

**Function:** `loadPLMspecs` / `parallel_load_spectra`
**Input:** `self.fnames` — a list of plain-text spectrum files, one file per `(x, y)` position.

Each file is expected to contain a header section, then tab-separated data lines of `wavelength \t counts`. The wavelength column is read once (from the first usable file) and assumed identical across all files. Loaded multi-threaded via `ThreadPoolExecutor`. This is the original/native project format and is unrelated to the array-based formats below.

---

## HDF5

**Function:** `loadHDF5`
**Package:** `h5py` (`pip install h5py`)
**Input:** `self.fnames` — exactly one `.h5` file.

Reads a 3D dataset plus 1D `X`, `Y`, and wavelength coordinate datasets, following the project's own HDF5 format spec (see `h5_format_spec.md`). Required structure:

```
example.h5
└── raw/
    ├── data          [3D dataset]  shape = (nX, nY, nLambda)  (or any permutation)
    ├── x_axis        [1D dataset]  shape = (nX,)
    ├── y_axis        [1D dataset]  shape = (nY,)
    └── wavelength    [1D dataset]  shape = (nLambda,)
```

**Required attributes** on the main dataset:
- `axis_order` — permutation of `("X", "Y", "Lambda")` describing the array's on-disk axis order
- `x_axis_path`, `y_axis_path`, `wavelength_path` — internal HDF5 paths to the 1D coordinate datasets

**Dataset selection:**
- If `self.hdf5_dataset_path` is set, that exact path is used.
- Else if a dataset exists at `raw/data`, it is used.
- Else the file is scanned for any 3D dataset; if exactly one is found it is used, if several are found the first is used (with a warning printed), if none are found an error is raised.

**Fallback behavior:** if `x_axis_path`/`y_axis_path`/`wavelength_path` attributes are missing, the loader falls back to plain index arrays (`0, 1, 2, ...`) for X/Y and wavelength.

See `h5_format_spec.md` for the full specification, a validation checklist, and the matching `save_3d_dataset()`/`load_3d_dataset()` reference functions.

---

## ENVI

**Function:** `loadENVI`
**Package:** `spectral` (`pip install spectral`)
**Input:** `self.fnames` — must include the ENVI header file (`.hdr`), typically alongside its paired raw binary data file (`.img`/`.dat`/`.raw`, same base name).

If `self.fnames` doesn't directly contain a `.hdr` path, the loader looks for `<first_filename_without_extension>.hdr` next to the first provided file.

Data is read via `spectral.open_image(hdr_path).load()`, which returns a `(rows, cols, bands)` = `(Y, X, Lambda)` array; this is transposed to the canonical `(X, Y, Lambda)`.

- **Wavelength axis:** taken from the header's band-center metadata (`img.bands.centers`) if present, otherwise a plain index array.
- **Spatial axes:** if the header's `map info` field is present, pixel size (fields 6 and 7) is used to scale `x_axis`/`y_axis`; otherwise unit index spacing is used.
- All other ENVI header fields are passed through as `metadata`.

---

## OME-TIFF

**Function:** `loadOMETIFF`
**Package:** `tifffile` (`pip install tifffile`)
**Input:** `self.fnames` — exactly one `.ome.tiff` / `.ome.tif` file.

The image is read via `tifffile.TiffFile(...).asarray()`. Any singleton `T` (time) or `Z` (depth) axes are squeezed out; the remaining `C` (channel) axis is treated as the spectral/`Lambda` axis, and `Y`/`X` as the spatial axes. If the axes cannot be reduced to exactly 3 dimensions (e.g. multiple time points or Z-slices present), loading fails with an explicit error — such stacks need to be split or reduced (e.g. max-projected) before loading.

- **Wavelength axis:** if OME-XML metadata is present and each `<Channel>` element has an `EmissionWavelength`, those values are used (must match the channel count exactly); otherwise a plain index array.
- **Spatial axes:** `PhysicalSizeX`/`PhysicalSizeY` from the OME-XML `<Pixels>` element are used to scale `x_axis`/`y_axis` if present, otherwise unit index spacing is used.
- Any OME-XML parsing failure is caught and reported as a warning; loading continues with index-based fallbacks.

---

## NetCDF

**Function:** `loadNetCDF`
**Packages:** `xarray`, `netCDF4` (`pip install xarray netCDF4`)
**Input:** `self.fnames` — exactly one `.nc` file.

Opened via `xarray.open_dataset()`.

- **Variable selection:** set `self.netcdf_variable` to pick an explicit data variable; otherwise the file is scanned for 3D variables (first match used, with a warning if several exist).
- **Axis role matching:** dimension names are matched case-insensitively against keyword sets — `x`/`col` → X, `y`/`row` → Y, `wavelength`/`lambda`/`band`/`wl`/`spectral` → Lambda. If this fails to uniquely resolve all three roles, the loader falls back to positional order `(Y, X, Lambda)`.
- **Coordinate values:** taken from the matching `xarray` coordinate variable if one exists (e.g. an actual `wavelength(wavelength)` coord with real values/units); otherwise a plain index array.
- Variable-level attributes (`da.attrs`) are passed through as `metadata`, along with which dimension was mapped to which role.

---

## Zarr

**Function:** `loadZarr`
**Package:** `zarr` (`pip install zarr`)
**Input:** `self.fnames` — exactly one path to a `.zarr` store (a directory, or any Zarr-supported store).

Follows **the same layout and attribute convention as the project's HDF5 format** (see the HDF5 section above and `h5_format_spec.md`), just stored via Zarr's group/array/attrs model instead of h5py:

```
example.zarr/
└── raw/
    ├── data          [3D array]  shape = (nX, nY, nLambda)  (or any permutation)
    ├── x_axis        [1D array]  shape = (nX,)
    ├── y_axis        [1D array]  shape = (nY,)
    └── wavelength    [1D array]  shape = (nLambda,)
```

Same required attributes on the main array (`axis_order`, `x_axis_path`, `y_axis_path`, `wavelength_path`), same dataset-selection logic (`self.zarr_dataset_path` → `raw/data` → auto-scan for a single 3D array), and the same index-array fallback if coordinate attributes are missing. A file written to satisfy the HDF5 spec can be re-created as a Zarr store with no format changes beyond the container type.

---

## Adding a new format

To add another loader:

1. Write a `loadXxx(self)` function that reads `self.fnames` and ends up with:
   - `cube` — a 3D array, canonical `(X, Y, Lambda)` order (or pass raw order + set `self.spectral_axis` and let `threaded_3D_array2SpectrumData` handle it — see its docstring)
   - `wavelength` — 1D array matching the Lambda dimension
   - `x_axis`, `y_axis` — 1D arrays matching the X/Y dimensions
   - `metadata` — a dict of anything worth keeping (units, acquisition settings, source format tag, etc.)
2. Call `threaded_3D_array2SpectrumData(self, cube, wl_array=wavelength, x_axis=x_axis, y_axis=y_axis, metadata=metadata)` at the end.
3. Add `'Xxx'` to the `loadingmethods` list at the top of the file and to `loadingmethodstofunctions` at the bottom.