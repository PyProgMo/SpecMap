# PLEM HDF5 Data Format Specification

This document defines the required layout for `.h5` files used in this project (hyperspectral PL mapping: X × Y × Wavelength data cubes). Files following this spec can be read with `load_3d_dataset()` without any manual configuration.

## Overview

An HDF5 file consists of:

1. **One 3D dataset** — the data cube (e.g. PL intensity or counts).
2. **Three 1D datasets** — the physical coordinate axes: `X`, `Y`, and `Lambda` (wavelength).
3. **Attributes** on the 3D dataset that declare the axis order and link to the coordinate datasets, so the array layout is self-describing rather than assumed by convention.

Canonical in-memory axis order after loading is always `(X, Y, Lambda)`, regardless of how the array happened to be stored on disk — the `axis_order` attribute lets the reader transpose correctly.

## Required Structure

```
example.h5
├── raw/
│   ├── data          [3D dataset]   shape = (nX, nY, nLambda)  (or any permutation — see axis_order)
│   ├── x_axis         [1D dataset]   shape = (nX,)
│   ├── y_axis          [1D dataset]   shape = (nY,)
│   └── wavelength    [1D dataset]   shape = (nLambda,)
```

Paths (`raw/data`, `raw/x_axis`, etc.) are the *default* convention used by `save_3d_dataset()`, but are not hardcoded requirements — the actual paths used are always recorded in attributes on the main dataset (see below), so a reader never has to guess.

## Required Attributes on the 3D Dataset

Stored on the dataset at `dataset_path` (default `raw/data`):

| Attribute | Type | Description |
|---|---|---|
| `axis_order` | tuple/list of 3 strings | Declares what each array axis represents, in order. Must be a permutation of `"X"`, `"Y"`, `"Lambda"`. Example: `("X", "Y", "Lambda")` means `cube.shape == (nX, nY, nLambda)`. |
| `x_axis_path` | string | Internal HDF5 path to the 1D dataset holding X coordinates. |
| `y_axis_path` | string | Internal HDF5 path to the 1D dataset holding Y coordinates. |
| `wavelength_path` | string | Internal HDF5 path to the 1D dataset holding wavelength values. |
| `axis_paths` | list of 3 strings | Convenience list `[x_axis_path, y_axis_path, wavelength_path]`. |

### Consistency rules

- `cube.shape[axis_order.index("X")] == len(x_axis)`
- `cube.shape[axis_order.index("Y")] == len(y_axis)`
- `cube.shape[axis_order.index("Lambda")] == len(wavelength)`

Files that violate these will fail validation in `load_3d_dataset()`.

## Recommended Attributes on the 1D Axis Datasets

Each coordinate dataset should declare its physical units:

| Dataset | Attribute | Example value |
|---|---|---|
| `x_axis` | `units` | `"um"` |
| `y_axis` | `units` | `"um"` |
| `wavelength` | `units` | `"nm"` (or `"eV"`) |

## Optional Metadata (`extra_attrs`)

Any additional acquisition metadata can be attached directly to the 3D dataset's attributes, e.g.:

- `integration_time_s`
- `stage_step_um`
- `laser_power_mW`
- `description`

These are preserved and returned in the `meta` dict by `load_3d_dataset()`, but are not required for the file to be valid.

## Minimal Example (Python)

```python
save_3d_dataset(
    filepath="example.h5",
    cube=cube,                        # shape (nX, nY, nLambda)
    x_axis=x_axis,                    # shape (nX,), e.g. stage position in um
    y_axis=y_axis,                    # shape (nY,)
    wavelength=wavelength,            # shape (nLambda,), e.g. nm
    axis_order=("X", "Y", "Lambda"),
    extra_attrs={"integration_time_s": 1.0, "stage_step_um": 0.5},
)
```

## Loading a Compliant File

```python
cube, x_axis, y_axis, wavelength, meta = load_3d_dataset("example.h5")
# cube.shape == (nX, nY, nLambda), regardless of on-disk axis_order
```

## Validation Checklist

A file is compliant if:

- [ ] Main dataset is 3D.
- [ ] Main dataset has an `axis_order` attribute that is a permutation of `X`, `Y`, `Lambda`.
- [ ] Main dataset has `x_axis_path`, `y_axis_path`, `wavelength_path` attributes pointing to valid, existing 1D datasets.
- [ ] The three 1D datasets' lengths match the corresponding dimensions of the main dataset (per `axis_order`).
- [ ] Axis datasets ideally declare `units`.

## Notes / Future Extensions

- Multiple data cubes (e.g. `raw` and `normalized`) can coexist in one file, each with its own `axis_order`/`axis_paths` attributes and its own (or shared) coordinate datasets — just point `dataset_path` at the one you want.
- A future revision could migrate from path-based attribute linking to HDF5's native **dimension scale** API (`h5py.Dataset.dims`), which is recognized natively by tools like HDFView and xarray. Not required for project compatibility at this time.
