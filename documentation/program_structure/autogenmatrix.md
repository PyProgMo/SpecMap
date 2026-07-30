# autogenmatrix — How Spectra Are Placed Into SpecDataMatrix

## Overview

`autogenmatrix` (in `lib9.py`) converts a flat list of `SpectrumData` objects (`self.specs`) into a 2-D grid called `self.SpecDataMatrix`.  
Each spectrum carries a physical `x-position` and `y-position` in its `data` dictionary.  
The method works out which row and column each spectrum belongs to and stores it at `self.SpecDataMatrix[row][col]` — **row = Y-index, column = X-index**.

---

## Step-by-Step Data Flow

### 1. Collect unique physical coordinates

```
_coord_tol = 1e-9   # tolerance for floating-point coordinate deduplication

for every spectrum i in self.specs:
    xpos = i.data['x-position']
    ypos = i.data['y-position']
    # add to self.mxcoords only if no existing value is within _coord_tol
    if not any(abs(xpos - cx) <= _coord_tol for cx in self.mxcoords):
        self.mxcoords.append(xpos)
    # same for self.mycoords
    if not any(abs(ypos - cy) <= _coord_tol for cy in self.mycoords):
        self.mycoords.append(ypos)

self.mxcoords = sorted(self.mxcoords)   # ascending
self.mycoords = sorted(self.mycoords)   # ascending
```

After this step `self.mxcoords` and `self.mycoords` are the **sorted lists of unique physical positions** found in the loaded files.  
Coordinates that differ only by floating-point noise (within 1e-9) are treated as identical, preventing phantom extra columns/rows.

---

### 2. Build the grid axes — `genmatgrid(xar, yar)`

`genmatgrid` constructs the regular-grid axes `self.PixAxX` and `self.PixAxY` and initialises the empty `SpecDataMatrix`.

#### Determining the grid step size

```
dxa = consecutive differences of self.mxcoords
dya = consecutive differences of self.mycoords

# Round all step values to 10 d.p. so near-identical floats collapse into one
dxa_rounded = [round(v, 10) for v in dxa]
dya_rounded = [round(v, 10) for v in dya]

if all rounded differences are equal:
    self.gdx = dxa_rounded[0]          # perfectly regular grid
else:
    self.gdx = most_frequent(dxa_rounded)   # irregular — fall back to majority step

(same logic for gdy)

self.gdx = round(self.gdx, 10)     # final guard against floating-point noise
self.gdy = round(self.gdy, 10)
```

#### Zero-step guard

```
if abs(self.gdx) < 1e-9 or abs(self.gdy) < 1e-9:
    log error via error_engine
    raise ValueError
```

If the computed step size is effectively zero (e.g., only one unique coordinate on an axis), `genmatgrid` logs an error and raises `ValueError` before any division-by-zero or infinite loop can occur.

#### Building the axis arrays

```
self.PixAxX = [round(i*gdx + x_min, 10)  for i in range(N_x)]
self.PixAxY = [round(i*gdy + y_min, 10)  for i in range(N_y)]

where N_x = int(round((x_max - x_min + gdx) / gdx))
      N_y = int(round((y_max - y_min + gdy) / gdy))
```

Each axis point is computed by multiplying the index by the step (`i * gdx + x_min`) to avoid accumulated rounding errors. The grid sizes use `int(round(...))` rather than plain `int(...)` to prevent silent truncation of values like `4.9999999…` to `4`, which would drop the rightmost/bottom column.

These are **evenly-spaced grids** that span from the minimum to the maximum observed coordinate.

#### Initialising the empty matrix

```
self.SpecDataMatrix = [[np.nan] * N_x  for _ in range(N_y)]
```

Every cell starts as `np.nan` (no spectrum assigned).

---

### 3. Populate the matrix — `SpecdataintoMatrix()`

For every spectrum `i` in `self.specs`:

```
x = i.data['x-position']
y = i.data['y-position']

xind = argmin(|self.PixAxX - x|)   # index in X-axis  (column)
yind = argmin(|self.PixAxY - y|)   # index in Y-axis  (row)

if self.SpecDataMatrix[yind][xind] is already a SpectrumData:
    if overwrite=True:
        overwrite the cell
    else:
        log a warning via error_engine.get_logger().warning(...)
        # warning includes grid index and both conflicting positions
else:
    self.SpecDataMatrix[yind][xind] = i
```

The mapping from physical coordinates to matrix indices is done by `deflib.closest_indices`, which performs a **nearest-neighbour lookup** on both axes independently.  
When a grid cell collision is detected and `overwrite=False` (the default), the event is written to the `.logfile` as a `WARNING` — no messagebox is raised and no spectrum is silently discarded without a record.

#### Index convention

| Matrix index | Axis array   | Physical direction |
|---|---|---|
| `[0]` (outer / row) | `self.PixAxY` | Y (vertical) |
| `[1]` (inner / column) | `self.PixAxX` | X (horizontal) |

So to retrieve the spectrum at physical position `(px, py)`:

```python
xind, yind = deflib.closest_indices(self.PixAxX, self.PixAxY, px, py)
spec = self.SpecDataMatrix[yind][xind]
```

---

## Bug Fixes Applied (commit 952a2cf)

The following bugs were identified and fixed in `autogenmatrix` / `genmatgrid` / `SpecdataintoMatrix`.

### Bug 1 — Floating-point coordinate mismatch  Fixed

**Scenario:** Two spectra share a physical coordinate in principle (e.g., both at `x = 10.0`) but their stored values differ by floating-point noise (e.g., `10.0` vs `10.000000001`).

**Former effect:**
- `autogenmatrix` treated them as two distinct X positions.
- The grid was over-populated with extra columns that are nearly identical in physical space.
- One spectrum ended up in a column shifted by one index relative to where it should be.

**Fix:** Replaced the exact `not in` equality check with a tolerance-based scan (`abs(xpos - cx) <= 1e-9`), so coordinates differing only by floating-point noise are treated as identical.

---

### Bug 2 — Irregular step sizes and the `most_frequent` heuristic  Fixed

**Scenario:** The scan has a non-uniform step size (e.g., most steps are 1 µm but one gap is 2 µm due to a missed measurement).

**Former effect:**
- Near-identical float steps such as `1.0` and `0.9999999999` were treated as different values.
- `most_freq_element` could pick the wrong majority step.
- Axis construction produced incorrect column positions, causing spectra to be mapped to wrong columns.

**Fix:** Step differences from `findif()` are now rounded to 10 decimal places before `remove_duplicates()` and `most_freq_element()`, so near-identical float steps collapse into one value.

---

### Bug 3 — Axis reconstruction accumulates rounding errors  Fixed

**Scenario:** `gdx` is an irrational or repeating decimal (e.g., `0.333…`).

**Former effect:**
- Each axis point accumulated a different rounding error from `i * gdx`, causing drift between the stored axis values and actual physical positions.
- `closest_indices` could map a spectrum to a neighbouring column due to this drift.

**Fix:** The existing `i * gdx + matstart` axis construction already avoids per-step accumulation. With Bug 2's step rounding now applied, axis values are consistent with the rounded step, ensuring `closest_indices` maps to the correct column.

---

### Bug 4 — Grid size under-count due to integer truncation  Fixed

**Scenario:** `(x_max - x_min + gdx) / gdx` yields something like `4.9999999…` due to floating-point arithmetic.

**Former effect:**
- `int(...)` truncated `4.9999…` to `4`, so the rightmost (or bottom) column was never created.
- The spectrum belonging there could not be placed in the matrix.
- `closest_indices` mapped it to the last existing column, colliding with the spectrum that legitimately lived there.

**Fix:** Changed to `int(round(...))` so the division result is rounded to the nearest integer before truncation.

```python
# Before — can silently drop last column
nx = int((self.matend[0] - self.matstart[0] + self.gdx) / self.gdx)

# After
nx = int(round((self.matend[0] - self.matstart[0] + self.gdx) / self.gdx))
```

---

### Bug 5 — Collision on the same grid cell (silent discard)  Fixed

**Scenario:** Two spectra have coordinates that both round to the same `(xind, yind)`.  
This can happen because of Bugs 1–4, or genuinely overlapping scan positions.

**Former effect:**
- `SpecdataintoMatrix` checked whether the cell already contained a `SpectrumData` object.
- If it did and `overwrite=False` (the default), the second spectrum was silently ignored with a bare `print()` message to the console.
- No exception was raised; the discarded spectrum was never stored anywhere and left no log entry.

**Fix:** The bare `print()` is replaced with `error_engine.get_logger().warning(...)`, so every collision is written to the `.logfile` with the grid index and both conflicting positions.

---

### Zero-step guard  Added

**Scenario:** All spectra share the same X or Y coordinate, producing zero consecutive differences and a computed `gdx` or `gdy` of zero.

**Former effect:** Division by zero or an infinite loop inside `genmatgrid`.

**Fix:** Added an `abs(gdx) < 1e-9` check immediately after step computation. If triggered, an error is logged via `error_engine` and a `ValueError` is raised before any division occurs.

---

### Bug 6 — `self.specs` is cleared before secondary calls to `SpecdataintoMatrix`

**Scenario:** At the end of `autogenmatrix`, `self.specs = []` is executed as a memory optimisation.  
If any downstream code (e.g., loading additional spectra without a full reinitialisation) calls `SpecdataintoMatrix()` a second time, `self.specs` is empty and nothing is written.

**Effect:** A silent no-op that leaves the matrix unchanged, potentially causing confusion when the user expects newly added spectra to appear.

**Root location:** End of `autogenmatrix` — `self.specs = []` after `gc.collect()`.

---

### Bug 7 — Single-spectrum edge case bypasses `genmatgrid`

**Scenario:** Exactly one spectrum is loaded.

**Effect:**
- `autogenmatrix` takes the `elif len(self.specs) == 1` branch and hardcodes  
  `self.PixAxX = [0]`, `self.PixAxY = [0]`.
- The axes no longer contain the actual physical coordinates; they contain `[0]`.
- Any code that later looks up a spectrum by physical coordinate via `closest_indices` will behave incorrectly because the axis no longer encodes real positions.

**Root location:** `autogenmatrix` `elif len(self.specs) == 1` branch.

---

## Quick Reference: Index Convention

```
self.SpecDataMatrix  shape: [N_y][N_x]

self.SpecDataMatrix[j][i]
                    │   └── column = X-index (self.PixAxX[i])
                    └────── row    = Y-index (self.PixAxY[j])

Physical → Matrix:
    i = argmin(|self.PixAxX - x_physical|)
    j = argmin(|self.PixAxY - y_physical|)
```
