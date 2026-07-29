# Performance Bottlenecks in lib9.py Fitting Routine

### 1. Per-Pixel Thread Creation Overhead (Critical)
In `fittoMatrixfitparams`, a new `threading.Thread` is created, started, and immediately `join()`'d for every single pixel in the hyperspectral image.
```python
fit_thread = thre.Thread(target=run_fit)
fit_thread.start()
fit_thread.join() # Wait for thread to complete
```
This forces totally sequential execution while introducing standard thread creation overhead (which is highly expensive). For an average 100x100 HSI, this creates and destroys 10,000 to 40,000+ threads sequentially, creating a massive artificial bottleneck.

### 2. Expensive O(N) Lookups in Inner Loops
List and dictionary lookups that remain constant are evaluated inside the double "for loop" + "while loop" logic.
```python
a = list(matl.fitkeys.keys()).index(self.selectwindowbox.get())
```
Additionally, `matl.addtofitparms.index('max_x') - len(...)` is repeatedly calculated ~12 times per pixel iteration. This should be cached once outside the loops.

### 3. Aggressive Sequential Retry Logic
If a fit fails, the code shifts the window slightly (`incmin`/`incmax`) and retries. Due to the `while tries < nmin+nmax` structure, a completely dead pixel will sequentially retry up to 40 times, each attempt creating a thread and failing, costing enormous time.

### 4. Code Duplication
The `if mode == 'fullHSI':` and `elif mode == 'roi':` branches have 160+ identical lines of code. This does not strictly affect performance execution speed, but drastically increases the complexity of optimizing the aforementioned points.

### 5. Remaining notes
The data are measured as unsigned 16-bit integers, but are converted to float64 for fitting. This is a waste of memory and CPU cycles. The data should be converted to float32 at the very least, or even better, float16 if the fitting routine can handle it. Depending on the grating the resolution is currently 2-4 pixels per nm. For the fit, the data must be converted into eV, which is a non-linear transformation. This means that the data are being oversampled by a factor of 2-4x for the fit, which is unnecessary and wasteful. The data should be resampled to 1 pixel per eV before fitting.
