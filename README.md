# SpecMap - Hyperspectral Data Analysis Tool

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Program Structure](#program-structure)
- [Core Features](#core-features)
- [Module Documentation](#module-documentation)
- [Usage Guide](#usage-guide)
- [Configuration](#configuration)
- [Advanced Documentation](#advanced-documentation)

## Overview

SpecMap is a comprehensive Python application for hyperspectral data analysis and visualization. The program provides a GUI-based interface built with Tkinter for loading, processing, analyzing, and exporting hyperspectral datasets. It supports various data formats and offers advanced features like spectral fitting, cosmic ray removal, background correction, and 2D Gaussian fitting for Clara imaging.

### Key Capabilities
- **Hyperspectral Data Processing**: Load and analyze hyperspectral image datasets
- **Spectral Analysis**: Plot, fit, and analyze individual spectra with various mathematical models
- **Image Processing**: Handle Clara imaging data with 2D Gaussian fitting capabilities
- **Data Export**: Export processed data and pixel matrices in CSV format
- **Interactive Visualization**: Click-to-analyze functionality with real-time plotting
- **Background Correction**: Multiple background correction methods including linear and cosmic ray removal
- **ROI Analysis**: Region of interest selection and analysis tools
- **Hyperspectral Standards Interface**: Convert SpecMap data to ENVI format and interface with industry-standard libraries (SPy, PySptools)


### Performance Optimizations 
SpecMap has been optimized for high-performance processing of large hyperspectral datasets:
- **30x faster** for typical workflows (Load → Process → Fit → Export)
- **40% less memory** usage for large datasets
- **Optimized derivative calculation**: 10-50x speedup using vectorized Savitzky-Golay filter
- **Efficient data handling**: Eliminated unnecessary array copies and memory allocations
- **Vectorized operations**: 200-300x speedup for spectrum averaging

For details, see [Optimization Documentation](documentation/OPTIMIZATION_README.md)

## Installation

### Requirements
The program requires the following Python packages (see `requirements.txt`):
- `tkinter` (GUI framework)
- `matplotlib` (plotting and visualization)
- `numpy` (numerical computations)
- `scipy` (scientific computing and curve fitting)
- `PIL/Pillow` (image processing)
- `pandas` (data manipulation)

### Setup guide for beginners
1. Clone the repository: `git clone [repository-url]`
2. create venv (windows cmd):
   ```
   python -m venv .venv
   .venv\Scripts\activate
   ```
3. make sure u are in the project directory and make sure all required packets are installed. Therefore run:
   ```
   pip install -r requirements.txt
   ```
4. Run the main application: `python main9.py`

### Hyperspectral Standard Interface (Optional)
SpecMap includes a Jupyter notebook for interfacing with industry-standard hyperspectral libraries:

 **[hyperspectral_standard_interface.ipynb](hyperspectral_standard_interface.ipynb)**

This notebook demonstrates:
- Converting SpecMap data to ENVI format (hyperspectral standard)
- Interfacing with Spectral Python (SPy) for advanced analysis
- Using PySptools for spectral unmixing and endmember extraction
- K-means clustering and other classification methods
- Data export in multiple formats (ENVI, NumPy, CSV)

**Quick Start Example:**
```python
python example_hyperspectral_conversion.py  # Standalone conversion script
```

**For Interactive Analysis:**
To use the full notebook, install optional dependencies:
```bash
pip install spectral pysptools  # Hyperspectral analysis libraries
pip install jupyter notebook    # If you need Jupyter
```

Then launch Jupyter:
```bash
jupyter notebook hyperspectral_standard_interface.ipynb
```

**Files:**
- `hyperspectral_standard_interface.ipynb` - Complete interactive tutorial
- `example_hyperspectral_conversion.py` - Standalone conversion script

### Troubleshooting
If you encounter TCL/Tkinter errors, you may need to:
1. Copy TCL folder from your Python installation to the .venv directory
2. Set environment variables:
   ```
   $projectdir = (Get-Location).Path
   $env:TCL_LIBRARY = "$projectdir\.venv\tcl\tcl8.6"
   $env:TK_LIBRARY  = "$projectdir\.venv\tcl\tk8.6"
   ```

## Program Structure

The application follows a modular architecture with clear separation of concerns:

```
SpecMap1/
├── main9.py              # Main application entry point
├── lib9.py               # Core hyperspectral analysis classes
├── deflib1.py            # Default settings and utility functions
├── claralib1.py          # Clara image processing module
├── mathlib3.py           # Mathematical fitting functions
├── PMclasslib1.py        # Pixel matrix and spectra data classes
├── newtonspeclib1.py     # Newton spectrum analysis
├── export1.py            # Data export functionality
├── defaults.txt          # Default configuration values
└── requirements.txt      # Python dependencies
```

## Core Features

### 1. Hyperspectral Data Loading and Processing
- **Multi-file Loading**: Process multiple spectral files from a directory
- **Format Support**: Handles various spectral data formats
- **Background Correction**: Multiple background subtraction methods
- **Cosmic Ray Removal**: Several algorithms for cosmic ray detection and removal
- **Power Correction**: Spectral power correction capabilities

### 2. Interactive Spectral Analysis
- **Pixel Selection**: Click-to-select pixels for detailed spectral analysis
- **Spectral Fitting**: Fit various mathematical models (Gaussian, Lorentzian, Voigt, etc.)
- **ROI Analysis**: Define and analyze regions of interest
- **Real-time Plotting**: Interactive matplotlib-based visualization

### 3. Clara Image Processing
- **2D Gaussian Fitting**: Advanced 2D Gaussian profile fitting for beam analysis
- **Scaling Support**: Multiple magnification factors (20x, 50x, 100x)
- **Area Calculations**: Beam area and parameter extraction

### 4. Data Export and Visualization
- **CSV Export**: Export pixel matrices and spectral data
- **Metadata Preservation**: Maintain acquisition parameters and settings
- **Interactive Colormaps**: Various colormap options for data visualization

## Module Documentation

### main9.py - Application Entry Point

**Class: FileProcessorApp**
- **Purpose**: Main application controller managing the GUI and coordinating between modules
- **Key Methods**:
  - `__init__(root, defaults)`: Initialize the main application window
  - `spec_loadfiles()`: Load and process hyperspectral data files
  - `cl_loadfiles()`: Load Clara imaging data
  - `newtonloadfiles()`: Load Newton spectrum data
  - `saveNanomap(filename)`: Save current analysis session
  - `loadhsisaved(filename)`: Load previously saved session

**GUI Structure**:
- **Tabbed Interface**: Multiple notebooks for different functionalities
  - Load Data: File loading and preprocessing options
  - Hyperspectra: Main analysis interface
  - Clara Image: Clara imaging tools
  - Export: Data export options
  - Newton Spectrum: Newton spectrum analysis
  - Settings: Configuration options

### lib9.py - Core Analysis Engine

**Class: SpectrumData**
- **Purpose**: Individual spectrum data container with metadata
- **Attributes**:
  - `WL`: Wavelength array
  - `PL`: Raw photoluminescence counts
  - `BG`: Background data
  - `PLB`: Background-corrected spectrum (PL - BG)
  - `data`: Metadata dictionary with acquisition parameters

**Class: XYMap**
- **Purpose**: Main hyperspectral dataset container and analysis engine
- **Key Features**:
  - Pixel matrix generation and management
  - Interactive plotting and pixel selection
  - Spectral fitting with multiple mathematical models
  - ROI selection and analysis
  - Colormap visualization

**Key Methods**:
- `loadfiles()`: Load and parse spectral data files
- `SpecdataintoMatrix()`: Organize spectra into spatial matrix
- `PlotSpectrum()`: Plot individual spectra
- `PlotPixMatrix()`: Generate pixel matrix visualizations
- `on_click()`: Handle mouse click events for pixel selection
- `fitSpectralData()`: Perform spectral fitting operations

**Class: Roihandler**
- **Purpose**: Region of interest selection and management
- **Features**: Polygon-based ROI definition with analysis capabilities

### deflib1.py - Utilities and Configuration

**Key Functions**:
- `initdefaults()`: Load and validate default configuration
- `save_defaults()`, `load_defaults()`: Configuration file management
- `remove_cosmics_*()`: Various cosmic ray removal algorithms
  - Linear interpolation method
  - Median filtering approach
  - Rolling mean technique
- `browse_folder()`, `select_file()`: File dialog utilities
- `power_correction()`: Spectral power correction

**Configuration System**:
- Persistent settings storage in `defaults.txt`
- Type validation and error handling
- Default value fallbacks

### mathlib3.py - Mathematical Models

**Fitting Functions**:
- `gaussianwind()`: Gaussian profile fitting
- `lorentzwind()`: Lorentzian profile fitting  
- `voigtwind()`: Voigt profile fitting (convolution of Gaussian and Lorentzian)
- `double_*()`: Double-peak fitting functions
- `linearwind()`: Linear baseline fitting

**Fit-Free Analysis**:
- **Stiffness**: Flank slope and asymmetry analysis
- **Derivative**: Savitzky-Golay based slope and inflection analysis
- **Moments**: Statistical moment analysis (COM, Variance, Skewness)
- **Decay/Rise**: Discrete derivative analysis for steepest slopes
- **Binning**: Resampling and difference analysis

**Parameter Estimation**:
- `estimate_voigt_params()`: Intelligent initial parameter estimation
- Peak detection and analysis utilities
- FWHM and area calculations

### claralib1.py - Clara Image Processing

**Class: imageprocessor**
- **Purpose**: Clara imaging data analysis and visualization
- **Features**:
  - 2D Gaussian profile fitting
  - Multi-magnification support (20x, 50x, 100x scaling)
  - Interactive plotting and analysis

**Key Functions**:
- `fit_gaussian_2d()`: Advanced 2D Gaussian fitting
- `gaussian_2d()`: 2D Gaussian mathematical model
- `plot2dfit()`: Visualization of fit results
- `area2dgaussian()`: Beam area calculations

### PMclasslib1.py - Data Structures

**Class: PMclass (Pixel Matrix Class)**
- **Purpose**: Container for 2D pixel matrix data with spatial information
- **Attributes**:
  - `PixMatrix`: 2D data array
  - `xax`, `yax`: Spatial axis arrays
  - `gdx`, `gdy`: Pixel spacing
  - `metadata`: Acquisition metadata

**Class: Spectra**
- **Purpose**: Individual spectrum data with export capabilities
- **Features**: Formatted data export with headers and metadata

### export1.py - Data Export

**Class: Exportframe**
- **Purpose**: GUI interface for data export operations
- **Features**:
  - CSV export with metadata headers
  - File dialog integration
  - Pixel matrix export functionality

### newtonspeclib1.py - Newton Spectrum Analysis

**Class: newtonspecopener**
- **Purpose**: Interface for Newton spectrum data analysis
- **Status**: Framework for future Newton-specific analysis features

## Usage Guide

### Basic Workflow

1. **Launch Application**:
   ```bash
   python main9.py
   ```

2. **Load Data**:
   - Select folder containing spectral files
   - Enter filename pattern and file extension
   - Configure preprocessing options (background correction, cosmic removal)
   - Click "Load data"

3. **Analyze Data**:
   - Use the interactive pixel matrix to select analysis points
   - View spectral plots in real-time
   - Apply mathematical fitting models
   - Define and analyze regions of interest

4. **Export Results**:
   - Navigate to Export tab
   - Choose export format and location
   - Export processed data and analysis results

### Advanced Features

**Spectral Fitting**:
- Select wavelength range for analysis
- Choose appropriate mathematical model (Gaussian, Lorentzian, Voigt)
- Adjust fit parameters and constraints
- View fit quality metrics

**ROI Analysis**:
- Define polygonal regions of interest
- Perform batch analysis on selected regions
- Generate averaged spectra from ROI areas

**Clara Image Analysis**:
- Load Clara imaging data
- Perform 2D Gaussian fitting for beam characterization
- Calculate beam parameters and areas

### Advanced Data Export Features 

#### Saving and Loading with ROI Masks
SpecMap now automatically saves and restores ROI (Region of Interest) masks when you save/load your analysis sessions. This means:
- **ROI masks persist** across sessions - no need to redraw them
- **Automatic validation** - ROI dimensions are checked against HSI data on load
- **Multiple ROI support** - save and restore any number of ROI masks

Example workflow:
```python
# Create ROI masks, then save
xymap.save_state('my_analysis.pkl')
# Output: Saved 2 ROI masks: ['bright_region', 'edge_region']

# Later, load and continue working
xymap.load_state('my_analysis.pkl')
# Output: Loaded 2 ROI masks (dimensions validated)
```

#### Spectral Averaging Options
When averaging HSI data to spectra, you can now work with multiple data types:
- **Spectrum (PL-BG)** - Background-corrected photoluminescence
- **First derivative** - Highlights spectral features and transitions
- **Second derivative** - Enhances fine spectral features
- **First derivative (normalized)** - Shape analysis independent of intensity
- **Second derivative (normalized)** - Enhanced feature detection

Use the new **"Export All Averaged Spectra"** button to:
1. Generate all spectral data types at once
2. Automatically calculate derivatives if needed
3. Export everything to a single CSV file

#### Comprehensive Data Export
New export capabilities include:
- **`exportHSIWithSpectra()`** - Export HSI matrix and all associated averaged spectra in separate files
- **`exportAllAveragedSpectra()`** - Generate and export all spectral data types in one operation
- **Metadata preservation** - All exports include wavelength ranges, acquisition parameters, and data type information

For detailed documentation, see [DATA_SAVE_LOAD_ENHANCEMENT.md](documentation/DATA_SAVE_LOAD_ENHANCEMENT.md)

## Configuration

### Default Settings
The application uses a configuration file (`defaults.txt`) to store user preferences:

- **File paths**: Default directories for data loading
- **Processing options**: Background correction, cosmic removal settings
- **Display preferences**: Font sizes, window dimensions, colormap choices
- **Analysis parameters**: Wavelength ranges, fitting constraints

### Customization
Users can modify defaults by:
1. Editing `defaults.txt` directly
2. Using the Settings tab in the application
3. Modifying the `defaults` dictionary in `deflib1.py`

### Cosmic Ray Removal Options
- **Linear Interpolation**: Fast, simple approach
- **Median Filtering**: Robust outlier removal
- **Rolling Mean**: Preserves overall spectral trends

### Background Correction
- **Multiple Backgrounds**: Use separate background for each spectrum
- **Linear Background**: Apply linear baseline correction
- **Power Correction**: Compensate for wavelength-dependent effects

## Technical Details

### Data Flow
1. **File Loading**: Raw spectral files → SpectrumData objects
2. **Matrix Organization**: SpectrumData → XYMap spatial matrix
3. **Processing**: Background correction, cosmic removal, power correction
4. **Analysis**: Interactive selection, fitting, ROI analysis
5. **Export**: Processed data → CSV files with metadata

### Memory Management
- Efficient numpy array operations
- On-demand data loading and processing
- Automatic cleanup of matplotlib figures and threads

### Error Handling
- Comprehensive exception handling throughout
- User-friendly error messages
- Graceful degradation for missing or corrupted data

## Advanced Documentation

### Mathematical Methods Documentation
For detailed mathematical descriptions of all cosmic ray removal, interpolation, and smoothing algorithms, see:

 **[COSMIC_REMOVAL_METHODS.md](COSMIC_REMOVAL_METHODS.md)**

This comprehensive document includes:
- Complete mathematical formulations for all 15+ methods
- Detection criteria and correction algorithms
- Statistical properties and computational complexity analysis
- Performance comparisons and use case recommendations
- References to academic literature

### Architecture Documentation
For detailed project structure and architecture:

 **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)**

### Software Documentation
For comprehensive developer documentation including program structure, GUI hierarchy, and detailed math library documentation:

 **[Software Documentation Summary](Software_Documentation/DOCUMENTATION_SUMMARY.md)**

Key documents include:
- **[MATHLIB_DOCUMENTATION.md](Software_Documentation/MATHLIB_DOCUMENTATION.md)**: Detailed mathematical formulations for all fit functions and analysis methods.
- **[PROGRAM_STRUCTURE.md](Software_Documentation/PROGRAM_STRUCTURE.md)**: Complete software tree and module hierarchy.
- **[GUI_HIERARCHY_AND_LIFECYCLE.md](Software_Documentation/GUI_HIERARCHY_AND_LIFECYCLE.md)**: GUI component structure and lifecycle management.

## Development and Version Control

The project uses Git for version control:
- **Main branch**: Protected, contains stable releases
- **Development**: Create feature branches for changes
- **Releases**: Tagged versions for distribution

For contributors:
1. Clone the repository: `git clone [repository-url]`
2. Create feature branch: `git checkout -b feature-name`
3. Make changes and commit: `git commit -m "description"`
4. Create pull request for review and merge

This documentation provides a complete overview of the SpecMap hyperspectral analysis application, covering all major features and implementation details.
