import numpy as np
from scipy.ndimage import median_filter
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
from scipy.ndimage import median_filter
import pandas as pd
from scipy.interpolate import UnivariateSpline
import os, copy
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

Notebooks = ['Load Data', 'Hyperspectra', 'Clara Image', 'HSI Plot', 'Newton Spectrum', 'TCSPC', 'HSI File Sorter', 'Plot Spectra', 'Cube2Image', 'Settings', 'Marwins Spekplotter']
# dir of this file
DIR = os.path.dirname(os.path.abspath(__file__))
# Default values for the application
DEFAULTS_FILE = f'{DIR}/defaults.txt'

# load defaults
def initdefaults():
    loadeddefaults = load_defaults()
    reqdefaults = defaults.copy()
    for i in loadeddefaults.keys():
        try:
            reqdefaults[i] = defaulttypes[i](loadeddefaults[i])
        except Exception as Error:
            if loadeddefaults[i] == 'None':
                pass
            else:
                print(f'Error: {Error} while loading Entries. Using default value: {i}={reqdefaults[i]}')
    return reqdefaults

def save_defaults(variables):
    """Save default values to a file."""
    with open(DEFAULTS_FILE, 'w') as file:
        for name, value in variables.items():
            file.write(f'{name} = {value}\n')

def load_defaults():
    """Load default values from a file."""
    variables = {}
    if os.path.exists(DEFAULTS_FILE):
        with open(DEFAULTS_FILE, 'r') as file:
            for line in file:
                # Split each line into name and value
                if '=' in line:
                    name, value = line.split('=', 1)
                    name = name.strip()
                    value = value.strip()
                    # Handle basic types
                    if value.isdigit():
                        value = int(value)
                    elif value.replace('.', '', 1).isdigit():
                        value = float(value)
                    elif value.lower() in ('true', 'false'):
                        value = value.lower() == 'true'
                    variables[name] = value
        # close the file
        del file
    return variables

# correlation methods for read data interpolation

def matrix_image_correction_Matrix(SpectrumDataMatrix, thresh, width):
    SpectrumDataMatrix_correlated = copy.copy(SpectrumDataMatrix)

    # Apply matrix image correction to the SpectrumData's XYMap data
    if SpectrumDataMatrix is not None:
        # if len(SpectrumData) < 3 or len(SpectrumData[0]) < 3: just pass
        if len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
            return SpectrumDataMatrix

        correctionweighting = gaussian_weight_matrix(dx=1.0, dy=1.0, sigma_x=1.0, sigma_y=1.0, size=3)
        print('Correctionweighting: ', correctionweighting)
        #return SpectrumDataMatrix_correlated # first testing

        # iterate over each spectrum in the XYMap and identify cosmics
        for i in range(len(SpectrumDataMatrix)-2):
            for j in range(len(SpectrumDataMatrix[i])-2):
                if hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                    existingmatrix = [[0,0,0],[0,0,0],[0,0,0]]
                    # check if surrounding pixels exist
                    for k in range(len(SpectrumDataMatrix[i][j].PLB)):
                        # fill existingmatrix
                        # iterate over surrounding pixels
                        for l in range(-1, 2):
                            for m in range(-1, 2):
                                # check if surrounding pixels exist
                                if SpectrumDataMatrix[i+l][j+m] is not None:
                                    existingmatrix[l+1][m+1] = 1
                    # if more than 80% of surrounding pixels exist, apply correction
                    if np.sum(existingmatrix) <= 0.8 * np.prod(np.array(existingmatrix).shape):
                        break
                        

                    for k in range(len(SpectrumDataMatrix[i][j].PLB)):
                        # correction, weight each pixel with the average of the surrounding pixels
                        local_matrix = np.array([[SpectrumDataMatrix[i-1][j-1].PLB[k], SpectrumDataMatrix[i-1][j].PLB[k], SpectrumDataMatrix[i-1][j+1].PLB[k]],
                                                  [SpectrumDataMatrix[i][j-1].PLB[k], SpectrumDataMatrix[i][j].PLB[k], SpectrumDataMatrix[i][j+1].PLB[k]],
                                                  [SpectrumDataMatrix[i+1][j-1].PLB[k], SpectrumDataMatrix[i+1][j].PLB[k], SpectrumDataMatrix[i+1][j+1].PLB[k]]])
                        # Apply a simple average filter
                        SpectrumDataMatrix_correlated[i][j].PLB[k] = np.mean(local_matrix)
        
    else:
        return SpectrumDataMatrix

    return SpectrumDataMatrix_correlated

def cosmic_correlation_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Identify and correct cosmic rays in a 2D spectrum matrix using Gaussian-weighted interpolation.
    
    This function detects pixels that differ from their surroundings more than expected 
    based on a Gaussian laser profile, then interpolates them using Gaussian-weighted 
    neighboring pixels.
    
    Parameters
    ----------
    SpectrumDataMatrix : 2D array of SpectrumData objects
        Matrix containing spectral data at each spatial position
    thresh : float
        Threshold multiplier for cosmic detection (e.g., 3.0 for 3-sigma)
    width : int
        Not used in this implementation, kept for compatibility
        
    Returns
    -------
    SpectrumDataMatrix_corrected : 2D array of SpectrumData objects
        Corrected spectrum matrix with cosmics interpolated
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    # Check if matrix is large enough
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    # Generate 3x3 Gaussian weighting kernel (normalized)
    gaussian_weights, _, _ = gaussian_weight_matrix(dx=1.0, dy=1.0, sigma_x=1.0, sigma_y=1.0, size=3)
    
    # Pre-compute weights without center
    weights_no_center = gaussian_weights.copy()
    weights_no_center[1, 1] = 0
    
    # Pre-compute neighbor mask (all except center)
    neighbor_mask = np.ones((3, 3), dtype=bool)
    neighbor_mask[1, 1] = False
    
    # Iterate over interior pixels (avoid edges)
    for i in range(1, len(SpectrumDataMatrix) - 1):
        for j in range(1, len(SpectrumDataMatrix[i]) - 1):
            # Check if center pixel has spectral data
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue
            
            # Check which surrounding pixels exist and build spectral cube
            valid_neighbors = np.ones((3, 3), dtype=bool)
            spectral_cube = []  # Will hold 3x3 arrays, one per wavelength
            
            # Build the spectral cube for all wavelengths at once
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    neighbor = SpectrumDataMatrix[i+di][j+dj]
                    if neighbor is None or not hasattr(neighbor, 'PLB'):
                        valid_neighbors[di+1, dj+1] = False
            
            # Need at least 80% of neighbors to be valid
            if np.sum(valid_neighbors) < 0.8 * 9:
                continue
            
            # Extract all spectra into a 3D array [3, 3, n_wavelengths]
            n_wl = len(SpectrumDataMatrix[i][j].PLB)
            spectral_cube = np.zeros((3, 3, n_wl))
            
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    if valid_neighbors[di+1, dj+1]:
                        spectral_cube[di+1, dj+1, :] = SpectrumDataMatrix[i+di][j+dj].PLB
            
            # Vectorized operations over all wavelengths
            center_spectrum = spectral_cube[1, 1, :]  # Shape: (n_wl,)
            
            # Calculate valid weights for neighbors (excluding center)
            valid_weights_neighbors = weights_no_center * valid_neighbors
            weight_sum = np.sum(valid_weights_neighbors)
            
            if weight_sum > 0:
                valid_weights_neighbors = valid_weights_neighbors / weight_sum
                
                # Expected values for all wavelengths at once: shape (n_wl,)
                expected_spectrum = np.sum(spectral_cube * valid_weights_neighbors[:, :, np.newaxis], axis=(0, 1))
                
                # Calculate std of neighbors across all wavelengths
                # Extract only valid neighbor values
                neighbor_valid_mask = neighbor_mask & valid_neighbors  # Exclude center and invalid
                neighbor_values = spectral_cube[neighbor_valid_mask, :]  # Shape: (n_neighbors, n_wl)
                
                if neighbor_values.shape[0] > 0:
                    neighbor_std = np.std(neighbor_values, axis=0)  # Shape: (n_wl,)
                    
                    # Detect cosmics: vectorized deviation check
                    deviations = np.abs(center_spectrum - expected_spectrum)
                    cosmic_mask = deviations > (thresh * neighbor_std)  # Boolean array
                    
                    # If any cosmics detected, interpolate them
                    if np.any(cosmic_mask):
                        # Calculate weights for all valid neighbors (including center position)
                        weights_all = gaussian_weights * valid_neighbors
                        weight_sum_all = np.sum(weights_all)
                        
                        if weight_sum_all > 0:
                            weights_all = weights_all / weight_sum_all
                            
                            # Corrected values for all wavelengths
                            corrected_spectrum = np.sum(spectral_cube * weights_all[:, :, np.newaxis], axis=(0, 1))
                            
                            # Replace only the cosmic-affected wavelengths
                            # Convert PLB to numpy array if it's a list, apply correction, convert back
                            plb_array = np.array(SpectrumDataMatrix_corrected[i][j].PLB)
                            plb_array[cosmic_mask] = corrected_spectrum[cosmic_mask]
                            SpectrumDataMatrix_corrected[i][j].PLB = plb_array.tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else plb_array
    
    return SpectrumDataMatrix_corrected

def adaptive_threshold_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Adaptive threshold cosmic removal using local spatial statistics.
    Calculates threshold based on local mean and std of surrounding pixels.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    # Iterate over interior pixels
    for i in range(1, len(SpectrumDataMatrix) - 1):
        for j in range(1, len(SpectrumDataMatrix[i]) - 1):
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue
            
            # Collect valid neighbor spectra
            neighbor_spectra = []
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    if (di != 0 or dj != 0):
                        neighbor = SpectrumDataMatrix[i+di][j+dj]
                        if neighbor is not None and hasattr(neighbor, 'PLB'):
                            neighbor_spectra.append(np.array(neighbor.PLB))
            
            if len(neighbor_spectra) < 4:  # Need sufficient neighbors
                continue
            
            # Stack neighbors: shape (n_neighbors, n_wavelengths)
            neighbor_stack = np.array(neighbor_spectra)
            
            # Calculate local statistics per wavelength
            local_mean = np.mean(neighbor_stack, axis=0)
            local_std = np.std(neighbor_stack, axis=0)
            
            # Get center spectrum
            center_spectrum = np.array(SpectrumDataMatrix[i][j].PLB)
            
            # Adaptive threshold based on local statistics
            deviation = np.abs(center_spectrum - local_mean)
            cosmic_mask = deviation > (thresh * local_std)
            
            # Replace cosmics with local mean
            if np.any(cosmic_mask):
                plb_array = center_spectrum.copy()
                plb_array[cosmic_mask] = local_mean[cosmic_mask]
                SpectrumDataMatrix_corrected[i][j].PLB = plb_array.tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else plb_array
    
    return SpectrumDataMatrix_corrected

def spectral_correlation_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Spectral correlation method - detects cosmics by checking spectral smoothness.
    Uses both spatial neighbors AND spectral continuity.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    for i in range(1, len(SpectrumDataMatrix) - 1):
        for j in range(1, len(SpectrumDataMatrix[i]) - 1):
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue
            
            center_spectrum = np.array(SpectrumDataMatrix[i][j].PLB)
            n_wl = len(center_spectrum)
            
            if n_wl < 3:
                continue
            
            # Calculate spectral second derivative (curvature)
            spectral_deriv2 = np.zeros(n_wl)
            spectral_deriv2[1:-1] = center_spectrum[:-2] - 2*center_spectrum[1:-1] + center_spectrum[2:]
            
            # Also check spatial consistency
            neighbor_spectra = []
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    if (di != 0 or dj != 0):
                        neighbor = SpectrumDataMatrix[i+di][j+dj]
                        if neighbor is not None and hasattr(neighbor, 'PLB'):
                            neighbor_spectra.append(np.array(neighbor.PLB))
            
            if len(neighbor_spectra) >= 4:
                neighbor_median = np.median(neighbor_spectra, axis=0)
                
                # Combined detection: spectral discontinuity AND spatial deviation
                std_deriv = np.std(spectral_deriv2)
                spatial_deviation = np.abs(center_spectrum - neighbor_median)
                spatial_std = np.std(spatial_deviation)
                
                spectral_anomaly = np.abs(spectral_deriv2) > thresh * std_deriv
                spatial_anomaly = spatial_deviation > thresh * spatial_std
                
                # Cosmic if BOTH spectral and spatial anomalies
                cosmic_mask = spectral_anomaly & spatial_anomaly
                
                if np.any(cosmic_mask):
                    plb_array = center_spectrum.copy()
                    plb_array[cosmic_mask] = neighbor_median[cosmic_mask]
                    SpectrumDataMatrix_corrected[i][j].PLB = plb_array.tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else plb_array
    
    return SpectrumDataMatrix_corrected

def bilateral_filter_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Bilateral filter - weights neighbors by both spatial distance AND spectral similarity.
    Preserves edges while removing cosmics.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    # Generate spatial Gaussian weights
    spatial_weights, _, _ = gaussian_weight_matrix(dx=1.0, dy=1.0, sigma_x=1.0, sigma_y=1.0, size=3)
    
    for i in range(1, len(SpectrumDataMatrix) - 1):
        for j in range(1, len(SpectrumDataMatrix[i]) - 1):
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue
            
            center_spectrum = np.array(SpectrumDataMatrix[i][j].PLB)
            n_wl = len(center_spectrum)
            
            # Build spectral cube
            spectral_cube = np.zeros((3, 3, n_wl))
            valid_mask = np.zeros((3, 3), dtype=bool)
            
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    neighbor = SpectrumDataMatrix[i+di][j+dj]
                    if neighbor is not None and hasattr(neighbor, 'PLB'):
                        spectral_cube[di+1, dj+1, :] = neighbor.PLB
                        valid_mask[di+1, dj+1] = True
            
            if np.sum(valid_mask) < 5:
                continue
            
            # Calculate intensity similarity weights (per wavelength)
            # Shape: (3, 3, n_wl)
            intensity_diff = np.abs(spectral_cube - center_spectrum[np.newaxis, np.newaxis, :])
            intensity_std = np.std(intensity_diff)
            
            if intensity_std > 0:
                intensity_weights = np.exp(-(intensity_diff**2) / (2 * (thresh * intensity_std)**2))
            else:
                intensity_weights = np.ones_like(intensity_diff)
            
            # Combined weights: spatial * intensity
            combined_weights = spatial_weights[:, :, np.newaxis] * intensity_weights * valid_mask[:, :, np.newaxis]
            
            # Normalize and apply
            weight_sum = np.sum(combined_weights, axis=(0, 1))
            weight_sum[weight_sum == 0] = 1  # Avoid division by zero
            
            filtered_spectrum = np.sum(spectral_cube * combined_weights, axis=(0, 1)) / weight_sum
            
            # Detect cosmics: large deviation from filtered result
            deviation = np.abs(center_spectrum - filtered_spectrum)
            cosmic_mask = deviation > thresh * np.std(deviation)
            
            if np.any(cosmic_mask):
                plb_array = center_spectrum.copy()
                plb_array[cosmic_mask] = filtered_spectrum[cosmic_mask]
                SpectrumDataMatrix_corrected[i][j].PLB = plb_array.tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else plb_array
    
    return SpectrumDataMatrix_corrected

def robust_median_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Robust median-based cosmic removal using spatial median and MAD.
    More robust to outliers than mean-based methods.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    for i in range(1, len(SpectrumDataMatrix) - 1):
        for j in range(1, len(SpectrumDataMatrix[i]) - 1):
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue
            
            # Collect neighbor spectra
            neighbor_spectra = []
            for di in range(-1, 2):
                for dj in range(-1, 2):
                    neighbor = SpectrumDataMatrix[i+di][j+dj]
                    if neighbor is not None and hasattr(neighbor, 'PLB'):
                        neighbor_spectra.append(np.array(neighbor.PLB))
            
            if len(neighbor_spectra) < 5:
                continue
            
            # Stack all neighbors including center
            all_spectra = np.array(neighbor_spectra)
            
            # Calculate median spectrum
            median_spectrum = np.median(all_spectra, axis=0)
            
            # Calculate MAD (Median Absolute Deviation) for robust threshold
            deviations = np.abs(all_spectra - median_spectrum[np.newaxis, :])
            mad = np.median(deviations, axis=0)
            
            # Get center spectrum
            center_spectrum = np.array(SpectrumDataMatrix[i][j].PLB)
            center_deviation = np.abs(center_spectrum - median_spectrum)
            
            # Detect cosmics using MAD-based threshold
            cosmic_mask = center_deviation > thresh * mad
            
            if np.any(cosmic_mask):
                plb_array = center_spectrum.copy()
                plb_array[cosmic_mask] = median_spectrum[cosmic_mask]
                SpectrumDataMatrix_corrected[i][j].PLB = plb_array.tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else plb_array
    
    return SpectrumDataMatrix_corrected

def iterative_cosmic_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Iterative cosmic removal - applies robust median method multiple times.
    Each iteration refines the result.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    # Apply robust median method twice
    for iteration in range(2):
        SpectrumDataMatrix_corrected = robust_median_Matrix(SpectrumDataMatrix_corrected, thresh, width)
    
    return SpectrumDataMatrix_corrected

def pca_anomaly(SpectrumData, thresh, width):
    pass

def pca_anomaly_Matrix(SpectrumDataMatrix, thresh, width):
    """
    PCA-based anomaly detection for cosmic rays.
    Identifies spectra that deviate from principal spectral components.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    # Collect all valid spectra for PCA
    all_spectra = []
    positions = []
    
    for i in range(len(SpectrumDataMatrix)):
        for j in range(len(SpectrumDataMatrix[i])):
            if hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                all_spectra.append(np.array(SpectrumDataMatrix[i][j].PLB))
                positions.append((i, j))
    
    if len(all_spectra) < 10:  # Need enough spectra for PCA
        return SpectrumDataMatrix
    
    # Convert to array: (n_spectra, n_wavelengths)
    spectra_matrix = np.array(all_spectra)
    
    # Perform PCA
    mean_spectrum = np.mean(spectra_matrix, axis=0)
    centered = spectra_matrix - mean_spectrum
    
    # SVD for PCA
    try:
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
        
        # Keep top components (explain 95% variance)
        variance_explained = (S**2) / np.sum(S**2)
        cumsum_var = np.cumsum(variance_explained)
        n_components = np.searchsorted(cumsum_var, 0.95) + 1
        n_components = min(n_components, len(S))
        
        # Project data onto principal components
        projected = U[:, :n_components] @ np.diag(S[:n_components])
        
        # Reconstruct spectra
        reconstructed = projected @ Vt[:n_components, :] + mean_spectrum
        
        # Calculate reconstruction error per spectrum
        reconstruction_error = np.sqrt(np.sum((spectra_matrix - reconstructed)**2, axis=1))
        
        # Threshold for anomaly detection
        error_median = np.median(reconstruction_error)
        error_mad = np.median(np.abs(reconstruction_error - error_median))
        
        # Identify anomalous spectra
        for idx, (i, j) in enumerate(positions):
            if reconstruction_error[idx] > error_median + thresh * error_mad:
                # Replace with reconstructed (cleaned) spectrum
                SpectrumDataMatrix_corrected[i][j].PLB = reconstructed[idx].tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else reconstructed[idx]
    
    except np.linalg.LinAlgError:
        # If PCA fails, return original
        pass
    
    return SpectrumDataMatrix_corrected

def gradient_based_Matrix(SpectrumDataMatrix, thresh, width):
    """
    Gradient-based cosmic detection using directional spatial gradients.
    Cosmics create high gradients in all directions, real features only in some.
    """
    SpectrumDataMatrix_corrected = copy.deepcopy(SpectrumDataMatrix)
    
    if SpectrumDataMatrix is None or len(SpectrumDataMatrix) < 3 or len(SpectrumDataMatrix[0]) < 3:
        return SpectrumDataMatrix
    
    for i in range(1, len(SpectrumDataMatrix) - 1):
        for j in range(1, len(SpectrumDataMatrix[i]) - 1):
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue
            
            center_spectrum = np.array(SpectrumDataMatrix[i][j].PLB)
            
            # Calculate gradients in 4 directions: horizontal, vertical, diagonal1, diagonal2
            gradients = []
            
            # Horizontal gradient
            if hasattr(SpectrumDataMatrix[i][j-1], 'PLB') and hasattr(SpectrumDataMatrix[i][j+1], 'PLB'):
                grad_h = np.abs(np.array(SpectrumDataMatrix[i][j+1].PLB) - np.array(SpectrumDataMatrix[i][j-1].PLB))
                gradients.append(grad_h)
            
            # Vertical gradient
            if hasattr(SpectrumDataMatrix[i-1][j], 'PLB') and hasattr(SpectrumDataMatrix[i+1][j], 'PLB'):
                grad_v = np.abs(np.array(SpectrumDataMatrix[i+1][j].PLB) - np.array(SpectrumDataMatrix[i-1][j].PLB))
                gradients.append(grad_v)
            
            # Diagonal gradients
            if hasattr(SpectrumDataMatrix[i-1][j-1], 'PLB') and hasattr(SpectrumDataMatrix[i+1][j+1], 'PLB'):
                grad_d1 = np.abs(np.array(SpectrumDataMatrix[i+1][j+1].PLB) - np.array(SpectrumDataMatrix[i-1][j-1].PLB))
                gradients.append(grad_d1)
            
            if hasattr(SpectrumDataMatrix[i-1][j+1], 'PLB') and hasattr(SpectrumDataMatrix[i+1][j-1], 'PLB'):
                grad_d2 = np.abs(np.array(SpectrumDataMatrix[i+1][j-1].PLB) - np.array(SpectrumDataMatrix[i-1][j+1].PLB))
                gradients.append(grad_d2)
            
            if len(gradients) < 2:
                continue
            
            # Stack gradients
            gradient_stack = np.array(gradients)
            
            # Cosmic rays have high gradients in ALL directions
            # Real features have high gradients in SOME directions
            mean_gradient = np.mean(gradient_stack, axis=0)
            std_gradient = np.std(gradient_stack, axis=0)
            
            # High mean gradient AND low std = uniform high gradient = cosmic
            # We want points where gradient is consistently high across directions
            high_gradient = mean_gradient > thresh * np.std(mean_gradient)
            uniform_gradient = std_gradient < 0.5 * mean_gradient  # Consistent across directions
            
            cosmic_mask = high_gradient & uniform_gradient
            
            if np.any(cosmic_mask):
                # Replace with median of valid neighbors
                neighbor_spectra = []
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        if (di != 0 or dj != 0):
                            neighbor = SpectrumDataMatrix[i+di][j+dj]
                            if neighbor is not None and hasattr(neighbor, 'PLB'):
                                neighbor_spectra.append(np.array(neighbor.PLB))
                
                if len(neighbor_spectra) > 0:
                    median_spectrum = np.median(neighbor_spectra, axis=0)
                    plb_array = center_spectrum.copy()
                    plb_array[cosmic_mask] = median_spectrum[cosmic_mask]
                    SpectrumDataMatrix_corrected[i][j].PLB = plb_array.tolist() if isinstance(SpectrumDataMatrix_corrected[i][j].PLB, list) else plb_array
    
    return SpectrumDataMatrix_corrected

# ========== SINGLE SPECTRUM METHODS (uncorrelated) ==========

# Matrix image correction method
def matrix_image_correction(data, thresh, width):
    return data # placeholder, no correction applied

def cosmic_correlation(data, thresh, width):
    return data # placeholder, no correction applied

def adaptive_threshold(data, thresh, width):
    """Single spectrum adaptive threshold - uses local statistics."""
    return data # Single spectrum version - no spatial correlation available

def spectral_correlation(data, thresh, width):
    """
    Detect and remove cosmics based on spectral continuity.
    Cosmic rays create sharp spikes that break spectral smoothness.
    """
    data = np.array(data)
    corrected = data.copy()
    
    # Calculate second derivative to find sharp spikes
    if len(data) < 3:
        return data
    
    # Use central differences for second derivative
    second_deriv = np.zeros_like(data)
    second_deriv[1:-1] = data[:-2] - 2*data[1:-1] + data[2:]
    
    # Identify spikes where second derivative is large
    std_deriv = np.std(second_deriv)
    cosmic_mask = np.abs(second_deriv) > thresh * std_deriv
    
    # Interpolate detected cosmics from neighbors
    for i in np.where(cosmic_mask)[0]:
        left = max(0, i - width)
        right = min(len(data), i + width + 1)
        # Use median of local window excluding the spike
        window = np.concatenate([data[left:i], data[i+1:right]])
        if len(window) > 0:
            corrected[i] = np.median(window)
    
    return corrected.tolist() if isinstance(data, list) else corrected

def bilateral_filter(data, thresh, width):
    """Single spectrum bilateral filter - spatial only (no intensity similarity)."""
    return data # Single spectrum version - simplified to median filter

def robust_median(data, thresh, width):
    """
    Robust median-based cosmic removal for single spectrum.
    Uses local median and MAD (Median Absolute Deviation) for detection.
    """
    data = np.array(data)
    corrected = data.copy()
    
    if len(data) < width:
        return data
    
    # Use rolling window median
    from scipy.ndimage import median_filter as medfilt
    median_filtered = medfilt(data, size=width)
    
    # Calculate deviation from median
    deviation = np.abs(data - median_filtered)
    
    # Use MAD for robust threshold
    mad = np.median(deviation)
    cosmic_mask = deviation > thresh * mad
    
    # Replace cosmics with local median
    corrected[cosmic_mask] = median_filtered[cosmic_mask]
    
    return corrected.tolist() if isinstance(data, list) else corrected

def iterative_cosmic(data, thresh, width):
    """
    Iterative cosmic removal - applies detection multiple times.
    """
    data = np.array(data)
    corrected = data.copy()
    
    # Apply median filter iteratively (2 passes)
    for iteration in range(2):
        from scipy.ndimage import median_filter as medfilt
        median_filtered = medfilt(corrected, size=width)
        deviation = np.abs(corrected - median_filtered)
        
        # Adaptive threshold: use local std
        std_dev = np.std(deviation)
        cosmic_mask = deviation > thresh * std_dev
        
        corrected[cosmic_mask] = median_filtered[cosmic_mask]
    
    return corrected.tolist() if isinstance(data, list) else corrected

def remove_cosmics_spacial(data, thresh, width):
    """
    Remove cosmics from a single spectrum using spatial information.
    This is a placeholder for a method that would require neighboring spectra.
    """
    return data # Single spectrum version - no spatial correlation available

def remove_cosmics_spacial_Matrix(SpectrumDataMatrix, threshold=5):
    X = len(SpectrumDataMatrix)
    Y = len(SpectrumDataMatrix[0])

    for i in range(X):
        for j in range(Y):
            if not hasattr(SpectrumDataMatrix[i][j], 'PLB'):
                continue

            center = SpectrumDataMatrix[i][j].PLB.copy()

            # collect 3x3 neighborhood
            neigh = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < X and 0 <= nj < Y:
                        if hasattr(SpectrumDataMatrix[ni][nj], 'PLB'):
                            neigh.append(SpectrumDataMatrix[ni][nj].PLB)

            if len(neigh) < 3:
                continue  # not enough stats

            neigh = np.array(neigh)  # shape: (Nneigh, M)

            # spatial median + std
            med = np.median(neigh, axis=0)
            std = np.std(neigh, axis=0)

            # avoid division problems
            std[std == 0] = 1e-12

            # detect cosmics
            mask = np.abs(center - med) > threshold * std

            # replace only cosmics
            center[mask] = med[mask]

            SpectrumDataMatrix[i][j].PLB = center

    return SpectrumDataMatrix

def whitaker_hayes_despike(flux, threshold=7, window=5):
    """Modified z-score despiking (Whitaker & Hayes 2018). https://doi.org/10.1016/j.chemolab.2018.06.009 """ 
    diff = np.diff(flux, prepend=flux[0])          # emphasize sharp spikes
    med = np.median(diff)
    mad = np.median(np.abs(diff - med)) or 1e-12
    mod_z = 0.6745 * (diff - med) / mad             # 0.6745 = MAD->sigma
    spikes = np.abs(mod_z) > threshold

    clean = flux.copy()
    idx = np.arange(len(flux))
    for i in np.where(spikes)[0]:
        lo, hi = max(0, i - window), min(len(flux), i + window + 1)
        good = ~spikes[lo:hi]
        if good.sum() >= 2:
            clean[i] = np.interp(i, idx[lo:hi][good], flux[lo:hi][good])
    return clean#, spikes

def gradient_based(data, thresh, width):
    """
    Gradient-based cosmic detection for single spectrum.
    Detects sharp changes in spectral gradient.
    """
    data = np.array(data)
    corrected = data.copy()
    
    if len(data) < 3:
        return data
    
    # Calculate gradient (first derivative)
    gradient = np.gradient(data)
    
    # Detect sudden changes in gradient
    gradient_change = np.abs(np.diff(gradient))
    
    if len(gradient_change) > 0:
        std_change = np.std(gradient_change)
        # Look for spikes in gradient change
        for i in range(1, len(data) - 1):
            if i < len(gradient_change) and gradient_change[i-1] > thresh * std_change:
                # Interpolate from neighbors
                corrected[i] = (data[i-1] + data[i+1]) / 2
    
    return corrected.tolist() if isinstance(data, list) else corrected

# Linear interpolation provides a simple and fast method for filling in missing values.
def remove_cosmics_linear(data, thresh, width):
    # Calculate the finite differences (first derivative)
    diff = np.diff(data)
    
    # Identify where the absolute value of the difference is greater than the threshold
    cosmic_starts = np.where(np.abs(diff) > thresh)[0]
    
    # Copy the data to avoid modifying the original array
    interpolated_data = np.copy(data)
    
    for start in cosmic_starts:
        # Define the range to look for the end of the cosmic
        end_range = min(start + width, len(data) - 1)
        
        # Identify the end of the cosmic within the given width
        for end in range(start + 1, end_range + 1):
            if np.abs(data[end] - data[start]) < thresh:
                # Perform linear interpolation between start and end
                interpolated_data[start:end+1] = np.linspace(data[start], data[end], end - start + 1)
                break
    
    return interpolated_data

#Median filtering is good for robust outlier removal.
def remove_cosmics_median_filter(data, thresh, width):
    # Apply a median filter with size `width`
    filtered_data = median_filter(data, size=width)
    
    # Replace cosmic rays with median values
    diff = np.abs(data - filtered_data)
    cosmic_indices = np.where(diff > thresh)[0]
    interpolated_data = np.copy(data)
    interpolated_data[cosmic_indices] = filtered_data[cosmic_indices]
    
    return interpolated_data

#Rolling mean is good for preserving the overall trend of the data.            
def remove_cosmics_rolling_mean(data, thresh, width):
    # Compute the rolling mean
    rolling_mean = pd.Series(data).rolling(window=width, center=True).mean().to_numpy()
    
    # Replace cosmic rays with rolling mean values
    diff = np.abs(data - rolling_mean)
    cosmic_indices = np.where(diff > thresh)[0]
    interpolated_data = np.copy(data)
    interpolated_data[cosmic_indices] = rolling_mean[cosmic_indices]
    
    return interpolated_data

#Spline interpolation is good for preserving the overall trend of the data.
def remove_cosmics_spline(data, thresh, width):
    # Create an array of indices
    x = np.arange(len(data))
    
    # Identify cosmic indices
    diff = np.diff(data)
    cosmic_indices = np.where(np.abs(diff) > thresh)[0]
    
    # Create a spline interpolation of the data
    spline = UnivariateSpline(x, data, s=width)
    
    # Replace cosmic rays with spline values
    interpolated_data = np.copy(data)
    for idx in cosmic_indices:
        end_range = min(idx + width, len(data) - 1)
        interpolated_data[idx:end_range+1] = spline(x[idx:end_range+1])
    
    return interpolated_data

#Nearest neighbor interpolation is good for preserving the overall trend of the data.
def remove_cosmics_nearest_neighbor(data, thresh, width):
    # Identify where the absolute value of the difference is greater than the threshold
    diff = np.diff(data)
    #cosmic_indices = np.where(np.abs(diff) > thresh)[0]

    cosmic_indices = np.where(np.abs(data) > thresh)[0]
    # Identify the start and end of cosmic rays
    start_indices = np.where(np.diff(cosmic_indices) != 1)[0]
    end_indices = np.append(start_indices[1:], len(cosmic_indices) - 1)
    cosmic_ranges = [(cosmic_indices[start], cosmic_indices[end]) for start, end in zip(start_indices, end_indices)]

    # Copy the data to avoid modifying the original array
    interpolated_data = np.copy(data)

    #print(cosmic_ranges)
    #print(interpolated_data)

    #print(len(interpolated_data), print(len(interpolated_data[0])))

    
    return interpolated_data

# Combined method: Linear interpolation followed by nearest neighbor averaging
def remove_cosmics_linear_then_neighbor(data, thresh, width):
    """
    Two-stage cosmic ray removal:
    1. First pass: Linear interpolation for cosmics detected by gradient changes
    2. Second pass: Nearest neighbor averaging for any remaining outliers
    
    Parameters:
    -----------
    data : array-like
        The spectrum data to clean
    thresh : float
        Threshold for cosmic detection
    width : int
        Width parameter for cosmic region detection
        
    Returns:
    --------
    corrected_data : array
        Spectrum with cosmics removed
    """
    data = np.array(data)
    
    # Stage 1: Linear interpolation for gradient-based cosmic detection
    # Calculate the finite differences (first derivative)
    diff = np.diff(data)
    
    # Identify where the absolute value of the difference is greater than the threshold
    cosmic_starts = np.where(np.abs(diff) > thresh)[0]
    
    # Copy the data to avoid modifying the original array
    corrected_data = np.copy(data)
    
    # Track which indices were corrected in stage 1
    corrected_indices = set()
    
    for start in cosmic_starts:
        # Define the range to look for the end of the cosmic
        end_range = min(start + width, len(data) - 1)
        
        # Identify the end of the cosmic within the given width
        for end in range(start + 1, end_range + 1):
            if np.abs(data[end] - data[start]) < thresh:
                # Perform linear interpolation between start and end
                corrected_data[start:end+1] = np.linspace(data[start], data[end], end - start + 1)
                corrected_indices.update(range(start, end + 1))
                break
    
    # Stage 2: Nearest neighbor averaging for remaining outliers
    # Apply median filter to get smoothed reference
    from scipy.ndimage import median_filter as medfilt
    
    if len(corrected_data) >= width:
        median_filtered = medfilt(corrected_data, size=width)
        
        # Calculate deviation from median
        deviation = np.abs(corrected_data - median_filtered)
        
        # Use MAD (Median Absolute Deviation) for robust threshold
        mad = np.median(deviation)
        
        # Detect remaining cosmics (excluding already corrected indices)
        if mad > 0:
            remaining_cosmic_mask = deviation > (thresh * 0.5)  # Use lower threshold for second pass
            
            # Only correct indices that weren't already fixed in stage 1
            for idx in np.where(remaining_cosmic_mask)[0]:
                if idx not in corrected_indices:
                    # Use local window average (nearest neighbors)
                    left = max(0, idx - width // 2)
                    right = min(len(corrected_data), idx + width // 2 + 1)
                    
                    # Get neighbors excluding the cosmic point itself
                    neighbors = np.concatenate([corrected_data[left:idx], corrected_data[idx+1:right]])
                    
                    if len(neighbors) > 0:
                        # Use median of neighbors for robust estimation
                        corrected_data[idx] = np.median(neighbors)
    
    return corrected_data.tolist() if isinstance(data, list) else corrected_data

# get grid dx and dy
def most_freq_element(arr):
    # Use a dictionary to count occurrences of each element
    count_dict = {}
    for value in arr:
        if value in count_dict:
            count_dict[value] += 1
        else:
            count_dict[value] = 1
    # Find the element with the highest count
    max_count = 0
    most_frequent = 1
    for key, count in count_dict.items():
        if count > max_count:
            max_count = count
            most_frequent = key
    return most_frequent

def remove_duplicates(arr):
    # Use a set to track unique elements
    seen = set()
    unique_arr = []
    for value in arr:
        if value not in seen:
            unique_arr.append(value)
            seen.add(value)
    return unique_arr

def findif(a):
    find = []
    for i in range(len(a)-1):
        find.append(a[i+1]-a[i])
    return find   

# Function stubs for menu commands
def file_new():
    print("New File")

def file_open():
    print("Open File")

def file_save():
    print("Save File")

def edit_undo():
    print("Undo")

def edit_redo():
    print("Redo")

def select_all():
    print("Select All")

def view_toggle():
    print("Toggle View")

def create_menu(root, menu_bar):
    # Create File menu
    file_menu = tk.Menu(menu_bar, tearoff=0)
    menu_bar.add_cascade(label="File", menu=file_menu)
    file_menu.add_command(label="New", command=file_new)
    file_menu.add_command(label="Open", command=file_open)
    file_menu.add_command(label="Save", command=file_save)
    file_menu.add_separator()
    file_menu.add_command(label="Exit", command=root.quit)

    # Create Edit menu
    edit_menu = tk.Menu(menu_bar, tearoff=0)
    menu_bar.add_cascade(label="Edit", menu=edit_menu)
    edit_menu.add_command(label="Undo", command=edit_undo)
    edit_menu.add_command(label="Redo", command=edit_redo)

    # Create Selection menu
    selection_menu = tk.Menu(menu_bar, tearoff=0)
    menu_bar.add_cascade(label="Selection", menu=selection_menu)
    selection_menu.add_command(label="Select All", command=select_all)

    # Create View menu
    view_menu = tk.Menu(menu_bar, tearoff=0)
    menu_bar.add_cascade(label="View", menu=view_menu)
    view_menu.add_command(label="Toggle View", command=view_toggle)

def on_click(event, PixMatrix):
    """
    Handles mouse click events on a matplotlib figure.
    Prints the coordinates and value at the clicked position in PixMatrix.
    """
    if event.inaxes:
        try:
            x = int(round(event.xdata))
            y = int(round(event.ydata))
            if 0 <= y < len(PixMatrix) and 0 <= x < len(PixMatrix[0]):
                value = PixMatrix[y][x]
                print(f"Clicked at (x={x}, y={y}), value: {value}")
            else:
                print(f"Clicked at (x={x}, y={y}), but out of bounds.")
        except Exception as e:
            print(f"Error handling click event: {e}")

def select_file(entry_var):
    file_path = filedialog.askopenfilename()
    if file_path:
        entry_var.set(file_path)

def select_folder(entry_var):
    folder_path = filedialog.askdirectory()
    if folder_path:
        entry_var.set(folder_path)

def browse_folder(folder_entry):
    folder_selected = filedialog.askdirectory()
    if folder_selected:
        folder_entry.delete(0, tk.END)
        folder_entry.insert(0, folder_selected)
    
def loadclaraimage(file):
    with open(file) as f:
        # Read header lines until we find 2 consecutive lines without ':'
        consecutive_no_colon = 0
        while consecutive_no_colon < 2:
            line = f.readline()
            if not line:  # End of file
                break
            if ':' not in line:
                consecutive_no_colon += 1
            else:
                consecutive_no_colon = 0
        
        # now: skip empty lines
        line = f.readline()
        while line.strip() == '':
            line = f.readline()

        # Read remaining data lines
        fload = f.readlines()
    
    x = []
    y = []
    data = []
    for i in fload:
        isplit = i.split('\n')[0].split('\t')
        x.append(float(isplit[0]))
        for j in range(1, len(isplit)):
            if isplit[j] == '':
                pass
            else:
                y.append(float(isplit[j]))
        data.append(y)
        y = []  
    return np.asarray(data)

# check if file ends with a number before the file extension
# increment the number by 1 and return the new filename
def increment_filename(file):
    fnoext = file.split('.')[-2]
    fext = file.split('.')[-1]
    num = 1
    # check if the file ends with a number and increment it by 1
    if fnoext[-1].isdigit():
        num = int(fnoext[-1]) + 1
        fnoext = fnoext[:-1]
    return fnoext + str(num) + '.' + fext

def closest_indices(X, Y, px, py):
    X = np.asarray(X)
    Y = np.asarray(Y)
    # Calculate the absolute differences between px and all elements in X
    diff_x = np.abs(X - px)
    # Calculate the absolute differences between py and all elements in Y
    diff_y = np.abs(Y - py)
    # Find the index of the minimum difference in X
    i = np.argmin(diff_x)
    # Find the index of the minimum difference in Y
    j = np.argmin(diff_y)
    return i, j

#roi = region of interest, returns a matrix with 1s in the region of interest and NaNs elsewhere
def highlight_roi(Mat, points):
    Mat = np.asarray(Mat)
    # Create a copy of the matrix initialized with NaN
    result = np.full_like(Mat, np.nan, dtype=float)
    # Get grid of all pixel coordinates in the matrix
    y, x = np.meshgrid(np.arange(Mat.shape[1]), np.arange(Mat.shape[0]))
    points_grid = np.vstack((x.ravel(), y.ravel())).T
    # Create a Path object from the points (closed polygon)
    polygon_path = Path(points)
    # Find which points are inside the polygon
    inside_mask = polygon_path.contains_points(points_grid)
    inside_mask = inside_mask.reshape(Mat.shape)
    # Set inside points to 1, leaving the rest as NaN
    result[inside_mask] = 1
    return np.transpose(result)

def correct_spectrum(spec, wl, speccorrect, wlcorrect):
    """
    Corrects a spectrum by interpolating a correction spectrum to match the original wavelength axis.

    Parameters:
    - spec (array-like): The original spectrum values.
    - wl (array-like): The wavelength axis of the original spectrum.
    - speccorrect (array-like): The correction spectrum values.
    - wlcorrect (array-like): The wavelength axis of the correction spectrum.

    Returns:
    - corrected_spec (np.ndarray): The corrected spectrum.
    """
    # check if the wavelength axes are sorted
    if not np.all(np.diff(wl) > 0):
        raise ValueError("Wavelength axis is not sorted.")
    if not np.all(np.diff(wlcorrect) > 0):
        raise ValueError("Correction wavelength axis is not sorted.")
    
    if wl[0] < wlcorrect[0]:
        # insert one element at the beginning of wlcorrect
        wlcorrect = np.insert(wlcorrect, 0, wlcorrect[0]-0.001)
        speccorrect = np.insert(speccorrect, 0, 0)
    if wl[-1] > wlcorrect[-1]:
        # append one element at the end of wlcorrect
        wlcorrect = np.append(wlcorrect, wlcorrect[-1]+0.001)
        speccorrect = np.append(speccorrect, 0)
    
    # Ensure wl and wlcorrect are numpy arrays
    wl = np.array(wl)
    wlcorrect = np.array(wlcorrect)
    dwlcorrect = np.diff(wlcorrect)
    dwl = np.diff(wl)
    # get the most frequent element of dwl
    dwlv = most_freq_element(dwl)
    
    # append dwl to dwlcorrect
    dwlcorrect = np.append(dwlcorrect, dwlv)

    # calculate smallest common divisor
    gcd = np.gcd.reduce(dwlcorrect.astype(int))
    # create a new wlcorrect with kgv
    wlcorrectinter = np.arange(wlcorrect[0], wlcorrect[-1], gcd)

    # interpolate speccorrect to match the new wlcorrect
    corrected_interp_spec = np.interp(wlcorrectinter, wlcorrect, speccorrect)
    # interpolate the original spectrum to match the new wlcorrect
    spec_wlc = np.interp(wlcorrectinter, wl, spec)
    # subtract the interpolated correction spectrum from the interpolated original spectrum
    corrected_spec = spec_wlc - corrected_interp_spec

    # interpolate the correction spectrum to match the original wavelength axis
    specc = np.interp(wl, wlcorrect, speccorrect)
    
    # Interpolate the correction spectrum to match the original wavelength axis
    #interpolated_correction = np.interp(wl, wlcorrect, speccorrect)

    # Subtract the interpolated correction spectrum from the original spectrum
    corrected_spec = np.array(spec) - specc #interpolated_correction

    return corrected_spec

def loadexpspec(fname):
    #Loads a spectrum from a file.
    wavelengths = []
    intensities = []
    
    with open(fname, 'r') as file:
        for line in file:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Parse data lines
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    wavelength = float(parts[0])
                    intensity = float(parts[1])
                    wavelengths.append(wavelength)
                    intensities.append(intensity)
                except ValueError:
                    continue  # Skip lines that can't be parsed
    return wavelengths, intensities

def testcorrect_spectrum():
    # spec in "C:\Users\volib\Documents\spectra\spect1.txt"
    wl, spec = loadexpspec("C:/Users/volib/Documents/spectra/spect1.txt")
    # corectspec in "C:\Users\volib\Documents\spectra\correctspec_t2.txt"
    wlcorrect, speccorrect = loadexpspec("C:/Users/volib/Documents/spectra/correctionspec_t2.txt")
    # Correct the spectrum
    corrected_spec = correct_spectrum(spec, wl, speccorrect, wlcorrect)
    # plot the corrected spectrum and the original spectrum
    plt.figure(figsize=(10, 6))
    plt.plot(wl, spec, label='Original Spectrum')
    plt.plot(wl, corrected_spec, label='Corrected Spectrum')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity')
    plt.title('Spectrum Correction')
    plt.legend()
    plt.show()

def remove_nan(arr):
    # Remove NaN values from the array and return it
    reta = []
    return reta

def return0():
    return 0

def fig_on_hoverevent(event, ax, fig, Z, x_range, y_range):
    def on_hover(event, ax, fig, Z, x_range, y_range):
        """
        Handles mouse hover event to display x, y, and z values dynamically for a 2D list or array.

        Parameters:
        - event: The mouse motion event.
        - ax: The axes of the plot.
        - fig: The matplotlib figure object.
        - Z: The 2D list (or array) of values.
        - x_range: Tuple (x_min, x_max) representing the x-coordinate range.
        - y_range: Tuple (y_min, y_max) representing the y-coordinate range.
        """
        if event.inaxes == ax:
            x_min, x_max = x_range
            y_min, y_max = y_range

            # Check if Z is a 2D list (list of lists)
            if not (isinstance(Z, list) and all(isinstance(row, list) for row in Z)):
                raise ValueError("Z must be a 2D list.")

            # Get the number of rows and columns
            num_rows = len(Z)
            num_cols = len(Z[0])

            # Calculate indices based on cursor position and plot extent
            x_idx = int((event.xdata - x_min) / (x_max - x_min) * (num_cols - 1))
            y_idx = int((event.ydata - y_min) / (y_max - y_min) * (num_rows - 1))

            if 0 <= x_idx < num_cols and 0 <= y_idx < num_rows:
                z_value = Z[y_idx][x_idx]  # Access the element in the 2D list
                ax.set_title(f"x: {event.xdata:.2f}, y: {event.ydata:.2f}, z: {z_value:.2f}")
                fig.canvas.draw_idle()

# todo: add a method that uses linear interpolation first and afterwards applies Nearest Neighbor average for remaining cosmics

cosmicfuncts = {
                'Linear Interpolation': remove_cosmics_linear, 
                'Median Filter': remove_cosmics_median_filter, 
                #'Rolling Mean': remove_cosmics_rolling_mean, 
                #'Spline Interpolation': remove_cosmics_spline,
                'Nearest Neighbor average': remove_cosmics_nearest_neighbor,
                'Linear + Neighbor (Combined)': remove_cosmics_linear_then_neighbor,
                'Matrix Image correction': matrix_image_correction, 
                'Cosmic Correlation': cosmic_correlation,
                'Adaptive Threshold': adaptive_threshold,
                'Spectral Correlation': spectral_correlation,
                'Bilateral Filter': bilateral_filter,
                'Robust Median': robust_median,
                'Iterative Cosmic': iterative_cosmic,
                'PCA Anomaly': pca_anomaly,
                'Gradient Based': gradient_based, 
                'spacial removal': remove_cosmics_spacial,
                'withaker-hayes despike': whitaker_hayes_despike
                } 

# Correlation functions for cosmic ray removal must be run by XYMap class instead of SpectrumData class only
correlationcosmicfuncts = {
                'Matrix Image correction': matrix_image_correction_Matrix,
                'Cosmic Correlation Matrix': cosmic_correlation_Matrix,
                'Adaptive Threshold Matrix': adaptive_threshold_Matrix,
                'Spectral Correlation Matrix': spectral_correlation_Matrix,
                'Bilateral Filter Matrix': bilateral_filter_Matrix,
                'Robust Median Matrix': robust_median_Matrix,
                'Iterative Cosmic Matrix': iterative_cosmic_Matrix,
                'PCA Anomaly Matrix': pca_anomaly_Matrix,
                'Gradient Based Matrix': gradient_based_Matrix,
                'spacial removal Matrix': remove_cosmics_spacial_Matrix,
                }
# Rolling Mean showed not to be good for cosmic removal
# Spline Interpolation showed not to be good for cosmic removal and took like forever to perform

speckeys = {'Wavelength axis': 'WL', 
            'Background (BG)': 'BG',
            'Counts (PL)': 'PL', 
            'Spectrum (PL-BG)': 'PLB'}

defaults={
    # General inits
    'Matrix_grid_dx': None,
    'Matrix_grid_dy': None,
    'windowsize_X': 900,
    'windowsize_Y': 900,
    # Load Data Notebook
    # Select Folder Frame
    'data_file': 'C:/Users/mol95ww/Desktop/Evaluation/data/VP/test7_MoS2_ML_2sec/test7_MoS2_ML_2sec',
    'filename': 'spectrum', 
    'file_extension': '.txt',
    # Background Subtraction Frame
    'multiple_Background': False,
    'linear_Background': True,
    'default_colormap': 'twilight', 
    # Cosmic Ray Removal Frame
    'remove_cosmics': True,
    'cosmic_threshold': 100,
    'cosmic_width': 10,
    'cosmic_method': 'Linear Interpolation',
    'calculate_first_derivative': False,
    'calculate_second_derivative': False, 
    'derivative_polynomial_order': 2,
    'derivative_Nfitpoints': 5,
    # Clara Image Frame
    'clara_image': 'C:/Users/mol95ww/Desktop/Evaluation/data/2024/qdot_100fach/Laser_in_zpos/145_0.asc',
    'newton_spectrum': 'C:/Users/mol95ww/Desktop/Evaluation/data/2024/Perovskite/N1_Sndoping/Pb-Sn_0_0625/PL_Sn_Image_0_0625_1250g_500lnm_470-1030nm.asc',
    # TCSPC Frame
    'tcspc_maindir': '',
    'tcspc_subdir': '',
    # Hyperspectra Notebook
    # cmap frame
    'lowest_wavelength': 500, 
    'highest_wavelength': 700,
    'colormap_threshold': 10000,
    'fontsize': 13,
    'maxfev': 2000,
    #speckeys
    'Wavelength axis': 'WL', 
    'Background (BG)': 'BG',
    'Counts (PL)': 'PL', 
    'Spectrum (PL-BG)': 'PLB',
    'data_set': 'Spectrum (PL-BG)',
    # buttonframe
    'selected_pixel_x': 0,
    'selected_pixel_y': 0,
    'selected_fit_function': 'gaussian',
    'seperate_fits': False,
    'save_hsi': "hsidata/hsi.txt", 
    'load_hsi_saved': "hsidata/hsi.txt",
    'enable_buttonmatrix': False, 
    'loadonstart': False, 
    'selected_WL_axis': 'Wavelength (nm)',
    'Fit_use_ROI_mask': True, 
    'roi_vis_type': 'cornerlines',
    'HSI_from_fitparam_useROI': True, 
    'power_correction': False, 
    'laser_spotsize_nm': 1000.0,
    # hsi specfilesorter defaults
    'hsifilesorter_maindir': '',
    'hsifilesorter_filename': '',
    'hsifilesorter_fileend': '.txt',
    'hsifilesorter_savedir': '',
    'hsifilesorter_processdir': False,
    'process_multiple_HSIs_bool': False,
    'save_multiple_HSIs_bool': True,
    # HSI Plot Options defaults
    'hsi_cmap': 'hot',
    'hsi_vmin': '',
    'hsi_vmax': '',
    'hsi_scalebar_length': 20.0,
    'hsi_scalebar_width': 2.0,
    'hsi_scalebar_pos_x': 2.0,
    'hsi_scalebar_pos_y': 2.0,
    'hsi_figsize_width': 7.0,
    'hsi_figsize_height': 6.0,
    'hsi_title': '',
    'hsi_cbar_unit': 'Kilo counts',
    'hsi_show_colorbar': True,
    'hsi_scalebar_fontsize': 12,
    'hsi_unit': '$\\mu m$',
    'do_normalize_HSI_var': False,
    'normalize_method': 'None',
    'normalize_wavelength': '600',
    'calc_norm_and_derive': False,
    'calc_norm_on_intensity': False, 
    # roi defaults
    'roi_plot_color': 'red',
    'roi_vis_type': 'cornerlines',
}

defaulttypes = {
    # General inits
    'windowsize_X': int,
    'windowsize_Y': int,
    # Load Data Notebook
    'data_file': str,
    'filename': str,
    'file_extension': str,
    'multiple_Background': bool,
    'linear_Background': bool,
    'default_colormap': str, 
    'remove_cosmics': bool,
    'cosmic_threshold': int,
    'cosmic_width': int,
    'cosmic_method': str,
    'calculate_first_derivative': bool,
    'calculate_second_derivative': bool,
    'derivative_polynomial_order': int,
    'derivative_Nfitpoints': int,
    # Clara Image Frame
    'clara_image': str,
    'newton_spectrum': str,
    'tcspc_maindir': str,
    'tcspc_subdir': str,
    # Hyperspectra Notebook
    'lowest_wavelength': float,
    'highest_wavelength': float,
    'colormap_threshold': int,
    'fontsize': int,
    'maxfev': int,
    'Wavelength axis': str,
    'Background (BG)': str,
    'Counts (PL)': str,
    'Spectrum (PL-BG)': str,
    'data_set': str,
    'selected_pixel_x': int,
    'selected_pixel_y': int,
    'selected_fit_function': str,
    'seperate_fits': bool,
    'save_hsi': str,
    'load_hsi_saved': str, 
    'Matrix_grid_dx': float,
    'Matrix_grid_dy': float,
    'enable_buttonmatrix': bool, 
    'loadonstart': bool,
    'selected_WL_axis': str,
    'Fit_use_ROI_mask': bool,
    'HSI_from_fitparam_useROI': bool, 
    'power_correction': bool,
    'laser_spotsize_nm': float,
    # hsi specfilesorter defaults
    'hsifilesorter_maindir': str,
    'hsifilesoter_filename': str,
    'hsifilesorter_fileend': str,
    'hsifilesorter_savedir': str,
    'hsifilesorter_processdir': str, 
    'process_multiple_HSIs_bool': bool,
    'save_multiple_HSIs_bool': bool,
    # HSI Plot Options types
    'hsi_cmap': str,
    'hsi_vmin': str,
    'hsi_vmax': str,
    'hsi_scalebar_length': float,
    'hsi_scalebar_width': float,
    'hsi_scalebar_pos_x': float,
    'hsi_scalebar_pos_y': float,
    'hsi_figsize_width': float,
    'hsi_figsize_height': float,
    'hsi_title': str,
    'hsi_cbar_unit': str,
    'hsi_show_colorbar': bool,
    'hsi_scalebar_fontsize': int,
    'hsi_unit': str,
    'do_normalize_HSI_var': bool,
    'normalize_method': str,
    'normalize_wavelength': str,
    'calc_norm_and_derive': bool,
    'calc_norm_on_intensity': bool, 
    'roi_plot_color': str,
    'roi_vis_type': str,
}

def testdefaults():
    for i in list(defaults.keys()):
        if i not in list(defaulttypes.keys()):
            print(f'{i} not in defaulttypes')
        if i not in list(defaults.keys()):
            print(f'{i} not in defaults')
        if defaulttypes[i] != type(defaults[i]):
            print(f'{i} not the same type')

def newton_loadfiles(file):
    data = None
    if not file:
        print("Error, Please select a file")
        return
    try:
        data = np.loadtxt(file, skiprows=34)
    except Exception as Error:
        print(f'Error: {Error}')
    
    return data

# WL Axis stuff
WL_values = ['Wavelength (nm)', 'Energy (eV)']

def wl_array_to_ev(wl):
    """
    Convert wavelength in nm to energy in eV.
    
    Parameters:
    - wl: Wavelength in nm.
    
    Returns:
    - Energy in eV.
    """
    h = 4.135667696e-15  # Planck's constant in eV·s
    c = 299792458e9      # Speed of light in nm/s
    for i in range(len(wl)): 
        if wl[i] != 0:
            wl[i] = h * c / wl[i]
        else:
            wl[i] = 0
    return wl  # Energy in eV

def gaussian_weight_matrix(dx=1.0, dy=1.0, sigma_x=None, sigma_y=None, size=3):
    """
    Generate a normalized 2D Gaussian weighting matrix.
    
    Parameters
    ----------
    dx, dy : float
        Sampling spacing along x and y axes (distance between spectra).
    sigma_x, sigma_y : float
        Standard deviations (beam waist) along x and y.
        If sigma_y is None, isotropic Gaussian is assumed.
    size : int
        Matrix dimension (must be odd). Default is 3 for 3x3.
        
    Returns
    -------
    w_norm : 2D numpy array
        Normalized Gaussian weights that sum to 1.
    x_grid, y_grid : 2D numpy arrays
        Coordinate grids (for reference).
    """
    if sigma_y is None:
        sigma_y = sigma_x

    # ensure odd size (so we have a center)
    if size % 2 == 0:
        raise ValueError("Matrix size must be odd (e.g. 3, 5, 7).")

    half = size // 2
    x = np.arange(-half, half + 1) * dx
    y = np.arange(-half, half + 1) * dy
    X, Y = np.meshgrid(x, y)

    # 2D Gaussian weights
    w = np.exp(-((X**2)/(2*sigma_x**2) + (Y**2)/(2*sigma_y**2))) # type: ignore

    # normalize to sum = 1
    w_norm = w / np.sum(w)
    return w_norm, X, Y        

# returnallfolders
def returnallfolders(maindir):
    allfolders = []
    for root, dirs, files in os.walk(maindir):
        for dir in dirs:
            allfolders.append(os.path.join(root, dir))
    return allfolders

class HSIPlotManager:
    """Manage creation and storage of HSI plots using the existing plot_HSI function.

    Usage:
      mgr = HSIPlotManager(tk_root, params)
      fig, ax = mgr.construct_plot(data, metadata)
      mgr.plot([data1, data2], [meta1, meta2])  # plots multiple and stores them

    The constructor accepts a Tk root element and a params object which may be either
    a dict of parameter names to values or a list/tuple of values (less recommended).
    Recognized params (defaults provided):
      cmap, vmin, vmax, scalebarlength, scalebarwidth, scalebarpos,
      figsize, title, xlabel, ylabel, show_colorbar, cbar_unit,
      scalebarfontsize, enable_drag, dx, unit
    """

    def __init__(self, tk_root=None, params=None):
        self.tk_root = tk_root
        # default parameters
        defaults = {
            'cmap': 'hot', 'vmin': None, 'vmax': None,
            'scalebarlength': 20, 'scalebarwidth': 2, 'scalebarpos': (2, 2),
            'figsize': (7, 6), 'title': '', 'xlabel': None, 'ylabel': None,
            'show_colorbar': True, 'cbar_unit': '', 'scalebarfontsize': 12,
            'enable_drag': True, 'dx': 1, 'unit': '$\\mu m$'
        }
        # normalize params into a dict
        self.params = defaults
        if params is not None:
            if isinstance(params, dict):
                for k, v in params.items():
                    if k in self.params:
                        self.params[k] = v
                    else:
                        # allow unknown params to be stored as well
                        self.params[k] = v
            elif isinstance(params, (list, tuple)):
                # If user passed a list, try to map in order of the defaults keys
                keys = list(defaults.keys())
                for i, v in enumerate(params):
                    if i < len(keys):
                        self.params[keys[i]] = v
            else:
                # leave defaults
                pass

        # storage for created plots
        # each entry: dict with keys 'fig', 'ax', 'data', 'metadata', 'params'
        self.plots = []

    def construct_plot(self, data, metadata=None, params_override=None):
        """Construct a single plot using plot_HSI and store the fig/ax.

        Args:
          data: 2D array-like image data
          metadata: optional metadata dict
          params_override: optional dict to override parameters for this plot

        Returns:
          fig, ax
        """
        # merge params
        call_params = dict(self.params)
        if params_override:
            call_params.update(params_override)

        # call the existing plot_HSI function; it returns fig, ax
        fig, ax = plot_HSI(data, metadata=metadata,
                           cmap=call_params.get('cmap'), #type: ignore
                           vmin=call_params.get('vmin'),
                           vmax=call_params.get('vmax'),
                           scalebarlength=call_params.get('scalebarlength'), #type: ignore
                           scalebarwidth=call_params.get('scalebarwidth'), #type: ignore
                           scalebarpos=call_params.get('scalebarpos'),
                           figsize=call_params.get('figsize'),
                           title=call_params.get('title'), #type: ignore
                           xlabel=call_params.get('xlabel'),
                           ylabel=call_params.get('ylabel'),
                           show_colorbar=call_params.get('show_colorbar'), #type: ignore
                           cbar_unit=call_params.get('cbar_unit'), #type: ignore
                           scalebarfontsize=call_params.get('scalebarfontsize'), #type: ignore
                           enable_drag=call_params.get('enable_drag')) #type: ignore
        self.plots.append({'fig': fig, 'ax': ax, 'data': data, 'metadata': metadata, 'params': call_params})
        return fig, ax

    def plot(self, data_list, metadata_list=None, params_override_list=None):
        """Plot multiple images and store them in self.plots.

        Args:
          data_list: iterable of 2D arrays
          metadata_list: optional iterable of metadata dicts (same length as data_list)
          params_override_list: optional iterable of params dicts (same length) to override defaults per-plot

        Returns:
          list of (fig, ax) tuples
        """
        if metadata_list is None:
            metadata_list = [None] * len(data_list)
        if params_override_list is None:
            params_override_list = [None] * len(data_list)

        results = []
        for d, m, p in zip(data_list, metadata_list, params_override_list):
            fig, ax = self.construct_plot(d, metadata=m, params_override=p)
            results.append((fig, ax))
        return results

    def clear(self):
        """Close and clear all stored figures."""
        for entry in self.plots:
            try:
                entry['fig'].clf()
                plt.close(entry['fig'])
            except Exception:
                pass
        self.plots = []

    def show_all(self):
        """Bring all stored figures to front/show them (useful when using interactive backends)."""
        for entry in self.plots:
            try:
                entry['fig'].canvas.manager.window.lift()
            except Exception:
                try:
                    entry['fig'].show()
                except Exception:
                    pass


# Example usage (put this in a new cell or run interactively):
# root = tk.Tk()
# mgr = HSIPlotManager(root, params={'cbar_unit': 'kHz', 'scalebarfontsize': 14})
# figs = mgr.plot([image1, image2], metadata_list=[meta1, meta2])
# each stored in mgr.plots

def plot_HSI(data, metadata=None, cmap='hot', vmin=None, vmax=None,
             scalebarlength=20, scalebarwidth=2, scalebarpos=(2, 2),
             figsize=(7, 6), title='', xlabel=None, ylabel=None,
             show_colorbar=True, cbar_unit='', scalebarfontsize=12,
             enable_drag=True, dx=1, unit='$\\mu m$'):
    """
    Plot HSI data with scale bar and colorbar similar to the reference image.
    
    Parameters:
    -----------
    data : 2D array
        The HSI image data to plot
    metadata : dict, optional
        Metadata dictionary (can contain 'dx', 'dy', 'unit', etc.)
    cmap : str
        Matplotlib colormap name (default 'hot')
    vmin, vmax : float, optional
        Color scale limits
    scalebarlength : float
        Length of scale bar in units (default 20)
    scalebarwidth : float
        Width of scale bar in units (default 2)
    scalebarpos : tuple
        Position of scale bar (x, y) in units from origin (default (2, 2))
    figsize : tuple
        Figure size in inches (default (7, 6))
    title : str
        Plot title
    xlabel, ylabel : str, optional
        Axis labels
    show_colorbar : bool
        Whether to show colorbar (default True)
    cbar_unit : str
        Unit label for colorbar (default '')
    scalebarfontsize : int
        Font size for scale bar text (default 12)
    enable_drag : bool
        Enable interactive dragging (not implemented yet)
    dx : float
        Pixel spacing (default 1)
    unit : str
        Unit string for scale bar (default '$\\mu m$')
    
    Returns:
    --------
    fig, ax : matplotlib figure and axes objects
    """
    
    # Extract dx from metadata if available
    if metadata and 'dx' in metadata:
        dx = metadata['dx']
    if metadata and 'unit' in metadata:
        unit = metadata['unit']
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot the image
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
    
    # Add colorbar
    if show_colorbar:
        cbar = fig.colorbar(im, ax=ax)
        if cbar_unit:
            cbar.set_label(cbar_unit, fontsize=scalebarfontsize)
        cbar.ax.tick_params(labelsize=scalebarfontsize-2)
    
    # Add scale bar
    # Convert scale bar length from real units to pixels
    scalebar_length_pixels = scalebarlength / dx
    scalebar_width_pixels = scalebarwidth / dx
    scalebar_x_pixels = scalebarpos[0] / dx
    scalebar_y_pixels = scalebarpos[1] / dx
    
    # Create scale bar rectangle (white)
    rect = Rectangle((scalebar_x_pixels, scalebar_y_pixels), 
                     scalebar_length_pixels, scalebar_width_pixels,
                     linewidth=1, edgecolor='white', facecolor='white')
    ax.add_patch(rect)
    
    # Add scale bar text
    text_x = scalebar_x_pixels
    text_y = scalebar_y_pixels + scalebar_width_pixels + (2 / dx)  # Slightly above the bar
    ax.text(text_x, text_y, f'{scalebarlength:.2f} {unit}', 
            color='white', fontsize=scalebarfontsize, 
            verticalalignment='bottom', fontweight='bold')
    
    # Set title and labels
    if title:
        ax.set_title(title, fontsize=scalebarfontsize+2)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=scalebarfontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=scalebarfontsize)
    
    # Remove axis ticks for cleaner look (like the reference image)
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.tight_layout()
    
    return fig, ax

### HSI normalize block: implement normalize HSI on a certain criteria


# datatype 2 unit decoding for integration: 
datatype2unit = {
    'Wavelength axis': 'nm',
    'Background (BG)': 'Counts',
    'Counts (PL)': 'Counts',
    'Spectrum (PL-BG)': 'Counts',
    'first_derivative': 'Counts',
    'second_derivative': 'Counts',
    } # rest: use 'A. U.' for arbitrary units for all other data types

# check definitions 
if __name__ == '__main__':
    try:
        testdefaults()
        print('All definitions are correct')
    except Exception as Error:
        print('Error while loading defaults')
        print(f'Error: {Error}')
    
    #testcorrect_spectrum()

