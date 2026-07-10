"""
HSI Normalization Module

This module provides normalization functionality for hyperspectral imaging (HSI) data.
It allows normalizing pixel intensities based on various criteria before plotting.
"""

import numpy as np

# Minimum normalization threshold to prevent numerical instability
MIN_NORMALIZATION_THRESHOLD = 1e-10


def normalize_derivatives_by_signal(spec_data_matrix, signal_key='PLB'):
    """
    Normalize first and second derivatives by the original signal at each wavelength.
    
    This function computes relative derivatives by dividing the derivative values
    by the corresponding signal values at each wavelength. This is useful because
    absolute derivatives depend on the signal magnitude:
    - A change from 10 to 100 counts gives absolute derivative 90
    - A change from 1 to 10 counts gives absolute derivative 9  
    - But relative derivatives (9x vs 9x) are comparable
    
    For each pixel (i,j) and wavelength k:
        Specdiff1_norm[k] = Specdiff1[k] / signal[k]
        Specdiff2_norm[k] = Specdiff2[k] / signal[k]
    
    Args:
        spec_data_matrix: 2D matrix of SpectrumData objects
        signal_key: Key for the signal to normalize by (default: 'PLB')
        
    Returns:
        None (modifies SpectrumData objects in place by adding 
              Specdiff1_norm and Specdiff2_norm attributes)
    """
    if spec_data_matrix is None or len(spec_data_matrix) == 0:
        return
    
    rows = len(spec_data_matrix)
    cols = len(spec_data_matrix[0]) if rows > 0 else 0
    
    for i in range(rows):
        for j in range(cols):
            spec_data = spec_data_matrix[i][j]
            
            # Check if valid SpectrumData object
            if spec_data is None or not hasattr(spec_data, 'dataokay'):
                continue
            
            if not spec_data.dataokay:
                continue
            
            # Get signal array
            if not hasattr(spec_data, signal_key):
                continue
            
            signal = getattr(spec_data, signal_key)
            if signal is None or len(signal) == 0:
                continue
            
            signal = np.asarray(signal, dtype=np.float32)
            
            # Normalize first derivative if available
            if hasattr(spec_data, 'Specdiff1') and spec_data.Specdiff1 is not None:
                specdiff1 = np.asarray(spec_data.Specdiff1, dtype=np.float32)
                if len(specdiff1) == len(signal):
                    # Create normalized derivative array (float32 to reduce RAM usage)
                    specdiff1_norm = np.zeros(len(specdiff1), dtype=np.float32)
                    # Only normalize where signal is above threshold
                    valid_mask = np.abs(signal) > MIN_NORMALIZATION_THRESHOLD
                    specdiff1_norm[valid_mask] = specdiff1[valid_mask] / signal[valid_mask]
                    spec_data.Specdiff1_norm = specdiff1_norm
            
            # Normalize second derivative if available
            if hasattr(spec_data, 'Specdiff2') and spec_data.Specdiff2 is not None:
                specdiff2 = np.asarray(spec_data.Specdiff2, dtype=np.float32)
                if len(specdiff2) == len(signal):
                    # Create normalized derivative array (float32 to reduce RAM usage)
                    specdiff2_norm = np.zeros(len(specdiff2), dtype=np.float32)
                    # Only normalize where signal is above threshold
                    valid_mask = np.abs(signal) > MIN_NORMALIZATION_THRESHOLD
                    specdiff2_norm[valid_mask] = specdiff2[valid_mask] / signal[valid_mask]
                    spec_data.Specdiff2_norm = specdiff2_norm


class HSINormalization:
    """
    Handles normalization of hyperspectral imaging data.
    
    Generates a normalization matrix that can be applied to HSI pixel data
    before plotting. The normalization is applied by multiplying each pixel's
    intensity by the corresponding normalization value.
    
    Normalization methods:
    - 'none': No normalization (matrix of 1s)
    - 'integrated_counts': Normalize by integrated counts in wavelength range
    - 'max_intensity': Normalize by maximum intensity in wavelength range
    - 'counts_at_wavelength': Normalize by counts at specific wavelength
    """
    
    def __init__(self, spec_data_matrix, wl_array):
        """
        Initialize the HSINormalization object.
        
        Args:
            spec_data_matrix: 2D matrix of SpectrumData objects
            wl_array: Wavelength array (1D numpy array)
        """
        self.spec_data_matrix = spec_data_matrix
        self.wl_array = wl_array
        self.normalization_matrix = None
        
        # Dictionary of normalization functions
        self.normalization_methods = {
            'none': self._normalize_none,
            'integrated_counts': self._normalize_integrated_counts,
            'max_intensity': self._normalize_max_intensity,
            'counts_at_wavelength': self._normalize_counts_at_wavelength,
            'normalize_intern': self._normalize_intern
        }
    
    def generate_normalization_matrix(self, method='none', params=None):
        """
        Generate normalization matrix based on selected method.
        
        Args:
            method: Normalization method name (str)
            params: Dictionary of parameters for the normalization method
                   Expected keys depend on method:
                   - 'integrated_counts': {'wl_start', 'wl_end', 'data_key'}
                   - 'max_intensity': {'wl_start', 'wl_end', 'data_key'}
                   - 'counts_at_wavelength': {'wavelength', 'data_key'}
        
        Returns:
            2D numpy array: Normalization matrix (same shape as spec_data_matrix)
        """
        if params is None:
            params = {}
        
        if method not in self.normalization_methods:
            print(f"Warning: Unknown normalization method '{method}'. Using 'none'.")
            method = 'none'
        
        # Call the appropriate normalization function
        self.normalization_matrix = self.normalization_methods[method](params)
        
        return self.normalization_matrix
    
    def _normalize_none(self, params):
        """
        No normalization - returns matrix of 1s.
        
        Args:
            params: Not used
            
        Returns:
            2D numpy array filled with 1.0
        """
        if self.spec_data_matrix is None or len(self.spec_data_matrix) == 0:
            return np.ones((1, 1))
        
        rows = len(self.spec_data_matrix)
        cols = len(self.spec_data_matrix[0]) if rows > 0 else 0
        return np.ones((rows, cols))
    
    def _normalize_integrated_counts(self, params):
        """
        Normalize by integrated counts in wavelength range.
        
        Each pixel is normalized by 1 / (integrated counts in range).
        
        Args:
            params: Dict with 'wl_start', 'wl_end', 'data_key' (e.g., 'PLB')
            
        Returns:
            2D numpy array: Normalization matrix
        """
        wl_start = params.get('wl_start', np.min(self.wl_array))
        wl_end = params.get('wl_end', np.max(self.wl_array))
        data_key = params.get('data_key', 'PLB')
        
        # Find wavelength indices
        wl_indices = np.where((self.wl_array >= wl_start) & (self.wl_array <= wl_end))[0]
        
        if len(wl_indices) == 0:
            print(f"Warning: No wavelength data in range [{wl_start}, {wl_end}]. Using no normalization.")
            return self._normalize_none(params)
        
        pix_start = wl_indices[0]
        pix_end = wl_indices[-1] + 1
        
        # Generate normalization matrix
        rows = len(self.spec_data_matrix)
        cols = len(self.spec_data_matrix[0]) if rows > 0 else 0
        norm_matrix = np.ones((rows, cols))
        
        for i in range(rows):
            for j in range(cols):
                spec_data = self.spec_data_matrix[i][j]
                
                # Check if valid SpectrumData object
                if spec_data is None or not hasattr(spec_data, 'dataokay'):
                    norm_matrix[i][j] = 1.0
                    continue
                
                if not spec_data.dataokay:
                    norm_matrix[i][j] = 1.0
                    continue
                
                # Get data array based on key
                if not hasattr(spec_data, data_key):
                    norm_matrix[i][j] = 1.0
                    continue
                
                data_array = getattr(spec_data, data_key)
                if data_array is None or len(data_array) == 0:
                    norm_matrix[i][j] = 1.0
                    continue
                
                # Calculate integrated counts
                integrated = np.sum(data_array[pix_start:pix_end])
                
                if integrated > MIN_NORMALIZATION_THRESHOLD:
                    norm_matrix[i][j] = 1.0 / integrated
                else:
                    norm_matrix[i][j] = 1.0
        
        return norm_matrix
    
    def _normalize_max_intensity(self, params):
        """
        Normalize by maximum intensity in wavelength range.
        
        Each pixel is normalized by 1 / (max intensity in range).
        
        Args:
            params: Dict with 'wl_start', 'wl_end', 'data_key' (e.g., 'PLB')
            
        Returns:
            2D numpy array: Normalization matrix
        """
        wl_start = params.get('wl_start', np.min(self.wl_array))
        wl_end = params.get('wl_end', np.max(self.wl_array))
        data_key = params.get('data_key', 'PLB')
        
        # Find wavelength indices
        wl_indices = np.where((self.wl_array >= wl_start) & (self.wl_array <= wl_end))[0]
        
        if len(wl_indices) == 0:
            print(f"Warning: No wavelength data in range [{wl_start}, {wl_end}]. Using no normalization.")
            return self._normalize_none(params)
        
        pix_start = wl_indices[0]
        pix_end = wl_indices[-1] + 1
        
        # Generate normalization matrix
        rows = len(self.spec_data_matrix)
        cols = len(self.spec_data_matrix[0]) if rows > 0 else 0
        norm_matrix = np.ones((rows, cols))
        
        for i in range(rows):
            for j in range(cols):
                spec_data = self.spec_data_matrix[i][j]
                
                # Check if valid SpectrumData object
                if spec_data is None or not hasattr(spec_data, 'dataokay'):
                    norm_matrix[i][j] = 1.0
                    continue
                
                if not spec_data.dataokay:
                    norm_matrix[i][j] = 1.0
                    continue
                
                # Get data array based on key
                if not hasattr(spec_data, data_key):
                    norm_matrix[i][j] = 1.0
                    continue
                
                data_array = getattr(spec_data, data_key)
                if data_array is None or len(data_array) == 0:
                    norm_matrix[i][j] = 1.0
                    continue
                
                # Calculate max intensity
                max_val = np.max(data_array[pix_start:pix_end])
                
                if max_val > MIN_NORMALIZATION_THRESHOLD:
                    norm_matrix[i][j] = 1.0 / max_val
                else:
                    norm_matrix[i][j] = 1.0
        
        return norm_matrix
    
    def _normalize_counts_at_wavelength(self, params):
        """
        Normalize by counts at specific wavelength.
        
        Each pixel is normalized by 1 / (counts at specified wavelength).
        
        Args:
            params: Dict with 'wavelength', 'data_key' (e.g., 'PLB')
            
        Returns:
            2D numpy array: Normalization matrix
        """
        target_wl = params.get('wavelength', np.mean(self.wl_array))
        data_key = params.get('data_key', 'PLB')
        
        # Find closest wavelength index
        wl_idx = np.argmin(np.abs(self.wl_array - target_wl))
        
        # Generate normalization matrix
        rows = len(self.spec_data_matrix)
        cols = len(self.spec_data_matrix[0]) if rows > 0 else 0
        norm_matrix = np.ones((rows, cols))
        
        for i in range(rows):
            for j in range(cols):
                spec_data = self.spec_data_matrix[i][j]
                
                # Check if valid SpectrumData object
                if spec_data is None or not hasattr(spec_data, 'dataokay'):
                    norm_matrix[i][j] = 1.0
                    continue
                
                if not spec_data.dataokay:
                    norm_matrix[i][j] = 1.0
                    continue
                
                # Get data array based on key
                if not hasattr(spec_data, data_key):
                    norm_matrix[i][j] = 1.0
                    continue
                
                data_array = getattr(spec_data, data_key)
                if data_array is None or len(data_array) == 0 or wl_idx >= len(data_array):
                    norm_matrix[i][j] = 1.0
                    continue
                
                # Get counts at wavelength
                counts = data_array[wl_idx]
                
                if counts > MIN_NORMALIZATION_THRESHOLD:
                    norm_matrix[i][j] = 1.0 / counts
                else:
                    norm_matrix[i][j] = 1.0
        
        return norm_matrix
    
    def _normalize_intern(self, params):
        """
        Marker for internal min-max normalization.
        
        This method returns a string marker because internal normalization
        requires access to the final pixel matrix, which is not available
        during normalization matrix generation. The actual normalization
        is performed in apply_normalization.
        
        Args:
            params: Not used
            
        Returns:
            str: 'normalize_intern'
        """
        return 'normalize_intern'
    
    @staticmethod
    def apply_normalization(pixel_matrix, normalization_matrix):
        """
        Apply normalization to pixel matrix.
        
        Multiplies each pixel value by the corresponding normalization value.
        
        Args:
            pixel_matrix: 2D numpy array of pixel values
            normalization_matrix: 2D numpy array of normalization values, or 'normalize_intern' string
            
        Returns:
            2D numpy array: Normalized pixel matrix
        """
        if isinstance(normalization_matrix, str) and normalization_matrix == 'normalize_intern':
            # Min-Max normalization to [0, 1] range
            min_val = np.nanmin(pixel_matrix)
            max_val = np.nanmax(pixel_matrix)
            
            if np.isclose(max_val, min_val):
                # Avoid division by zero if all values are the same
                return np.zeros_like(pixel_matrix)
            
            return (pixel_matrix - min_val) / (max_val - min_val)

        if normalization_matrix is None:
            return pixel_matrix
        
        # Element-wise multiplication
        return pixel_matrix * normalization_matrix
