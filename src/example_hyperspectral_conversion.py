#!/usr/bin/env python3
"""
Example script demonstrating SpecMap to Hyperspectral Standard conversion.

This script shows how to:
1. Load SpecMap data
2. Convert to ENVI format (hyperspectral standard)
3. Export to NumPy format

For full interactive examples, see: hyperspectral_standard_interface.ipynb
"""

import numpy as np
import os
import sys

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def specmap_to_envi(xymap, output_filename):
    """
    Convert SpecMap hyperspectral data to ENVI format.
    
    ENVI format is the industry standard for hyperspectral data.
    It consists of two files:
    - .img: Binary data file
    - .hdr: ASCII header file with metadata
    
    Parameters:
    -----------
    xymap : lib9.XYMap
        SpecMap hyperspectral data object with HSI cube and wavelengths
    output_filename : str
        Output filename (without extension)
        
    Returns:
    --------
    tuple : (data_file, header_file) paths
    """
    if not hasattr(xymap, 'HSI') or xymap.HSI is None:
        raise ValueError("No hyperspectral data found in XYMap object")
    
    hsi_data = xymap.HSI
    
    # Create ENVI metadata
    metadata = {
        'lines': hsi_data.shape[0],      # Number of rows
        'samples': hsi_data.shape[1],    # Number of columns
        'bands': hsi_data.shape[2],      # Number of spectral bands
        'header offset': 0,
        'file type': 'ENVI Standard',
        'data type': 4,                   # 4 = 32-bit float
        'interleave': 'bip',             # Band Interleaved by Pixel
        'byte order': 0,                  # Little endian
        'wavelength': xymap.WL.tolist(),
        'wavelength units': 'nm',
        'description': 'SpecMap hyperspectral data cube',
    }
    
    img_file = output_filename + '.img'
    hdr_file = output_filename + '.hdr'
    
    # Write binary data file
    hsi_data.astype(np.float32).tofile(img_file)
    
    # Write header file
    with open(hdr_file, 'w') as f:
        f.write('ENVI\n')
        for key, value in metadata.items():
            if isinstance(value, list):
                # Format list values
                f.write(f'{key} = {{\n')
                for i, v in enumerate(value):
                    f.write(f' {v}')
                    if i < len(value) - 1:
                        f.write(',')
                    if (i + 1) % 10 == 0:  # Line break every 10 values
                        f.write('\n')
                f.write('}\n')
            else:
                f.write(f'{key} = {value}\n')
    
    print(f"ENVI files created:")
    print(f"  Data: {img_file}")
    print(f"  Header: {hdr_file}")
    
    return img_file, hdr_file


def export_specmap_to_numpy(xymap, output_file):
    """
    Export SpecMap data as NumPy archive (.npz).
    
    This format is ideal for sharing data with Python users.
    
    Parameters:
    -----------
    xymap : lib9.XYMap
        SpecMap hyperspectral data object
    output_file : str
        Output filename (with or without .npz extension)
    """
    if not hasattr(xymap, 'HSI') or xymap.HSI is None:
        raise ValueError("No hyperspectral data found")
    
    # Ensure .npz extension
    if not output_file.endswith('.npz'):
        output_file += '.npz'
    
    np.savez_compressed(
        output_file,
        hsi_cube=xymap.HSI,
        wavelengths=xymap.WL,
        x_axis=xymap.xax,
        y_axis=xymap.yax,
    )
    
    print(f"Saved to {output_file}")
    print(f"  To load: data = np.load('{output_file}')")
    print(f"  Access with: data['hsi_cube'], data['wavelengths'], etc.")


def main():
    """
    Example usage - demonstrates loading and converting SpecMap data.
    
    Note: This requires GUI libraries. For headless environments,
    modify the imports or use the data structures directly.
    """
    print("=" * 70)
    print("SpecMap to Hyperspectral Standard Conversion Example")
    print("=" * 70)
    
    # Check for test data
    data_folder = "testdatasets/HSI20240903_I01"
    if not os.path.exists(data_folder):
        print("\nExample data folder not found:", data_folder)
        print("This is a template script. To use it:")
        print("1. Load your SpecMap data using lib9.XYMap")
        print("2. Call specmap_to_envi() or export_specmap_to_numpy()")
        print("\nSee hyperspectral_standard_interface.ipynb for full examples.")
        return
    
    # For actual usage with GUI libraries available:
    try:
        import lib9
        import deflib1 as deflib
        
        print(f"\nSpecMap modules loaded")
        print(f"Loading data from: {data_folder}")
        
        # Initialize and load data
        filename_pattern = "HSI20240903_I01_"
        file_extension = ".txt"
        
        xymap = lib9.XYMap(
            filename_pattern=filename_pattern,
            data_folder=data_folder,
            file_extension=file_extension,
            loadeachbg=False,
            linearbg=False,
            removecosmics=False,
        )
        
        xymap.loadfiles()
        xymap.SpecdataintoMatrix()
        
        if xymap.HSI is not None:
            print(f"\nHyperspectral cube loaded")
            print(f"  Shape: {xymap.HSI.shape} (rows × cols × wavelengths)")
            print(f"  Wavelength range: {xymap.WL[0]:.2f} - {xymap.WL[-1]:.2f} nm")
            
            # Convert to ENVI
            print("\nConverting to ENVI format...")
            specmap_to_envi(xymap, "specmap_export")
            
            # Export to NumPy
            print("\nExporting to NumPy format...")
            export_specmap_to_numpy(xymap, "specmap_export.npz")
            
            print("\n" + "=" * 70)
            print("Conversion complete!")
            print("=" * 70)
            print("\nGenerated files:")
            print("  - specmap_export.img (ENVI binary data)")
            print("  - specmap_export.hdr (ENVI header)")
            print("  - specmap_export.npz (NumPy archive)")
            print("\nThese files can be used with:")
            print("  - Spectral Python (SPy): spectral.io.envi.open('specmap_export.hdr')")
            print("  - PySptools: for spectral unmixing and analysis")
            print("  - Any ENVI-compatible software")
            
    except ImportError as e:
        print(f"\nCannot load SpecMap modules: {e}")
        print("This example requires GUI libraries (tkinter).")
        print("For headless environments, see hyperspectral_standard_interface.ipynb")


if __name__ == '__main__':
    main()
