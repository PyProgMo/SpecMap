# Hypersperctra Notebook Documentation
Notebook Frame Design:

![Hyperspectra Notebook Frame](ProjectImages/Hyperspectra_notebook_frame.png)

This notebook is the heart of the SpecMap software. It is used to visualize and analyze hyperspectral data. The notebook provides various tools for data manipulation, visualization, and analysis. The notebook allows users to perform tasks such as:
- Visualizing hyperspectral images as 2D Images
- Extracting spectra from selected regions
- Performing spectral analysis and classification
- Display single spectra
- Average spectra within regions of interest
- Fit spectra with predefined models and display the results of those models as 2D images. 
- export processed data and results in various formats for further analysis or reporting.

The notebook is designed to be user-friendly and intuitive, allowing users to easily navigate through the different functionalities and perform complex analyses with minimal effort. The notebook also provides options for customizing the visualization settings, such as color maps, scaling, and display options, to enhance the visual representation of the hyperspectral data.

# Description of the Notebook: What do I see?

![WL Range and Colormap Threshold](ProjectImages/wl_threshold_selection.png)

- "Lowest Wavelength" and "Highest Wavelength" are the minimum and maximum wavelengths of the hyperspectral data, respectively. These values are used to define the range of wavelengths that will be displayed in the hyperspectral image. It defines the limits for the integraton and for the Fit of the spectra. Each HSI gets the start and end wavelength entered on the creation into its metadata. 
- "Colormap threshold \ Counts" uses the integration thresholds as fit thresholds. If a pixel-wise fit (all fits are pixel-wise) is below the threshold, no fit is carried out and the pixel obtains a np.nan value. Pixels with np.nan values are displayed white in the HSI. 


![Create HSI column](ProjectImages/create_HSI_column.png)

- "Select Data Set": The user can select the data set to be used for creating the hyperspectral image (HSI). The available data sets are listed in a dropdown menu, and the user can choose the desired data set for analysis. Available datasets are: the raw data, background, background subtracted data (default) and the derivatives calculated in the loading tab. 
- "Create intensity colormap" creates a HSI by integration of the selected Data Set from Lowest Wavelength to Highest Wavelength. All available HSIs are listed in the "Select HSI Image combobox". 
- "Create spectral maximum colormap" use the "Select Fit function" to fit each spectrum (respecting the Colormap threshold) and create a HSI of the spectral shift usually in nm. If the "Fit: use Selected ROI mask is True, the fit is only carried out for pixels within the selected ROI mask. If the "Fit: use Selected ROI mask is False, the fit is carried out for all pixels in the HSI. The resulting HSI will contain the fitted parameter values for each pixel, which can be used for further analysis or visualization. 
Note: this function does only update the fitparameters and displays them afterwards. If a Pixel is excluded by the ROI but a fitparameter exists (is not np.nan), it will be displayed in the HSI. If a Pixel is excluded by the ROI and no fitparameter exists (is np.nan), it will be displayed white in the HSI. To focus on the required area, multiply a ROI mask with the HSI to display only the area of interest. 
- "Update spectral fit maxima" does the same as Create spectral maximum colormap but does not carry out new fits. Just display the existing fitparameters as a HSI. 

# Tasks and Functionalities
1. Create a HSI by integration: 
![Create HSI by Integration](ProjectImages/Hyperspectra_usage1.png)
Create an Image that Integrates all spectra from a start wavelength to an end wavelength and displays the result as an Image. 
