import numpy as np
import os, sys
from scipy.optimize import curve_fit
from scipy.special import jv
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import tkinter as tk

cl_scalings = {'20x': 0.568, '50x': 0.2272, '100x': 0.0568}

class imageprocessor():
    def __init__(self, Notebook, imagefile, loadfunct, metadata, dx, dy):
        self.Notebook = Notebook
        self.imagefile = imagefile
        self.loadfunct = loadfunct
        self.metadata = metadata
        self.dx = dx
        self.dy = dy
        self.g2dpopt = None
        self.buildnotebook()
    
    def fit2dgaussian(self):
        self.g2dpopt = fit_gaussian_2d(self.imagedata, self.dx, self.dy)
    
    def buildnotebook(self):
        # load the image 
        self.imagedata = self.loadfunct(self.imagefile)
        # create a new frame for the image processing
        self.image_frame = tk.Frame(self.Notebook, borderwidth=5, relief="ridge")
        self.image_frame.grid(row=0, column=0, sticky='nsew')

        # create a new frame for the image processing
        self.image_frame = tk.Frame(self.Notebook, borderwidth=5, relief="ridge")
        self.image_frame.grid(row=0, column=0, sticky='nsew')
        self.plotimageb = tk.Button(self.image_frame, text='Plot Image', command=self.plotimage)
        self.plotimageb.grid(row=0, column=0)
        self.fitg2Dbutton = tk.Button(self.image_frame, text='Fit 2D Gaussian', command=lambda: self.fit2dgaussian())
        self.fitg2Dbutton.grid(row=0, column=1)
        self.plotfitbutton = tk.Button(self.image_frame, text='Plot Fit', command=lambda: plot2dfit(self.imagedata, self.g2dpopt, self.dx, self.dy))
        self.plotfitbutton.grid(row=0, column=2)
        self.area = tk.Button(self.image_frame, text='Area', command=lambda: area2dgaussian(self.imagefile, self.g2dpopt, np.exp(-2), self.dx, self.dy))
        self.area.grid(row=0, column=3)

        # plot: add entry for the image label
        tk.Label(self.image_frame, text='Image Label:').grid(row=1, column=0)
        self.imagelabel = tk.Entry(self.image_frame, width=60)
        self.imagelabel.grid(row=1, column=1, columnspan=7)
        tk.Label(self.image_frame, text='Legend Label:').grid(row=2, column=0)
        self.legendlabel = tk.Entry(self.image_frame, width=60)
        self.legendlabel.grid(row=2, column=1, columnspan=7)

    def plotimage(self):
        fig, ax = plt.subplots()
        cim = ax.imshow(self.imagedata, cmap='viridis')
        ''''''
        #fig.colorbar()
        # Get current ticks
        current_xticks = np.arange(self.imagedata.shape[1])
        current_yticks = np.arange(self.imagedata.shape[0])
        # Multiply ticks by constants
        new_xticks = np.multiply(current_xticks, self.dx).round().astype(int)
        new_yticks = np.multiply(current_yticks, self.dy).round().astype(int)
        # Set new ticks
        ax.set_xticks(current_xticks)  # Set the ticks to be at the indices of the current ticks
        ax.set_yticks(current_yticks)
        # Set tick labels to the new values
        ax.set_xticklabels(new_xticks)
        ax.set_yticklabels(new_yticks)
        # Set the axis labels auto adjust
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  # Adjust the number of bins to fit the plot size

        # Add axis labels
        ax.set_xlabel('x axis in $\\mu m$')
        ax.set_ylabel('y axis in $\\mu m$')

        # add colorbar with Legend label
        cbar = fig.colorbar(cim)
        legend_label_text = self.legendlabel.get()
        if legend_label_text:
            cbar.set_label(legend_label_text)
        else:
            cbar.set_label('Counts')
        ''''''
        # Show the plot
        plt.tight_layout()  # Adjust layout to avoid overlapping
        plt.show()
    
def gaussian_2d(coords, x0, y0, sigma_x, sigma_y, amplitude):
    """
    Compute a 2D Gaussian function.
    """
    x, y = coords
    return amplitude * np.exp(
        -(((x - x0)**2) / (2 * sigma_x**2) + ((y - y0)**2) / (2 * sigma_y**2))
    ).ravel()

def fit_gaussian_2d(data, dx, dy):
    """
    Fit a 2D Gaussian function to the data.
    
    Parameters:
    - data: 2D numpy array containing the data to fit.
    - dx: Spacing along the x-axis.
    - dy: Spacing along the y-axis.
    
    Returns:
    - popt: Optimal parameters of the 2D Gaussian (x0, y0, sigma_x, sigma_y, amplitude).
    """
    # Create a meshgrid for the x and y values
    x = np.arange(data.shape[1]) * dx
    y = np.arange(data.shape[0]) * dy
    X, Y = np.meshgrid(x, y)

    # Initial guess for the parameters
    initial_guess = (
        x[np.argmax(data) % data.shape[1]],  # x0 guess
        y[np.argmax(data) // data.shape[1]],  # y0 guess
        1,  # sigma_x guess
        1,  # sigma_y guess
        np.max(data)  # amplitude guess
    )

    # Fit the 2D Gaussian function to the data
    popt, _ = curve_fit(
        gaussian_2d, (X.ravel(), Y.ravel()), data.ravel(), p0=initial_guess
    )
    return popt

def popt2fwhm(popt):
    """
    Convert the parameters of a 2D Gaussian fit to Full Width at Half Maximum (FWHM).
    """
    sigma_x, sigma_y = popt[2], popt[3]
    fwhm_x = 2 * np.sqrt(2 * np.log(2)) * sigma_x
    fwhm_y = 2 * np.sqrt(2 * np.log(2)) * sigma_y
    return fwhm_x, fwhm_y

def area2dgaussian(data, popt, thresh, dx, dy):
    # return the size of the area of the gaussian
    x0, y0, sigma_x, sigma_y, amplitude = popt
    xbelowthresh = find_x_thresh(x0, sigma_x, amplitude, thresh)
    ybelowthresh = find_x_thresh(y0, sigma_y, amplitude, thresh)
    xsize = abs(x0 - xbelowthresh)/2
    ysize = abs(y0 - ybelowthresh)/2
    print('fit widht x and y in mum:', xsize, ysize, 'dx and dy:', dx, dy, 'x0, y0:', x0, y0, 'xsize, ysize:', xsize, ysize, 'xbelowthresh, ybelowthresh:', xbelowthresh, ybelowthresh)
    # Calculate the area of ellipse 
    print('Area of ellipse:', round(np.pi * xsize * ysize, 3), 'mum^2')
    return np.pi * sigma_x * sigma_y

def plot2dfit(data, popt, dx, dy):
    """
    Plot the original 2D data and the 2D Gaussian fit.
    
    Parameters:
    - data: 2D numpy array containing the original data.
    - popt: Optimal parameters of the 2D Gaussian (x0, y0, sigma_x, sigma_y, amplitude).
    - dx: Spacing along the x-axis.
    - dy: Spacing along the y-axis.
    """
    # Extract fit parameters
    x0, y0, sigma_x, sigma_y, amplitude = popt
    
    # Create a meshgrid for the fit
    x = np.arange(data.shape[1]) * dx
    y = np.arange(data.shape[0]) * dy
    X, Y = np.meshgrid(x, y)

    # Compute the 2D Gaussian fit
    fit = amplitude * np.exp(
        -(((X - x0)**2) / (2 * sigma_x**2) + ((Y - y0)**2) / (2 * sigma_y**2))
    )

    # Plot the original data
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].imshow(data, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis')
    ax[0].set_title('Original Data')
    ax[0].set_xlabel('X in mum')
    ax[0].set_ylabel('Y in mum')
    
    # Plot the fitted data
    ax[1].imshow(fit, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis')
    ax[1].set_title('2D Gaussian Fit')
    ax[1].set_xlabel('X in mum')
    ax[1].set_ylabel('Y in mum')

    # add colorbars
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes((0.87, 0.15, 0.05, 0.7))
    fig.colorbar(ax[0].imshow(data, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis'), cax=cbar_ax)
    cbar_ax.set_ylabel('Counts')
    
    plt.tight_layout()
    plt.show()

def find_x_thresh(x0, sigma_x, amplitude, thresh):
    """
    Find x_thresh, the x-coordinate where the Gaussian amplitude falls to thresh.

    Parameters:
    - x0: Center of the Gaussian.
    - sigma_x: Standard deviation of the Gaussian along the x-axis.
    - amplitude: Peak amplitude of the Gaussian.
    - thresh: Threshold value to find x_thresh.

    Returns:
    - x_thresh: The x-coordinate where the amplitude equals thresh.
    # np.sqrt(-2 * sigma_x**2 * np.log(thresh / amplitude) is the formula for x_thresh
    # This formula is derived from the Gaussian function.
    # We can use this formula to find x_thresh without iterating over the Gaussian function.
    """
    if thresh > amplitude:
        raise ValueError("Threshold cannot exceed the Gaussian amplitude.")
    
    # Solve for x_thresh using the Gaussian formula
    x_thresh = x0 + np.sqrt(-2 * sigma_x**2 * np.log(thresh / amplitude))
    return x_thresh

