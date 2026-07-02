import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.special import wofz
from scipy.optimize import fminbound
from scipy.special import jv
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import hilbert


# window functions
def voigt_pseudo(x, cen, fwhm, eta):
    sigma = fwhm / 2.354820045
    gamma = fwhm / 2.0
    G = np.exp(-0.5 * ((x - cen) / sigma)**2)
    L = gamma**2 / ((x - cen)**2 + gamma**2)
    return eta * L + (1 - eta) * G

def doublevoigt5param(x, A, EX, fwhm, r, B, delta=30.0, eta=0.5):
    """
    A     : overall amplitude
    EX    : exciton (center of main peak)
    fwhm  : shared width
    r     : relative amplitude (trion/exciton)
    B     : background
    delta : fixed splitting (meV), default ~30 meV for MoSe2
    eta   : mixing parameter (0=Gaussian, 1=Lorentzian)
    """
    peak1 = voigt_pseudo(x, EX, fwhm, eta)
    peak2 = voigt_pseudo(x, EX - delta, fwhm, eta)
    
    return A * (peak1 + r * peak2) + B

def doublevoigt5paramwind(x, A, EX, fwhm, r, B):
    return doublevoigt5param(x, A, EX, fwhm, r, B)

def double_voigtwind(x, amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2):
    return voigtwind(x, amp1, cen1, wid1, gamma1) + voigtwind(x, amp2, cen2, wid2, gamma2)

def double_gaussianwind(x, amp1, cen1, wid1, amp2, cen2, wid2):
    return gaussianwind(x, amp1, cen1, wid1) + gaussianwind(x, amp2, cen2, wid2)

def double_lorentzwind(x, amp1, cen1, wid1, amp2, cen2, wid2):
    return lorentzwind(x, amp1, cen1, wid1) + lorentzwind(x, amp2, cen2, wid2)

def gaussianwind(x, amp, cen, wid):
    return amp * np.exp(-(x - cen)**2 / (2 * wid**2))

def lorentzwind(x, amp, cen, wid):
    return amp / (1 + (x - cen)**2 / wid**2)

def voigtwindapprox(x, anp, cen, sig, gam):
    # returns the pseudo-voigt provile as a sum of gaussian and lorentzian profiles
    return anp * (1 - gam / (2 * sig) * np.sqrt(2 * np.pi) * np.exp(-2 * np.log(2) * ((x - cen) / sig)**2) + (1 - gam) / np.pi / sig / (1 + ((x - cen) / sig)**2))

# def voigtwind(x, amp, cen, wid, gamma): old version
#    return amp * np.real(wofz(((x - cen) + 1j*gamma) / wid / np.sqrt(2))) / wid / np.sqrt(2*np.pi)

def voigtwind_complex(x, amp, cen, wid, gamma):
    sigma = wid / np.sqrt(2 * np.log(2))  # Convert FWHM to standard deviation
    z = ((x - cen) + 1j * gamma) / (sigma * np.sqrt(2))
    return amp * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))

# pseudo-voigt function for faster fitting
def voigtwind(x, amp, cen, fwhm, eta):
    sigma = fwhm / 2.354820045
    gamma = fwhm / 2.0
    G = np.exp(-0.5*((x-cen)/sigma)**2)
    L = gamma**2 / ((x-cen)**2 + gamma**2)
    return amp * (eta*L + (1-eta)*G)

def linearwind(x, a, b):
    return np.multiply(a, x) + b

# estimate initial parameters for fitting functions

def doublevoigt_5param_approx(x, y, delta=30.0):
    """
    Estimate starting parameters for:
    I(E) = A * [V(EX) + r * V(EX - delta)] + B

    Returns:
    A, EX, fwhm, r, B
    """

    # --- 1. Background (robust) ---
    B = np.percentile(y, 5)  # lower baseline estimate
    y_corr = y - B

    # --- 2. Exciton peak (dominant) ---
    idx_max = np.argmax(y_corr)
    EX = x[idx_max]
    A = y_corr[idx_max]

    # --- 3. FWHM (merged peak) ---
    half_max = A / 2
    above = np.where(y_corr > half_max)[0]

    if len(above) > 2:
        fwhm = x[above[-1]] - x[above[0]]
    else:
        fwhm = (x[-1] - x[0]) / 10  # fallback

    # --- 4. Estimate trion contribution r ---
    # look at intensity at EX - delta
    target_E = EX - delta
    idx_trion = np.argmin(np.abs(x - target_E))

    I_trion = y_corr[idx_trion]

    # crude ratio estimate
    r = I_trion / A

    # constrain to physical range
    r = np.clip(r, 0.05, 1.0)

    # --- 5. Final amplitude scaling ---
    # A should represent total prefactor
    # since model is A*(1 + r*...), reduce slightly
    A = A / (1 + r)
    return A, EX, fwhm, r, B

def estimate_voigt_params(x, y):
    """
    Estimate good starting values for fitting a Voigt profile.
    
    Parameters:
    - x: Array of x values.
    - y: Array of y values.
    
    Returns:
    - amp_init: Estimated amplitude of the Voigt peak.
    - cen_init: Estimated center (location) of the Voigt peak.
    - wid_init: Estimated width (FWHM) of the Voigt peak.
    - gamma_init: Estimated Lorentzian width (gamma), set to a small value initially.
    """
    # Estimate amplitude as the max value of y
    amp_init = np.max(y)
    
    # Estimate center as the x value where y is maximum
    cen_init = x[np.argmax(y)]
    
    # Estimate FWHM (Full Width at Half Maximum)
    half_max = amp_init / 2
    
    # Find where y crosses half_max to calculate the width
    left_idx = np.where(y > half_max)[0][0]  # First index where y > half max
    right_idx = np.where(y > half_max)[0][-1]  # Last index where y > half max
    fwhm = x[right_idx] - x[left_idx]  # Approximate FWHM as width
    
    # Set initial width (wid) as the FWHM
    wid_init = fwhm
    
    # Set an initial small value for gamma
    gamma_init = 0.1  # You can adjust this based on data
    
    return amp_init, cen_init, wid_init, gamma_init

def estimate_double_voigt_params(x, y):
    """
    Estimate good starting values for fitting a double Voigt profile.
    
    Parameters:
    - x: Array of x values.
    - y: Array of y values.
    
    Returns:
    - amp_init: List of estimated amplitudes for the two Voigt peaks.
    - cen_init: List of estimated centers (locations) for the two Voigt peaks.
    - wid_init: List of estimated widths (FWHM) for the two Voigt peaks.
    - gamma_init: List of estimated Lorentzian widths (gamma) for both peaks.
    """
    # Find two prominent peaks in the data
    peaks, _ = find_peaks(y, height=0.1, distance=10)  # Adjust height and distance if needed
    
    if len(peaks) < 2:
        raise ValueError("Less than two peaks found in the data. Make sure the data contains two distinct peaks.")

    # Sort the peaks by prominence (highest peaks first)
    sorted_peaks = peaks[np.argsort(y[peaks])][::-1]  # Sort peaks by height, highest first
    
    amp_init = []
    cen_init = []
    wid_init = []
    gamma_init = []

    # Estimate parameters for each peak
    for peak_idx in sorted_peaks[:2]:  # Only consider the top 2 peaks
        amp = y[peak_idx]  # Amplitude is the peak height
        cen = x[peak_idx]  # Center is the x-value at the peak

        # Estimate FWHM (Full Width at Half Maximum)
        half_max = amp / 2
        left_idx = np.where(y > half_max)[0][0]  # First index where y > half max
        right_idx = np.where(y > half_max)[0][-1]  # Last index where y > half max
        fwhm = x[right_idx] - x[left_idx]  # Approximate FWHM as width

        # Set initial values for this peak
        amp_init.append(amp)
        cen_init.append(cen)
        wid_init.append(fwhm)
        gamma_init.append(0.1)  # Small initial value for gamma
    
    return amp_init[0], cen_init[0], wid_init[0], gamma_init[0], amp_init[1], cen_init[1], wid_init[1], gamma_init[1]

def estimate_double_gaussian_params(x, y):
    """
    Estimate good starting values for fitting a double Gaussian function.
    
    Parameters:
    - x: Array of x values.
    - y: Array of y values.
    
    Returns:
    - A1_init, mu1_init, sigma1_init: Estimated parameters for the first Gaussian.
    - A2_init, mu2_init, sigma2_init: Estimated parameters for the second Gaussian.
    """
    # Find two prominent peaks in the data
    peaks, _ = find_peaks(y, height=0.1, distance=len(x) // 5)  # Adjust parameters as necessary

    if len(peaks) < 2:
        raise ValueError("Less than two peaks found in the data. Ensure the data has two distinct peaks.")

    # Sort the peaks by height (prominence), highest peaks first
    sorted_peaks = peaks[np.argsort(y[peaks])][::-1]

    # Initialize arrays for the parameters
    amp_init = []
    cen_init = []
    sigma_init = []

    # Estimate parameters for each peak (top 2 peaks)
    for peak_idx in sorted_peaks[:2]:
        # Amplitude is the height at the peak
        A = y[peak_idx]
        
        # Center is the x-value at the peak
        mu = x[peak_idx]

        # Estimate FWHM by finding where the peak crosses half its maximum
        half_max = A / 2
        left_idx = np.where(y[:peak_idx] < half_max)[0][-1]  # Left crossing
        right_idx = np.where(y[peak_idx:] < half_max)[0][0] + peak_idx  # Right crossing

        fwhm = x[right_idx] - x[left_idx]
        
        # Estimate sigma from FWHM (Gaussian relationship: FWHM ≈ 2.355 * sigma)
        sigma = fwhm / 2.355

        # Append initial estimates
        amp_init.append(A)
        cen_init.append(mu)
        sigma_init.append(sigma)

    # Return initial parameters for both Gaussians
    A1_init, mu1_init, sigma1_init = amp_init[0], cen_init[0], sigma_init[0]
    A2_init, mu2_init, sigma2_init = amp_init[1], cen_init[1], sigma_init[1]
    
    return A1_init, mu1_init, sigma1_init, A2_init, mu2_init, sigma2_init

def estimate_double_lorentz_params(x, y):
    """
    Estimate good starting values for fitting a double Lorentzian function.
    
    Parameters:
    - x: Array of x values.
    - y: Array of y values.
    
    Returns:
    - A1_init, mu1_init, gamma1_init: Estimated parameters for the first Lorentzian.
    - A2_init, mu2_init, gamma2_init: Estimated parameters for the second Lorentzian.
    """
    # Find two prominent peaks in the data
    peaks, _ = find_peaks(y, height=0.1, distance=len(x) // 5)  # Adjust parameters if necessary

    if len(peaks) < 2:
        raise ValueError("Less than two peaks found in the data. Ensure the data has two distinct peaks.")

    # Sort peaks by height (prominence), highest peaks first
    sorted_peaks = peaks[np.argsort(y[peaks])][::-1]

    # Initialize arrays for the parameters
    amp_init = []
    cen_init = []
    gamma_init = []

    # Estimate parameters for each peak (top 2 peaks)
    for peak_idx in sorted_peaks[:2]:
        # Amplitude is the height at the peak
        A = y[peak_idx]
        
        # Center is the x-value at the peak
        mu = x[peak_idx]

        # Estimate HWHM by finding where the peak crosses half its maximum
        half_max = A / 2
        left_idx = np.where(y[:peak_idx] < half_max)[0][-1]  # Left crossing
        right_idx = np.where(y[peak_idx:] < half_max)[0][0] + peak_idx  # Right crossing

        fwhm = x[right_idx] - x[left_idx]
        
        # For Lorentzian, the HWHM is directly related to gamma
        gamma = fwhm / 2

        # Append initial estimates
        amp_init.append(A)
        cen_init.append(mu)
        gamma_init.append(gamma)

    # Return initial parameters for both Lorentzians
    A1_init, mu1_init, gamma1_init = amp_init[0], cen_init[0], gamma_init[0]
    A2_init, mu2_init, gamma2_init = amp_init[1], cen_init[1], gamma_init[1]
    
    return A1_init, mu1_init, gamma1_init, A2_init, mu2_init, gamma2_init

# fit window functions to data  
def fitdoublegaussiantospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    if guess is None:
        initialguess = estimate_double_gaussian_params(x, y)
    else:
        initialguess = guess[0:6]
    #[np.max(y), x[np.argmax(y)], np.std(x), np.max(y), x[np.argmax(y)], np.std(x)] #old fit estimate
    fitdata, pcov = curve_fit(double_gaussianwind, x, y, p0=initialguess, maxfev=maxfev)
    amp1_fit, cen1_fit, wid1_fit, amp2_fit, cen2_fit, wid2_fit = fitdata
    return amp1_fit, cen1_fit, wid1_fit, amp2_fit, cen2_fit, wid2_fit, pcov

def fitdoublelorentztospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    if guess is None:
        initialguess = estimate_double_lorentz_params(x, y)
    else:
        initialguess = guess[0:6]
    #[np.max(y), x[np.argmax(y)], np.std(x), np.max(y), x[np.argmax(y)], np.std(x)] old fit estimate
    fitdata, pcov = curve_fit(double_lorentzwind, x, y, p0=initialguess, maxfev=maxfev)
    amp1_fit, cen1_fit, wid1_fit, amp2_fit, cen2_fit, wid2_fit = fitdata
    return amp1_fit, cen1_fit, wid1_fit, amp2_fit, cen2_fit, wid2_fit, pcov

def fitdoublevoigttospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    if guess is None:
        initialguess = estimate_double_voigt_params(x, y)
    else:
        initialguess = guess[0:8]
    #[np.max(y), x[np.argmax(y)], np.std(x), np.std(x), np.max(y), x[np.argmax(y)], np.std(x), np.std(x)] old fit estimate
    fitdata, pcov = curve_fit(double_voigtwind, x, y, p0=initialguess, maxfev=maxfev,
                              # tighten fith convergence by setting bounds
                              ftol = 1e-10, xtol=1e-10, gtol=1e-10,
                              )
    amp1_fit, cen1_fit, wid1_fit, gamma1_fit, amp2_fit, cen2_fit, wid2_fit, gamma2_fit = fitdata
    return amp1_fit, cen1_fit, wid1_fit, gamma1_fit, amp2_fit, cen2_fit, wid2_fit, gamma2_fit, pcov

def fitvoigttospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    #initialguess = [np.max(y), x[np.argmax(y)], np.std(x)*2, np.std(x)]
    if guess is None:
        initialguess = estimate_voigt_params(x, y)
    else:
        initialguess = guess[0:4]
    fitdata, pcov = curve_fit(voigtwind, x, y, p0=initialguess, maxfev=maxfev)
    amp_fit, cen_fit, wid_fit, gamma_fit = fitdata
    return amp_fit, cen_fit, wid_fit, gamma_fit, pcov

def fitlorentztospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    if guess is None:
        initialguess = [np.max(y), x[np.argmax(y)], np.std(x)*2.4]
    else:
        initialguess = guess[0:3]
    fitdata, pcov = curve_fit(lorentzwind, x, y, p0=initialguess, maxfev=maxfev)
    amp_fit, cen_fit, wid_fit = fitdata
    return amp_fit, cen_fit, wid_fit, pcov

def fitgaussiantospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    if guess is None:
        initialguess = [np.max(y), x[np.argmax(y)], np.std(x)]
    else:
        initialguess = guess[0:3]
    #print('initialguess:', initialguess)
    #print('generated guess:', [np.max(y), x[np.argmax(y)], np.std(x)])
    fitdata, pcov = curve_fit(gaussianwind, x, y, p0=initialguess, maxfev=maxfev)
    amp_fit, cen_fit, wid_fit = fitdata
    return amp_fit, cen_fit, wid_fit, pcov
    
def fitlinetospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    a_guess = (np.max(y) - np.min(y)) / (np.max(x) - np.min(x))
    b_guess = np.mean(y) - a_guess * np.mean(x)
    initial_guess = [a_guess, b_guess]
    fitdata, pconf = curve_fit(linearwind, x, y, p0=initial_guess, maxfev=maxfev)
    a, b = fitdata
    return a, b, pconf

def fitdoublevoigt5paramtospec(start, end, WL, PLB, maxfev=10000, guess=None):
    x = WL[start: end]
    y = PLB[start: end]
    if guess is None:
        initialguess = doublevoigt_5param_approx(x, y)
    else:
        initialguess = guess[0:5]
    fitdata, pcov = curve_fit(doublevoigt5param, x, y, p0=initialguess, maxfev=maxfev)
    A_fit, EX_fit, fwhm_fit, r_fit, B_fit = fitdata
    return A_fit, EX_fit, fwhm_fit, r_fit, B_fit, pcov

# get maxima of fitted functions
def getmaxdoublegaussian(xmin, xmax, amp1, cen1, wid1, amp2, cen2, wid2):
    # use newton method to find max
    #x = Newtonmax(lambda x: -double_gaussianwind(x, amp1, cen1, wid1, amp2, cen2, wid2), cen1, tol=1e-6, maxiter=10000, xmin=xmin, xmax=xmax)
    #y = double_gaussianwind(xmax, amp1, cen1, wid1, amp2, cen2, wid2)
    x, y = Maxbyinsert(lambda x: -double_gaussianwind(x, amp1, cen1, wid1, amp2, cen2, wid2), [xmin, xmax], 10000)
    #success, x, fun = find_max_of_fit(lambda x: -double_gaussianwind(x, amp1, cen1, wid1, amp2, cen2, wid2), xmin=xmin, xmax=xmax)
    return x, y # -fun

def getmaxdoublelorentz(xmin, xmax, amp1, cen1, wid1, amp2, cen2, wid2):
    #x = Newtonmax(lambda x: -double_lorentzwind(x, amp1, cen1, wid1, amp2, cen2, wid2), cen1, tol=1e-6, maxiter=10000, xmin=xmin, xmax=xmax)
    #y = double_lorentzwind(x, amp1, cen1, wid1, amp2, cen2, wid2)
    x, y = Maxbyinsert(lambda x: -double_lorentzwind(x, amp1, cen1, wid1, amp2, cen2, wid2), [xmin, xmax], 10000)
    #success, x, fun = find_max_of_fit(lambda x: -double_lorentzwind(x, amp1, cen1, wid1, amp2, cen2, wid2), xmin=xmin, xmax=xmax)
    return x, y #-fun

def getmaxdoublevoigt(xmin, xmax, amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2):
    #x = Newtonmax(lambda x: -double_voigtwind(x, amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2), cen1, tol=1e-6, maxiter=10000, xmin=xmin, xmax=xmax)
    #y = double_voigtwind(x, amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2)
    x, y = Maxbyinsert(lambda x: -double_voigtwind(x, amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2), [xmin, xmax], 10000)
    #success, x, fun = find_max_of_fit(lambda x: -double_voigtwind(x, amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2), xmin=xmin, xmax=xmax)
    return x, y
    
# max of 1D Fit functions can be obtained by reading their parameters
def getmaxgaussian(xmin, xmax, amp, cen, wid):
    return cen, amp

def getmaxlorentz(xmin, xmax, amp, cen, wid):
    return cen, amp

def getmaxvoigt(xmin, xmax, amp, cen, wid, gamma):
    return cen, amp

def getmaxlinear(xmin, xmax, a, b):
    pass

def getmaxdoublevoigt5param(xmin, xmax, A, EX, fwhm, r, B):
    #x = Newtonmax(lambda x: -doublevoigt_5param(x, A, EX, fwhm, r, B), EX, tol=1e-6, maxiter=10000, xmin=xmin, xmax=xmax)
    #y = doublevoigt_5param(x, A, EX, fwhm, r, B)
    x, y = Maxbyinsert(lambda x: -doublevoigt5param(x, A, EX, fwhm, r, B), [xmin, xmax], 10000)
    #success, x, fun = find_max_of_fit(lambda x: -doublevoigt_5param(x, A, EX, fwhm, r, B), xmin=xmin, xmax=xmax)
    return x, y #-fun

# 2D Pixmatrix Fit functions
# 2d gaussian function
def gaussian2drot(coords, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x, y = coords
    x, y = x - xo, y - yo
    a = (np.cos(theta)**2) / (2 * sigma_x**2) + (np.sin(theta)**2) / (2 * sigma_y**2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x**2) + (np.sin(2 * theta)) / (4 * sigma_y**2)
    c = (np.sin(theta)**2) / (2 * sigma_x**2) + (np.cos(theta)**2) / (2 * sigma_y**2)
    g = amplitude * np.exp(-(a * x**2 + 2 * b * x * y + c * y**2)) + offset
    return g.ravel()

# 2d gaussian function
def gaussian2d(coords, amplitude, xo, yo, sigma_x, sigma_y, offset):
    x, y = coords
    x, y = x - xo, y - yo
    scale = amplitude / (2 * np.pi * sigma_x * sigma_y)
    g = scale * np.exp(-(x**2 / (2 * sigma_x**2) + y**2 / (2 * sigma_y**2))) + offset
    return g.ravel()

def fitgaussiand2dtomatrixrot(inpdata, plotfit, gdx, gdy, colormap, maxfev=10000):
    data = np.array(inpdata)
    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])
    x, y = np.meshgrid(x, y)
    initialguess = [np.max(data), np.argmax(data) % data.shape[1], np.argmax(data) // data.shape[1], 1, 1, 0, 0]
    popt, pcov = curve_fit(gaussian2drot, (x, y), data.ravel(), p0=initialguess, maxfev=maxfev)
    #fwhmx = 2 * np.sqrt(2 * np.log(2)) * fitdata[3]
    #fwhmy = 2 * np.sqrt(2 * np.log(2)) * fitdata[4]
    data_fited = gaussian2drot((x, y), *popt).reshape(data.shape)
    fwhmx = np.sqrt(2 * np.log(2)) * popt[3] * gdx * 2 
    fwhmy = np.sqrt(2 * np.log(2)) * popt[4] * gdy * 2 
    print('FWHM X = {} mum, FWHM Y = {} mum'.format(round(fwhmx, 2), round(fwhmy, 2)))
    # print when the fited function falls below 1/e of the maximum
    #beamx = np.where(data_fited > np.max(data_fited) / np.e)[1]*gdx
    #beamy = np.where(data_fited > np.max(data_fited) / np.e)[0]*gdx
    #print('Beam X = {} mum, Beam Y = {} mum'.format(np.amax(beamx)-np.amin(beamx), np.amax(beamy)-np.amin(beamy)))

    #calculate beam waist
    x0 = popt[1] * gdx
    sigma_x = popt[3] * gdx
    amplitude = popt[0]
    y0 = popt[2] * gdy
    sigma_y = popt[4] * gdy
    bwx_start = x0 - sigma_x*np.sqrt(2)
    bwx_end = x0 + sigma_x*np.sqrt(2)
    bwy_start = y0 - sigma_y*np.sqrt(2)
    bwy_end = y0 + sigma_y*np.sqrt(2)
    print('Beam Waist X = {} mum, Beam Waist Y = {} mum'.format(abs(bwx_end-bwx_start), abs(bwy_end-bwy_start)))

    # plot data_fited
    if plotfit == True:
        fig, ax = plt.subplots()
        plt.imshow(data_fited, cmap=colormap)

        # Get current ticks
        current_xticks = np.arange(data.shape[1])
        current_yticks = np.arange(data.shape[0])
        # Multiply ticks by constants
        new_xticks = np.multiply(current_xticks, gdx).round(3)
        new_yticks = np.multiply(current_yticks, gdy).round(3)
        # Set new ticks
        ax.set_xticks(current_xticks)  # Set the ticks to be at the indices of the current ticks
        ax.set_yticks(current_yticks)
        # Set tick labels to the new values
        ax.set_xticklabels(new_xticks)
        ax.set_yticklabels(new_yticks)
        # Set the axis labels auto adjust
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  # Adjust the number of bins to fit the plot size
        # Set the title and labels
        ax.set_title('Fitted Gaussian 2D')
        ax.set_xlabel('X-axis in \u03bcm')
        ax.set_ylabel('Y-axis in \u03bcm')
        # label the colormap with the units counts
        cbar = plt.colorbar()
        cbar.set_label('Counts')

        plt.show()

def fitgaussian2dtomatrix(inpdata, plotfit, gdx, gdy, colormap, maxfev=10000):
    data = np.array(inpdata)
    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])
    x, y = np.meshgrid(x, y)
    initialguess = [np.max(data), np.argmax(data) % data.shape[1], np.argmax(data) // data.shape[1], 1, 1, 0]
    popt, pcov = curve_fit(gaussian2d, (x, y), data.ravel(), p0=initialguess, maxfev=maxfev)

    data_fited = gaussian2d((x, y), *popt).reshape(data.shape)
    fwhmx = np.sqrt(2 * np.log(2)) * popt[3] * gdx * 2 
    fwhmy = np.sqrt(2 * np.log(2)) * popt[4] * gdy * 2 
    print('FWHM X = {} mum, FWHM Y = {} mum'.format(round(fwhmx, 2), round(fwhmy, 2)))

    if plotfit == True:
        fig, ax = plt.subplots()  
        plt.imshow(data_fited, cmap=colormap)
        # Get current ticks
        current_xticks = np.arange(data.shape[1])
        current_yticks = np.arange(data.shape[0])
        # Multiply ticks by constants
        new_xticks = np.multiply(current_xticks, gdx).round(3)
        new_yticks = np.multiply(current_yticks, gdy).round(3)
        # Set new ticks
        ax.set_xticks(current_xticks)  # Set the ticks to be at the indices of the current ticks
        ax.set_yticks(current_yticks)
        # Set tick labels to the new values
        ax.set_xticklabels(new_xticks)
        ax.set_yticklabels(new_yticks)
        # Set the axis labels auto adjust
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  # Adjust the number of bins to fit the plot size
        # Set the title and labels
        ax.set_title('Fitted Gaussian 2D')
        ax.set_xlabel('X-axis in \u03bcm')
        ax.set_ylabel('Y-axis in \u03bcm')
        # label the colormap with the units counts
        cbar = plt.colorbar()
        cbar.set_label('Counts')
        plt.show()

    #calculate beam waist
    x0 = popt[1] * gdx
    sigma_x = popt[3] * gdx
    amplitude = popt[0]
    y0 = popt[2] * gdy
    sigma_y = popt[4] * gdy
    bwx_start = x0 - sigma_x*np.sqrt(2)
    bwx_end = x0 + sigma_x*np.sqrt(2)
    bwy_start = y0 - sigma_y*np.sqrt(2)
    bwy_end = y0 + sigma_y*np.sqrt(2)
    print('Beam Waist X = {} mum, Beam Waist Y = {} mum'.format(abs(bwx_end-bwx_start), abs(bwy_end-bwy_start)))

# FWHM functions
def getlorentzfwhm(fitparams):
    # returns the FWHM of a lorentzian fit
    # fitparams: amp, cen, wid
    amp, cen, wid = fitparams[0:3]
    fwhm = 2 * wid
    return fwhm

def getgaussianfwhm(fitparams):
    # returns the FWHM of a gaussian fit
    # fitparams: amp, cen, wid
    amp, cen, wid = fitparams[0:3]
    fwhm = 2 * np.sqrt(2 * np.log(2)) * wid
    return fwhm

def getvoigtfwhm(fitparams):
    # returns the FWHM of a voigt fit
    # fitparams: amp, cen, wid, gamma
    amp, cen, wid, gamma = fitparams[0:4]
    # FWHM of a Voigt profile is not straightforward, so it is calculated on a grid using the fwhmbygrid function
    # 
    fwhm = fwhmbygrid(wid, cen - 5 * wid, cen + 5 * wid, npoints=10000)
    return fwhm

def getdoublegaussianfwhm(fitparams):
    # returns the FWHM of a double gaussian fit
    # fitparams: amp1, cen1, wid1, amp2, cen2, wid2
    amp1, cen1, wid1, amp2, cen2, wid2 = fitparams[0:6]
    fwhm1 = 2 * np.sqrt(2 * np.log(2)) * wid1
    fwhm2 = 2 * np.sqrt(2 * np.log(2)) * wid2
    return fwhm1 + fwhm2

def getdoublelorentzfwhm(fitparams):
    # returns the FWHM of a double lorentzian fit
    # fitparams: amp1, cen1, wid1, amp2, cen2, wid2
    amp1, cen1, wid1, amp2, cen2, wid2 = fitparams[0:6]
    fwhm1 = 2 * wid1
    fwhm2 = 2 * wid2
    return fwhm1 + fwhm2

def getdoublevoigtfwhm(xgrid, fitparams):
    amp1, cen1, wid1, gamma1, amp2, cen2, wid2, gamma2 = fitparams

    y = (
        voigtwind(xgrid, amp1, cen1, wid1, gamma1) +
        voigtwind(xgrid, amp2, cen2, wid2, gamma2)
    )

    ymax = np.max(y)
    half = ymax / 2

    idx = np.where(y >= half)[0]
    if len(idx) < 2:
        return np.nan

    i1, i2 = idx[0], idx[-1]

    # linear interpolation (left crossing)
    x1 = np.interp(half, [y[i1-1], y[i1]], [xgrid[i1-1], xgrid[i1]])

    # right crossing
    x2 = np.interp(half, [y[i2], y[i2+1]], [xgrid[i2], xgrid[i2+1]])

    return x2 - x1

def getdoublevoigt5paramfwhm(fitparam): #x, A, EX, fwhm, r, B, delta=30, eta=0.5) # let do this without xgrid as parameter, we can use the fitparam to generate the xgrid
    A, EX, fwhm, r, B = fitparam # older versions also had delta and eta here, but we will ignore them for now
    # delta and eta 
    x = np.linspace(EX - 5 * fwhm, EX + 5 * fwhm, 10000)

    y = doublevoigt5param(x, A, EX, fwhm, r, B)#, delta, eta) 

    ymax = np.max(y)
    half = ymax / 2

    # find crossings
    indices = np.where(y >= half)[0]
    if len(indices) < 2:
        return np.nan

    i1, i2 = indices[0], indices[-1]

    # use x to find the corresponding x values for the half maximum
    x1 = np.interp(half, [y[i1-1], y[i1]], [x[i1-1], x[i1]])
    x2 = np.interp(half, [y[i2], y[i2+1]], [x[i2], x[i2+1]])

    return x2 - x1

def fwhmbygrid(f, wlstart, wlend, npoints=100000):
    # Generate x values
    x = np.linspace(wlstart, wlend, npoints)
    # Calculate y values
    y = f(x)
    # Find the maximum y value
    ymax = np.max(y)
    # Find the half maximum value
    half_max = ymax / 2
    # Find where y crosses half_max
    left_idx = np.where(y > half_max)[0][0]  # First index where y > half max
    right_idx = np.where(y > half_max)[0][-1]  # Last index where y > half max
    # Calculate FWHM
    fwhm = x[right_idx] - x[left_idx]
    return fwhm

# newton's method to find the maximum of a function
def Newtonmax(f, x0, tol=1e-8, maxiter=10000, xmin=0, xmax=1000):
    # Initialize the iteration counter
    iter = 0
    # Initialize the error
    error = tol + 1
    # Iterate until the error is less than the tolerance or the maximum number of iterations is reached
    while error > tol and iter < maxiter and x0 > xmin and x0 < xmax:
        # Calculate the derivative of the function at the current point
        f_prime = (f(x0 + tol) - f(x0 - tol)) / (2 * tol)
        # Calculate the next point using Newton's method
        x1 = x0 - f(x0) / f_prime
        # Calculate the error
        error = np.abs(x1 - x0)
        # Update the current point
        x0 = x1
        # Increment the iteration counter
        iter += 1
    return x0

def Maxbyinsert(f, wl, npoints):
    # Find the maximum of a function by inserting points
    # f: function to be maximized
    # wl: wavelength range
    # npoints: number of points to be inserted
    x = np.linspace(wl[0], wl[1], npoints)
    y = f(x)
    max_index = np.argmin(y)
    return x[max_index], y[max_index]

# another method to find the maximum of a function in an interval
def find_max_of_fit(fitfunc, tol=1e-6, maxiter=10000, xmin=500, xmax=600):
    #Find the maximum of a fitted function within the interval [xmin, xmax].
    notfound = True
    if xmax < xmin+10:
        xmax = xmin+100
    x0 = xmin
    while notfound == True: # if the maximum is not found, increase the initial guess
        result = minimize(lambda x: -fitfunc(x), x0=550, bounds=[(xmin, xmax)])
        if result.success:
            notfound = False
        else:
            if x0 > xmax:
                return False, x0, -fitfunc(x0)
            x0 += 10
    #print(result.x, -result.fun)
    #return result.success, result.x, -result.fun
    return result.success, result.x, fitfunc(result.x) # type: ignore   

# calculate the r_squared between fit and data
# r_squared = 1 - (ss_res / ss_tot) 
def calc_r_squared(fit, data): # args: data(observations), y_fit(model)
    # Note: lib9.py passes (observations, model) as (fit, data)
    # So 'fit' argument contains observations, 'data' argument contains model
    
    # Convert to numpy arrays to handle lists and ensure correct math
    obs = np.array(fit)
    model = np.array(data)
    
    # Calculate sums of squares
    # ss_res = sum((y_obs - y_model)^2)
    ss_res = np.sum((obs - model)**2)
    
    # ss_tot = sum((y_obs - mean(y_obs))^2)
    # We must calculate variance on OBSERVATIONS, not model
    ss_tot = np.sum((obs - np.mean(obs))**2)
    
    if ss_tot == 0:
        r_squared = 0.0
    else:
        r_squared = 1 - (ss_res / ss_tot)
        
    return r_squared, ss_res, ss_tot
# ss_res is the sum of the squared residuals
# ss_tot is the total sum of squares

# + is for ss_res, ss_tot, r_squared, pixstart, pixend, wlstart, wlend, fwhm, max_x, max_y
addtofitparms = ['ss_res', 'ss_tot', 'r_squared', 'fwhm', 'pixstart', 'pixend', 'wlstart', 'wlend', 'max_x', 'max_y', 'fit_status'] # Note: this one is essential to exist
unitstoaddfit = ['Counts', 'Counts', '', 'nm', 'nm', 'nm', 'nm', 'nm', 'nm', ''] # add further units to the array after r_squared
# add further parameters to the array after r_squared
def buildfitparas():
    fa = []
    for i in fitkeys.keys():
        fa.append([])#[np.nan])
        for j in range(0, fitkeys[i][4]+len(addtofitparms)):
            fa[-1].append(np.nan)
    return fa

# return a list of all fit parameters by their names
def getlistofallFitparameters():
    fl = []
    for i in fitkeys.keys():
        fl.append([])
        for j in range(0, fitkeys[i][4]):
            fl[-1].append(i + ' ' + str(j))
        for j in range(0, len(addtofitparms)):
            fl[-1].append(i + ' ' + addtofitparms[j])
    return fl

def getlistofallFitparaminone():
    fl = []
    for i in fitkeys.keys():
        for j in range(0, fitkeys[i][4]):
            fl.append(i + ' ' + str(j))
        for j in range(0, len(addtofitparms)):
            fl.append(i + ' ' + addtofitparms[j])
    return fl
    
def getindexofFitparameter(fl, fitname):# fl is a list of all fit parameters, returns the index of the fit parameter in the list
    for i in range(0, len(fl)):
        if fitname in fl[i]:
            return i, fl[i].index(fitname)
    return -1, -1

# obtain oscillations
# if a the signal is overlayed by a faster oscillation, the oscillation period is being isolated

def smooth_background(signal, window_length=51, polyorder=3):
    """
    Smooth the signal to extract the background using Savitzky-Golay filter.
    
    Parameters:
    - signal: Array of signal values
    - window_length: Length of the filter window (must be odd)
    - polyorder: Order of the polynomial used to fit the samples
    
    Returns:
    - Smoothed background signal
    """    
    # Ensure window_length is odd and less than signal length
    if window_length % 2 == 0:
        window_length += 1
    if window_length > len(signal):
        window_length = len(signal) if len(signal) % 2 == 1 else len(signal) - 1
    if window_length < polyorder + 2:
        window_length = polyorder + 2
        if window_length % 2 == 0:
            window_length += 1
    
    # Apply Savitzky-Golay filter for smoothing
    background = savgol_filter(signal, window_length, polyorder)
    return background

def extract_phase_evolution(oscillations, wl, maxima_indices, minima_indices):
    """
    Extract the phase evolution (instantaneous frequency) from oscillations.
    
    The phase changes as the oscillations speed up or slow down over the energy/wavelength range.
    
    Parameters:
    - oscillations: Array of isolated oscillations
    - wl: Wavelength/energy array
    - maxima_indices: Indices of maxima in oscillations
    - minima_indices: Indices of minima in oscillations
    
    Returns:
    - instantaneous_freq: Array of instantaneous frequency at each peak
    - freq_centers: Energy/wavelength values where frequency is calculated
    - phase_trend_slope: Linear fit slope of frequency evolution (chirp rate)
    - phase_trend_intercept: Linear fit intercept
    - mean_period: Mean oscillation period
    - period_std: Standard deviation of period (measure of chirp)
    """
    
    # Method 1: Calculate instantaneous frequency from peak spacing
    # Combine maxima and minima for better frequency estimation
    all_peaks = np.sort(np.concatenate([maxima_indices, minima_indices]))
    
    if len(all_peaks) < 3:
        # Not enough peaks to estimate phase evolution
        return (np.array([]), np.array([]), 0.0, 0.0, 0.0, 0.0)
    
    # Calculate period (spacing) between consecutive peaks
    peak_wl = wl[all_peaks]
    periods = np.diff(peak_wl)  # Period in energy/wavelength units
    
    # Instantaneous frequency is 1/period
    instantaneous_freq = 1.0 / periods
    freq_centers = (peak_wl[:-1] + peak_wl[1:]) / 2  # Midpoint between peaks
    
    # Fit linear trend to frequency evolution: freq(E) = slope * E + intercept
    # This captures the "chirp" or phase evolution
    if len(freq_centers) > 1:
        freq_center_mean = np.mean(freq_centers)
        coeffs = np.polyfit(freq_centers - freq_center_mean, instantaneous_freq, 1)
        phase_trend_slope = coeffs[0]  # How fast frequency changes with energy
        phase_trend_intercept = coeffs[1] - coeffs[0] * freq_center_mean  # Initial frequency
    else:
        phase_trend_slope = 0.0
        phase_trend_intercept = instantaneous_freq[0] if len(instantaneous_freq) > 0 else 0.0
    
    # Calculate mean and std of period
    mean_period = np.mean(periods)
    period_std = np.std(periods)
    
    # Method 2: Use Hilbert transform for continuous phase (stored for potential future use)
    try:
        analytic_signal = hilbert(oscillations)
        instantaneous_phase = np.unwrap(np.angle(analytic_signal))
        # Instantaneous frequency from phase derivative
        hilbert_freq = np.gradient(instantaneous_phase) / np.gradient(wl) / (2 * np.pi)
    except:
        hilbert_freq = None
    
    return (instantaneous_freq, freq_centers, phase_trend_slope, 
            phase_trend_intercept, mean_period, period_std)

def fitoscillationtospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Extract oscillations from spectrum and return fit parameters including phase evolution.
    
    Parameters:
    - start: Start index of the spectrum region
    - end: End index of the spectrum region
    - WL: Wavelength/energy array
    - PLB: Signal intensity array
    - maxfev: Not used, kept for compatibility
    - guess: Optional dict with 'window_length', 'polyorder', 'prominence', 'distance'
    
    Returns:
    - Tuple of fit parameters compatible with fitkeys structure:
      Scalar parameters (12 values):
        0: background_mean, 1: osc_amplitude, 2: n_maxima, 3: n_minima, 
        4: mean_max, 5: mean_min, 6: freq_est, 7: osc_range,
        8: phase_chirp (frequency evolution slope), 9: initial_freq (frequency at start),
        10: mean_period, 11: period_std (variability)
      Array parameters:
        background, oscillations, maxima_indices, minima_indices,
        maxima_wl, minima_wl, maxima_values, minima_values,
        instantaneous_freq, freq_centers
    """
    x = WL[start:end]
    y = PLB[start:end]
    
    # Parse optional parameters from guess
    window_length = 51
    polyorder = 3
    prominence = None
    distance = 5
    
    if guess is not None and isinstance(guess, dict):
        window_length = guess.get('window_length', 51)
        polyorder = guess.get('polyorder', 3)
        prominence = guess.get('prominence', None)
        distance = guess.get('distance', 5)
    
    # Call isolate_oscillation
    (background, oscillations, 
     maxima_indices, minima_indices,
     maxima_wl, minima_wl,
     maxima_values, minima_values) = isolate_oscillation(
        y, x, window_length, polyorder, prominence, distance
    )
    
    # Calculate summary statistics for fitkeys compatibility
    background_mean = np.mean(background)
    osc_amplitude = np.std(oscillations)
    n_maxima = len(maxima_indices)
    n_minima = len(minima_indices)
    mean_max = np.mean(maxima_values) if len(maxima_values) > 0 else 0
    mean_min = np.mean(minima_values) if len(minima_values) > 0 else 0
    
    # Estimate oscillation frequency (peaks per unit energy/wavelength)
    if n_maxima > 1:
        freq_est = n_maxima / (x[-1] - x[0])
    else:
        freq_est = 0
    
    osc_range = np.max(oscillations) - np.min(oscillations) if len(oscillations) > 0 else 0
    
    # Extract phase evolution (NEW!)
    (instantaneous_freq, freq_centers, phase_chirp, 
     initial_freq, mean_period, period_std) = extract_phase_evolution(
        oscillations, x, maxima_indices, minima_indices
    )
    
    # Return scalar parameters first (12 values for fitkeys compatibility), then arrays
    return (background_mean, osc_amplitude, n_maxima, n_minima, 
            mean_max, mean_min, freq_est, osc_range,
            phase_chirp, initial_freq, mean_period, period_std,  # NEW phase parameters
            background, oscillations, 
            maxima_indices, minima_indices,
            maxima_wl, minima_wl,
            maxima_values, minima_values,
            instantaneous_freq, freq_centers)  # NEW phase arrays

def getoscillationdata(xmin, xmax, *fitparams):
    """
    Extract oscillation data from fit parameters.
    This is a placeholder to maintain compatibility with fitkeys structure.
    
    Returns:
    - None (oscillation data is already in the fit parameters)
    """
    return None

def isolate_oscillation(signal, wl, window_length=51, polyorder=3, prominence=None, distance=5): 
    """
    Isolate oscillations from a signal that consists of a peak overlayed by faster oscillations.
    
    Parameters:
    - signal: Array of signal values (DR Signal)
    - wl: Array of wavelength/energy values corresponding to the signal
    - window_length: Length of the smoothing filter window (default: 51, must be odd)
    - polyorder: Order of polynomial for smoothing (default: 3)
    - prominence: Required prominence of peaks (default: auto-calculated)
    - distance: Minimum distance between peaks (default: 5)
    
    Returns:
    - background: Smoothed background signal
    - oscillations: Isolated oscillations (signal - background)
    - maxima_indices: Indices of maxima in the oscillations
    - minima_indices: Indices of minima in the oscillations
    - maxima_wl: Wavelength/energy values at maxima
    - minima_wl: Wavelength/energy values at minima
    - maxima_values: Signal values at maxima
    - minima_values: Signal values at minima
    """
    
    # Step 1: Smooth the signal to get the background
    background = smooth_background(signal, window_length, polyorder)
    
    # Step 2: Extract oscillations by subtracting background
    oscillations = signal - background
    
    # Step 3: Find maxima in the oscillations
    if prominence is None:
        # Auto-calculate prominence as a fraction of the oscillation range
        osc_range = np.max(oscillations) - np.min(oscillations)
        prominence = osc_range * 0.1  # 10% of the range
    
    maxima_indices, _ = find_peaks(oscillations, prominence=prominence, distance=distance)
    
    # Step 4: Find minima by inverting the oscillations
    minima_indices, _ = find_peaks(-oscillations, prominence=prominence, distance=distance)
    
    # Step 5: Extract wavelength and signal values at peaks
    maxima_wl = wl[maxima_indices] if len(maxima_indices) > 0 else np.array([])
    minima_wl = wl[minima_indices] if len(minima_indices) > 0 else np.array([])
    maxima_values = oscillations[maxima_indices] if len(maxima_indices) > 0 else np.array([])
    minima_values = oscillations[minima_indices] if len(minima_indices) > 0 else np.array([])
    
    return (background, oscillations, 
            maxima_indices, minima_indices,
            maxima_wl, minima_wl,
            maxima_values, minima_values)


# dictionary of window functions and their corresponding fit functions
# and functions to get the maxima of the fitted functions
# key: window function name
# value: list of window function, fit function, and function to get maxima of the fit function
# parnum: number of parameters for the fit function

'''
Note about the fit
fitparamunits array contains the following:

'''

# ============================================================================
# Stiffness Analysis Functions
# ============================================================================
# These functions analyze peak shape through flank slopes and top curvature
# without assuming specific line shapes (Voigt, Gaussian, etc.)

def find_baseline_region(intensity, noise_threshold=0.1):
    """
    Find baseline regions at start and end of spectrum.
    Returns indices where intensity is below noise_threshold * max(intensity).
    
    Parameters:
    -----------
    intensity : array
        Intensity values
    noise_threshold : float
        Fraction of max intensity to consider as baseline (default 0.1)
    
    Returns:
    --------
    left_base_idx : int
        Index of left baseline point
    right_base_idx : int
        Index of right baseline point
    """
    max_intensity = np.max(intensity)
    threshold = noise_threshold * max_intensity
    
    # Find left baseline (first point where intensity rises above threshold)
    left_base_idx = 0
    for i in range(len(intensity)):
        if intensity[i] > threshold:
            left_base_idx = max(0, i - 1)
            break
    
    # Find right baseline (last point where intensity drops below threshold)
    right_base_idx = len(intensity) - 1
    for i in range(len(intensity) - 1, -1, -1):
        if intensity[i] > threshold:
            right_base_idx = min(len(intensity) - 1, i + 1)
            break
    
    return left_base_idx, right_base_idx

def calculate_flank_slopes(energy, intensity, peak_fraction=0.5, smooth=False):
    """
    Calculate left and right flank slopes using linear regression.
    
    Parameters:
    -----------
    energy : array
        Energy values (x-axis)
    intensity : array
        Intensity values (y-axis)
    peak_fraction : float
        Fraction of peak height to use for flank definition (default 0.5 = 50%)
    smooth : bool
        Apply smoothing before calculation (default False)
    
    Returns:
    --------
    a_L : float
        Left flank slope (counts/eV)
    a_R : float
        Right flank slope (counts/eV)
    E_left_base : float
        Left baseline energy
    E_right_base : float
        Right baseline energy
    """
    # Smooth if requested
    if smooth and len(intensity) > 5:
        intensity = gaussian_filter1d(intensity, sigma=1.0)
    
    # Find peak position
    peak_idx = np.argmax(intensity)
    peak_intensity = intensity[peak_idx]
    
    # Find baseline points
    left_base_idx, right_base_idx = find_baseline_region(intensity, noise_threshold=0.1)
    
    # Define flank regions (from baseline to peak_fraction of peak height)
    flank_threshold = peak_fraction * peak_intensity
    
    # Left flank: from left baseline to where intensity reaches flank_threshold
    left_flank_end = peak_idx
    for i in range(left_base_idx, peak_idx):
        if intensity[i] >= flank_threshold:
            left_flank_end = i
            break
    
    # Right flank: from peak to where intensity drops below flank_threshold
    right_flank_start = peak_idx
    for i in range(peak_idx, right_base_idx):
        if intensity[i] <= flank_threshold:
            right_flank_start = i
            break
    
    # Linear fit for left flank
    if left_flank_end > left_base_idx:
        x_left = energy[left_base_idx:left_flank_end+1]
        y_left = intensity[left_base_idx:left_flank_end+1]
        if len(x_left) > 1:
            poly_coeffs = np.polyfit(x_left - np.mean(x_left), y_left, 1)
            a_L = float(poly_coeffs[0])
        else:
            a_L = 0.0
    else:
        a_L = 0.0
    
    # Linear fit for right flank
    if right_flank_start < right_base_idx:
        x_right = energy[right_flank_start:right_base_idx+1]
        y_right = intensity[right_flank_start:right_base_idx+1]
        if len(x_right) > 1:
            poly_coeffs = np.polyfit(x_right - np.mean(x_right), y_right, 1)
            a_R = float(poly_coeffs[0])
        else:
            a_R = 0.0
    else:
        a_R = 0.0
    
    E_left_base = float(energy[left_base_idx])
    E_right_base = float(energy[right_base_idx])
    
    return a_L, a_R, E_left_base, E_right_base

def calculate_peak_curvature(energy, intensity, window_fraction=0.1):
    """
    Calculate peak-top curvature using quadratic fit.
    
    Parameters:
    -----------
    energy : array
        Energy values (x-axis)
    intensity : array
        Intensity values (y-axis)
    window_fraction : float
        Fraction of energy range to use for peak-top fit (default 0.1 = 10%)
    
    Returns:
    --------
    curvature : float
        Second derivative coefficient c_2 from quadratic fit (counts/eV²)
    E_max : float
        Energy at peak maximum
    I_max : float
        Intensity at peak maximum
    """
    # Find peak position
    peak_idx = np.argmax(intensity)
    E_max = float(energy[peak_idx])
    I_max = float(intensity[peak_idx])
    
    # Define window around peak
    energy_range = float(energy[-1] - energy[0])
    window_size = window_fraction * energy_range
    
    # Find indices within window
    window_mask = np.abs(energy - E_max) <= window_size / 2
    x_window = energy[window_mask]
    y_window = intensity[window_mask]
    
    # Quadratic fit: I(E) = c_0 + c_1*(E - E_max) + c_2*(E - E_max)²
    if len(x_window) >= 3:
        # Shift to center at E_max for better numerical stability
        x_shifted = x_window - E_max
        coeffs = np.polyfit(x_shifted, y_window, 2)
        curvature = float(coeffs[0])  # c_2 coefficient (second derivative / 2)
    else:
        curvature = 0.0
    
    return curvature, E_max, I_max

def calculate_asymmetry(a_L, a_R):
    """
    Calculate asymmetry metric from flank slopes.
    
    Parameters:
    -----------
    a_L : float
        Left flank slope (positive for rising edge)
    a_R : float
        Right flank slope (negative for falling edge)
    
    Returns:
    --------
    S_asym : float
        Asymmetry metric: (|a_R| - a_L) / (|a_R| + a_L)
        Range: -1 to +1
        0 = symmetric
        > 0 = steeper right side
        < 0 = steeper left side
    """
    abs_aR = float(np.abs(a_R))
    abs_aL = float(np.abs(a_L))
    denominator = abs_aR + abs_aL
    
    if denominator > 0:
        S_asym = float((abs_aR - abs_aL) / denominator)
    else:
        S_asym = 0.0
    
    return S_asym

def fitstiffnesstospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Perform stiffness analysis on spectrum.
    
    Parameters:
    -----------
    start : int
        Start index in arrays
    end : int
        End index in arrays
    WL : array
        Energy/wavelength array
    PLB : array
        Intensity array
    maxfev : int
        Not used (for compatibility with other fit functions)
    guess : array or None
        Not used (for compatibility with other fit functions)
    
    Returns:
    --------
    E_max : float
        Energy at peak maximum (eV)
    I_max : float
        Intensity at peak maximum (counts)
    a_L : float
        Left flank slope (counts/eV)
    a_R : float
        Right flank slope (counts/eV)
    S_asym : float
        Asymmetry metric (dimensionless)
    curvature : float
        Peak-top curvature (counts/eV²)
    pcov : None
        Placeholder for covariance matrix (not applicable here)
    """
    # Extract spectrum region - convert to numpy arrays explicitly
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    # Handle edge cases
    if len(x) < 3 or np.max(y) <= 0:
        # Return zeros if spectrum is too short or empty
        return float(0.0), float(0.0), float(0.0), float(0.0), float(0.0), float(0.0), None
    
    # Find peak position first
    peak_idx = np.argmax(y)
    E_max = float(x[peak_idx])
    I_max = float(y[peak_idx])
    
    # Calculate flank slopes
    try:
        a_L, a_R, _, _ = calculate_flank_slopes(x, y, peak_fraction=0.5, smooth=False)
        a_L = float(a_L)
        a_R = float(a_R)
    except:
        a_L = float(0.0)
        a_R = float(0.0)
    
    # Calculate peak curvature
    try:
        curvature, _, _ = calculate_peak_curvature(x, y, window_fraction=0.1)
        curvature = float(curvature)
    except:
        curvature = float(0.0)
    
    # Calculate asymmetry
    try:
        S_asym = calculate_asymmetry(a_L, a_R)
        S_asym = float(S_asym)
    except:
        S_asym = float(0.0)
    
    return E_max, I_max, a_L, a_R, S_asym, curvature, None

def stiffness_window(x, E_max, I_max, a_L, a_R, S_asym, curvature):
    """
    Construct a split Gaussian peak from stiffness parameters for visualization.
    Widths are estimated from flank slopes to provide a visual check of the analysis.
    """
    # Ensure x is a numpy array for element-wise operations
    x = np.array(x)
    
    # Factor relating max slope to width for Gaussian: max_slope = I_max / (sigma * sqrt(e))
    # sigma = I_max / (max_slope * sqrt(e))
    # sqrt(e) approx 1.6487
    factor = 1.6487
    
    # Estimate widths from slopes (handle potential zeros)
    if abs(a_L) > 1e-9:
        sigma_L = I_max / (abs(a_L) * factor)
    else:
        sigma_L = 1.0 # Fallback
        
    if abs(a_R) > 1e-9:
        sigma_R = I_max / (abs(a_R) * factor)
    else:
        sigma_R = 1.0 # Fallback
    
    # Construct split Gaussian
    y = np.zeros_like(x)
    
    # Left side (x < E_max)
    mask_L = x < E_max
    if np.any(mask_L):
        y[mask_L] = I_max * np.exp(-(x[mask_L] - E_max)**2 / (2 * sigma_L**2))
    
    # Right side (x >= E_max)
    mask_R = x >= E_max
    if np.any(mask_R):
        y[mask_R] = I_max * np.exp(-(x[mask_R] - E_max)**2 / (2 * sigma_R**2))
        
    return y

def getstiffnessfwhm(params):
    """
    Estimate FWHM from stiffness parameters (using split Gaussian approximation).
    params: [E_max, I_max, a_L, a_R, S_asym, curvature]
    """
    try:
        I_max = params[1]
        a_L = params[2]
        a_R = params[3]
        
        factor = 1.6487 # sqrt(e)
        
        if abs(a_L) > 1e-9:
            sigma_L = I_max / (abs(a_L) * factor)
        else:
            sigma_L = 0.0
            
        if abs(a_R) > 1e-9:
            sigma_R = I_max / (abs(a_R) * factor)
        else:
            sigma_R = 0.0
            
        # FWHM approx 2.355 * sigma (average sigma)
        fwhm = 2.355 * (sigma_L + sigma_R) / 2
        return fwhm
    except:
        return 0.0

def getmaxstiffness(xmin, xmax, E_max, I_max, a_L, a_R, S_asym, curvature):
    """
    Get maximum position from stiffness analysis parameters.
    
    Parameters:
    -----------
    xmin, xmax : float
        Energy range (not used, for compatibility)
    E_max : float
        Energy at peak maximum
    I_max : float
        Intensity at peak maximum
    a_L, a_R, S_asym, curvature : float
        Other stiffness parameters (not used here)
    
    Returns:
    --------
    E_max : float
        Energy at peak maximum
    I_max : float
        Intensity at peak maximum
    """
    return E_max, I_max

# ============================================================================
# Derivative Analysis Functions (Savitzky-Golay)
# ============================================================================

def fitderivativestospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Perform derivative-based analysis using Savitzky-Golay filter.
    Calculates max slope (1st derivative) and peak curvature (2nd derivative).
    Extended to include inflection point width and asymmetry.
    
    Parameters:
    -----------
    start, end, WL, PLB, maxfev, guess: Standard fit function parameters
    
    Returns:
    --------
    E_max : float
        Energy at peak maximum
    I_max : float
        Intensity at peak maximum
    Max_Slope_Rise : float
        Maximum positive slope (rising edge)
    Max_Slope_Fall : float
        Maximum negative slope (falling edge)
    Curvature_Top : float
        Second derivative at peak maximum (curvature)
    S_asym_deriv : float
        Asymmetry based on max slopes: (|Fall| - Rise) / (|Fall| + Rise)
    Inflection_Width : float
        Distance between max slope positions (approx width)
    Inflection_Asym : float
        Asymmetry based on inflection point positions
    pcov : None
    """
    # Extract spectrum region
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    # Handle edge cases
    if len(x) < 5 or np.max(y) <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
        
    # Parameters for Savitzky-Golay
    # Window length must be odd and <= len(x)
    window_length = min(11, len(x))
    if window_length % 2 == 0:
        window_length -= 1
    if window_length < 5:
        window_length = len(x) if len(x) % 2 != 0 else len(x) - 1
        
    polyorder = 3
    
    try:
        # Calculate derivatives
        # delta x is needed for correct scaling of derivatives
        dx = np.mean(np.diff(x))
        
        # 0th derivative (smoothed intensity)
        y_smooth = savgol_filter(y, window_length, polyorder, deriv=0)
        
        # 1st derivative (slope)
        d1 = savgol_filter(y, window_length, polyorder, deriv=1, delta=dx)
        
        # 2nd derivative (curvature)
        d2 = savgol_filter(y, window_length, polyorder, deriv=2, delta=dx)
        
        # Find peak position from smoothed data
        peak_idx = np.argmax(y_smooth)
        E_max = float(x[peak_idx])
        I_max = float(y_smooth[peak_idx])
        
        # Max slope on rising edge (before peak)
        if peak_idx > 0:
            idx_rise = np.argmax(d1[:peak_idx])
            max_slope_rise = float(d1[idx_rise])
            E_rise = float(x[idx_rise])
        else:
            max_slope_rise = 0.0
            E_rise = E_max
            
        # Max slope on falling edge (after peak) - should be negative
        if peak_idx < len(d1) - 1:
            # search in [peak_idx:]
            idx_fall_local = np.argmin(d1[peak_idx:])
            idx_fall = peak_idx + idx_fall_local
            max_slope_fall = float(d1[idx_fall])
            E_fall = float(x[idx_fall])
        else:
            max_slope_fall = 0.0
            E_fall = E_max
            
        # Curvature at peak (2nd derivative at peak index)
        curvature_top = float(d2[peak_idx])
        
        # Asymmetry metric based on max slopes
        abs_fall = abs(max_slope_fall)
        abs_rise = abs(max_slope_rise)
        denom = abs_fall + abs_rise
        
        if denom > 0:
            s_asym = float((abs_fall - abs_rise) / denom)
        else:
            s_asym = 0.0
            
        # Inflection Width and Asymmetry
        inflection_width = float(E_fall - E_rise)
        
        width_L = E_max - E_rise
        width_R = E_fall - E_max
        
        if inflection_width > 1e-9:
            inflection_asym = float((width_R - width_L) / inflection_width)
        else:
            inflection_asym = 0.0
            
        return E_max, I_max, max_slope_rise, max_slope_fall, curvature_top, s_asym, inflection_width, inflection_asym, None
        
    except Exception as e:
        print(f"Derivative analysis failed: {e}")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None

def derivative_window(x, E_max, I_max, slope_rise, slope_fall, curvature, s_asym, infl_width, infl_asym):
    """
    Construct a visualization for derivative analysis.
    Uses a split Gaussian approximated from the inflection point positions.
    This provides a better reconstruction for asymmetric peaks than slope magnitude alone.
    """
    x = np.array(x)
    
    # If inflection width is available, use it to determine sigmas
    if infl_width > 1e-9:
        # W_L = W/2 * (1 - A)
        # W_R = W/2 * (1 + A)
        # For Gaussian, inflection point is at sigma. So sigma_L = W_L, sigma_R = W_R
        sigma_L = (infl_width / 2.0) * (1.0 - infl_asym)
        sigma_R = (infl_width / 2.0) * (1.0 + infl_asym)
        
        # Safety check
        if sigma_L <= 0: sigma_L = 1.0
        if sigma_R <= 0: sigma_R = 1.0
        
    else:
        # Fallback to slope magnitude method
        factor = 0.6065
        if abs(slope_rise) > 1e-9:
            sigma_L = (I_max * factor) / abs(slope_rise)
        else:
            sigma_L = 1.0
            
        if abs(slope_fall) > 1e-9:
            sigma_R = (I_max * factor) / abs(slope_fall)
        else:
            sigma_R = 1.0
        
    y = np.zeros_like(x)
    
    # Left side
    mask_L = x < E_max
    if np.any(mask_L):
        y[mask_L] = I_max * np.exp(-(x[mask_L] - E_max)**2 / (2 * sigma_L**2))
        
    # Right side
    mask_R = x >= E_max
    if np.any(mask_R):
        y[mask_R] = I_max * np.exp(-(x[mask_R] - E_max)**2 / (2 * sigma_R**2))
        
    return y

def getderivativefwhm(params):
    """
    Estimate FWHM from derivative parameters.
    params: [E_max, I_max, slope_rise, slope_fall, curvature, s_asym, infl_width, infl_asym]
    """
    try:
        # If inflection width is available (index 6), use it
        if len(params) > 6:
            infl_width = params[6]
            if infl_width > 1e-9:
                # For Gaussian, inflection points are at +/- sigma
                # Distance is 2*sigma
                # FWHM = 2.355 * sigma = 1.177 * (2*sigma) = 1.177 * infl_width
                return 1.177 * infl_width
        
        # Fallback to slope method
        I_max = params[1]
        slope_rise = params[2]
        slope_fall = params[3]
        
        factor = 0.6065 # exp(-0.5)
        
        if abs(slope_rise) > 1e-9:
            sigma_L = (I_max * factor) / abs(slope_rise)
        else:
            sigma_L = 0.0
            
        if abs(slope_fall) > 1e-9:
            sigma_R = (I_max * factor) / abs(slope_fall)
        else:
            sigma_R = 0.0
            
        return 2.355 * (sigma_L + sigma_R) / 2
    except:
        return 0.0

def getmaxderivative(xmin, xmax, E_max, I_max, *args):
    return E_max, I_max

# ============================================================================
# Moment Analysis Functions (Variance & Quantiles)
# ============================================================================

def fitmomentstospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Perform moment-based analysis (fit-free).
    Calculates Center of Mass, Variance (Sigma), Skewness, and Quantile Widths.
    Robust for asymmetric peaks.
    
    Parameters:
    -----------
    start, end, WL, PLB, maxfev, guess: Standard fit function parameters
    
    Returns:
    --------
    COM : float
        Center of Mass (First Moment)
    Total_Int : float
        Total Integrated Intensity (Sum)
    Sigma : float
        Standard Deviation (Square root of Second Moment)
    Skewness : float
        Third Standardized Moment
    Quantile_Width : float
        Width containing 68% of data (16% to 84% percentiles)
    Quantile_Asym : float
        Asymmetry based on quantiles: ( (x84-x50) - (x50-x16) ) / (x84-x16)
    pcov : None
    """
    # Extract spectrum region
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    # Handle edge cases
    if len(x) < 3 or np.sum(y) <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
        
    try:
        # Ensure y is positive (baseline subtraction might be needed in real usage, 
        # but here we assume PLB is already baseline corrected or we take abs/clip)
        y_proc = np.maximum(y, 0)
        total_intensity = np.sum(y_proc)
        
        if total_intensity == 0:
             return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
             
        # Normalized probability distribution
        prob = y_proc / total_intensity
        
        # 1. Center of Mass (First Moment)
        com = np.sum(x * prob)
        
        # 2. Variance & Sigma (Second Moment)
        variance = np.sum(((x - com)**2) * prob)
        sigma = np.sqrt(variance)
        
        # 3. Skewness (Third Moment)
        if sigma > 0:
            skewness = np.sum(((x - com)**3) * prob) / (sigma**3)
        else:
            skewness = 0.0
            
        # 4. Quantile Width (16% - 84%)
        cumsum = np.cumsum(prob)
        # Interpolate to find x values for specific probabilities
        # We need strictly increasing cumsum for interp, but it might have flats.
        # Usually fine for noisy data.
        
        x16 = np.interp(0.16, cumsum, x)
        x50 = np.interp(0.50, cumsum, x)
        x84 = np.interp(0.84, cumsum, x)
        
        quantile_width = x84 - x16
        
        # 5. Quantile Asymmetry
        # (Right_Half_Width - Left_Half_Width) / Total_Width
        if quantile_width > 0:
            quantile_asym = ((x84 - x50) - (x50 - x16)) / quantile_width
        else:
            quantile_asym = 0.0
            
        return float(com), float(total_intensity), float(sigma), float(skewness), float(quantile_width), float(quantile_asym), None
        
    except Exception as e:
        print(f"Moment analysis failed: {e}")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None

def moment_window(x, com, total_int, sigma, skewness, q_width, q_asym):
    """
    Construct a visualization for moment analysis.
    Uses a Gaussian defined by COM and Sigma (Second Moment).
    """
    x = np.array(x)
    
    if sigma <= 0:
        return np.zeros_like(x)
        
    # Gaussian area = A = total_int * dx (approx)
    # But total_int here is sum(y).
    # If x is uniformly spaced with spacing dx: sum(y) * dx = Integral.
    # Gaussian integral = Amplitude * sigma * sqrt(2*pi).
    # So Amplitude = (sum(y) * dx) / (sigma * sqrt(2*pi)).
    
    # Estimate dx
    if len(x) > 1:
        dx = np.mean(np.diff(x))
    else:
        dx = 1.0
        
    amplitude = (total_int * dx) / (sigma * np.sqrt(2 * np.pi))
    
    # Simple Gaussian reconstruction
    y = amplitude * np.exp(-(x - com)**2 / (2 * sigma**2))
    
    return y

def getmomentfwhm(params):
    """
    Estimate FWHM from moment parameters.
    params: [COM, Total_Int, Sigma, Skewness, Quantile_Width, Quantile_Asym]
    """
    try:
        # Use Quantile Width (contains 68% of data)
        # For Gaussian, 16%-84% is exactly 2*sigma.
        # FWHM = 2.355 * sigma.
        # So FWHM approx 1.177 * Quantile_Width
        q_width = params[4]
        return 1.177 * q_width
    except:
        return 0.0

def getmaxmoments(xmin, xmax, com, total_int, sigma, *args):
    """
    Return max position and intensity of the reconstructed Gaussian.
    """
    # Reconstruct max intensity
    # Amplitude = (total_int * dx) / (sigma * sqrt(2*pi))
    # We don't have dx here easily without x array.
    # But we can approximate or just return total_int for now?
    # fitkeys expects (x_max, y_max).
    # Let's try to estimate y_max assuming some average dx or just return 0 for y_max if not critical.
    # Actually, for plotting "fitmaxX" and "fitmaxY", we need reasonable values.
    # Let's assume dx is small or just return total_int as a proxy for "magnitude".
    # OR, better: The user might want to map COM as the "position".
    
    # For intensity, let's return total_int (Integrated Intensity) as it's a robust measure of signal strength.
    return com, total_int

# ============================================================================
# Center of Mass Analysis
# ============================================================================

def fitcomtospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Calculate Center of Mass (First Moment) of the spectrum.
    Fit-free method. Respects asymmetric peaks by integrating over the full shape.
    Formula: COM = Sum(lambda_i * I_i) / Sum(I_i)
    """
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    if len(x) < 2 or np.sum(y) <= 0:
        return 0.0, 0.0, None
        
    try:
        # Ensure y is positive
        y_proc = np.maximum(y, 0)
        total_intensity = np.sum(y_proc)
        
        if total_intensity == 0:
            return 0.0, 0.0, None
            
        com = np.sum(x * y_proc) / total_intensity
        
        return float(com), float(total_intensity), None
    except Exception as e:
        return 0.0, 0.0, None

def com_window(x, com, total_int):
    """
    Window function for Center of Mass.
    Returns zeros as there is no shape fit to display.
    """
    return np.zeros_like(np.array(x))

def getmaxcom(xmin, xmax, com, total_int, *args):
    return com, total_int

def getcomfwhm(params):
    return 0.0

# ============================================================================
# Max Decay Slope Analysis
# ============================================================================

def fitdecaytospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Calculate Maximum Decay Slope on the right flank.
    Finds peak maximum, then computes discrete derivatives on the right side.
    Returns the steepest fall (minimum slope value).
    Robustness: Uses 5th percentile of slopes if enough points are available.
    """
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    if len(x) < 3 or np.max(y) <= 0:
        return 0.0, 0.0, 0.0, 0.0, None
        
    try:
        # Find peak index
        k_max = np.argmax(y)
        E_max = float(x[k_max])
        I_max = float(y[k_max])
        
        # Check if there are points to the right
        if k_max >= len(x) - 2:
            return E_max, I_max, 0.0, E_max, None
            
        # Extract right side
        x_right = x[k_max:]
        y_right = y[k_max:]
        
        # Calculate discrete derivatives
        # s_i = (I_{i+1} - I_i) / (lambda_{i+1} - lambda_i)
        dy = np.diff(y_right)
        dx = np.diff(x_right)
        
        # Avoid division by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            slopes = dy / dx
            
        # Filter out NaNs or Infs
        valid_mask = np.isfinite(slopes)
        slopes = slopes[valid_mask]
        
        if len(slopes) == 0:
            return E_max, I_max, 0.0, E_max, None
            
        # We are looking for the steepest fall, which is the minimum (most negative) slope.
        # Robustness: Ignore last few points if array is long enough
        if len(slopes) > 5:
            # Ignore last 2 points to avoid tail noise
            slopes_robust = slopes[:-2]
        else:
            slopes_robust = slopes
            
        if len(slopes_robust) == 0:
             slopes_robust = slopes
             
        # Robustness: Use 5th percentile instead of absolute min if enough points
        if len(slopes_robust) >= 10:
            s_min = np.percentile(slopes_robust, 5)
        else:
            s_min = np.min(slopes_robust)
            
        # Find the energy where this slope occurs (approximate)
        # We need to map back to original x array
        # This is a bit tricky with percentile, so let's just find the index of the value closest to s_min
        idx_min = np.argmin(np.abs(slopes - s_min))
        # x_right has length N, slopes has length N-1. 
        # Slope i corresponds to interval between x_right[i] and x_right[i+1].
        # Let's assign it to x_right[i] (or midpoint)
        slope_energy = float(x_right[idx_min])
        
        # Return absolute value as requested? 
        # "Use |s_min| as a flank-steepness measure."
        # But usually fit functions return the actual parameter. 
        # Let's return the signed slope, the user can take abs in their head or we label it "Max Decay Slope".
        # Actually, for "Max Decay Slope", a large negative number is the "Max Decay".
        # Let's return the signed value s_min.
        
        # --- Left Flank (Rise) ---
        # Check if there are points to the left
        if k_max < 2:
            s_max = 0.0
            rise_slope_energy = E_max
        else:
            x_left = x[:k_max+1] # Include peak for continuity
            y_left = y[:k_max+1]
            
            dy_left = np.diff(y_left)
            dx_left = np.diff(x_left)
            
            with np.errstate(divide='ignore', invalid='ignore'):
                slopes_left = dy_left / dx_left
                
            valid_mask_left = np.isfinite(slopes_left)
            slopes_left = slopes_left[valid_mask_left]
            
            if len(slopes_left) == 0:
                s_max = 0.0
                rise_slope_energy = E_max
            else:
                # Robustness: Ignore first few points?
                if len(slopes_left) > 5:
                    slopes_left_robust = slopes_left[2:] # Ignore first 2 points
                else:
                    slopes_left_robust = slopes_left
                
                if len(slopes_left_robust) == 0:
                    slopes_left_robust = slopes_left

                # We want max positive slope.
                # Robustness: 95th percentile
                if len(slopes_left_robust) >= 10:
                    s_max = np.percentile(slopes_left_robust, 95)
                else:
                    s_max = np.max(slopes_left_robust)
                    
                # Find energy
                idx_max = np.argmin(np.abs(slopes_left - s_max))
                rise_slope_energy = float(x_left[idx_max])

        return E_max, I_max, float(s_min), slope_energy, float(s_max), rise_slope_energy, None
        
    except Exception as e:
        print(f"Decay analysis failed: {e}")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None

def decay_window(x, E_max, I_max, s_min, slope_energy, s_max, rise_slope_energy):
    """
    Visualization for Decay Slope.
    Draws a tangent line at the point of steepest descent.
    """
    x = np.array(x)
    y = np.zeros_like(x)
    
    # We want to draw a line segment representing the slope.
    # Line equation: y - y0 = m(x - x0) => y = m(x - x0) + y0
    # We know m = s_min and x0 = slope_energy.
    # We don't strictly know y0 (intensity at slope_energy) passed in params.
    # But we can estimate it or just draw a line that looks "slope-like".
    # Better: The window function is usually used to show the "fit". 
    # Since we don't have a full fit, maybe we just return zeros?
    # Or we can try to reconstruct a line segment of length, say, 10% of the range.
    
    # Let's just return zeros to keep it clean, as "fit-free" implies no model to show.
    # Or, if we want to be fancy, we could try to show the slope.
    # But without y0, it's hard to place it vertically correct.
    
    return y

def getmaxdecay(xmin, xmax, E_max, I_max, *args):
    return E_max, I_max

def getdecayfwhm(params):
    return 0.0

# ============================================================================
# Binning Decay Analysis
# ============================================================================

def fitbinningtospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Bin spectrum to 11 points (10 intervals) and calculate intensity changes.
    Returns Start Intensity and 10 differences (slopes).
    User requested 10 parameters, so we use 11 points to get 10 intervals.
    """
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    if len(x) < 2:
        return tuple([0.0]*11) + (None,)
        
    try:
        # Resample to 11 points
        x_bins = np.linspace(x[0], x[-1], 11)
        y_bins = np.interp(x_bins, x, y)
        
        # Calculate differences (Counts)
        diffs = np.diff(y_bins)
        
        # Return Start Intensity + 10 diffs
        params = [float(y_bins[0])] + [float(d) for d in diffs]
        return tuple(params) + (None,)
    except Exception as e:
        print(f"Binning analysis failed: {e}")
        return tuple([0.0]*11) + (None,)

def binning_window(x, start_int, *diffs):
    """
    Reconstruct binned spectrum.
    """
    x = np.array(x)
    if len(x) == 0: return x
    
    # Reconstruct y_bins points
    y_points = [start_int]
    current_val = start_int
    for d in diffs:
        current_val += d
        y_points.append(current_val)
        
    y_points = np.array(y_points)
    
    # x coordinates for these points
    x_bins = np.linspace(x[0], x[-1], 11)
    
    # Interpolate to original x
    y_recon = np.interp(x, x_bins, y_points)
    
    return y_recon

def getmaxbinning(xmin, xmax, start_int, *diffs):
    # Reconstruct to find max
    y_points = [start_int]
    current_val = start_int
    for d in diffs:
        current_val += d
        y_points.append(current_val)
    
    max_idx = np.argmax(y_points)
    max_val = y_points[max_idx]
    
    # x position
    x_bins = np.linspace(xmin, xmax, 11)
    max_x = x_bins[max_idx]
    
    return max_x, max_val

def getbinningfwhm(params):
    return 0.0

# ============================================================================
# First Derivative Analysis (Finite Differences)
# ============================================================================

def fitderivative1tospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Calculate First Derivative using finite differences.
    Extracts edge sharpness and peak position from derivative features.
    
    Parameters:
    -----------
    start, end, WL, PLB, maxfev, guess: Standard fit function parameters
    
    Returns:
    --------
    max_pos_deriv : float
        Maximum positive derivative (steepest rise)
    energy_max_pos : float
        Energy where max positive derivative occurs
    max_neg_deriv : float
        Maximum negative derivative (steepest fall)
    energy_max_neg : float
        Energy where max negative derivative occurs
    deriv_range : float
        Total derivative span (max - min)
    zero_cross_energy : float
        Energy where derivative crosses zero (peak position estimate)
    pcov : None
    """
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    # Handle edge cases
    if len(x) < 3 or np.max(y) <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
    
    try:
        # Calculate first derivative using finite differences
        dy = np.diff(y)
        dx = np.diff(x)
        
        # Avoid division by zero
        with np.errstate(divide='ignore', invalid='ignore'):
            deriv = dy / dx
        
        # Filter out NaNs and Infs
        valid_mask = np.isfinite(deriv)
        deriv_valid = deriv[valid_mask]
        
        if len(deriv_valid) < 2:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
        
        # x coordinates for derivatives (midpoints)
        x_deriv = (x[:-1] + x[1:]) / 2
        x_deriv_valid = x_deriv[valid_mask]
        
        # Trim edges to avoid boundary effects
        if len(deriv_valid) > 4:
            trim_start = 1
            trim_end = -1
            deriv_trimmed = deriv_valid[trim_start:trim_end]
            x_deriv_trimmed = x_deriv_valid[trim_start:trim_end]
        else:
            deriv_trimmed = deriv_valid
            x_deriv_trimmed = x_deriv_valid
        
        if len(deriv_trimmed) == 0:
            deriv_trimmed = deriv_valid
            x_deriv_trimmed = x_deriv_valid
        
        # Find max positive derivative (95th percentile for robustness)
        pos_derivs = deriv_trimmed[deriv_trimmed > 0]
        if len(pos_derivs) >= 10:
            max_pos_deriv = float(np.percentile(pos_derivs, 95))
        elif len(pos_derivs) > 0:
            max_pos_deriv = float(np.max(pos_derivs))
        else:
            max_pos_deriv = 0.0
        
        # Find energy at max positive derivative
        if max_pos_deriv > 0:
            idx_max_pos = np.argmin(np.abs(deriv_trimmed - max_pos_deriv))
            energy_max_pos = float(x_deriv_trimmed[idx_max_pos])
        else:
            energy_max_pos = float(x[0])
        
        # Find max negative derivative (5th percentile for robustness)
        neg_derivs = deriv_trimmed[deriv_trimmed < 0]
        if len(neg_derivs) >= 10:
            max_neg_deriv = float(np.percentile(neg_derivs, 5))
        elif len(neg_derivs) > 0:
            max_neg_deriv = float(np.min(neg_derivs))
        else:
            max_neg_deriv = 0.0
        
        # Find energy at max negative derivative
        if max_neg_deriv < 0:
            idx_max_neg = np.argmin(np.abs(deriv_trimmed - max_neg_deriv))
            energy_max_neg = float(x_deriv_trimmed[idx_max_neg])
        else:
            energy_max_neg = float(x[-1])
        
        # Calculate derivative range
        deriv_range = float(max_pos_deriv - max_neg_deriv)
        
        # Find zero-crossing (peak position)
        # Look for sign changes in derivative
        sign_changes = np.where(np.diff(np.sign(deriv_valid)))[0]
        
        if len(sign_changes) > 0:
            # Find peak in original spectrum
            peak_idx = np.argmax(y)
            peak_energy = float(x[peak_idx])
            
            # Find zero-crossing closest to peak
            zero_crossings = []
            for idx in sign_changes:
                # Interpolate to find exact zero-crossing
                if idx < len(deriv_valid) - 1:
                    x1, x2 = x_deriv_valid[idx], x_deriv_valid[idx + 1]
                    d1, d2 = deriv_valid[idx], deriv_valid[idx + 1]
                    
                    if d2 != d1:
                        # Linear interpolation
                        zero_x = x1 - d1 * (x2 - x1) / (d2 - d1)
                        zero_crossings.append(zero_x)
            
            if len(zero_crossings) > 0:
                # Choose zero-crossing closest to peak
                zero_crossings = np.array(zero_crossings)
                closest_idx = np.argmin(np.abs(zero_crossings - peak_energy))
                zero_cross_energy = float(zero_crossings[closest_idx])
            else:
                zero_cross_energy = peak_energy
        else:
            # No zero-crossing found, use peak position
            zero_cross_energy = float(x[np.argmax(y)])
        
        return max_pos_deriv, energy_max_pos, max_neg_deriv, energy_max_neg, deriv_range, zero_cross_energy, None
        
    except Exception as e:
        print(f"First derivative analysis failed: {e}")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None

def derivative1_window(x, max_pos_deriv, energy_max_pos, max_neg_deriv, energy_max_neg, deriv_range, zero_cross_energy):
    """
    Window function for first derivative visualization.
    Returns zeros (fit-free method, no reconstruction).
    """
    return np.zeros_like(np.array(x))

def getmaxderivative1(xmin, xmax, max_pos_deriv, energy_max_pos, max_neg_deriv, energy_max_neg, deriv_range, zero_cross_energy):
    """
    Return peak position and max derivative magnitude.
    """
    return zero_cross_energy, max_pos_deriv

def getderivative1fwhm(params):
    """
    No FWHM estimate from first derivative alone.
    """
    return 0.0

# ============================================================================
# Second Derivative Analysis (Finite Differences)
# ============================================================================

def fitderivative2tospec(start, end, WL, PLB, maxfev=10000, guess=None):
    """
    Calculate Second Derivative (curvature) using finite differences.
    Extracts peak sharpness, inflection points, and effective peak width.
    
    Parameters:
    -----------
    start, end, WL, PLB, maxfev, guess: Standard fit function parameters
    
    Returns:
    --------
    max_neg_curv : float
        Maximum negative curvature (peak sharpness)
    energy_max_neg_curv : float
        Energy where max negative curvature occurs (peak center)
    max_pos_curv : float
        Maximum positive curvature (valley/shoulder sharpness)
    curv_range : float
        Total curvature span
    left_infl_energy : float
        Left inflection point energy
    right_infl_energy : float
        Right inflection point energy
    infl_width : float
        Distance between inflection points
    pcov : None
    """
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    # Handle edge cases
    if len(x) < 4 or np.max(y) <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
    
    try:
        # Calculate first derivative
        dy = np.diff(y)
        dx = np.diff(x)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            deriv1 = dy / dx
        
        # Filter first derivative
        valid_mask1 = np.isfinite(deriv1)
        deriv1_valid = deriv1[valid_mask1]
        x_deriv1 = (x[:-1] + x[1:]) / 2
        x_deriv1_valid = x_deriv1[valid_mask1]
        
        if len(deriv1_valid) < 3:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
        
        # Calculate second derivative from first derivative
        dy2 = np.diff(deriv1_valid)
        dx2 = np.diff(x_deriv1_valid)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            deriv2 = dy2 / dx2
        
        # Filter second derivative
        valid_mask2 = np.isfinite(deriv2)
        deriv2_valid = deriv2[valid_mask2]
        x_deriv2 = (x_deriv1_valid[:-1] + x_deriv1_valid[1:]) / 2
        x_deriv2_valid = x_deriv2[valid_mask2]
        
        if len(deriv2_valid) < 2:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None
        
        # Trim edges
        if len(deriv2_valid) > 4:
            trim_start = 1
            trim_end = -1
            deriv2_trimmed = deriv2_valid[trim_start:trim_end]
            x_deriv2_trimmed = x_deriv2_valid[trim_start:trim_end]
        else:
            deriv2_trimmed = deriv2_valid
            x_deriv2_trimmed = x_deriv2_valid
        
        if len(deriv2_trimmed) == 0:
            deriv2_trimmed = deriv2_valid
            x_deriv2_trimmed = x_deriv2_valid
        
        # Find max negative curvature (5th percentile for robustness)
        neg_curvs = deriv2_trimmed[deriv2_trimmed < 0]
        if len(neg_curvs) >= 10:
            max_neg_curv = float(np.percentile(neg_curvs, 5))
        elif len(neg_curvs) > 0:
            max_neg_curv = float(np.min(neg_curvs))
        else:
            max_neg_curv = 0.0
        
        # Find energy at max negative curvature
        if max_neg_curv < 0:
            idx_max_neg_curv = np.argmin(np.abs(deriv2_trimmed - max_neg_curv))
            energy_max_neg_curv = float(x_deriv2_trimmed[idx_max_neg_curv])
        else:
            # Fallback to peak position
            energy_max_neg_curv = float(x[np.argmax(y)])
        
        # Find max positive curvature (95th percentile)
        pos_curvs = deriv2_trimmed[deriv2_trimmed > 0]
        if len(pos_curvs) >= 10:
            max_pos_curv = float(np.percentile(pos_curvs, 95))
        elif len(pos_curvs) > 0:
            max_pos_curv = float(np.max(pos_curvs))
        else:
            max_pos_curv = 0.0
        
        # Calculate curvature range
        curv_range = float(max_pos_curv - max_neg_curv)
        
        # Find inflection points (zero-crossings of second derivative)
        sign_changes = np.where(np.diff(np.sign(deriv2_valid)))[0]
        
        inflection_points = []
        for idx in sign_changes:
            if idx < len(deriv2_valid) - 1:
                x1, x2 = x_deriv2_valid[idx], x_deriv2_valid[idx + 1]
                d1, d2 = deriv2_valid[idx], deriv2_valid[idx + 1]
                
                if d2 != d1:
                    # Linear interpolation
                    infl_x = x1 - d1 * (x2 - x1) / (d2 - d1)
                    inflection_points.append(infl_x)
        
        # Identify left and right inflections relative to peak
        peak_energy = float(x[np.argmax(y)])
        
        if len(inflection_points) >= 2:
            inflection_points = np.array(inflection_points)
            inflection_points_sorted = np.sort(inflection_points)
            
            # Find inflections closest to peak on left and right
            left_infls = inflection_points_sorted[inflection_points_sorted < peak_energy]
            right_infls = inflection_points_sorted[inflection_points_sorted > peak_energy]
            
            if len(left_infls) > 0:
                left_infl_energy = float(left_infls[-1])  # Closest to peak on left
            else:
                # Fallback: 10% of range to the left
                left_infl_energy = float(peak_energy - 0.1 * (x[-1] - x[0]))
            
            if len(right_infls) > 0:
                right_infl_energy = float(right_infls[0])  # Closest to peak on right
            else:
                # Fallback: 10% of range to the right
                right_infl_energy = float(peak_energy + 0.1 * (x[-1] - x[0]))
        
        elif len(inflection_points) == 1:
            # Only one inflection found
            infl = inflection_points[0]
            if infl < peak_energy:
                left_infl_energy = float(infl)
                right_infl_energy = float(peak_energy + 0.1 * (x[-1] - x[0]))
            else:
                left_infl_energy = float(peak_energy - 0.1 * (x[-1] - x[0]))
                right_infl_energy = float(infl)
        
        else:
            # No inflections found, use defaults
            left_infl_energy = float(peak_energy - 0.1 * (x[-1] - x[0]))
            right_infl_energy = float(peak_energy + 0.1 * (x[-1] - x[0]))
        
        # Calculate inflection width
        infl_width = float(right_infl_energy - left_infl_energy)
        
        return max_neg_curv, energy_max_neg_curv, max_pos_curv, curv_range, left_infl_energy, right_infl_energy, infl_width, None
        
    except Exception as e:
        print(f"Second derivative analysis failed: {e}")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None

def derivative2_window(x, max_neg_curv, energy_max_neg_curv, max_pos_curv, curv_range, left_infl_energy, right_infl_energy, infl_width):
    """
    Window function for second derivative visualization.
    Returns zeros (fit-free method, no reconstruction).
    """
    return np.zeros_like(np.array(x))

def getmaxderivative2(xmin, xmax, max_neg_curv, energy_max_neg_curv, max_pos_curv, curv_range, left_infl_energy, right_infl_energy, infl_width):
    """
    Return peak position (from curvature) and magnitude.
    """
    # Return energy at max negative curvature as peak position
    # Magnitude is just the curvature value (negative)
    return energy_max_neg_curv, abs(max_neg_curv)

def getderivative2fwhm(params):
    """
    Return inflection width as FWHM estimate.
    params: [max_neg_curv, energy_max_neg_curv, max_pos_curv, curv_range, 
             left_infl_energy, right_infl_energy, infl_width]
    """
    try:
        return params[6]  # inflection width
    except:
        return 0.0

# ============================================================================
# Derivative Points Analysis (Zero Points of D1 and Maxima of D2)
# ============================================================================

def fitderivativepoints(start, end, WL, PLB, maxfev=10000, guess=None, deriv1=None, deriv2=None):
    """
    Find up to 5 zero points of the 1st derivative and up to 10 maxima of the 2nd derivative.
    Uses pre-calculated derivatives passed as arguments.
    Returns 30 parameters + None.
    """
    x = np.array(WL[start:end])
    y = np.array(PLB[start:end])
    
    # Initialize result array: 30 zeros (15 pairs of x,y)
    results = [0.0] * 30
    
    if len(x) < 4:
        return *results, None
        
    try:
        # --- Handle First Derivative (Zero Crossings) ---
        # We need deriv1 to correspond to the ROI or be sliceable
        if deriv1 is not None and len(deriv1) > 0:
            if len(deriv1) == len(WL): # Full spectrum derivative
                 d1_roi = np.array(deriv1[start:end])
                 x_d1 = x
            elif len(deriv1) == len(x): # Already sliced
                 d1_roi = np.array(deriv1)
                 x_d1 = x
            else:
                 d1_roi = []
        else:
            d1_roi = [] # Do not calculate locally
            
        if len(d1_roi) > 1:
            # Find zero crossings (sign changes)
            sign_changes = np.where(np.diff(np.sign(d1_roi)))[0]
            
            d1_points = []
            for idx in sign_changes:
                # Interpolate x where D1 is zero
                y1 = d1_roi[idx]
                y2 = d1_roi[idx+1]
                if y2 != y1:
                     x_cross = x_d1[idx] - y1 * (x_d1[idx+1] - x_d1[idx]) / (y2 - y1)
                else:
                     x_cross = x_d1[idx]
                
                # Intensity at x_cross from original spectrum
                int_cross = np.interp(x_cross, x, y)
                d1_points.append((x_cross, int_cross))
            
            # Sort by Intensity descending
            d1_points.sort(key=lambda p: p[1], reverse=True)
            
            for i in range(min(5, len(d1_points))):
                results[2*i] = float(d1_points[i][0])
                results[2*i+1] = float(d1_points[i][1])

        # --- Handle Second Derivative (Maxima) ---
        if deriv2 is not None and len(deriv2) > 0:
            if len(deriv2) == len(WL):
                 d2_roi = np.array(deriv2[start:end])
                 x_d2 = x
            elif len(deriv2) == len(x):
                 d2_roi = np.array(deriv2)
                 x_d2 = x
            else:
                 d2_roi = []
        else:
             d2_roi = []

        if len(d2_roi) > 2:
            # Find local maxima 
            d2_points = []
            for i in range(1, len(d2_roi)-1):
                if d2_roi[i] > d2_roi[i-1] and d2_roi[i] > d2_roi[i+1]:
                    curr_x = x_d2[i]
                    curr_d2_val = d2_roi[i]
                    # Corresponding Intensity in original spectrum
                    curr_int = np.interp(curr_x, x, y)
                    d2_points.append((curr_x, curr_int, curr_d2_val))
            
            # Sort by curvature magnitude (d2_val) descending
            d2_points.sort(key=lambda p: p[2], reverse=True)
            
            for i in range(min(10, len(d2_points))):
                offset = 10 + 2*i
                results[offset] = float(d2_points[i][0])
                results[offset+1] = float(d2_points[i][1])

        return *results, None

    except Exception as e:
        print(f"Derivative points analysis failed: {e}")
        return *([0.0]*30), None

def derivativepoints_window(x, *params):
    return np.zeros_like(np.array(x))

def getmaxderivativepoints(xmin, xmax, *params):
    # Return first point WL and Int
    return params[0], params[1]

def getderivativepointsfwhm(params):
    return 0.0

# ============================================================================
# End of Derivative Analysis Functions
# ============================================================================


# todo: add derivative points methods to fitkeys that use the calculated derivative if it has been calculated. 
# if not not calculated, init with zeros.

# explained fitkeys structure:
#            'name':[
#               window_function, 
#               fit_function, 
#               get_max_function, 
#               'Label for GUI', 
#               number_of_output_parameters, 
#               ['list of parameter names'], 
#               ['list of parameter units'], 
#               fwhm_function, 
#               fit_state
#               ]     

fitkeys = {'lorentz':[lorentzwind, fitlorentztospec, getmaxlorentz, 'Lorentz fit', 3, ['Lorentzian amplitude', 'Lorentzian center', 'Lorentzian width'], ['Counts', 'nm', 'nm'], getlorentzfwhm, 0],
           'gaussian':[gaussianwind, fitgaussiantospec, getmaxgaussian, 'Gaussian fit', 3, ['Gaussian amplitude', 'Gaussian center', 'Gaussian width'], ['Counts', 'nm', 'nm'], getgaussianfwhm, 0],
           'voigt':[voigtwind, fitvoigttospec, getmaxvoigt, 'Voigt fit', 4, ['Voigt amplitude', 'Voigt center', 'Voigt width', 'Voigt gamma'], ['Counts', 'nm', 'nm', 'nm'], getvoigtfwhm, 0], 
           'linear':[linearwind, fitlinetospec, getmaxlinear, 'Linear fit', 2, ['Linear slope', 'Linear offset'], ['nm', 'Counts'], None, 0],
           'double lorentz':[double_lorentzwind, fitdoublelorentztospec, getmaxdoublelorentz, 'Double Lorentz fit', 6, ['Double Lorentzian amplitude 1', 'Double Lorentzian center 1', 'Double Lorentzian width 1', 'Double Lorentzian amplitude 2', 'Double Lorentzian center 2', 'Double Lorentzian width 2'], ['Counts', 'nm', 'nm', 'Counts', 'nm', 'nm'], getdoublelorentzfwhm, 0], 
           'double gaussian':[double_gaussianwind, fitdoublegaussiantospec, getmaxdoublegaussian, 'Double Gaussian fit', 6, ['Double Gaussian amplitude 1', 'Double Gaussian center 1', 'Double Gaussian width 1', 'Double Gaussian amplitude 2', 'Double Gaussian center 2', 'Double Gaussian width 2'], ['Counts', 'nm', 'nm', 'Counts', 'nm', 'nm'], getdoublegaussianfwhm, 0], 
           'double voigt':[
               double_voigtwind, fitdoublevoigttospec, getmaxdoublevoigt, 'Double Voigt fit', 8, ['Double Voigt amplitude 1', 'Double Voigt center 1', 'Double Voigt width 1', 'Double Voigt gamma 1', 'Double Voigt amplitude 2', 'Double Voigt center 2', 'Double Voigt width 2', 'Double Voigt gamma 2'], ['Counts', 'nm', 'nm', 'nm', 'Counts', 'nm', 'nm', 'nm'], getdoublevoigtfwhm, 0
               ],
           'oscillations': [
               None,  # No window function for oscillations (uses isolate_oscillation internally)
               fitoscillationtospec,  # Fitting function
               getoscillationdata,  # Placeholder for compatibility
               'Oscillation Extraction with Phase', 
               12,  # Number of scalar output parameters (including phase evolution)
               ['Background mean', 'Oscillation amplitude', 'Number of maxima', 'Number of minima', 
                'Mean maxima value', 'Mean minima value', 'Oscillation frequency estimate', 'Oscillation range',
                'Phase chirp rate', 'Initial frequency', 'Mean period', 'Period std dev'],
               ['Counts', 'Counts', '', '', 'Counts', 'Counts', 'eV^-1', 'Counts',
                'eV^-2', 'eV^-1', 'eV', 'eV'],
               None, # No FWHM function
               0 # fit state: 0= not fitted, 1=fitted, 2=failed
           ],
           'stiffness': [
               stiffness_window,  # Window function (returns zeros)
               fitstiffnesstospec,  # Analysis function
               getmaxstiffness,  # Get peak position
               'Stiffness Analysis',  # Label for GUI
               6,  # Number of output parameters
               ['Peak Energy', 'Peak Intensity', 'Left Stiffness', 'Right Stiffness', 'Asymmetry', 'Curvature'],
               ['eV', 'Counts', 'Counts/eV', 'Counts/eV', '', 'Counts/eV²'],
               getstiffnessfwhm,  # FWHM function (returns 0.0)
               0  # fit state: 0= not fitted, 1=fitted, 2=failed
           ],
           'derivative': [
               derivative_window,
               fitderivativestospec,
               getmaxderivative,
               'Derivative Analysis (SG)',
               8,
               ['Peak Energy', 'Peak Intensity', 'Max Rise Slope', 'Max Fall Slope', 'Peak Curvature', 'Slope Asymmetry', 'Inflection Width', 'Inflection Asymmetry'],
               ['eV', 'Counts', 'Counts/eV', 'Counts/eV', 'Counts/eV²', '', 'eV', ''],
               getderivativefwhm,
               0
           ],
           'moments': [
               moment_window,
               fitmomentstospec,
               getmaxmoments,
               'Moment Analysis (Fit-Free)',
               6,
               ['Center of Mass', 'Total Intensity', 'Sigma (Variance)', 'Skewness', 'Quantile Width (16-84%)', 'Quantile Asymmetry'],
               ['eV', 'Counts', 'eV', '', 'eV', ''],
               getmomentfwhm,
               0
           ],
           'com': [
               com_window,
               fitcomtospec,
               getmaxcom,
               'Center of Mass',
               2,
               ['Center of Mass', 'Integrated Intensity'],
               ['eV', 'Counts'],
               getcomfwhm,
               0
           ],
           'decay': [
               decay_window,
               fitdecaytospec,
               getmaxdecay,
               'Max Decay/Rise Slope',
               6,
               ['Peak Energy', 'Peak Intensity', 'Max Decay Slope', 'Decay Slope Energy', 'Max Rise Slope', 'Rise Slope Energy'],
               ['eV', 'Counts', 'Counts/eV', 'eV', 'Counts/eV', 'eV'],
               getdecayfwhm,
               0
           ],
           'binning': [
               binning_window,
               fitbinningtospec,
               getmaxbinning,
               'Binning Decay (10 bins)',
               11,
               ['Start Intensity'] + [f'Bin Diff {i+1}' for i in range(10)],
               ['Counts'] + ['Counts']*10,
               getbinningfwhm,
               0
           ],
           'derivative1': [
               derivative1_window,
               fitderivative1tospec,
               getmaxderivative1,
               'First Derivative (FD)',
               6,
               ['Max Positive Derivative', 'Energy at Max Pos. Deriv.', 'Max Negative Derivative', 
                'Energy at Max Neg. Deriv.', 'Derivative Range', 'Zero-Crossing Energy'],
               ['Counts/eV', 'eV', 'Counts/eV', 'eV', 'Counts/eV', 'eV'],
               getderivative1fwhm,
               0
           ],
           'derivative2': [
               derivative2_window,
               fitderivative2tospec,
               getmaxderivative2,
               'Second Derivative (FD)',
               7,
               ['Max Negative Curvature', 'Energy at Max Neg. Curv.', 'Max Positive Curvature', 
                'Curvature Range', 'Left Inflection Energy', 'Right Inflection Energy', 'Inflection Width'],
               ['Counts/eV²', 'eV', 'Counts/eV²', 'Counts/eV²', 'eV', 'eV', 'eV'],
               getderivative2fwhm,
               0
           ],
           'derivative_points': [
               derivativepoints_window,
               fitderivativepoints,
               getmaxderivativepoints,
               'Derivative Points',
               30,
               ['D1 Zero 1 WL', 'D1 Zero 1 Int', 'D1 Zero 2 WL', 'D1 Zero 2 Int', 'D1 Zero 3 WL', 'D1 Zero 3 Int', 'D1 Zero 4 WL', 'D1 Zero 4 Int', 'D1 Zero 5 WL', 'D1 Zero 5 Int',
                'D2 Max 1 WL', 'D2 Max 1 Int', 'D2 Max 2 WL', 'D2 Max 2 Int', 'D2 Max 3 WL', 'D2 Max 3 Int', 'D2 Max 4 WL', 'D2 Max 4 Int', 'D2 Max 5 WL', 'D2 Max 5 Int',
                'D2 Max 6 WL', 'D2 Max 6 Int', 'D2 Max 7 WL', 'D2 Max 7 Int', 'D2 Max 8 WL', 'D2 Max 8 Int', 'D2 Max 9 WL', 'D2 Max 9 Int', 'D2 Max 10 WL', 'D2 Max 10 Int'],
               ['nm', 'Counts'] * 15,
               getderivativepointsfwhm,
               0
           ], 
# explained fitkeys structure:
#            'name':[
#               window_function, 
#               fit_function, 
#               get_max_function, 
#               'Label for GUI', 
#               number_of_output_parameters, 
#               ['list of parameter names'], 
#               ['list of parameter units'], 
#               fwhm_function, 
#               fit_state
#               ]  
# example:'double voigt':[
#               double_voigtwind, 
#               fitdoublevoigttospec, 
#               getmaxdoublevoigt, 
#               'Double Voigt fit', 
#               8, 
#               ['Double Voigt amplitude 1', 'Double Voigt center 1', 'Double Voigt width 1', 'Double Voigt gamma 1', 'Double Voigt amplitude 2', 'Double Voigt center 2', 'Double Voigt width 2', 'Double Voigt gamma 2'], 
#               ['Counts', 'nm', 'nm', 'nm', 'Counts', 'nm', 'nm', 'nm'], 
#               getdoublevoigtfwhm, 
#               0
#               ] 

           'doublevoigt5params': [
               doublevoigt5paramwind, 
               fitdoublevoigt5paramtospec, 
               getmaxdoublevoigt5param, 
                'Double Voigt 5-Parameter Fit', 
                5, 
                ['Amplitude 1', 'Center 1', 'Width 1', 'Amplitude 2', 'Center 2'], 
                ['Counts', 'nm', 'nm', 'Counts', 'nm'], 
                getdoublevoigt5paramfwhm,
                0
           ]
           }

# units to add fit is not defined in context so I rely on surrounding code to handle units
fitunits = {'lorentz': fitkeys['lorentz'][6][:]+ unitstoaddfit,
            'gaussian': fitkeys['gaussian'][6][:] + unitstoaddfit,
            'voigt': fitkeys['voigt'][6][:] + unitstoaddfit,
            'linear': fitkeys['linear'][6][:] + unitstoaddfit,
            'double lorentz': fitkeys['double lorentz'][6][:] + unitstoaddfit,
            'double gaussian': fitkeys['double gaussian'][6][:] + unitstoaddfit,
            'double voigt': fitkeys['double voigt'][6][:] + unitstoaddfit, 
            'oscillations': fitkeys['oscillations'][6][:]+ unitstoaddfit,
            'stiffness': fitkeys['stiffness'][6][:] + unitstoaddfit,
            'derivative': fitkeys['derivative'][6][:] + unitstoaddfit,
            'moments': fitkeys['moments'][6][:] + unitstoaddfit,
            'com': fitkeys['com'][6][:] + unitstoaddfit,
            'decay': fitkeys['decay'][6][:] + unitstoaddfit,
            'binning': fitkeys['binning'][6][:] + unitstoaddfit,
            'derivative1': fitkeys['derivative1'][6][:] + unitstoaddfit,
            'derivative2': fitkeys['derivative2'][6][:] + unitstoaddfit,
            'derivative_points': fitkeys['derivative_points'][6][:] + unitstoaddfit,
            'doublevoigt5params': fitkeys['doublevoigt5params'][6][:] + unitstoaddfit
            }

# fitparametersparis: dict of the fit parameters and their units. Key of getlistofallFitparameters() is the key of the fitkeys dictionary
# fitunitparis: dict of the fit parameters and their units. Key of getlist
fitunitparis = {}
all_fit_params = getlistofallFitparameters()
for i in range(len(all_fit_params)):
    fit_type = list(fitkeys.keys())[i]
    #print(i, fit_type, all_fit_params[i])
    #for j, param_name in enumerate(all_fit_params[i]):
    #    if j < fitkeys[fit_type][4]:  # Regular fit parameters
    #        fitunitparis[param_name] = fitunits[fit_type][j]
    #    else:  # Additional parameters (ss_res, ss_tot, etc.)
    #        fitunitparis[param_name] = fitunits[fit_type][j]

if __name__ == '__main__':
    print('This is a library of window functions and their corresponding fit functions.')
    print('Use the fitkeys dictionary to access the functions.')

    ''' print fit units:'''
    print('Fit units:')
    for key in fitunits.keys():
        print(f"{key}: {fitunits[key]}")

    print('Fit unit paris:')
    for key in fitunitparis.keys():
        print(f"{key}: {fitunitparis[key]}")
    
    ''' print fit parameters:
    print('Fit parameters:')
    for key in fitkeys.keys():
        print(f"{key}: {fitkeys[key][5]}")

    #print fit keys:
    print('Fit keys:')
    for key in fitkeys.keys():
        print(key)
    '''

    #print fit parameters:
    #print(getlistofallFitparameters())
    #print(len(getlistofallFitparameters()))
    