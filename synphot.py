"""

ariel@ufsc
May 2017

Provides tools for synthetic photometry

"""


import numpy as np
from scipy.interpolate import interp1d


def resampler(x_old, y_old, x_new):

    interp = interp1d(x_old, y_old, bounds_error=False, fill_value=(0., 0.))

    y_new = interp(x_new)
    
    return y_new
    

def synflux(wl, flux, filter_curve):
    
    # Reading filter if needed:
    if type(filter_curve) is str: 
        wl_filter, transmitance = np.genfromtxt(filter_curve).transpose()
    else:
        wl_filter, transmitance = filter_curve[0], filter_curve[1]
    
    # Resampling filter and spectrum to 1\AA intervals:
    wl_new = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5)
    
    T = resampler(wl_filter, transmitance, wl_new)
    flux = resampler(wl, flux, wl_new)
    
    # Convolution:
    synflux = np.trapz(flux * transmitance * wl_new, dx=1)
    synflux /= np.trapz(wl_new * transmitance, dx=1)

    return synflux


def synmag(wl, flux, filter_curve):

    # Reading filter if needed:
    if type(filter_curve) is str: 
        wl_filter, transmitance = np.genfromtxt(filter_curve).transpose()
    else:
        wl_filter, transmitance = filter_curve[0], filter_curve[1]

    # Resampling filter and spectrum to 1 Angstrom intervals:
    wl_new = np.arange(np.round(wl_filter[0]) - 5, np.round(wl_filter[-1]) + 5)

    transmitance = resampler(wl_filter, transmitance, wl_new)
    flux = resampler(wl, flux, wl_new)

    # Calculate magnitudes and errors:
    m_ab = -2.5 * np.log10(np.trapz(flux * transmitance * wl_new, dx=1) / np.trapz(transmitance / wl_new, dx=1)) - 2.41

    return m_ab


def pivot_wavelength(filter_curve):

    # Reading filter if needed:
    if type(filter_curve) is str: 
        wl_filter, transmitance = np.genfromtxt(filter_curve).transpose()
    else:
        wl_filter, transmitance = filter_curve[0], filter_curve[1]
    
    # Calculating pivot_wavelength    
    pivot_wl = np.trapz(transmitance * wl_filter, dx=1) / np.trapz(transmitance * (wl_filter**-1), dx=1)
    pivot_wl = np.sqrt(pivot_wl)
    
    return pivot_wl


def effective_wavelength(wl, spectrum, filter_curve):
    """

    This is defined as the mean wavelength of the filter weighted by transmission of the filter
    and spectrum of the source


    :param wl:
    :param spectrum:
    :param filter_curve:
    :return:
    """

    # Reading filter if needed:
    if type(filter_curve) is str: 
        wl_filter, transmitance = np.genfromtxt(filter_file).transpose()
    else:
        wl_filter, transmitance = filter_curve[0], filter_curve[1]
    
    # Resampling filter and spectrum to 1\AA intervals:
    wl_filter = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5,1)
    
    transmitance = resampler(wl_filter, transmitance, wl_filter)
    spectrum = resampler(wl, spectrum, wl_filter)
    
    # Convolution time:
    effective_wl = np.trapz(wl_filter**2 * spectrum * transmitance, dx=1)
    effective_wl /= np.trapz(spectrum * transmitance * wl_filter, dx=1)
    
    return effective_wl
