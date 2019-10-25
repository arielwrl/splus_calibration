"""

ariel@iag
10/2019

Fits stellar models to Gaia color-magnitude diagrams

"""

import numpy as np
import os
from astropy.io import fits
from synphot import synmag


def read_models_coelho2014(model_dir):
    """

    Reads the Coelho (2014) models to a dictionary

    Parameters
    ------------
    model_dir: str
        directory on which models are stored

    ------------
    return: dict
        Dictionary with models

    """

    file_list = os.listdir(model_dir)

    # Read a reference model to get wavelngths
    model = fits.open(model_dir+file_list[0])

    log_wl = np.array([model[0].header['CRVAL1'] + i * model[0].header['CDELT1']
                       for i in range(model[0].header['NAXIS1'])])

    lin_wl = 10**log_wl

    # Create and fill dictionary

    models_dict = {'file': np.array(file_list),
                   'log_wl': log_wl,
                   'lin_wl': lin_wl,
                   'flux': np.zeros((len(file_list), len(log_wl)))
                   }

    for i in range(len(file_list)):
        print('>>> Reading model ', i+1, ' of ', len(file_list)+1, ':', file_list[i])

        model = fits.open(model_dir+file_list[i])

        flux = model[0].data

        models_dict['flux'][i] = flux

    return models_dict


def normalize_gaia_data(gaia_data):
    """

    Normalizing stellar data by radius to compare with theoretical models

    :param gaia_data:
    :return:
    """

    gaia_norm_mags = {'G': gaia_data['Gmag'] - 2.5 * np.log10(4 * np.pi * gaia_data['Rad'] ** 2),
                      'RP': gaia_data['RPmag'] - 2.5 * np.log10(4 * np.pi * gaia_data['Rad'] ** 2),
                      'BP': gaia_data['BPmag'] - 2.5 * np.log10(4 * np.pi * gaia_data['Rad'] ** 2)
                      }

    return gaia_norm_mags


def calc_syn_mags(model_dict):
    """

    Calculates synthetic magnitudes normalized by stellar radius

    :param model_dict:
    :return:
    """

    rp_syn = np.array([synmag(model_dict['lin_wl'], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.Grp.dat')
                       for i in range(len(model_dict['flux']))])

    bp_syn = np.array([synmag(model_dict['lin_wl'], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.Gbp.dat')
                       for i in range(len(model_dict['flux']))])

    g_syn = np.array([synmag(model_dict['lin_wl'], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.G.dat')
                      for i in range(len(model_dict['flux']))])

    synmags = {'G': g_syn,
               'BP': bp_syn,
               'RP': rp_syn
               }

    return synmags


def fit_models_coelho2014(model_dict, gaia_data):
    """

    Parameters
    ------------

    model_dict: dict
        Output from read_models_coelho2014

    gaia_data: astropy table
        Output from matching.find_gaia_stars

    ------------
    return: array-like
        File names of the models corresponding to Gaia observations

    """

    gaia_norm_mags = normalize_gaia_data(gaia_data)
    syn_mags = calc_syn_mags(model_dict)

    for i in range(len(gaia_data)):
        model_norm_mags = {'G': g_syn - 2.5 * np.log10(4 * np.pi * gaia_data['Rad'][i] ** 2),
                           'RP': rp_syn - 2.5 * np.log10(4 * np.pi * gaia_data['Rad'][i] ** 2),
                           'BP': bp_syn - 2.5 * np.log10(4 * np.pi * gaia_data['Rad'][i] ** 2)
                           }



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matching import find_gaia_stars
    from astropy.table import Table

    model_dir = 'data/s_coelho14_sed/'

    models = read_models_coelho2014(model_dir)

    for i in [0, 10, 100, 1000, 1500, 2000, 2500, 3000]:
        plt.plot(models['lin_wl'], models['flux'][i])

    tile_coords = Table.read('data/tiles_new20190701_clean.csv', format='ascii.csv')
    test_field = tile_coords['NAME'] == 'STRIPE82_0001'
    test_ra, test_dec = tile_coords['RA'][test_field], tile_coords['DEC'][test_field]

    gaia_data = find_gaia_stars(test_ra, test_dec)

