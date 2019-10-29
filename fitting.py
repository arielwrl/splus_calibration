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


def read_models_ngsl(model_dir, table_fname):
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

    model_table = Table.read(table_fname)

    model_dict = {'name': np.array([model_table['Target Name'][i].lower() for i in range(len(model_table))]),
                  'wl': np.zeros((len(model_table), 2885)),
                  'flux': np.zeros((len(model_table), 2885))
                  }

    for i in range(len(model_table)):
        print('>>> Reading model ', i + 1, model_dict['name'][i])

        model = fits.open(model_dir + 'h_stis_ngsl_' + model_dict['name'][i].strip() + '_v2.fits')

        wl = np.array([model[1].data[j][0] for j in range(2885)])
        flux = np.array([model[1].data[j][1] for j in range(2885)])

        model_dict['wl'][i] = wl
        model_dict['flux'][i] = flux

    return model_dict


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


def calc_syn_mags(model_dict, lib='ngsl'):
    """

    Calculates synthetic magnitudes normalized by stellar radius

    :param model_dict:
    :return:
    """

    if lib == 'coelho2014':

        rp_syn = np.array([synmag(model_dict['lin_wl'], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.Grp.dat')
                           for i in range(len(model_dict['flux']))])

        bp_syn = np.array([synmag(model_dict['lin_wl'], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.Gbp.dat')
                           for i in range(len(model_dict['flux']))])

        g_syn = np.array([synmag(model_dict['lin_wl'], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.G.dat')
                          for i in range(len(model_dict['flux']))])

    if lib == 'ngsl':
        rp_syn = np.array([synmag(model_dict['wl'][i], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.Grp.dat')
                           for i in range(len(model_dict['flux']))])

        bp_syn = np.array([synmag(model_dict['wl'][i], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.Gbp.dat')
                           for i in range(len(model_dict['flux']))])

        g_syn = np.array([synmag(model_dict['wl'][i], model_dict['flux'][i], 'data/filters/GAIA_GAIA2.G.dat')
                          for i in range(len(model_dict['flux']))])

    synmags = {'G': g_syn,
               'BP': bp_syn,
               'RP': rp_syn
               }

    return synmags


def fit_models_coelho2014(model_dict, gaia_data):
    """

    FIXME: Ongoing work

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


def fit_models_ngsl(ngsl_table, gaia_matched_catalog):
    """

    FIXME: Using gaia observed magnitudes for ngsl stars, models are slightly bluer than gaia data.

    Fits NGSL models to gaia Data

    Parameters
    ------------
    model_dict:
    ngsl_table:
    gaia_catalog:

    ------------
    return:

    """

    model_colors = ngsl_table['bp_rp']
    obs_colors = gaia_matched_catalog['BP-RP']

    model_mags = 5 + 5 * np.log10(ngsl_table['parallax']) - 15 + ngsl_table['phot_g_mean_mag']
    obs_mags = 5 + 5 * np.log10(gaia_matched_catalog['Plx']) - 15 + gaia_matched_catalog['Gmag']

    gaia_matched_catalog['best_fit'] = np.zeros(len(gaia_matched_catalog), dtype=int)

    for i in range(len(gaia_matched_catalog)):
        print('>>> Fitting ', i+1, 'of ', len(gaia_matched_catalog))

        chi2s = np.array([np.sqrt((obs_colors[i]-model_colors[j])**2 + (obs_mags[i]-model_mags[j])**2)
                          for j in range(len(model_colors))])

        best_chi2 = np.min(chi2s[~np.isnan(chi2s)])
        gaia_matched_catalog['best_fit'][i] = np.argwhere(chi2s == best_chi2)[0][0]

    # Let's store the distance to the templates and to each star, which will be necessary to find zero points

    d_star = 3.085677e18 * 1000 / gaia_matched_catalog['Plx']
    d_model = 3.085677e18 * 1000 / ngsl_table['parallax'][gaia_matched_catalog['best_fit']]

    gaia_matched_catalog['d_star'] = d_star
    gaia_matched_catalog['d_model'] = d_model

    gaia_matched_catalog['fitted_template_wl'] = model_dict['wl'][fitted_catalog['best_fit']]
    gaia_matched_catalog['fitted_template_flux'] = model_dict['flux'][fitted_catalog['best_fit']]

    normalization = 4 * np.pi * (d_model/d_star)**2

    flux_obsframe = gaia_matched_catalog['fitted_template_flux'] * normalization[:, np.newaxis]

    gaia_matched_catalog['fitted_template_obsframe'] = flux_obsframe

    # Synphot in models

    filter_list = ['SPLUS_F378', 'SPLUS_F395', 'SPLUS_F410', 'SPLUS_F430', 'SPLUS_F515', 'SPLUS_F660', 'SPLUS_F861',
                   'SPLUS_U', 'SPLUS_G', 'SPLUS_R', 'SPLUS_I', 'SPLUS_Z']
    filter_dir = 'data/filters/'

    filter_dict = {}
    for filt in filter_list:
        filter_dict[filt] = np.genfromtxt(filter_dir + filt + '.dat').transpose()

    gaia_matched_catalog['syn_mags_splus'] = np.zeros((len(gaia_matched_catalog), 12))
    for i in range(len(gaia_matched_catalog)):

        wl, flux = gaia_matched_catalog['fitted_template_wl'][i], gaia_matched_catalog['fitted_template_obsframe'][i]

        m_syn = np.array([synmag(wl, flux, filter_curve=filter_dict[filter_list[j]]) for j in range(12)])

        gaia_matched_catalog['syn_mags_splus'][i] = m_syn

    return gaia_matched_catalog


def splus_synmags(fitted_catalog, models_dict):
    """

    Calculates splus magnitudes for a set of fitted models

    :param fitted_catalog:
    :param models_dict:
    :return:
    """

    # Put models in the distance of the observed_stars





if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matching import *
    from astropy.table import Table

    # For Coelho 2014 models
    # FIXME: Coelho (2014) fitting procedure to be done
    #
    # model_dir = 'data/s_coelho14_sed/'
    #
    # models = read_models_coelho2014(model_dir)
    #
    # for i in [0, 10, 100, 1000, 1500, 2000, 2500, 3000]:
    #     plt.plot(models['lin_wl'], models['flux'][i])
    #
    # tile_coords = Table.read('data/tiles_new20190701_clean.csv', format='ascii.csv')
    # test_field = tile_coords['NAME'] == 'STRIPE82_0001'
    # test_ra, test_dec = tile_coords['RA'][test_field], tile_coords['DEC'][test_field]
    #
    # gaia_data = find_gaia_stars(test_ra, test_dec)

    # For NGSL models
    model_dir = '/home/ariel/Workspace/S-PLUS/splus_calibration/data/NGSL/stis_ngsl_v2/'
    table_fname = 'data/NGSL/ngsl_gaia.fits'

    ngsl_table = Table.read(table_fname)

    models = read_models_ngsl(model_dir, table_fname)

    for i in [0, 10, 50, 60, 70, 100, 150, 200, 250, 300, 350]:
        plt.plot(models['wl'][i], models['flux'][i])

    model_mags = calc_syn_mags(models)

    tile_coords = Table.read('data/tiles_new20190701_clean.csv', format='ascii.csv')
    test_field = tile_coords['NAME'] == 'STRIPE82_0001'
    test_ra, test_dec = tile_coords['RA'][test_field], tile_coords['DEC'][test_field]

    catalog = read_splus_catalog('data/STRIPE82-0001/', 'STRIPE82-0001', 'R')
    catalog_stars = find_splus_stars(catalog)
    gaia_data = find_gaia_stars(test_ra, test_dec)
    matched_data = match_splus_gaia(catalog_stars, gaia_data)

    fitted_catalog = fit_models_ngsl(ngsl_table, matched_data)

    plot_hr(matched_data)
    plt.scatter(ngsl_table['bp_rp'], 5 + 5 * np.log10(ngsl_table['parallax']) - 15 + ngsl_table['phot_g_mean_mag'], s=15)
    plt.scatter(model_mags['BP']-model_mags['RP'], 5 + 5 * np.log10(ngsl_table['parallax']) - 15 + model_mags['G'], s=5)

    best_fits = fitted_catalog['best_fit']
    plt.scatter(ngsl_table['bp_rp'][best_fits], 5 + 5 * np.log10(ngsl_table['parallax'])[best_fits]
                - 15 + ngsl_table['phot_g_mean_mag'][best_fits], s=2)





