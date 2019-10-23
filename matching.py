"""

ariel@iag
09/2019

Find Gaia sources in S-PLUS catalog to perform calibration

"""

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astropy.io import fits
from astropy.table import Table, join, vstack

Vizier.ROW_LIMIT = 5000


def read_splus_catalog(data_path, field, band):
    """

    Reads splus catalogs provided by the pipeline

    Parameters:
    ------------
    data_path: str
        Path to data directory
    field: str
        S-PLUS field
    band: str
        Photometric band in capital letters

    ------------
    return: astropy table
        Catalog of detected sources
    """

    data = fits.open(data_path+field+'_'+band+'_dual.catalog')[2].data

    catalog = Table()

    for field in data.dtype.names:
        catalog[field] = data[field]

    return catalog


def find_splus_stars(catalog, criterion=0.9):
    """

    Finds stars using Source Extractor CLASS_STAR

    Parameters:
    ------------

    catalog: astropy table
        Catalog of sources
    criterion: float
        Criterion to identify stars

    ------------
    return: structured array
        Catalog of stars in the field
    """

    flag_star = (catalog['CLASS_STAR'] > criterion) & (catalog['FLAGS'] == 0)

    print('>>> We found ', flag_star.sum(), ' stars in the field')

    catalog_stars = catalog[flag_star]

    return catalog_stars


def find_gaia_stars(ra, dec):
    """

    FIXME: This is a preliminary implementation that searches gaia sources around the average RA and Dec of the field
           and finds matches with the input catalog from S-PLUS

    Finds Gaia data for stars in an S-PLUS field

    Parameters:
    ------------

    catalog_stars: structured array
        Catalog of stars in an S-PLUS field

    ------------
    return:
        Catalog of Gaia sources

    """

    splus_coords = SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg), frame='icrs')

    gaia_data = Vizier.query_region(splus_coords, radius=1 * u.deg, catalog=['gaia'])['I/345/gaia2']

    # Only keep objects with detected proper motions to exclude contamination by galaxies
    # FIXME: We are selecting objects with pre-calculated radius
    pm_flag = (gaia_data['pmRA'] > 0) & (gaia_data['Rad'] > 0.5) & (gaia_data['Rad'] < 1.5)
    gaia_data = gaia_data[pm_flag]

    print('>>> Found', len(gaia_data), ' matches with Gaia')

    return gaia_data


def plot_hr(gaia_data):
    """

    Plot HR diagram for quality control

    Parameters:
    ------------

    gaia_data: output of find_gaia_stars

    """

    plt.scatter(gaia_data['BP-RP'], 5 + 5 * np.log10(gaia_data['Plx']) - 15 +  gaia_data['Gmag'], s=5)
    plt.gca().invert_yaxis()

    plt.show()


# This is just for testing purposes:
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()

    tile_coords = Table.read('data/tiles_new20190701_clean.csv', format='ascii.csv')

    test_field = tile_coords['NAME'] == 'STRIPE82_0001'

    test_ra, test_dec = tile_coords['RA'][test_field], tile_coords['DEC'][test_field]

    catalog = read_splus_catalog('data/STRIPE82-0001/', 'STRIPE82-0001', 'R')
    catalog_stars = find_splus_stars(catalog)
    gaia_data = find_gaia_stars(catalog_stars, test_ra, test_dec)

    plot_hr(matched_catalog)

    # This is just some thoughts about models:
    #
    # model = fits.open('data/s_coelho14_sed/t09000_g+2.0_m05p00_sed.fits')
    #
    # log_wl = np.array([model[0].header['CRVAL1'] + i * model[0].header['CDELT1']
    #                    for i in range(model[0].header['NAXIS1'])])
    #
    # flux = model[0].data



