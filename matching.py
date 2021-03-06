"""

ariel@iag
09/2019

Find Gaia sources in S-PLUS catalog to perform calibration

"""

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.table import Table, hstack
import matplotlib.pyplot as plt


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

    splus_coords = SkyCoord(ra=ra, dec=dec, unit=(u.hour, u.deg), frame='icrs', equinox='J2000')

    vizier_search = Vizier(columns=['*', 'RAJ2000', 'DEJ2000'], catalog='I/345')
    vizier_search.ROW_LIMIT = 999999999

    gaia_data = vizier_search.query_region(splus_coords, radius=Angle(1.0, "deg"))[0]

    # Only keep objects with detected proper motions to exclude contamination by galaxies
    # FIXME: We are selecting objects with pre-calculated radius
    flag = (gaia_data['pmRA'] > 0) & (gaia_data['Rad'] > 0.5) & (gaia_data['Rad'] < 100.0)
    gaia_data = gaia_data[flag]

    print('>>> Found', len(gaia_data), ' matches with Gaia')

    return gaia_data


def plot_hr(gaia_data):
    """

    Plot HR diagram for quality control

    Parameters:
    ------------

    gaia_data: output of find_gaia_stars

    """

    plt.scatter(gaia_data['BP-RP'], 5 + 5 * np.log10(gaia_data['Plx']) - 15 + gaia_data['Gmag'], s=5)
    plt.gca().invert_yaxis()

    plt.show()


def match_splus_gaia(splus_catalog, gaia_catalog):
    """

    Match the splus and gaia catalogs and return the matched table

    """

    print('matching SPLUS with Gaia...')

    a = SkyCoord(ra=splus_catalog['ALPHA_J2000'], dec=splus_catalog['DELTA_J2000'], unit=(u.deg, u.deg))
    b = SkyCoord(ra=gaia_catalog['RAJ2000'], dec=gaia_catalog['DEJ2000'], unit=(u.deg, u.deg))
    idx, d2d, d3d = a.match_to_catalog_3d(b)

    mask = d2d < 1.0 * u.arcsec
    matched_catalog = hstack([splus_catalog[mask], gaia_catalog[idx[mask]]])

    return matched_catalog


# This is just for testing purposes:
if __name__ == '__main__':
    plt.ion()

    data_path = 'data/STRIPE82-0001/'

    tile_coords = Table.read('data/tiles_new20190701_clean.csv', format='ascii.csv')

    test_field = tile_coords['NAME'] == 'STRIPE82_0001'

    test_ra, test_dec = tile_coords['RA'][test_field], tile_coords['DEC'][test_field]

    catalog = read_splus_catalog(data_path, 'STRIPE82-0001', 'R')
    catalog_stars = find_splus_stars(catalog)
    gaia_data = find_gaia_stars(test_ra, test_dec)
    match_data = match_splus_gaia(catalog_stars, gaia_data)

    plot_hr(gaia_data)

    # print('saving fig', datapath + 'dif_coords.png')
    # plt.scatter(matched_catalog['ALPHA_J2000'] - matched_catalog['RAJ2000'],
    #             matched_catalog['DELTA_J2000'] - matched_catalog['DEJ2000'],
    #             s=15, color='b')
    # plt.xlabel(r'$\Delta \mathrm{RA\ [SPLUS - Gaia]}$')
    # plt.ylabel(r'$\Delta \mathrm{DEC\ [SPLUS - Gaia]}$')
    # plt.grid()
    # plt.savefig(datapath + 'dif_coords.png', format='png', dpi=120)
    # plt.close()
