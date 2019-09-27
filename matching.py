"""

ariel@iag
09/2019

Find Gaia sources in S-PLUS catalog to perform calibration

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from astropy.io import fits
from astropy.table import Table, join


def read_splus_catalog(data_path, field, band):
    """

    Reads splus catalogs provided by the pipeline

    :param data_path: str
        Path to data directory
    :param field: str
        S-PLUS field
    :param band: str
        Photometric band in capital letters
    :return: astropy table
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

    :param catalog: astropy table
        Catalog of sources
    :param criterion: float
        Criterion to identify stars
    :return: structured array
        Catalog of stars in the field
    """

    flag_star = catalog['CLASS_STAR'] > criterion

    print('>>> We found ', flag_star.sum(), ' stars in the field')

    catalog_stars = catalog[flag_star]

    return catalog_stars


def find_gaia_stars(catalog_stars):
    """

    FIXME: This is a preliminary implementation that searches gaia sources around the average RA and Dec of the field
           and finds matches with the input catalog from S-PLUS

    Finds Gaia data for stars in an S-PLUS field

    :param catalog_stars: structured array
        Catalog of stars in an S-PLUS field
    :return:
        Catalog of Gaia sources

    """

    coord = SkyCoord(ra=np.mean(catalog_stars['ALPHA_J2000']), dec=np.mean(catalog_stars['DELTA_J2000'])
                     , unit=(u.degree, u.degree), frame='icrs')
    gaia_data = Gaia.cone_search(coordinate=coord, radius=u.Quantity(2, u.degree))

    matched_catalog = join(catalog_stars, gaia_data.get_data(), join_type='left')

    return matched_catalog


# This is just for testing purposes:
if __name__ == '__main__':

    catalog = read_splus_catalog('data/STRIPE82-0001/', 'STRIPE82-0001', 'R')
    catalog_stars = find_splus_stars(catalog)
    matched_catalog = find_gaia_stars(catalog_stars)

