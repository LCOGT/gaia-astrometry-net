from glob import glob
from astropy.io import fits
import os
import logging
import multiprocessing as mp
from lcogt_logging import LCOGTFormatter
import datetime
import sys
import gzip
from astropy.io import ascii
from astrometry.util.util import *
import sqlite3
from astropy.table import Table

from gaia_astrometry_index_files.dbs import make_gaia_db
from gaia_astrometry_index_files.utils import nside_to_healpixels, great_circle_distance

logger = logging.getLogger('GaiaCatalogs')
logger.setLevel(logging.DEBUG)
logging.StreamHandler(sys.stdout)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = LCOGTFormatter()
handler.setFormatter(formatter)
logger.addHandler(handler)

fovs = {'p6': 60.0,
        'p5': 42.426406871192846,
        'p4': 24.0,
        'p3': 16.97056274847714,
        'p2': 12.0,
        'p1': 8.48528137423857,
        'p0': 6.0,
        'm1': 4.242640687119285,
        'm2': 3.0,
        'm3': 2.1213203435596424}


def create_index_files(catalog_directory='/net/fsfs.lco.gtn/data/AstroCatalogs/GAIA-DR2/gaia_source/csv'):
    """
    Make the gaia astrometry.net index files
    :param catalog_directory: directory with the catalog files
    :param n_sources_for_all_sky: number of bright sources that should be taken from each
                                  individual catalog to put in the all sky catalog
    :return:
    """
    # First let's make the data a little more useable by putting the catalog into an sqlite db.
    initial_catalogs = glob(os.path.join(catalog_directory, '*.csv.gz'))

    if not os.path.exists('gaia.db') or os.path.getsize('gaia.db') < 1e6:
        # Now we can make a database of catalog. This is going to be one heck of a db.
        make_gaia_db(initial_catalogs)

    healpixels = nside_to_healpixels(16)
    healpix_catalogs = make_gaia_healpix_catalogs(healpixels, db_address='sqlite:///gaia.db', ncpu=6)

    # Now we need to make index files down to scale=0
    # This is currently what is provided by the astrometry.net people and is
    # probably about 30 GB
    index_scales = np.arange(0, 7)[::-1]
    for scale in index_scales:
        make_individual_index_files(healpix_catalogs, scale, ncpu=6)

    index_scales = np.arange(-3, 0)[::-1]
    for scale in index_scales:
        make_individual_index_files(healpix_catalogs, scale, ncpu=6)

    save_index_file_meta_data('gaia-dr2-index-files.dat', healpixels, range(-3, 7), 16)
    # New index files! Woot!


def save_index_file_meta_data(meta_data_file_name, healpixels, scales, nside):
    meta_data = {'filename': [], 'ra': [], 'dec': [], 'radius': [], 'fov': []}
    for healpixel in healpixels:
        for scale in scales:
            meta_data['filename'].append(make_index_file_name(nside, scale, healpixel['index']))
            meta_data['ra'].append(healpixel['ra'])
            meta_data['dec'].append(healpixel['dec'])
            meta_data['radius'].append(healpixel['radius'])
            meta_data['fov'].append(fovs[scale_to_filename_string(scale)])
    Table(meta_data).write(meta_data_file_name, format='ascii')


def save_catalog(catalog, counter):
    hdu = fits.BinTableHDU(catalog)
    output_name = 'gaia-fullcatalog-{counter}.fits'.format(counter=counter)
    hdu.writeto(output_name, clobber=True)
    return output_name


# For each heal pix
def make_individual_healpix_catalog(healpixel):
    db_address = 'sqlite:///gaia.db'
    if os.path.exists(get_healpix_catalog_name(healpix_id=healpixel['index'],
                                               nside=healpixel['nside'], allsky=False)):
        return get_healpix_catalog_name(healpix_id=healpixel['index'],
                                        nside=healpixel['nside'], allsky=False)
    else:
        # For each healpixel fits file get the sources in this heal pixel
        stars_in_healpixel = get_sources_in_healpixel(healpixel, db_address)
        stars_in_healpixel = Table(stars_in_healpixel)
        # Sort the healpix catalog descending by flux
        stars_in_healpixel.sort('flux')
        stars_in_healpixel.reverse()
        # Write out the healpix catalog to a fits file
        return write_out_healpix_catalog(stars_in_healpixel, healpixel['index'], healpixel['nside'])


def make_gaia_healpix_catalogs(healpixels, db_address='sqlite:///gaia.db', ncpu=6):
    p = mp.Pool(ncpu)
    catalog_names = p.map(make_individual_healpix_catalog, healpixels)
    p.close()
    return catalog_names


def parse_catalog(filename):
    with gzip.open(filename, 'rt') as csv_file:
        csv_data = csv_file.read()
    catalog_data = ascii.read(csv_data, 'fast_csv')
    return catalog_data


def query_sources(ramin, ramax, decmin, decmax, db_address='gaia.db'):
    connection = sqlite3.connect(db_address)
    cursor = connection.cursor()
    sql_command = 'select sources.g_flux, sources.ra, sources.dec from sources, positions ' \
                  'where positions.ramin >= {ramin} and positions.ramax <= {ramax} ' \
                  'and positions.decmin >= {decmin} and positions.decmax <= {decmax} ' \
                  'and positions.id = sources.id'
    sql_command = sql_command.format(ramin=ramin, ramax=ramax, decmin=decmin, decmax=decmax)
    cursor.execute(sql_command)
    rows = cursor.fetchall()
    connection.close()
    sources_in_healpix = {'flux': [], 'ra': [], 'dec':[]}
    for row in rows:
        sources_in_healpix['flux'].append(row[0])
        sources_in_healpix['ra'].append(row[1])
        sources_in_healpix['dec'].append(row[2])
    return sources_in_healpix


def get_sources_in_healpixel(healpixel, db_address):

    logger.info('Getting sources in healpixel: {hp_id}'.format(hp_id=healpixel['index']))
    # Open the file

    decmin = healpixel['dec'] - healpixel['radius']
    decmax = healpixel['dec'] + healpixel['radius']
    # See http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates#RefBronstein
    delta_ra = np.rad2deg(np.arcsin(np.sin(np.deg2rad(healpixel['radius'])) / np.cos(np.deg2rad(healpixel['dec']))))
    ramin = healpixel['ra'] - delta_ra
    ramax = healpixel['ra'] + delta_ra

    # If one of the poles is included, just do all sources in the polar cap
    if decmin < -90.0 or decmax > 90.0:
        ramin = 0.0
        ramax = 360.0
        sources_in_healpix = query_sources(ramin, ramax, decmin, decmax)
    elif ramin < 0.0:
        sources_in_healpix = query_sources(ramin, ramax, decmin, decmax)
        ramin = 360.0 - delta_ra
        ramax = 360.0
        sources_from_wrap = query_sources(ramin, ramax, decmin, decmax)
        for key in sources_in_healpix:
            sources_in_healpix[key] += sources_from_wrap[key]

    elif ramax > 360.0:
        sources_in_healpix = query_sources(ramin, ramax, decmin, decmax)
        ramin = 0.0
        ramax = delta_ra
        sources_from_wrap = query_sources(ramin, ramax, decmin, decmax)
        for key in sources_in_healpix:
            sources_in_healpix[key] += sources_from_wrap[key]

    else:
        sources_in_healpix = query_sources(ramin, ramax, decmin, decmax)

    return sources_in_healpix


def make_all_sky_catalog(catalog_names, n_sources_per_catalog):
    logger.info('Making all sky catalog')
    initial_catalog = np.array(fits.getdata(catalog_names[0], 1))
    allsky_catalog = np.zeros(n_sources_per_catalog * len(catalog_names),
                              dtype=initial_catalog.dtype)
    del initial_catalog
    # For each healpix catalog
    counter = 0
    for filename in catalog_names:
        logger.info('Adding {filename} to all sky catalog'.format(filename=filename))
        catalog_data = np.array(fits.getdata(filename, 1))
        if len(catalog_data) >= n_sources_per_catalog:
            # Save the brightest stars in the all sky catalog
            allsky_catalog[counter:counter+n_sources_per_catalog] = catalog_data[:n_sources_per_catalog]
            counter += n_sources_per_catalog
        else:
            allsky_catalog[counter:counter+len(catalog_data)] = catalog_data
            counter += len(catalog_data)

    allsky_catalog = allsky_catalog[:counter]
    # Remove duplicates
    logger.info('Removing duplicates from all sky catalog')
    allsky_catalog = remove_duplicate_sources(allsky_catalog, len(catalog_names))
    logger.info('Sorting all sky catalog')
    # Sort the all sky catalog by flux
    sorted_indices = np.argsort(allsky_catalog['phot_g_mean_flux'])[::-1]
    # Save the all sky catalog
    return write_out_healpix_catalog(allsky_catalog[sorted_indices], allsky=True)


def remove_duplicate_sources(catalog, n_healpix):
    # Sort the catalog by ra, dec
    sorted_indices = np.lexsort((catalog['dec'], catalog['ra']))
    catalog = catalog[sorted_indices]

    is_duplicate = np.zeros(len(catalog), dtype=bool)
    # Go through each source
    for i, source in enumerate(catalog):
        # If the source is already flagged as a duplicate, move on
        if not is_duplicate[i]:
            # Check the next sources to see if they have the same ra and dec
            # We only need to check for the number of healpixels because that
            # will be the maximum number of duplicates
            # otherwise the computation gets out of hand
            offsets = great_circle_distance(catalog[i + 1:i + n_healpix]['ra'],
                                            catalog[i+1:i+n_healpix]['dec'],
                                            source['ra'], source['dec'])
            is_duplicate[i+1:i+n_healpix] = offsets < 1e-6

    # Remove the duplicates
    return catalog[~is_duplicate]


def make_all_sky_index_files(catalog, scale):
    # Run the index making program from astrometry.net on the all sky catalog
    logger.info('Making all sky index file for scale: {scale}'.format(scale=scale))
    scale_filename_str = scale_to_filename_string(scale)
    index_name = 'gaia-index-{scale}.fits'.format(scale=scale_filename_str)
    if not os.path.exists(index_name):
        os.system('build-astrometry-index -i {catalog} -o {index_name} -P {scale} '
                  '-A ra -D dec -I {unique_id} -j 0.1 -M &> {log_name}'
                  ''.format(catalog=catalog, index_name=index_name,
                            scale=scale, unique_id=make_unique_id(scale),
                            log_name=index_name.replace('.fits', '.log')))


def make_individual_index_files(catalogs, scale, margin=1, ncpu=6):
    p = mp.Pool(ncpu)
    p.map(make_single_index_file, [(catalog, scale, margin) for catalog in catalogs])
    p.close()


def make_index_file_name(nside, scale, healpixel_id):
    return 'gaia-index-{scale}-{nside}-{hp_id}.fits'.format(nside=nside, hp_id=healpixel_id,
                                                            scale=scale_to_filename_string(scale))


def make_single_index_file(args):
    catalog, scale, margin = args
    nside, healpix_id = catalog[:-5].split('-')[2:4]
    logger.info('Making index file. nside: {nside}, healpix: {hp_id}, scale: {scale} '.format(nside=nside,
                                                                                              hp_id=healpix_id,
                                                                                              scale=scale))
    index_name = make_index_file_name(nside, scale, healpix_id)

    if not os.path.exists(index_name):
        # For the smaller indexes run the astrometry.net index maker on the healpix files
        os.system('build-astrometry-index -i {catalog} -o {index_name} -P {scale} '
                  '-H {healpix_id} -s {nside} -m {margin} -A ra -D dec -I {unique_id} -j 0.1 '
                  '&> {log_name}'.format(catalog=catalog, index_name=index_name, scale=scale,
                                         margin=margin, nside=nside, healpix_id=healpix_id,
                                         unique_id=make_unique_id(scale),
                                         log_name=index_name.replace('.fits', '.log')))


def scale_to_filename_string(scale):
    return '{0:+d}'.format(scale).replace('+', 'p').replace('-', 'm')


def write_out_healpix_catalog(catalog, healpix_id=0, nside=0, allsky=False):
    table_hdu = fits.BinTableHDU(catalog)
    catalog_name = get_healpix_catalog_name(healpix_id=healpix_id, nside=nside, allsky=allsky)
    table_hdu.writeto(catalog_name, overwrite=True)
    return catalog_name


def get_healpix_catalog_name(healpix_id=0, nside=0, allsky=False):
    if allsky:
        catalog_name = 'gaia-catalog-allsky.fits'
    else:
        catalog_name = 'gaia-catalog-{nside}-{hp_id}.fits'.format(hp_id=healpix_id, nside=nside)
    return catalog_name


def make_unique_id(scale):
    now = datetime.datetime.now()
    # 0 for minus 1 for plus
    unique_id = '{yr}{m}{d}{s:+03d}'.format(yr=now.year % 2000, m=now.month, d=now.day, s=scale)
    unique_id = unique_id.replace('+', '1')
    unique_id = unique_id.replace('-', '0')
    return unique_id
