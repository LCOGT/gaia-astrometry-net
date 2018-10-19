import numpy as np
from glob import glob
from astropy.io import fits
from sqlalchemy import create_engine, pool
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, Integer, Float, Index
from sqlalchemy.ext.declarative import declarative_base
import os
import logging
import multiprocessing as mp
from lcogt_logging import LCOGTFormatter
import datetime
import sys
from astrometry.util.util import *

Base = declarative_base()

logger = logging.getLogger('GaiaCatalogs')
logger.setLevel(logging.DEBUG)
logging.StreamHandler(sys.stdout)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = LCOGTFormatter()
handler.setFormatter(formatter)
logger.addHandler(handler)


def create_index_files(catalog_directory='/home/cmccully/gaia_20180201/raw/',
                       n_sources_for_all_sky=100000):
    """
    Make the gaia astrometry.net index files
    :param catalog_directory: directory with the catalog files
    :param n_sources_for_all_sky: number of bright sources that should be taken from each
                                  individual catalog to put in the all sky catalog
    :return:
    """
    # First let's make the data a little more useable by putting the catalog into fewer,
    # bigger files
    initial_catalogs = glob(os.path.join(catalog_directory, 'GaiaSource*.fits'))
    # Check to see if the munged catalogs already exist
    munged_catalogs = glob('gaia-fullcatalog*.fits')
    if len(munged_catalogs) != 53:
        munged_catalogs = munge_initial_small_catalogs(initial_catalogs)

    healpixels = nside_to_healpixels(8)
    healpix_catalogs = make_gaia_healpix_catalogs(healpixels, munged_catalogs, ncpu=6)

    if not os.path.exists('gaia-catalog-allsky.fits'):
        all_sky_catalog = make_all_sky_catalog(healpix_catalogs, n_sources_for_all_sky)
    else:
        all_sky_catalog = 'gaia-catalog-allsky.fits'

    # Start from the large scales and make the index files from the all sky catalog
    index_scales = np.arange(7, 20)[::-1]
    for scale in index_scales:
        make_all_sky_index_files(all_sky_catalog, scale)

    # Now we need to make index files down to scale=0
    # This is currently what is provided by the astrometry.net people and is
    # probably about 30 GB
    index_scales = np.arange(0, 7)[::-1]
    for scale in index_scales:
        make_individual_index_files(healpix_catalogs, scale, ncpu=6)

    if not os.path.exists('gaia.db'):
        # Now we can make a database of catalog. This is going to be one heck of a db.
        make_gaia_db('sqlite:///gaia.db', munged_catalogs)


    index_scales = np.arange(-3, 0)[::-1]
    for scale in index_scales:
        make_individual_index_files(healpix_catalogs, scale, ncpu=6)

    # New index files! Woot!



def munge_initial_small_catalogs(initial_catalogs):
    output_catalogs = []
    first_catalog = parse_catalog(initial_catalogs[0])
    counter = 0
    useful_columns = ['ra', 'ra_error', 'dec', 'dec_error',
                      'phot_g_mean_flux', 'phot_g_mean_flux_error']
    munged_data = np.zeros(100 * len(first_catalog), dtype=first_catalog.dtype)
    munged_data = munged_data[useful_columns]
    sources_per_catalog = len(first_catalog)
    file_counter = 0
    # Go through each of the catalog files
    for catalog in initial_catalogs[:-1]:
        logger.info('Munging data from {filename}'.format(filename=catalog))
        # Save the ra, dec, flux, and corresponding errors in batches of 100
        catalog_data = parse_catalog(catalog)
        munged_catalog_location = slice(counter * sources_per_catalog,
                                        (counter + 1) * sources_per_catalog)
        for column in useful_columns:
            munged_data[column][munged_catalog_location] = catalog_data[column]
        counter += 1
        if counter == 100:
            output_catalogs.append(save_catalog(munged_data, file_counter))
            counter = 0
            file_counter += 1

    # Write out the last file
    last_catalog = parse_catalog(initial_catalogs[-1])
    munged_catalog_location = slice(counter * sources_per_catalog,
                                    counter * sources_per_catalog + len(last_catalog))
    for column in useful_columns:
        munged_data[column][munged_catalog_location] = last_catalog[column]
    munged_data = munged_data[:counter * sources_per_catalog + len(last_catalog)]
    output_catalogs.append(save_catalog(munged_data, file_counter))
    return output_catalogs


def save_catalog(catalog, counter):
    hdu = fits.BinTableHDU(catalog)
    output_name = 'gaia-fullcatalog-{counter}.fits'.format(counter=counter)
    hdu.writeto(output_name, clobber=True)
    return output_name


def nside_to_healpixels(nside):
    healpixel_indexes = np.arange(12 * nside**2)

    # Get the pixel positions of each healpix
    healpix_ras = [healpix_to_radec(i, nside, 0.5, 0.5)[0] * 180.0 / np.pi for i in healpixel_indexes]
    healpix_decs = [healpix_to_radec(i, nside, 0.5, 0.5)[1] * 180.0 / np.pi for i in healpixel_indexes]

    # Get the radias+1 degree overlap for nside=2
    ras_along_center_ring = [healpix_ras[i] for i in healpixel_indexes if healpix_decs[i] == 0.0]

    healpix_radius = np.abs(ras_along_center_ring[0] - ras_along_center_ring[1]) / 2.0
    return [{'index': i, 'ra': healpix_ras[i], 'dec': healpix_decs[i],
             'radius': healpix_radius, 'nside': nside} for i in healpixel_indexes]


class Star(Base):
    __tablename__ = 'sources'
    id = Column(Integer, primary_key=True, autoincrement=True)
    ra = Column(Float) #, index=True
    dec = Column(Float) # , index=True
    ra_error = Column(Float)
    dec_error = Column(Float)
    g_flux = Column(Float)
    g_flux_error = Column(Float)


def great_circle_distance(ra1, dec1, ra2, dec2):
    ra1_rad, dec1_rad = np.deg2rad([ra1, dec1])
    ra2_rad, dec2_rad = np.deg2rad([ra2, dec2])
    distance_rad = np.arccos(np.sin(dec1_rad) * np.sin(dec2_rad) + np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra2_rad - ra1_rad))
    return np.rad2deg(distance_rad)


def test_great_circle_distance():
    np.testing.assert_allclose(great_circle_distance(0, 0, 100, 0), 100.0, atol=1e-5)
    np.testing.assert_allclose(great_circle_distance(0, 0, 181, 0), 179.0, atol=1e-5)
    np.testing.assert_allclose(great_circle_distance(0, 0, 0, 10), 10, atol=1e-5)
    np.testing.assert_allclose(great_circle_distance(0, 75, 180, 75), 30.0, atol=1e-5)


def get_session(db_address):
    """
    Get a connection to the database.

    Returns
    -------
    session: SQLAlchemy Database Session
    """
    # Build a new engine for each session. This makes things thread safe.
    engine = create_engine(db_address, poolclass=pool.NullPool)
    Base.metadata.bind = engine

    # We don't use autoflush typically. I have run into issues where SQLAlchemy would try to flush
    # incomplete records causing a crash. None of the queries here are large, so it should be ok.
    db_session = sessionmaker(bind=engine, autoflush=False, expire_on_commit=False)
    session = db_session()

    return session


def create_db(db_address):
    # Create an engine for the database
    engine = create_engine(db_address)

    # Create all tables in the engine
    # This only needs to be run once on initialization.
    Base.metadata.create_all(engine)


def make_gaia_db(db_address, catalogs):
    create_db(db_address)
    session = get_session(db_address)
    for catalog_file in catalogs:
        logger.info('Adding {filename} to db.'.format(filename=catalog_file))
        catalog_data = parse_catalog(catalog_file)

        for chunk in range(0, len(catalog_data), 100000):
            logger.info('Adding records starting at {chunk} for {filename}'.format(filename=catalog_file, chunk=chunk))
            session.bulk_insert_mappings(Star, [{'ra': row['ra'], 'dec': row['dec'],
                                                 'ra_error': row['ra_error'],
                                                 'dec_error': row['dec_error'],
                                                 'g_flux': row['phot_g_mean_flux'],
                                                 'g_flux_error': row['phot_g_mean_flux_error']}
                                                for row in catalog_data[chunk: chunk + 100000]])
            session.commit()
    session.close()
    ra_index = Index('ra_index', Star.ra)
    engine = create_engine(db_address, poolclass=pool.NullPool)
    Base.metadata.bind = engine
    ra_index.create(bind=engine)
    dec_index = Index('dec_index', Star.dec)
    dec_index.create(bind=engine)


def make_gaia_healpix_catalogs(healpixels, catalogs, ncpu=6):
    catalog_names = []
    # For each heal pix
    for healpixel in healpixels:
        if os.path.exists(get_healpix_catalog_name(healpix_id=healpixel['index'],
                                                   nside=healpixel['nside'], allsky=False)):
            catalog_names.append(get_healpix_catalog_name(healpix_id=healpixel['index'],
                                                          nside=healpixel['nside'], allsky=False))
        else:
            # For each catalog fits file get the sources in this heal pixel
            # (Use multiple processes for this)
            p = mp.Pool(ncpu)
            sources_parameters = [(filename, healpixel) for filename in catalogs]
            stars_in_healpixel = p.map(get_sources_in_healpixel, sources_parameters)
            p.close()
            stars_in_healpixel = flatten_catalog_arrays(stars_in_healpixel)

            # Sort the healpix catalog descending by flux
            sorted_indices = np.argsort(stars_in_healpixel['phot_g_mean_flux'])[::-1]
            # Write out the healpix catalog to a fits file
            catalog_names.append(write_out_healpix_catalog(stars_in_healpixel[sorted_indices],
                                                           healpixel['index'], healpixel['nside']))

    return catalog_names


def parse_catalog(filename):
    catalog_data = np.array(fits.getdata(filename, 1))
    return catalog_data


def get_sources_in_healpixel(args):
    filename, healpixel = args

    logger.info('Getting sources in healpixel: {hp_id} from file: {filename}'.format(hp_id=healpixel['index'],
                                                                                     filename=filename))
    # Open the file
    catalog_data = parse_catalog(filename)
    offsets = great_circle_distance(catalog_data['ra'], catalog_data['dec'], healpixel['ra'], healpixel['dec'])
    # Return the sources within the distance threshold
    # We could also now run the catalog through hpsplit and only take the largest catalog
    sources_in_healpix = catalog_data[offsets <= healpixel['radius']].copy()
    del catalog_data
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
            offsets = great_circle_distance(catalog[i+1:i+n_healpix]['ra'],
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
                  '-A ra -D dec -I {unique_id} -j 0.1 &> {log_name}'
                  ''.format(catalog=catalog, index_name=index_name,
                            scale=scale, unique_id=make_unique_id(scale),
                            log_name=index_name.replace('.fits', '.log')))


def make_individual_index_files(catalogs, scale, margin=1, ncpu=6):
    p = mp.Pool(ncpu)
    p.map(make_single_index_file, [(catalog, scale, margin) for catalog in catalogs])
    p.close()


def make_single_index_file(args):
    catalog, scale, margin = args
    nside, healpix_id = catalog[:-5].split('-')[2:4]
    scale_name_str = scale_to_filename_string(scale)
    logger.info('Making index file. nside: {nside}, healpix: {hp_id}, scale: {scale} '.format(nside=nside,
                                                                                              hp_id=healpix_id,
                                                                                              scale=scale))
    index_name = 'gaia-index-{scale}-{nside}-{hp_id}.fits'.format(nside=nside,
                                                                  hp_id=healpix_id,
                                                                  scale=scale_name_str)

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


def flatten_catalog_arrays(list_of_catalogs):
    # Find the total number of stars in the catalog
    n_total_sources = np.sum([len(catalog) for catalog in list_of_catalogs])
    # Create an array of corresponding size
    merged_catalog = np.zeros(n_total_sources, dtype=list_of_catalogs[0].dtype)
    # copy the results from each individual array to the new array
    counter = 0
    for catalog in list_of_catalogs:
        merged_catalog[counter:counter + len(catalog)] = catalog
        counter += len(catalog)
    return merged_catalog


def write_out_healpix_catalog(catalog, healpix_id=0, nside=0, allsky=False):
    table_hdu = fits.BinTableHDU(catalog)
    catalog_name = get_healpix_catalog_name(healpix_id=healpix_id, nside=nside, allsky=allsky)
    table_hdu.writeto(catalog_name, clobber=True)
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


def plot_source_density():
    source_density = np.zeros((1000, 1000))
    merged_catalogs = glob('gaia-fullcatalog*.fits')
    for catalog in merged_catalogs:
        print(catalog)
        # Read in the file
        data = fits.getdata(catalog, 1)
        # Convert the decs to thetas: theta = 90 - dec
        theta = 90 - data['dec']
        # calculate z: z = cos(theta)
        z = np.cos(np.deg2rad(theta))
        # if abs(z) < 2/3:
        # xs = ra (in radians)
        xs = np.deg2rad(data['ra'].copy())
        # ys = (3 pi / 8) z
        ys = 3.0 * np.pi / 8.0 * z
        # else:
        # sigma = 2 - (3 (1 - abs(z)))^1/2
        sigma = 2 - (3 * (1.0 - np.abs(z)))**0.5
        # sigma(-z) = -sigma(z)
        sigma[z < 0] = -sigma[z < 0]
        # phi_l = ra % pi/2
        phi_l = np.deg2rad(data['ra']) % (np.pi / 2.0)
        # xs = ra - (abs(sigma(z)) - 1)(phi_l - pi / 4)
        polar_cap = np.abs(z) > (2.0 / 3.0)
        xs[polar_cap] = np.deg2rad(data['ra'][polar_cap]) - (np.abs(sigma[polar_cap]) - 1.0)*(phi_l[polar_cap] - (np.pi / 4.0))
        # ys = pi / 4 * sigma(z)
        ys[polar_cap] = np.pi / 4.0 * sigma[polar_cap]
        source_density += np.histogram2d(xs, ys, bins=(1000, 1000),
                                         range=[[0.0, 2.0 * np.pi], [-np.pi / 2.0, np.pi / 2.0]])[0]
    return source_density


if __name__ == "__main__":
    run()
