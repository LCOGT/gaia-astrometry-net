from glob import glob
from astropy.io import fits
from sqlalchemy import create_engine, pool
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, Integer, Float
from sqlalchemy.ext.declarative import declarative_base
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

Base = declarative_base()

logger = logging.getLogger('GaiaCatalogs')
logger.setLevel(logging.DEBUG)
logging.StreamHandler(sys.stdout)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = LCOGTFormatter()
handler.setFormatter(formatter)
logger.addHandler(handler)


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

    # New index files! Woot!


def save_catalog(catalog, counter):
    hdu = fits.BinTableHDU(catalog)
    output_name = 'gaia-fullcatalog-{counter}.fits'.format(counter=counter)
    hdu.writeto(output_name, clobber=True)
    return output_name


def nside_to_healpixels(nside):
    healpixel_indexes = np.arange(12 * nside**2, dtype=int)

    # Get the pixel positions of each healpix
    healpix_ras = [healpix_to_radec(int(i), nside, 0.5, 0.5)[0] * 180.0 / np.pi for i in healpixel_indexes]
    healpix_decs = [healpix_to_radec(int(i), nside, 0.5, 0.5)[1] * 180.0 / np.pi for i in healpixel_indexes]

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


class Position(Base):
    __tablename__ = 'positions'
    id = Column(Integer, primary_key=True, autoincrement=True)
    ramin = Column(Float)
    ramax = Column(Float)
    decmin = Column(Float)
    decmax = Column(Float)


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


def parse_row(line):
    values = line.split(',')
    return values[5], values[6], values[7], values[8], values[47], values[48]


def make_gaia_db(catalogs):
    connection = sqlite3.connect('gaia.db')
    cursor = connection.cursor()
    cursor.execute('PRAGMA LOCKING_MODE=EXCLUSIVE;')
    cursor.execute('PRAGMA SYNCHRONOUS=OFF;')
    cursor.execute('PRAGMA journal_mode=memory;')
    cursor.execute("PRAGMA count_changes=OFF;")
    cursor.execute('PRAGMA TEMP_STORE=memory;')
    cursor.execute('PRAGMA auto_vacuum = 0;')
    cursor.execute('PRAGMA foreign_keys=OFF;')
    cursor.execute('begin transaction')
    cursor.execute('CREATE TABLE IF NOT EXISTS sources (id INTEGER NOT NULL PRIMARY KEY, ra FLOAT, dec FLOAT, ra_error FLOAT, dec_error FLOAT, g_flux FLOAT, g_flux_error FLOAT);')
    cursor.execute('COMMIT')
    for catalog_file in catalogs:
        logger.info('Adding {filename} to db.'.format(filename=catalog_file))
        cursor.execute('begin transaction')
        with gzip.open(catalog_file, 'rt') as csv_file:
            columns = csv_file.readline()
            logger.info('Making rows variable')
            rows = [parse_row(line) for line in csv_file]
        logger.info('Adding records to db')
        cursor.executemany('INSERT INTO sources (ra, ra_error, dec, dec_error, g_flux, g_flux_error) values(?, ?, ?, ?, ?, ?)', rows)
        cursor.execute('COMMIT')
    cursor.execute('CREATE VIRTUAL TABLE if not exists positions using rtree(id, ramin, ramax, decmin, decmax);')
    cursor.execute('insert into positions (id, ramin, ramax, decmin, decmax) select id, ra, ra, dec, dec from sources;')
    cursor.close()


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

        # Sort the healpix catalog descending by flux
        stars_in_healpixel.sort('g_flux').reverse()
        # Write out the healpix catalog to a fits file
        return write_out_healpix_catalog(stars_in_healpixel, healpixel['index'], healpixel['nside'])


def make_gaia_healpix_catalogs(healpixels, db_address='sqlite:////gaia.db', ncpu=6):
    p = mp.Pool(ncpu)
    catalog_names = p.map(make_individual_healpix_catalog, healpixels)
    p.close()
    return catalog_names


def parse_catalog(filename):
    with gzip.open(filename, 'rt') as csv_file:
        csv_data = csv_file.read()
    catalog_data = ascii.read(csv_data, 'fast_csv')
    return catalog_data


def query_sources(db_address, ramin, ramax, decmin, decmax):
    db_session = get_session(db_address)
    sources_in_healpix = db_session.query(Position, Star).filter(Position.ramin >= ramin).filter(Position.ramax <= ramax)
    sources_in_healpix = sources_in_healpix.filter(Position.decmin >= decmin).filter(Position.decmax <= decmax)
    sources_in_healpix = sources_in_healpix.all()
    db_session.close()
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
        sources_in_healpix = query_sources(db_address, ramin, ramax, decmin, decmax)
    elif ramin < 0.0:
        sources_in_healpix = query_sources(db_address, ramin, ramax, decmin, decmax)
        ramin = 360.0 - delta_ra
        ramax = 360.0
        sources_in_healpix += query_sources(db_address, ramin, ramax, decmin, decmax)
    elif ramax > 360.0:
        sources_in_healpix = query_sources(db_address, ramin, ramax, decmin, decmax)
        ramin = 0.0
        ramax = delta_ra
        sources_in_healpix += query_sources(db_address, ramin, ramax, decmin, decmax)
    else:
        sources_in_healpix = query_sources(db_address, ramin, ramax, decmin, decmax)

    # Return the sources within the distance threshold
    source_catalog = {'flux': [source.g_flux for source, position in sources_in_healpix],
                      'ra': [source.ra for source, position in sources_in_healpix],
                      'dec': [source.dec for source, position in sources_in_healpix]}
    return source_catalog


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
                  '-A ra -D dec -I {unique_id} -j 0.1 -M &> {log_name}'
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
