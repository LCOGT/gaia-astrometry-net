import os
import subprocess
import shlex
import argparse
import logging

from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units
import numpy as np
from gaia_astrometry_index_files.catalog import make_catalog
from gaia_astrometry_index_files.utils import great_circle_distance

logger = logging.getLogger(__name__)

ASTROMETRY_NET_COMMAND = 'solve-field --crpix-center --no-verify --no-tweak ' \
      ' --radius {radius} --ra {ra} --dec {dec} --guess-scale ' \
      '--scale-units arcsecperpix --scale-low {scale_low} --scale-high {scale_high} ' \
      '--no-plots -N none --no-remove-lines ' \
      '--code-tolerance 0.003 --pixel-error 1 -d 1-200 ' \
      '--solved none --match none --rdls {catalog_source_file} --wcs {wcs_name} --corr none --overwrite ' \
      '-X X -Y Y -s FLUX --width {nx} --height {ny} -b {config_name} {catalog_name}'


def get_relevant_index_files(ra, dec, radius, index_file_path):
    index_file_meta_data = ascii.read(os.path.join(index_file_path, 'gaia-dr2-index-files.dat'))
    offsets = great_circle_distance(index_file_meta_data['ra'], index_file_meta_data['dec'], ra, dec)
    relevant_healpixels = offsets <= (index_file_meta_data['radius'] + radius)
    return list(index_file_meta_data['filename'][relevant_healpixels])


def save_config(filename, index_files, index_file_path):
    lines= ['inparallel\n', 'cpulimit 300\n', 'add_path {path}\n'.format(index_file_path)]
    for index_file in index_files:
        lines.append('index {index}\n'.format(index=index_file))
    with open(filename, 'w') as f:
        f.writelines(lines)


def solve_frame():
    parser = argparse.ArgumentParser(description='Run astrometry.net on an LCO fits file using the GAIA DR2 Index files')
    parser.add_argument('--ra', default=None)
    parser.add_argument('--dec', default=None)
    parser.add_argument('--radius', default=2.0)
    parser.add_argument('--index-file-path', dest='index_file_path')
    parser.add_argument('--filename')
    args = parser.parse_args()
    data, header = fits.getdata(args.filename)
    if args.ra is None or args.dec is None:
        ra, dec = parse_ra_dec(header)

    if args.ra is not None:
        ra = args.ra
    if args.dec is not None:
        dec = args.dec

    pixel_scale = header['PIXSCALE']
    # Skip the image if we don't have some kind of initial RA and Dec guess
    if np.isnan(ra) or np.isnan(dec):
        logger.error('Skipping WCS solution. No initial pointing guess from header.')
        return

    basename = os.path.basename(args.filename)
    catalog_name = os.path.join(os.getcwd(), basename.replace('.fits', '.cat.fits'))
    try:
        catalog = make_catalog(data, header)
        catalog[:40].write(catalog_name)
    except:
        logger.error('Could not produce source catalog')
        return

    index_files = get_relevant_index_files(ra, dec, args.radius, args.index_file_path)
    config_file_name = 'astrometry.cfg'
    save_config(config_file_name, index_files, args.index_file_path)
    # Run astrometry.net
    wcs_name = os.path.join(os.getcwd(), basename.replace('.fits', '.wcs.fits'))
    catalog_source_name = basename.replace('.fits', '.rdls.fits')
    command = ASTROMETRY_NET_COMMAND.format(ra=ra, dec=dec, scale_low=0.9 * pixel_scale, radius=args.radius,
                                            scale_high=1.1 * pixel_scale, wcs_name=wcs_name,
                                            catalog_name=catalog_name, nx=data.shape[1], ny=data.shape[0],
                                            config_name=config_file_name, catalog_source_name=catalog_source_name)

    try:
        console_output = subprocess.check_output(shlex.split(command))
    except subprocess.CalledProcessError:
        logger.error('Astrometry.net threw an error.')
        return

    # Copy the WCS keywords into original image
    new_header = fits.getheader(wcs_name)

    header_keywords_to_update = ['CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2', 'CRVAL1',
                                 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    for keyword in header_keywords_to_update:
        header[keyword] = new_header[keyword]

    # Update the RA and Dec header keywords
    header['RA'], header['DEC'] = get_ra_dec_in_sexagesimal(header['CRVAL1'], header['CRVAL2'])

    # Add the RA and Dec values to the catalog
    add_ra_dec_to_catalog(header, catalog)

    save_solved_image(basename.replace('.fits', 'solved.fits'), header, data, catalog)
    save_astrometry_catalog_regions(basename.replace('.fits', '.gaia.reg'), catalog_source_name)


def save_astrometry_catalog_regions(filename, catalog_fits_file):
    hdu = fits.open(catalog_fits_file)
    lines = ['{ra} {dec}\n'.format(ra=row['RA'], dec=row['DEC']) for row in hdu[1].data]
    with open(filename, 'w') as f:
        f.writelines(lines)


def save_solved_image(filename, header, data, catalog):
    primary_hdu = fits.PrimaryHDU(data=data, header=header)
    catalog_hdu = fits.BinTableHDU(catalog, name='CAT')
    fits.HDUList([primary_hdu, catalog_hdu]).writeto(filename, overwrite=True)


def add_ra_dec_to_catalog(header, catalog):
    image_wcs = WCS(header)
    ras, decs = image_wcs.all_pix2world(catalog['x'], catalog['y'], 1)
    catalog['ra'] = ras
    catalog['dec'] = decs
    catalog['ra'].unit = 'degree'
    catalog['dec'].unit = 'degree'
    catalog['ra'].description = 'Right Ascension'
    catalog['dec'].description = 'Declination'


def get_ra_dec_in_sexagesimal(ra, dec):
    """
    Convert a decimal RA and Dec to sexagesimal

    Parameters
    ----------
    ra : float
         Right Ascension in decimal form
    dec : float
         Declination in decimal form

    Returns
    -------
    tuple of str : RA, Dec converted to a string

    """
    coord = SkyCoord(ra, dec, unit=(units.deg, units.deg))
    coord_str = coord.to_string('hmsdms', precision=4, pad=True)
    ra_str, dec_str = coord_str.split()
    ra_str = ra_str.replace('h', ':').replace('m', ':').replace('s', '')
    dec_str = dec_str.replace('d', ':').replace('m', ':').replace('s', '')
    # Return one less digit of precision for the dec
    dec_str = dec_str[:-1]
    return ra_str, dec_str


def parse_ra_dec(header):
    try:
        coord = SkyCoord(header.get('RA'), header.get('DEC'), unit=(units.hourangle, units.degree))
        ra = coord.ra.deg
        dec = coord.dec.deg
    except (ValueError, TypeError):
        # Fallback to CRVAL1 and CRVAL2
        try:
            coord = SkyCoord(header.get('CRVAl1'), header.get('CRVAL2'), unit=(units.degree, units.degree))
            ra = coord.ra.deg
            dec = coord.dec.deg
        except (ValueError, TypeError):
            # Fallback to Cat-RA and CAT-DEC
            try:
                coord = SkyCoord(header.get('CAT-RA'), header.get('CAT-DEC'), unit=(units.hourangle, units.degree))
                ra = coord.ra.deg
                dec = coord.dec.deg
            except (ValueError, TypeError) as e:
                logger.error('Could not get initial pointing guess. {0}'.format(e),
                             extra_tags={'filename': header.get('ORIGNAME')})
                ra, dec = np.nan, np.nan
    return ra, dec
