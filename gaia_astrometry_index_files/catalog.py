import sep
import numpy as np
from astropy.table import Table

# Note that threshold is number of sigma, not an absolute number because we provide the error
# array to SEP.
THRESHOLD = 10.0
MIN_AREA = 9


def prune_nans_from_table(table):
    nan_in_row = np.zeros(len(table), dtype=bool)
    for col in table.colnames:
        nan_in_row |= np.isnan(table[col])
    return table[~nan_in_row]


def save_catalog_meta_data(sources):
    catalog = sources['x', 'y', 'xpeak', 'ypeak',
                      'flux', 'fluxerr', 'peak', 'fwhm',
                      'a', 'b', 'theta', 'kronrad', 'ellipticity',
                      'x2', 'y2', 'xy', 'flag']

    # Add the units and description to the catalogs
    catalog['x'].unit = 'pixel'
    catalog['x'].description = 'X coordinate of the object'
    catalog['y'].unit = 'pixel'
    catalog['y'].description = 'Y coordinate of the object'
    catalog['xpeak'].unit = 'pixel'
    catalog['xpeak'].description = 'X coordinate of the peak'
    catalog['ypeak'].unit = 'pixel'
    catalog['ypeak'].description = 'Windowed Y coordinate of the peak'
    catalog['flux'].unit = 'count'
    catalog['flux'].description = 'Flux within a Kron-like elliptical aperture'
    catalog['fluxerr'].unit = 'count'
    catalog['fluxerr'].description = 'Error on the flux within Kron aperture'
    catalog['peak'].unit = 'count'
    catalog['peak'].description = 'Peak flux (flux at xpeak, ypeak)'
    catalog['fwhm'].unit = 'pixel'
    catalog['fwhm'].description = 'FWHM of the object'
    catalog['a'].unit = 'pixel'
    catalog['a'].description = 'Semi-major axis of the object'
    catalog['b'].unit = 'pixel'
    catalog['b'].description = 'Semi-minor axis of the object'
    catalog['theta'].unit = 'degree'
    catalog['theta'].description = 'Position angle of the object'
    catalog['kronrad'].unit = 'pixel'
    catalog['kronrad'].description = 'Kron radius used for extraction'
    catalog['ellipticity'].description = 'Ellipticity'
    catalog['x2'].unit = 'pixel^2'
    catalog['x2'].description = 'Variance on X coordinate of the object'
    catalog['y2'].unit = 'pixel^2'
    catalog['y2'].description = 'Variance on Y coordinate of the object'
    catalog['xy'].unit = 'pixel^2'
    catalog['xy'].description = 'XY covariance of the object'
    catalog['flag'].description = 'Bit mask of extraction/photometry flags'

    catalog.sort('flux')
    catalog.reverse()
    return catalog


def make_catalog(data, header):
    # Set the number of source pixels to be 5% of the total. This keeps us safe from
    # satellites and airplanes.
    sep.set_extract_pixstack(int(data.shape[1] * data.shape[0] * 0.05))

    data = data.copy()
    error = (np.abs(data) + header['RDNOISE'] ** 2.0) ** 0.5
    mask = data > 0.9 * header['SATURATE']

    # Fits can be backwards byte order, so fix that if need be and subtract
    # the background
    try:
        bkg = sep.Background(data, mask=mask, bw=32, bh=32, fw=3, fh=3)
    except ValueError:
        data = data.byteswap(True).newbyteorder()
        bkg = sep.Background(data, mask=mask, bw=32, bh=32, fw=3, fh=3)
    bkg.subfrom(data)

    # Do an initial source detection
    sources = sep.extract(data, THRESHOLD, mask=mask, minarea=MIN_AREA, err=error, deblend_cont=0.005)

    # Convert the detections into a table
    sources = Table(sources)

    # We remove anything with a detection flag >= 8
    # This includes memory overflows and objects that are too close the edge
    sources = sources[sources['flag'] < 8]

    sources = prune_nans_from_table(sources)

    # Calculate the ellipticity
    sources['ellipticity'] = 1.0 - (sources['b'] / sources['a'])

    # Fix any value of theta that are invalid due to floating point rounding
    # -pi / 2 < theta < pi / 2
    sources['theta'][sources['theta'] > (np.pi / 2.0)] -= np.pi
    sources['theta'][sources['theta'] < (-np.pi / 2.0)] += np.pi

    # Calculate the kron radius
    kronrad, krflag = sep.kron_radius(data, sources['x'], sources['y'],
                                      sources['a'], sources['b'],
                                      sources['theta'], 6.0)
    sources['flag'] |= krflag
    sources['kronrad'] = kronrad

    # Calcuate the equivilent of flux_auto
    flux, fluxerr, flag = sep.sum_ellipse(data, sources['x'], sources['y'],
                                          sources['a'], sources['b'],
                                          np.pi / 2.0, 2.5 * kronrad,
                                          subpix=1, err=error)
    sources['flux'] = flux
    sources['fluxerr'] = fluxerr
    sources['flag'] |= flag

    # Calculate the FWHMs of the stars:
    fwhm = 2.0 * (np.log(2) * (sources['a'] ** 2.0 + sources['b'] ** 2.0)) ** 0.5
    sources['fwhm'] = fwhm

    # Cut individual bright pixels. Often cosmic rays
    sources = sources[fwhm > 1.0]

    # Update the catalog to match fits convention instead of python array convention
    sources['x'] += 1.0
    sources['y'] += 1.0

    sources['xpeak'] += 1
    sources['ypeak'] += 1

    sources['theta'] = np.degrees(sources['theta'])

    return save_catalog_meta_data(sources)
