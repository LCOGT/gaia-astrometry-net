from glob import glob

import numpy as np
from astropy.io import fits


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
