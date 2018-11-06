import numpy as np
from astrometry.util.util import healpix_to_radec


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


def test_great_circle_distance():
    np.testing.assert_allclose(great_circle_distance(0, 0, 100, 0), 100.0, atol=1e-5)
    np.testing.assert_allclose(great_circle_distance(0, 0, 181, 0), 179.0, atol=1e-5)
    np.testing.assert_allclose(great_circle_distance(0, 0, 0, 10), 10, atol=1e-5)
    np.testing.assert_allclose(great_circle_distance(0, 75, 180, 75), 30.0, atol=1e-5)


def great_circle_distance(ra1, dec1, ra2, dec2):
    ra1_rad, dec1_rad = np.deg2rad([ra1, dec1])
    ra2_rad, dec2_rad = np.deg2rad([ra2, dec2])
    distance_rad = np.arccos(np.sin(dec1_rad) * np.sin(dec2_rad) + np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra2_rad - ra1_rad))
    return np.rad2deg(distance_rad)
