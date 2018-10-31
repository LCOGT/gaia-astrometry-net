"""
gaia-astrometry-net
Create astrometry.net index files from the GAIA catalogs.

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0
October 2018
"""
from setuptools import setup, find_packages


setup(name='gaia-astrometry-net',
      author=['Curtis McCully'],
      author_email=['cmccully@lco.global'],
      version='0.0.1',
      packages=find_packages(),
      package_dir={'gaia_astrometry_index_files': 'gaia_astrometry_index_files'},
      include_package_data=True,
      install_requires=['numpy', 'sqlalchemy', 'astrometry', 'astropy', 'lcogt_logging'],
      entry_points={'console_scripts': ['create_gaia_index_files=gaia_astrometry_index_files.main:create_index_files',
                                        'solve-astrometry-gaia=gaia_astrometry_index_files.solve:solve_frame']})
