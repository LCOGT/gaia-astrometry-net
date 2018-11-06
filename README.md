# gaia-astrometry-net
*Code to generate astrometry.net index files from the GAIA DR2 catalog*

## Prerequisites
- numpy
- sqlalchemy
- astrometry
- astropy
- lcogt_logging
- sep
- pyyaml


## Installation
All of the dependencies below install automatically when you run python setup.py
install (or via pip) except *astrometry*. The *astrometry* package is not 
currently on PyPi and is provided directly by astrometry.net. This can be installed
by cloning the astrometry.net repo at https://github.com/dstndstn/astrometry.net
and running python setup.py install in the top directory.

Once the *astrometry* package is installed, the code can be installed in the 
usual way by running python setup.py install.

## Usage
There are two main purposes for the code in this repo. The first is to generate
astrometry.net index files from the GAIA DR2 catalog. The second is wrap 
astrometry.net to use the new index files to solve for the astrometry of an
image.

The index file creation logic is mostly in index_file_creator.py. To create the index files
go into an empty directory and run 
```
create_gaia_index_files
``` 
from the command line.
This will then use the network mounted of the GAIA DR2 csv files to create astrometry.net 
index files. The first issue is that csv files do not cover consistent patches of the sky.
As such, the first thing that we do is build an SQLite database with all of the positions 
and fluxes of the GAIA sources. This DB includes two tables: sources and positions.
The positions table is indexed using an r-tree index so that cone searches are
fast (even for the 1.5 billion GAIA sources). Generating the DB can take a few days
of cpu time to complete. The final DB will be ~200 GB. After we have the queryable 
catalog (the DB), we split the catalog into healpixels (see https://healpix.jpl.nasa.gov/). 
The sources in each healpixel are stored in a binary table extension of a fits file
that astrometry.net can read. Finally, we use astrometry.net's internal tools
to create index files that astrometry.net can use for solves. We create a range
of index files that have quads for many different image scales
(large field of view to small field of view images). The meta data for the index
files is stored in an ascii file called `gaia-dr2-index-files.dat`.

Once the index files are created, we need to wrap astrometry.net to use them. 
The code to do this is mostly in `solve.py`. 

## Examples
To solve an image using the GAIA astrometry.net index files, you can run
```
solve-astrometry-gaia --filename {filename} --index-file-path {path} --ra {ra} --dec {dec} --radius {radius}
```
`filename` is the name of the image to solve.
The `index-file-path` argument should point to the directory with the GAIA astrometry index files.
By default, the code will try to extra RA and DEC from the header, but the 
user can manually override those values using `--ra` and `--dec`. The default
search radius is 2 degrees but this can also be overridden by the user.

## Support
[API documentation]()  
[Create an issue](https://issues.lco.global/)
