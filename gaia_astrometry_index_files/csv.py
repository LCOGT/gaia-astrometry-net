import gzip
import logging
logger = logging.getLogger(__name__)


def parse_csv(catalog_file):
    with gzip.open(catalog_file, 'rt') as csv_file:
        columns = csv_file.readline().split(',')
        indices = [columns.index(i) for i in ['ra', 'ra_error', 'dec', 'dec_error', 'phot_g_mean_flux', 'phot_g_mean_flux_error']]
        logger.info('Making rows variable')
        rows = [parse_row(line, indices) for line in csv_file]
    return rows


def parse_row(line, indicies):
    # This only works for GAIA DR2. We could get the headers from
    values = line.split(',')
    return (values[index] for index in indicies)
