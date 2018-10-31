import gzip
import sqlite3
from sqlalchemy import Column, Integer, Float, create_engine, pool
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from gaia_astrometry_index_files.main import logger

Base = declarative_base()


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


def parse_row(line):
    values = line.split(',')
    return values[5], values[6], values[7], values[8], values[47], values[48]