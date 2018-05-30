import os

from astropy.table import Table, vstack, join
from astropy.io.votable import from_table
from astropy.io.votable.tree import VOTableFile, Resource

from .logger import get_logger

def load_table(tabname=None, table_id=None, append=True):
    """Load an xml table.

    Parameters:
        tabname (str, optional): table file name to load.
        table_id (str, optional): table ID. Required if tabname given.
        append (bool, default=True): true if new data are appended.
    """
    logger = get_logger(__name__)

    if table_id:
        meta = {'ID':table_id, 'name':table_id+'.xml'}
    try:
        if append:
            logger.debug('Trying to open %s', tabname)
            table = Table.read(tabname, table_id=table_id)
            logger.info('Table loaded: %s', os.path.basename(tabname))
        else:
            logger.info('Creating new table')
            table = Table(meta=meta)
    except IOError:
        logger.info('Creating new table')
        table = Table(meta=meta)

    return table

def update_table(table, table_new, table_id):
    try:
        table_up = Table([table_new], names=table_new.keys(),
                        meta={'ID':table_id, 'name':table_id+'.xml'})
    except ValueError:
        table_up = table_new

    try:
        table = vstack([table, table_up])
    except TypeError:
        table = table_up
    
    return table

def save_table(table, tabname):
    results = VOTableFile()
    resource = Resource()
    results.resources.append(resource)
    resource.tables.append(from_table(table).get_first_table())
    results.to_xml(tabname)


