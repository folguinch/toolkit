import os

from astropy.table import Table as aspyTable
from astropy.table import vstack, join
from astropy.io.votable import from_table
from astropy.io.votable.tree import VOTableFile, Resource

from .logger import get_logger

def load_table(tabname=None, table_id=None, append=True,
        logger=get_logger(__name__)):
    """Load an xml table.

    Parameters:
        tabname (str, optional): table file name to load.
        table_id (str, optional): table ID. Required if tabname given.
        append (bool, default=True): true if new data are appended.
    """
    if table_id:
        meta = {'ID':table_id, 'name':table_id+'.xml'}
    try:
        if append:
            logger.debug('Trying to open %s', tabname)
            table = aspyTable.read(tabname, table_id=table_id)
            logger.info('Table loaded: %s', os.path.basename(tabname))
        else:
            logger.info('Creating new table')
            table = aspyTable(meta=meta)
    except IOError:
        logger.info('Creating new table')
        table = aspyTable(meta=meta)

    return table

def update_table(table, table_new, table_id):
    try:
        table_up = aspyTable([table_new], names=table_new.keys(),
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

class Table:
    """Define a table object.

    Attributes:
        name (str): file name.
        table_id (str): table ID within the file.
        data (astropy.Table): table.
        logger (logging.logger): system logger.
    """

    def __init__(self, filename=None, table_id=None, append=True):
        """Initialize a new table object.

        Parameters:
            filename (str, default=None): data filename.
            table_id (str, default=None): table ID.
            append (bool, default=True): append data if file exist.
        """
        self.logger = get_logger(__name__)

        if filename:
            self.name = os.path.expanduser(filename)
        else:
            self.name = None

        if table_id:
            self.table_id = table_id
        elif filename:
            self.table_id = os.path.splitext(os.path.basename(self.name))[0]
            self.logger.warn('Using table ID from file name: %s', self.table_id)
        else:
            self.logger.warn('Using deafult table ID name')
            self.table_id = 'table'

        self.data = load_table(tabname=self.name, table_id=self.table_id,
                append=append, logger=self.logger)

    def add_row(self, row):
        """Insert a new row to the table.

        Parameters:
            row (OrderedDict): new row.
        """
        self.data = update_table(self.data, row, self.table_id)

    def write(self, tablename=None):
        """Save table to disk.

        Parameters:
            tablename (str, optional): update the table file name.
        """
        if tablename:
            self.name = tablename
        save_table(self.data, self.name)

