import sqlite3 as sq

class myDatabase(object):
    """Object to manage databases.

    Attributes:
        db (sqlite3 database): database
        cursor (sqlite3 cursor): database cursor
        fmt (str): parameter format
    """

    def __init__(self, filename):
        self.db = sq.connect(filename)
        self.cursor = self.db.cursor()
        #self.fmt = []

    def commit(self):
        """Commit database"""
        self.db.commit()

    def execute(self, cmd):
        """Execute operation"""
        self.cursor.execute(cmd)

    def close(self):
        """Close database"""
        self.commit()
        self.db.close()

    def create_table(self, table, keys, fmts):
        """Create table
        
        Parameters:
            table (str): table name
            keys (list): column names
            fmts (list): data types
        """
        cols = ['%s %s,' % (key,fmt) for key, fmt in zip(keys, fmts)]
        cols = ' '.join(cols)
        cmd = 'CREATE TABLE IF NOT EXIST %s (%s)' % (table, cols)
        self.execute(cmd)

        #self.fmt = dbfmt_to_strfmt(fmts)

    def update(self, table, values):
        """Insert new database entries
        
        Parameters:
            table (str): table name
            values (list): values to insert
        """
        #vals = ', '.join(self.fmt)
        #vals = vals % values
        #self.execute.execute("INSERT INTO %s VALUES (%s)" % (table, vals))
        cmd = ','.join(['?'] * len(values))
        cmd = "INSERT INTO %s VALUES (%s)" % (table, cmd)
        self.execute(cmd, values)
        self.commit()

    def query(self, table, idcol, condition, values):
        """Select column satisfying all the conditions

        Parameters:
            table (str): table name
            idcol (str): column to select or coma separated columns
            condition (str): condition(s) to satisfy
            values (dict): key and values to fill the condition
        """
        cmd = 'SELECT %s FROM %s WHERE %s' % (idcol, table, condition)
        self.execute(cmd, values)

def load_database(filename, table, keys, fmts):
    """Load a database and create table.

    Parameters:
        filename (str): database file name
        table (str): table name
        keys (list): column names
        fmts (list): value formats
    """
    db = myDatabase(filename)
    db.create_table(table, keys, fmts)

    return db

def dbfmt_to_strfmt(fmts):
    """Convert database format to string format"""
    newfmts = []

    for fmt in fmts:
        if 'TEXT' in fmt:
            newfmts += ["'%s'"]
        elif 'REAL' in fmt:
            newfmts += ['%f']
        elif 'INTEGER' in fmt:
            newfmt += ['%i']
        else:
            raise NotImplementedError('Unrecognized format: %s' % fmt)

    return newfmts

