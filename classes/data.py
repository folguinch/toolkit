import os
from abc import ABCMeta, abstractmethod

from ..logger import get_logger

class Data(object):
    """Data ABC.

    Attributes:
        address: file name.
        data: the data in the file.
        logger: logging manager.
    """

    __metaclass__ = ABCMeta

    def __init__(self, address=None, data=None):
        """Defines a new data object.

        Parameters:
            address: filename
        """
        self.data = data
        self.logger = get_logger(__name__)

        try:
            self.address = os.path.realpath(os.path.expanduser(address))
        except AttributeError:
            self.address = address

        if self.address and os.path.isfile(self.address):
            self.logger.debug('Load file: %s', address)
            self.logger.info('Load file: %s', os.path.basename(address))
            self.load()
        else:
            self.logger.info('Nothing to load from: %s', self.address)

    @abstractmethod
    def load(self):
        """Open file in *address*."""
        pass

    @abstractmethod
    def save(self, file_name=None):
        """Saves the file in *address* or in a new address if provided.
        
        Parameters:
            file_name (default=None): new file name.
        """
        pass
