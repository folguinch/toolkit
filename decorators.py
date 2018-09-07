import os
from functools import wraps
from datetime import datetime

from .logger import get_logger

"""Decorator classes and functions"""

REGISTERED_CLASSES = {}
REGISTERED_FUNCTIONS = {}

def register_class(cls):
    REGISTERED_CLASSES[cls.__name__.lower()] = cls
    return cls

def register_function(fn):
    REGISTERED_FUNCTIONS[fn.__name__.lower()] = fn
    return fn

def logged_class(name, filename):

    class ClassWrapper(object):
        __metaclass__=type

        def __init__(self, cls, name=name, filename=filename):
            self.cls = cls
            print name, filename
            self.logger = get_logger(name, filename)

        def __call__(self, *args, **kwargs):
            kwargs.setdefault('logger', self.logger)
            self.logger.debug('Initializing class')
            cls = self.cls(*args, **kwargs)
            return cls

    return ClassWrapper

def timed(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        ti = datetime.now()
        results = fn(*args, **kwargs)
        tf = datetime.now()
        logger = get_logger(__name__)
        logger.info('Function: %s executed in %s', fn.__name__, tf-ti)
        return results
    return wrapper

class checkPaths(object):

    """Normalize the files or directories paths, and file(s) existance if 
    requested"""

    def __init__(self, fn):
        """Initialize decorator"""
        self.fn = fn

    def __call__(self, *args, **kwargs):
        """Use os.path to normalize path of files/directories in args"""
        paths = []
        for arg in args:
            try:
                path = os.path.realpath(os.path.expanduser(arg))
                if os.path.isfile(path):
                    if kwargs.set_default('overwrite', False):
                        print('File %s already exists' % os.path.basename(path))
                        exit()
                    else:
                        paths += [path]
                else:
                    path += [path]
            except AttributeError:
                continue

        self.fn(*args, **kwargs)

