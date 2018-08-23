import os

"""Decorator classes and functions"""

DECORATED_CLASSES = {}
DECORATED_FUNCTIONS = {}

def register_class(cls):
    REGISTERED_CLASSES[cls.__name__.lower()] = cls
    return cls

def register_function(fn):
    REGISTERED_FUNCTIONS[fn.__name__.lower()] = fn
    return fn

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
