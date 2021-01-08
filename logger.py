from datetime import datetime
from typing import Optional, Union
import logging
import logging.handlers as handlers
import os
import pathlib

def get_level():
    numeric_level = getattr(logging, args.loglevel[0].upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.loglevel[0])
    return numeric_level

def get_logger(name: str, 
               filename: Union[str, pathlib.Path] = 'debug.log', 
               verbose: Optional[str] = None, 
               timestamp: bool = False,
               stdoutlevel: int = logging.INFO, 
               filelevel: int = logging.DEBUG,
               maxBytes: int = 5242880, 
               backupCount: int = 5) -> logging.Logger:
    """Creates a new logger.

    Verbose levels are:
      - v: basic logging with INFO level for stdout and DEBUG level for file
        logging. Equivalent to verbose=None.
      - vv: looging with DEBUG level for stdout and file logging.
      - vvv: same as verbose vv but add a timestamp to the file name and add
        time to the stdout log.
    Additional appereances of the character v will be ignored.

    Args:
      name: name of the logger.
      filename: optional; file name of the log.
      verbose: optional; verbose level.
      timestamp: optional; add timestamp to log file name.
      stdoutlevel: optional; logging level fot std output logging.
      filelevel: optional; logging level for file logging.
      maxBytes: optional; maximum size of logging file in bytes.
      backupCount: optional; maximum number of log files to rotate.
    """
    # Filter verbose
    if verbose is None:
        pass
    elif verbose == 'v':
        stdoutlevel = logging.INFO
        filelevel = logging.DEBUG
    elif verbose.startswith('vv'):
        stdoutlevel = filelevel = logging.DEBUG
        if verbose.startswith('vvv'):
            timestamp = True
    else:
        pass

    # Create logger
    logger = logging.getLogger(name)
    if not len(logger.handlers):
        logger.setLevel(filelevel)

        # File handler
        filename = pathlib.Path(filename).expanduser().resolve()
        if timestamp:
            tstamp = '.' + datetime.now().isformat(timespec='milliseconds')
            filename = filename.with_suffix(tstamp + filename.suffix)
        filefmt = ('%(asctime)s [%(levelname)s] - %(filename)s '
                   '(%(funcName)s:%(lineno)s): %(message)s')
        fh = handlers.RotatingFileHandler(filename, maxBytes=maxBytes, 
                                          backupCount=backupCount)
        fh.setLevel(filelevel)
        fh.setFormatter(logging.Formatter(filefmt))

        # Stream handler
        sh = logging.StreamHandler()
        sh.setLevel(stdoutlevel)
        if stdoutlevel == logging.DEBUG:
            streamfmt = '- %(filename)s (%(funcName)s:%(lineno)s): %(message)s'
            if timestamp:
                streamfmt = '%(asctime)s [%(levelname)s] ' + streamfmt
            else:
                streamfmt = '%(levelname)s ' + streamfmt
        else:
            streamfmt = '%(levelname)s: %(message)s'
        sh.setFormatter(logging.Formatter(streamfmt))

        # Register handlers
        logger.addHandler(fh)
        logger.addHandler(sh)

    return logger

