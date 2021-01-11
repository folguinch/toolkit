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

def get_stdout_format(level: int = logging.INFO, 
                      timestamp: bool = False) -> str:
    """Determines the stdout format.

    Args:
      level: optional; logging level.
      timestamp: optional; include timestamp?
    """
    if level == logging.DEBUG:
        fmt = '- %(filename)s (%(funcName)s:%(lineno)s): %(message)s'
        if timestamp:
            fmt = '%(asctime)s [%(levelname)s] ' + fmt
        else:
            fmt = '%(levelname)s ' + fmt
    elif level == logging.INFO:
        fmt = '%(levelname)s: %(message)s'
    else:
        raise NotImplementedError(f'fmt for level {level} not implemented')

    return fmt


def get_stdout_handler(level: int, timestamp: bool = False, 
                       fmt: Optional[str] = None) -> logging.StreamHandler:
    """Create an standard output handler with the given log level.

    Args:
      level: logging level.
      timestamp: optional; add time stamp to default format.
      fmt: optional; text format.
    """
    sh = logging.StreamHandler()
    sh.setLevel(level)
    if fmt is None:
        fmt = get_stdout_format(level=level, timestamp=timestamp)
    sh.setFormatter(logging.Formatter(fmt))

    return sh

def get_file_handler(level: int, filename: Union[str, pathlib.Path], 
                     timestamp: bool = False, maxBytes: int = 5242880, 
                     backupCount: int = 5, 
                     fmt: Optional = None) -> handlers.RotatingFileHandler:
    """Create a file logging handler.

    Args:
      level: logging level.
      filename: logging file name
    """
    filename = pathlib.Path(filename).expanduser().resolve()
    if timestamp:
        tstamp = '.' + datetime.now().isoformat(timespec='milliseconds')
        filename = filename.with_suffix(tstamp + filename.suffix)
    if fmt is None:
        fmt = ('%(asctime)s [%(levelname)s] - %(filename)s '
               '(%(funcName)s:%(lineno)s): %(message)s')
    fh = handlers.RotatingFileHandler(filename, maxBytes=maxBytes,
                                      backupCount=backupCount)
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter(fmt))

    return fh

def _levels_from_verbose(verbose: str, 
                         stdoutlevel: int = logging.INFO, 
                         filelevel: int = logging.DEBUG, 
                         timestamp: bool = False) -> dict:
    """Return logging levels depending on the verbose value.

    Accepted verbose values are:
      - None: use the defaults in stdoutlevel and filelevel.
      - v: basic logging with INFO level for stdout and DEBUG level for file
        logging. 
      - vv: looging with DEBUG level for stdout and file logging.
      - vvv: same as verbose vv but add a timestamp to the file name and add
        time to the stdout log.
    Additional appereances of the character v will be ignored.

    Args:
      verbose: verbose level.
      stdoutlevel: optional; default value for stdout level.
      filelevel: optional; default value for file level.
      timestamp: optional; include a log timestamp.
    Returns:
      Dictionary with the values of stdoutlevel and filelevel
    """
    levels = {'stdout': stdoutlevel, 'file': filelevel, 'timestamp': timestamp}
    if verbose is None:
        pass
    elif verbose == 'v':
        levels['stdout'] = logging.INFO
        levels['file'] = logging.DEBUG
    elif verbose.startswith('vv'):
        levels['stdout'] = levels['file'] = logging.DEBUG
        if verbose.startswith('vvv'):
            levels['timestamp'] = True
    else:
        pass

    return levels

def get_logger(name: str, 
               filename: Union[None, str, pathlib.Path] = None,
               verbose: Optional[str] = None, 
               timestamp: bool = False,
               stdoutlevel: int = logging.INFO, 
               filelevel: int = logging.DEBUG,
               maxBytes: int = 5242880, 
               backupCount: int = 5) -> logging.Logger:
    """Creates a new logger.

    Verbose levels overwrite the other parameter values. Accepted verbose
    values are:
      - None: use the other keyword arguments.
      - v: basic logging with INFO level for stdout and DEBUG level for file
        logging. 
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
    levels = _levels_from_verbose(verbose, stdoutlevel=stdoutlevel,
                                  filelevel=filelevel, timestamp=timestamp)

    # Create logger
    logger = logging.getLogger(name)
    if not len(logger.handlers):
        logger.setLevel(levels['file'])

        # File handler
        if filename is not None:
            fh = get_file_handler(levels['file'], filename, 
                                  timestamp=levels['timestamp'], 
                                  maxBytes=maxBytes, 
                                  backupCount=backupCount)
        else:
            fh = None

        # Stream handler
        sh = get_stdout_handler(levels['stdout'], 
                                timestamp=levels['timestamp'])

        # Register handlers
        if fh is not None:
            logger.addHandler(fh)
        logger.addHandler(sh)

    return logger

def get_stdout_logger(name: str, 
                      verbose: Optional[str] = None, 
                      timestamp: bool = False,
                      level: int = logging.INFO) -> logging.Logger:
    """Creates a new standard output logger.

    Args:
      name: name of the logger.
      verbose: optional; verbose level.
      timestamp: optional; add timestamp to log file name.
      level: optional; logging level fot std output logging.
    """
    return get_logger(name, filename=None, verbose=verbose,
                      timestamp=timestamp, stdoutlevel=level)

def update_logger(logger: logging.Logger, 
                  filename: Union[None, str, pathlib.Path] = None,
                  verbose: Optional[str] = None,
                  timestamp: bool = False,
                  stdoutlevel: int = logging.INFO, 
                  filelevel: int = logging.DEBUG,
                  maxBytes: int = 5242880, 
                  backupCount: int = 5) -> logging.Logger:
    """Update logger handlers.

    Args:
      logger: logger instance.
      filename: optional; file name for new file handler.
      verbose: optional; update logging levels.
      timestamp: optional; default timestamp values for new/updated handlers.
      stdoutlevel: optional; change standard output level.
      filelevel: optional; change file level.
      maxBytes: optional; maximum size of logging file in bytes.
      backupCount: optional; maximum number of log files to rotate.
    Returns:
      An updated logger instance.
    """
    # Updated levels
    levels = _levels_from_verbose(verbose, stdoutlevel=stdoutlevel,
                                  filelevel=filelevel, timestamp=timestamp)
    # Loop through the handlers
    nfilehandlers = 0
    for handler in logger.handlers:
        if hasattr(handler, 'baseFilename'):
            # Update level
            handler.setLevel(levels['file'])
            nfilehandlers += 1
        else:
            # Update level
            handler.setLevel(levels['stdout'])

            # Update format
            fmt = get_stdout_format(level=levels['stdout'], 
                                    timestamp=levels['timestamp'])
            handler.setFormatter(logging.Formatter(fmt))

    # Create file handler if needed
    if nfilehandlers == 0 and filename is not None:
        fh = get_file_handler(levels['file'], filename, 
                              timestamp=levels['timestamp'], 
                              maxBytes=maxBytes, 
                              backupCount=backupCount)
        logger.addHandler(fh)

    return logger

