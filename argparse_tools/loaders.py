from typing import Optional

from configparseradv import configparser
from spectral_cube import SpectralCube

def load_spectral_cube(args, cubename: Optional[str] = None):
    """Read a spectral_cube object and store it in args.

    Args:
      args: NameSpace to store the cube.
      cubename: optional; spectral cube file name.
    """
    if cubename is None:
        cubename = args.cubename
    try:
        args.cube = SpectralCube.read(cubename)
    except ValueError:
        args.cube = SpectralCube.read(cubename[0])

def load_config(args, config: Optional[str] = None):
    """Read a configuration file and store it in args.

    Args:
      args: NameSpace to store the configparser.
      config: optional; configuration file name.
    """
    if config is None:
        config = args.config
    args.config = configparser.ConfigParserAdv()
    args.config.read(config)

