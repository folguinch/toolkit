from typing import Optional

from configparseradv import configparser
from spectral_cube import SpectralCube

def load_spectral_cube(args, cubename: Optional[str] = None, 
                       use_dask: bool = False) -> None:
    """Read a spectral_cube object and store it in args.

    Args:
      args: NameSpace to store the cube.
      cubename: optional; spectral cube file name.
      use_dask: optional; use dask?
    """
    if cubename is None:
        cubename = args.cubename
    try:
        args.cube = SpectralCube.read(cubename, use_dask=use_dask)
    except ValueError:
        args.cube = SpectralCube.read(cubename[0], use_dask=use_dask)

    if use_dask:
        args.cube.use_dask_scheduler('threads', num_workers=10)

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

