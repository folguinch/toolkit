from spectral_cube import SpectralCube

from ..myconfigparser import myConfigParser

def load_spectral_cube(args, cubename=None):
    if cubename is None:
        cubename = args.cubename
    try:
        args.cube = SpectralCube.read(cubename)
    except ValueError:
        args.cube = SpectralCube.read(cubename[0])

def load_config(args, config=None):
    if config is None:
        config = args.config
    args.config = myConfigParser()
    args.config.read(config)

