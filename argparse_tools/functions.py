"""Functions for helping the processing of arguments."""
from typing import Optional

def positions_to_pos(args: 'argparse.ArgumentParser',
                     wcs: Optional['astropy.WCS'] = None) -> None:
    """Store the positions in the `args` object.

    Args:
      args: argument parser.
      wcs: optional; WCS to convert coordinates to pixels.
    """
    # Get the positions
    if args.position:
        if len(args.position)%2 != 0:
            raise ValueError('Odd number of position values')
        for xy in zip(args.position[::2],args.position[1::2]):
            args.pos += [xy]
    elif args.reference:
        args.pos = args.reference
    elif args.coordinate:
        args.pos = args.coordinate
    else:
        for src in args.source:
            args.pos += [src.position]

    # Convert SkyCoords to pixels if wcs
    if wcs is not None:
        try:
            args.pos = [[p.ra.degree, p.dec.degree] for p in args.pos]
            args.pos = wcs.all_world2pix(args.pos, 0)
        except AttributeError:
            pass
