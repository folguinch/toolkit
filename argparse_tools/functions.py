"""Functions for helping the processing of arguments."""
from typing import Optional

def positions_to_pixels(args: 'argparse.ArgumentParser',
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
            pos += [xy]
    elif args.reference:
        pos = args.reference
    elif args.coordinate:
        pos = args.coordinate
    else:
        try:
            if args.source:
                for src in args.source:
                    pos += [src.position]
            else:
                pos = None
        except AttributeError:
            pos = None

    # Convert SkyCoords to pixels if wcs
    if wcs is not None and pos is not None:
        try:
            pos = [[p.ra.degree, p.dec.degree] for p in pos]
            pos = wcs.all_world2pix(pos, 0)
        except AttributeError:
            pass

    # Store the value
    try:
        # This should be deleted in the future
        args.pos = pos
    except AttributeError:
        pass
    args.position = args.coordinate = args.reference = pos
