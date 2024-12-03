"""Functions for helping the processing of arguments."""
from typing import Optional, Dict, Sequence

from astropy.wcs.utils import pixel_to_skycoord

def positions_to_pixels(args: 'argparse.Namespace',
                        wcs: Optional['astropy.wcs.WCS'] = None) -> None:
    """Store the positions in the `args` object as pixels.

    Args:
      args: Argument parser.
      wcs: Optional. WCS to convert coordinates to pixels.
    """
    # Get the positions
    pos = []
    if args.position:
        if len(args.position)%2 != 0:
            raise ValueError('Odd number of position values')
        for xy in zip(args.position[::2],args.position[1::2]):
            pos.append(xy)
    elif args.reference:
        pos = args.reference
    elif args.coordinate:
        pos = args.coordinate
    else:
        try:
            if args.source:
                for src in args.source:
                    pos.append(src.position)
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
    #try:
    #    # This should be deleted in the future
    #    args.pos = pos
    #except AttributeError:
    #    pass
    args.position = args.coordinate = args.reference = pos

def pixels_to_positions(args: 'argparse.Namespace',
                        wcs: Optional['astropy.wcs.WCS'] = None) -> None:
    """Store the positions in `args` as `SkyCoords`.

    Args:
      args: Argument parser.
      wcs: Optional. WCS to convert coordinates to pixels.
    """
    # Get the positions
    pos = []
    if args.position:
        if len(args.position)%2 != 0:
            raise ValueError('Odd number of position values')
        for x, y in zip(args.position[::2], args.position[1::2]):
            pos.append(pixel_to_skycoord(x, y, wcs))
    elif args.reference:
        pos = args.reference
    elif args.coordinate:
        pos = args.coordinate
    else:
        try:
            if args.source:
                for src in args.source:
                    pos.append(src.position)
            else:
                pos = None
        except ValueError:
            pos.append(args.source.position)
        except AttributeError:
            pos = None

    # Store the value
    args.position = args.coordinate = args.reference = pos

def source_properties(args: 'argparse.Namespace',
                      properties: Sequence) -> Dict:
    """Read source properties.

    The source `properties` can be specified through command line arguments or
    `astro_source.source.Source` configuration. Command line options take
    precedence.

    Args:
      args: Command line arguments.
      properties: Properties to read.

    Returns:
      A dictionary with the property values.
    """
    args_dict = vars(args)
    try:
        source = args.source[0]
    except AttributeError:
        source = None
    props_values = {}
    for prop in properties:
        if prop in args_dict and (val := args_dict[prop]) is not None:
            props_values[prop] = val
        elif source is not None:
            props_values[prop] = source.config.getquantity('INFO', prop,
                                                           fallback=None)
        else:
            props_values[prop] = None

    return props_values
