"""Tool for working with regions."""
from typing import Optional, Sequence
import argparse
import sys

from astropy.coordinates import SkyCoord
from regions import Regions, PolygonSkyRegion
from shapely import LineString
import astropy.units as u

from ..argparse_tools import actions, parents

def _load_region(args: argparse.Namespace):
    """Load region from line args."""
    args.reg = Regions.read(str(args.region[0]), format=args.format[0]).pop()

def _line_to_poly(args:argparse.Namespace):
    """Convert a line region to a polygon by dilation."""
    poly_region = sky_line_to_poly(args.reg, args.dilation)

    # Save file
    filename = args.output[0]
    poly_region.write(str(filename), format=args.format[0])

def region_to_buffer(region: Regions, dilation: u.Quantity):
    """Convert a line region to a `shapely` polygon by dilation."""
    # Define shapely line
    vertices = list(zip(region.vertices.ra.deg, region.vertices.dec.deg))
    line = LineString(vertices)

    # Dilate and convert to region
    buff = line.buffer(dilation.to(u.deg).value)

    return buff

def sky_line_to_poly(region: Regions, dilation: u.Quantity,
                     frame: str = 'icrs'):
    """Convert a line region to a polygon by dilating around vertices.

    Args:
      region: region object.
      dilation: amount of dilation.
      frame: optional; coordinate system.
    """
    buff = region_to_buffer(region, dilation)
    poly = SkyCoord(list(buff.exterior.coords), unit='deg', frame=frame)
    poly_reg = PolygonSkyRegion(vertices=poly)

    return poly_reg

def transform_region(args: Optional[Sequence] = None):
    """Convert a region into another shape.

    Currently it converts a set of points (a line) to a polygon by dilation.

    Args:
      args: arguments.
    """
    # Steps
    pipe = [_load_region, _line_to_poly]
    # Setup
    args_parents = [
        parents.logger('debug_transform_region.log'),
        parents.verify_files('region',
                             region={'help': 'Region file', 'nargs': 1}),
    ]
    parser = argparse.ArgumentParser(
        description='Convert between compatible region shapes',
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=args_parents,
        conflict_handler='resolve',
    )
    parser.add_argument('--format', nargs=1, default=['crtf'],
                        help='Region format')
    parser.add_argument('output', action=actions.NormalizePath, nargs=1,
                        help='Output region')
    parser.add_argument('dilation', nargs=2, action=actions.ReadQuantity,
                        help='Dilation amount')
    parser.set_defaults(reg=None)

    # Read arguments
    if args is None:
        args = sys.argv[1:]
    args = parser.parse_args(args)

    # Process data
    for step in pipe:
        step(args)

if __name__ == '__main__':
    transform_region(sys.argv[1:])
