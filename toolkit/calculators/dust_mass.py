#!/bin/python3
import os, argparse
from itertools import product

from astropy.analytic_functions import blackbody_nu
#from astroSource.source import LoadSourcefromConfig
import astropy.units as u
import numpy as np
import toolkit.argparse_tools.actions as actions 

def dust_mass(nu: u.Quantity[u.GHz],
              flux: u.Quantity[u.Jy],
              temp: u.Quantity[u.K],
              distance: u.Quantity[u.kpc],
              kappa: u.Quantity[u.cm**2/u.g] = 1.0*u.cm**2/u.g,
              Rdg: float = 0.01) -> u.Quantity[u.M_sun]:
    """Calculate dust mass.
    
    Args:
      nu: Frequency.
      flux: Flux density.
      temp: Dust temperature.
      distance: Distance to the source.
      kappa: Optional. Dust opacity.
      Rdg: Optional. Dust-to-gas ratio.
    """
    nu = nu.to(u.GHz, equivalencies=u.spectral())
    d = distance.to(u.cm)
    bnu = blackbody_nu(nu, temp).to(u.W / u.m**2 / u.Hz,
                                    u.dimensionless_angles())
    bnu = bnu.cgs.decompose()
    mass = flux.decompose() * d**2 / (Rdg * kappa * bnu)

    return mass.to(u.M_sun)

def preprocess(args):
    # Put units
    args.rdg = args.rdg[0]
    args.kappa = args.kappa[0] * u.cm**2/u.g
    args.temp = args.temp * u.K
    
    # Fluxes
    if args.fluxfile is not None:
        unit = args.fluxfile[1]['F']
        args.flux = args.fluxfile[0]['F'] * unit
    else:
        args.flux = np.array([args.flux.value]) * args.flux.unit

    # Frequency
    if args.wav is not None:
        args.nu = args.wav.to(u.GHz, equivalencies=u.spectral())

    # Distance
    if args.source is not None:
        args.distance = args.source.distance

def get_dust_mass(args):
    for fnu,t in product(args.flux, args.temp):
        Md = dust_mass(args.nu, fnu, t, args.distance, kappa=args.kappa,
                Rdg=args.rdg)
        args.results += [[fnu, t, Md]]

def generate_text(args):
    lines = ["Dust mass results"]
    lines += ["Frequency: {0}".format(args.nu.to(u.GHz))]
    lines += ["Distance: {0}".format(args.distance.to(u.kpc))]
    lines += ["Dust opacity: {0}".format(args.kappa.to(u.cm**2/u.g))]
    lines += ["Dust to gas ratio: %f" % args.rdg]
    lines += ["Flux\tTemp\tMass"]
    lines += ["{0.unit:s}\t{1.unit:s}\t{2.unit:s}".format(*tuple(args.results[0]))]
    for r in args.results:
        lines += ["{0.value:.1f}\t{1.value:.1f}\t{2.value:.1f}".format(*tuple(r))]

    return '\n'.join(lines)

def postprocess(args):
    # Print results
    text = generate_text(args)
    print(text)

    if args.output:
        with open(args.output[0], 'w') as out:
            out.write(text)

def dust_mass_calculator(args: Optional[Sequence[str]] = None) -> None:
    # Command line arguments
    parser = argparse.ArgumentParser(
        description='Dust mass calculator',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-k', '--kappa', nargs=2, metavar=('VALUE', 'UNIT'),
                        action=actions.ReadQuantity,
                        help='Dust opacity')
    parser.add_argument('-r', '--rdg', nargs=1, type=float, default=[0.01],
                        help='Dust-to-gas ratio')
    #parser.add_argument('--output', nargs=1, 
    #        help='Output file name')
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('-f', '--flux', metavar=('FLUX', 'UNIT'),
                        nargs=2, action=actions.ReadQuantity,
                        help='Source flux')
    group1.add_argument('--fluxfile', metavar='FILENAME', 
                        action=actions.LoadMixedStructArray,
                        help='File name to read fluxes')
    #group2 = parser.add_mutually_exclusive_group(required=True)
    #group2.add_argument('-s', '--source', action=LoadSourcefromConfig,
    #        help='Load source to get distance')
    parser.add_argument('-d', '--distance', metavar=('DIST', 'UNIT'), 
                        action=actions.ReadQuantity,
                        help='Distance to source')
    group3 = parser.add_mutually_exclusive_group(required=True)
    group3.add_argument('--nu', metavar=('FREQ', 'UNIT'),
                        nargs=2, action=actions.ReadQuantity,
                        help='Frequency of observations')
    group3.add_argument('--wav', metavar=('WAVE', 'UNIT'),
                        nargs=2, action=actions.ReadQuantity,
                        help='Wavelength of observations')
    group4 = parser.add_mutually_exclusive_group(required=True)
    group4.add_argument('-t', '--temp', metavar='TEMP', nargs='+', 
                        action=actions.ReadQuantity,
                        help='Dust temperature list')
    group4.add_argument('--trange', dest='temp', metavar=('LOW', 'HIGH', 'N'), nargs=3, 
                        action=actions.ArrayFromRange,
                        help='Dust temperature range in Kelvin')
    args = parser.parse_args()
    for fn in args.pipe:
        fn(args)

if __name__=='__main__':
    dust_mass_calculator(sys.argv[1:])

