"""Tools for extracting information from 2-D images."""
from typing import Optional, Union, Sequence, List, Callable, Tuple

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from radio_beam import Beam
import astropy.units as u
import numpy as np

from .masking import emission_mask, mask_structures
from ..maths import distance_array
from ..converters import quantity_from_hdu, array_to_hdu

def squeeze_image(image: fits.PrimaryHDU) -> fits.PrimaryHDU:
    """Reduce number of axes not used in input image."""
    header = WCS(image, naxis=['longitude', 'latitude']).to_header()

    return fits.PrimaryHDU(data=np.squeeze(image.data), header=header)

def pixels_per_beam(image: fits.PrimaryHDU):
    """Number of pixels in the beam area."""
    # Beam
    beam = Beam.from_fits_header(image.header)

    # Get pixel area
    pixarea = WCS(image, naxis=['longitude', 'latitude'])
    pixarea = pixarea.proj_plane_pixel_area()

    # Number of pixels
    npix = beam.sr / pixarea
    npix = npix.to(u.one).value

    return npix

def stats_at_position(image: fits.PrimaryHDU,
                      position: SkyCoord,
                      radius: u.Quantity,
                      stats: Sequence[Callable] = (np.mean, np.std)
                      ) -> Sequence[u.Quantity]:
    """Calculate statistics in a circle centered at a given position.

    Args:
      image: the map.
      position: center of the circular region.
      radius: radius of the circular region.
      stats: optional; statistical functions.
    """
    # Get distance map
    data = quantity_from_hdu(image)
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    pixsize = np.sqrt(wcs.proj_plane_pixel_area())
    x, y = skycoord_to_pixel(position, wcs)
    dist = distance_array(data.shape, (x, y))
    dist = dist * pixsize

    # Stats
    stats_val = []
    for stat in stats:
        stats_val.append(stat(data[dist < radius]))

    return stats_val

def stats_in_beam(image: fits.PrimaryHDU,
                  position: SkyCoord,
                  stats: Sequence[Callable] = (np.mean, np.std),
                  beam_radius_factor: float = 1,
                  ) -> Sequence[u.Quantity]:
    """Calculate statistics in a beam sized circle centered at a position.

    If `beam_radius_factor` is given, then the radius of the circle will be
    multiplied by this number.

    Args:
      image: the map.
      position: center of the circular region.
      stats: optional; statistical functions.
      beam_radius_factor: optional; multiplicative factor for beam raidus.
    """
    # Beam
    beam = Beam.from_fits_header(image.header)
    radius = np.sqrt(beam.sr / np.pi).to(u.deg) * beam_radius_factor

    return stats_at_position(image, position, radius, stats=stats)

def image_cutout(image: fits.PrimaryHDU, position: SkyCoord,
                 size: Union[u.Quantity, Sequence],
                 filename: Optional['pathlib.Path'] = None) -> fits.PrimaryHDU:
    """Generate a cutout from imput data.

    Args:
      image: data to cut from.
      position: central position of the cutout.
      size: size of the cutout.
      filename: optional; file name to save the cutout.
    """
    # Check 2D image
    aux = squeeze_image(image)

    # Get cutout
    cutout = Cutout2D(aux.data, position=position, size=size, wcs=WCS(aux))

    # Save
    header = cutout.wcs.to_header()
    header['BUNIT'] = image.header['BUNIT']
    cutout = fits.PrimaryHDU(data=cutout.data, header=header)
    if filename is not None:
        cutout.writeto(filename, overwrite=True)

    return cutout

def emission_peaks(image: fits.PrimaryHDU,
                   mask: Optional[np.array] = None,
                   threshold: Optional[u.Quantity] = None,
                   nsigma: float = 5,
                   min_area: float = 9,
                   log: Callable = print
                   ) -> Union[List[SkyCoord], List[u.Quantity]]:
    """Find emission peaks in input image.

    First it creates a threshold mask from `threshold` or `nsigma*rms` with
    `rms` estimated from the median absolute deviation of the data. Then it
    identifies connected structures in the mask and filters out structures with
    small areas. Finally it looks at the position of the peaks at each
    structure.

    If the input image has beam parameters in the header, then `min_area` is
    the minimum number of beam areas for mask structures. Otherwise it is the
    minimum number of pixels for mask structures.

    Args:
      images: input image.
      mask: optional; peak mask.
      threshold: optional; emission threshold for mask.
      nsigma: optional; number of rms levels for threshold mask.
      min_area: optional; number of beam areas or pixels.
      log: optional; logging function.

    Returns:
      A list of coordinates of the peaks in the input image.
      The radius containing the mask structure
    """
    # Create mask
    if mask is None:
        mask = emission_mask(image, threshold=[threshold], nsigma=nsigma,
                             log=log)

    # Fiter small masks
    if ('BMIN' in image.header and
        'BMAJ' in image.header and
        'BPA' in image.header):
        min_area = min_area * pixels_per_beam(image)
    mask, labels, nlabels = mask_structures(mask, min_area=min_area)
    log(f'{nlabels} structures detected in mask')

    # Find max per each structure
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    pixsize = np.sqrt(wcs.proj_plane_pixel_area())
    positions = []
    radii = []
    for label in range(1, nlabels+1):
        # Masked data
        masked_data = np.ma.array(np.squeeze(image.data), mask=labels != label)

        # Get max
        ymax, xmax = np.unravel_index(np.nanargmax(masked_data),
                                      masked_data.shape)
        positions.append(SkyCoord.from_pixel(xmax, ymax, wcs))
        log(f'Found peak at pixel: {xmax}, {ymax}')
        log(f'Peak position: {positions[-1]}')

        # Get radius
        distance = np.ma.array(distance_array(masked_data.shape, (xmax, ymax)),
                               mask=masked_data.mask)
        radii.append(np.nanmax(distance) * pixsize)

    return positions, radii

def intensity_gradient(image: fits.PrimaryHDU,
                       distance: Optional[u.Quantity] = None
                       ) -> Tuple[fits.PrimaryHDU]:
    """Calculate a gradient map of the input image.

    If distance is given the spatial unit of the gradient modulus will be
    converted from arcsec to au.

    Args:
      image: input image.
      distance: optional; distance to the source.

    Returns:
      The modulus and direction arrays of the emission gradient.
    """
    # Calculate modulus and direction
    data = quantity_from_hdu(image)
    grow, gcol = np.gradient(np.squeeze(data.value))
    mod = np.sqrt(grow**2 + gcol**2) * data.unit / u.pixel
    drc = np.degrees(np.arctan2(grow, gcol)) * u.deg

    # Convert pixels
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    pixsize = np.sqrt(wcs.proj_plane_pixel_area())
    print(pixsize.to(u.arcec))
    pixel_scale = u.pixel_scale(pixsize / u.pixel)
    mod = mod.to(data.unit / u.arcsec, equivalencies=pixel_scale)
    if distance is not None:
        mod = mod.to(data.unit / u.au,
                     equivalencies=[(u.arcsec, u.au,
                                     lambda x: x * distance.to(u.pc).value,
                                     lambda x: x / distance.to(u.pc).value)])

    return array_to_hdu(mod, image), array_to_hdu(drc, image)

def minimal_radius(image: fits.PrimaryHDU, position: SkyCoord):
    """The radius where a circle centered at `position` contains valid data."""
    # Convert to pixel
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    ypos, xpos = skycoord_to_pixel(position, wcs)

    # Radius
    pixsize = np.sqrt(wcs.proj_plane_pixel_area())
    distance = distance_array(np.squeeze(image.data).shape, (xpos, ypos),
                              mask=np.isnan(np.squeeze(image.data)))
    return np.nanmax(distance) * pixsize
