"""Tools for extracting information from 2-D images."""
from typing import Optional, Union, Sequence, List, Callable, Tuple, Dict

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from radio_beam import Beam
from scipy import ndimage
import astropy.units as u
import numpy as np

from .masking import emission_mask, mask_structures, plot_mask
from ..maths import distance_array
from ..converters import quantity_from_hdu, array_to_hdu

def copy_header_keys(ref_header: Dict,
                     new_header: Dict,
                     keys: Sequence = ('BUNIT', 'BMIN', 'BMAJ', 'BPA')
                     ) -> Dict:
    """Copy selected keys from one header into another header."""
    for key in keys:
        new_header[key] = ref_header[key]

    return new_header

def squeeze_image(image: fits.PrimaryHDU) -> fits.PrimaryHDU:
    """Reduce number of axes not used in input image."""
    header = WCS(image, naxis=['longitude', 'latitude']).to_header()

    return fits.PrimaryHDU(data=np.squeeze(image.data), header=header)

def get_coord_axes(image: fits.PrimaryHDU) -> Tuple[u.Quantity]:
    """Get the coordinate axes from header."""
    wcs = WCS(img.header)

    nx = image.data.shape[1]
    ny = image.data.shape[0]

    xaxis = wcs.wcs_pix2world(np.arange(nx), np.zeros(nx), 0)[0]
    yaxis = wcs.wcs_pix2world(np.zeros(ny), np.arange(ny), 0)[1]

    xunit = u.Unit(image.header.get('CUNIT1'))
    yunit = u.Unit(image.header.get('CUNIT2'))

    return xaxis*xunit, yaxis*yunit

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

def get_peak(image: fits.PrimaryHDU) -> Tuple[SkyCoord, u.Quantity]:
    """Get the coordinate of the peak and peak value."""
    data = quantity_from_hdu(image)
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    ymax, xmax = np.unravel_index(np.nanargmax(data), data.shape)
    position = SkyCoord.from_pixel(xmax, ymax, wcs)

    return position, data[ymax, xmax]

def position_in_image(position: SkyCoord, image: fits.PrimaryHDU) -> bool:
    """Is the `position` in the `image`?"""
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    x, y = skycoord_to_pixel(position, wcs)

    try:
        return image.data[y, x] is not None
    except IndexError:
        return False

def positions_in_image(positions: Sequence[SkyCoord],
                       image: fits.PrimaryHDU) -> List[SkyCoord]:
    """Filter positions present in an image."""
    filtered = []
    for position in positions:
        if position_in_image(position, image):
            filtered.append(position)

    return filtered

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
    header = copy_header_keys(image.header, header)
    cutout = fits.PrimaryHDU(data=cutout.data, header=header)
    if filename is not None:
        cutout.writeto(filename, overwrite=True)

    return cutout

def identify_structures(
    image: fits.PrimaryHDU,
    mask: Optional[np.array] = None,
    calcmask: bool = False,
    threshold: Optional[u.Quantity] = None,
    nsigma: float = 5,
    min_area: float = 9,
    plot: Optional['pathlib.Path'] = None,
    log: Callable = print,
) -> Tuple[List[SkyCoord], List[SkyCoord], List[u.Quantity]]:
    """Find valid emission structures in input image.

    If `mask` is not given, it first it creates a mask with valid (non-NaN)
    data. This mask can be combined with a threshold mask if `calcmask` is set
    to `True`. This mask is generated from `threshold` or `nsigma*rms` with
    `rms` estimated from the median absolute deviation of the data. Then it
    identifies connected structures in the mask and filters out structures with
    small areas. Finally it looks for the centroid and lenghts along each axes
    for each structure.

    If the input image has beam parameters in the header, then `min_area` is
    the minimum number of beam areas for mask structures. Otherwise it is the
    minimum number of pixels for mask structures.

    Args:
      images: input image.
      mask: optional; peak mask.
      calcmask: optional; calculate a threshold mask?
      threshold: optional; emission threshold for mask.
      nsigma: optional; number of rms levels for threshold mask.
      min_area: optional; number of beam areas or pixels.
      log: optional; logging function.

    Returns:
      A list of coordinates of the centers of structures in the input image.
      The length along x and y axes of each structure.
    """
    # Create mask
    if mask is None:
        mask = ~np.isnan(image.data)
        if calcmask:
            mask = mask & emission_mask(image, threshold=[threshold],
                                        nsigma=nsigma, log=log)

    # Label mask structures
    if ('BMIN' in image.header and
        'BMAJ' in image.header and
        'BPA' in image.header):
        min_area = min_area * pixels_per_beam(image)
    mask, labels, nlabels = mask_structures(mask, min_area=min_area)
    log(f'{nlabels} structures detected in mask')

    # Find centroid
    wcs = WCS(image, naxis=['longitude', 'latitude'])
    pixsize = np.sqrt(wcs.proj_plane_pixel_area()).to(u.arcsec)

    # Find objects
    objects = ndimage.find_objects(labels)

    # Convert to physical quantities
    centroids_coord = []
    lengths = []
    if plot is not None:
        fig, ax = plot_mask(mask, figsize=(15, 15), layout='tight')
    for slcy, slcx in objects:
        cenx = (slcx.start + slcx.stop) / 2
        ceny = (slcy.start + slcy.stop) / 2
        centroids_coord.append(SkyCoord.from_pixel(cenx, ceny, wcs))
        lengths.append((abs(slcy.start - slcy.stop) * pixsize,
                        abs(slcx.start - slcx.stop) * pixsize))
        log(f'Structure centroid: {cenx}, {ceny}')
        log(f'Centroid coordinate: {centroids_coord[-1]}')
        log(f'Structure size: {lengths[-1][1]} x {lengths[-1][0]}')
        if plot is not None:
            ax.scatter(cenx, ceny, s=25, c='m', marker='o')
            ax.text(cenx, ceny,
                    (f'{centroids_coord[-1].ra.deg:.6f}, '
                     f'{centroids_coord[-1].dec.deg:.6f}'),
                    color='c', ha='center', va='bottom')

    # Save plot
    if plot is not None:
        fig.savefig(plot)

    return centroids_coord, lengths

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
    pixel_scale = u.pixel_scale(pixsize / u.pixel)
    mod = mod / (1.*u.pixel).to(u.arcsec, equivalencies=pixel_scale) * u.pixel
    if distance is not None:
        mod = (mod * u.arcsec / (1.*u.arcsec).to(
            u.au, equivalencies=[(u.arcsec, u.au,
                                  lambda x: x * distance.to(u.pc).value,
                                  lambda x: x / distance.to(u.pc).value)]))

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
