"""Tools for generating masks from image data."""
from typing import Optional, Callable, Sequence, Tuple, Dict, List

from astropy.wcs import WCS
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import scipy.ndimage as ndimg

from ..maths import quick_rms

def emission_mask(*images: fits.PrimaryHDU,
                  initial_mask: Optional[npt.ArrayLike] = None,
                  threshold: Optional[Sequence[u.Quantity]] = None,
                  nsigma: float = 5,
                  log: Callable = print) -> npt.ArrayLike:
    """Generate a mask from images with emission over a threshold.

    An inital mask with the same dimessions as the data in args can be given.
    The logic `and` operator is used to combine all the masks.

    If threshold is not given, an rms value is estimated for each image in
    `args` and the threshold will be `nsigma*rms`.

    Args:
      images: Input images.
      initial_mask: Optional. Initial mask.
      threshold: Optional. Emission threshold per image.
      nsigma: Optional. Number of rms levels for threshold.
      log: Optional. Logging function.

    Returns:
      A mask array.
    """
    # Initial mask
    if initial_mask is not None:
        mask = np.squeeze(initial_mask)
    else:
        mask = np.ones(np.squeeze(images[0].data).shape, dtype=bool)

    # Check threshold
    if threshold is not None:
        if len(threshold) == len(images):
            pass
        elif len(threshold) == 1:
            threshold = threshold * len(images)
        else:
            raise ValueError('Threshold length does not match images')

    # Iterate over images
    for i, image in enumerate(images):
        # Check shapes
        if np.squeeze(image.data).shape != mask.shape:
            log('Image cannot be combined into mask')
            continue

        # Threshold
        if threshold is None or threshold[i] is None:
            rms = quick_rms(image.data)
            thresh = rms * nsigma
        else:
            data_unit = u.Unit(image.header['BUNIT'])
            thresh = threshold[i].to(data_unit).value

        # Update mask
        mask = mask & (np.squeeze(image.data) > thresh)

    return mask

def mask_structures(mask: npt.ArrayLike, min_area: Optional = None) -> Tuple:
    """Identify mask structures and index them.

    Args:
      maks: Input mask.
      min_area: Optional. Minimum area in pixels of mask structures.

    Returns:
      A mask array with small structures filtered out.
      An array with the labeled structures.
      The number of structures identified.
    """
    # Label mask structures
    structure = ndimg.generate_binary_structure(mask.ndim, 1)
    labels, nlabels = ndimg.label(mask, structure=structure)

    # Filter small structures
    if min_area is not None:
        component_sizes = np.bincount(labels.ravel())
        small = component_sizes < min_area
        small_mask = small[labels]
        mask[small_mask] = False

        # Re-label
        labels, nlabels = ndimg.label(mask, structure=structure)

    return mask, labels, nlabels

def split_mask_structures(mask: npt.ArrayLike,
                          min_area: Optional = None,
                          padding: float = 0.,
                          ) -> Tuple[npt.ArrayLike, List[npt.ArrayLike]]:
    """Identify mask structures and generate a mask for each.

    The generated masks have the same dimensions as the input mask.

    Args:
      maks: Input mask.
      min_area: Optional. Minimum area in pixels of mask structures.
      padding: Optional. Border as a fraction of length of each axis range.

    Returns:
      A mask array with small structures filtered out.
      A list with all the sub-masks.
    """
    # Get final mask and structures
    mask, labels, _ = mask_structures(mask, min_area=min_area)

    # Identify mask structures
    objects = ndimg.find_objects(labels)

    # Iterate over objects and generate individual masks
    shape = mask.shape
    submasks = []
    for slcy, slcx in objects:
        # Padding
        padx = abs(slcx.start - slcx.stop) * padding
        pady = abs(slcy.start - slcy.stop) * padding

        # Valid indices
        slc = (slice(max(0, slcy.start - pady),
                     min(shape[0], slcy.stop + pady)),
               slice(max(0, slcx.start - padx),
                     min(shape[1], slcx.stop + padx)))

        # New mask
        newmask = np.zeros(shape, dtype=bool)
        newmask[slc] = True
        submasks.append(newmask)

    return mask, submasks

def position_in_mask(position: 'astropy.coordinates.SkyCoord',
                     mask: np.array,
                     wcs: WCS) -> bool:
    """Returns value of position in mask."""
    x, y = skycoord_to_pixel(position, wcs)

    return mask[int(y), int(x)]

def plot_mask(
    mask: npt.ArrayLike,
    scatter: Optional[Sequence['astropy.coordinates.SkyCoord']] = None,
    wcs: Optional[WCS] = None,
    scatter_kwds: Optional[Dict] = None,
    **kwargs,
) -> Tuple['matplotlib.figure.Figure', 'matplotlib.axes.Axes']:
    """Create a plot of `mask`.

    Args:
      mask: Mask array.
      scatter: Optional. Sky coordinates of scatter points.
      wcs: Optional. Coordiante system to interpret the scatter points.
      scatter_kwds: Optional. Keywords for `plt.scatter`.
      kwargs: Otional. Keywords for `plt.subplots`.

    Returns:
      A tuple with the figure and axis.
    """
    # Create figure
    figkwds = {'figsize': (15, 15), 'layout': 'tight'}
    figkwds.update(kwargs)
    fig, ax = plt.subplots(1, 1, **figkwds)
    ax.imshow(mask.astype(int), vmin=0, vmax=1, cmap='inferno',
              origin='lower')
    ax.set_xlabel('x (pix)')
    ax.set_ylabel('y (pix)')

    # Scatters
    if scatter is not None and wcs is not None:
        scatter_kwds_defaults = {'s': 25, 'c': 'm', 'marker': 'x'}
        if scatter_kwds is not None:
            scatter_kwds_defaults.update(scatter_kwds)
        for scat in scatter:
            x, y = skycoord_to_pixel(scat, wcs)
            ax.scatter(x, y, **scatter_kwds_defaults)
            ax.text(x, y,
                    (f'{scat.ra.deg:.6f}, '
                     f'{scat.dec.deg:.6f}'),
                    color='c', ha='center', va='bottom')

    return fig, ax
