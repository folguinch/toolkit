"""Tools for generating masks from image data."""
from typing import Optional, Callable, Sequence, Tuple, Dict

from astropy.wcs import WCS
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimg

from ..maths import quick_rms

def emission_mask(*images: fits.PrimaryHDU,
                  initial_mask: Optional[np.array] = None,
                  threshold: Optional[Sequence[u.Quantity]] = None,
                  nsigma: float = 5,
                  log: Callable = print) -> np.array:
    """Generate a mask from images with emission over a threshold.

    An inital mask with the same dimessions as the data in args can be given.
    The logic `and` operator is used to combine all the masks.

    If threshold is not given, an rms value is estimated for each image in
    `args` and the threshold will be `nsigma*rms`.

    Args:
      images: input images.
      initial_mask: optional; initial mask.
      threshold: optional; emission threshold per image.
      nsigma: optional; number of rms levels for threshold.
      log: optional; logging function.

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

def mask_structures(mask: np.array, min_area: Optional = None) -> Tuple:
    """Identify mask structures and index them.

    Args:
      maks: input mask.
      min_area: optional; minimum area in pixels of mask structures.

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

def position_in_mask(position: 'astropy.coordinates.SkyCoord',
                     mask: np.array,
                     wcs: WCS) -> bool:
    """Returns value of position in mask."""
    x, y = skycoord_to_pixel(position, wcs)

    return mask[int(y), int(x)]

def plot_mask(
    mask: np.array,
    scatter: Optional[Sequence['astropy.coordinates.SkyCoord']] = None,
    wcs: Optional[WCS] = None,
    scatter_kwds: Dict = {},
    **kwargs,
) -> Tuple['matplotlib.Figure', 'matplotlib.Axes']:
    """Create a plot of `mask`.
    
    The `kwargs` are passed to `plt.subplot`.
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
        scatter_kwds_defaults = {'s': 25, 'c': 'm', 'marker': 'o'}
        scatter_kwds_defaults.update(scatter_kwds)
        for scat in scatter:
            x, y = skycoord_to_pixel(scat, wcs)
            ax.scatter(x, y, **scatter_kwds_defaults)
            ax.text(x, y,
                    (f'{scat.ra.deg:.6f}, '
                     f'{scat.dec.deg:.6f}'),
                    color='c', ha='center', va='bottom')

    return fig, ax
