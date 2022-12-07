"""Tools for working with `SpectralCube` objects."""
from typing import Callable, List, Optional, Sequence, TypeVar, Tuple, Union
import pathlib

from spectral_cube import SpectralCube
import astropy
import astropy.wcs.utils as aputils
import astropy.units as u
import numpy as np

from ..maths import quick_rms

ConfigParserAdv = TypeVar('ConfigParserAdv')
Projection = TypeVar('Projection')
Path = TypeVar('Path', pathlib.Path, str)
Position = TypeVar('Position', astropy.coordinates.SkyCoord, Tuple[int, int])

def get_restfreq(cube: SpectralCube) -> u.Quantity:
    """Get rest frequency from cube header.

    Args:
      cube: `SpectralCube` object.

    Returns:
      The rest frequency from cube header as a `Quantity` object.
    """
    try:
        restfreq = cube.header['RESTFRQ'] * u.Hz
    except KeyError:
        restfreq = cube.header['RESTFREQ'] * u.Hz
    return restfreq

def get_cube_rms(cube: SpectralCube, use_header: bool = False,
                 sampled: bool = False, log: Callable = print) -> u.Quantity:
    """Do a quick calculation of the rms.

    The `sampled` calculation is recommended for big data cubes.

    Args:
      cube: spectral cube.
      use_header: use value stored in header (if any).
      sampled: optional; determine the rms from a sample of channels.
      log: optional; logging function.
    """
    # From header
    if use_header and 'RMS' in cube.header:
        log('Using cube rms in header')
        return cube.header['RMS'] * cube.unit

    # From cube
    try:
        if not sampled:
            log('Calculating rms from cube')
            rms = quick_rms(cube.unmasked_data[:])
        else:
            log('Calculating rms from 20% of channels')
            nchans = len(cube.spectral_axis)
            fraction = int(0.2 * nchans)
            chans = np.linspace(1, 9, fraction)
            chans = (np.floor(chans / 10 * nchans)).astype(int)
            rms = quick_rms(cube.unmasked_data[chans])
    except TypeError:
        log('Calculating rms from cube')
        rms = quick_rms(cube.unmasked_data[:])
    return rms

def vel_axis(cube: SpectralCube,
             restfreq: Optional[u.Quantity] = None,
             vel_units: u.Unit = u.km/u.s) -> u.Quantity:
    """Get the spectral axis in velocity units.

    If the cube spectral axis units are equivalent to vel_units but restfreq is
    not None, it recalculates the velocity at the given restfreq.

    Args:
      cube: spectral cube to extract the spectral axis.
      restfreq: optional; restfrequency of the spectral axis.
      vel_units: optional; the velocity units.
    """
    # Current spectral axis
    spaxis = cube.spectral_axis

    # Cases if input is already in correct units
    if spaxis.unit.is_equivalent(vel_units) and restfreq is None:
        return spaxis.to(vel_units)
    elif spaxis.unit.is_equivalent(vel_units):
        # Convert to frequency
        datarestfreq = get_restfreq(cube)
        spaxis = cube.with_spectral_unit(u.GHz, velocity_convention='radio',
                rest_value=datarestfreq).spectral_axis

    # Convert frequency to velocity
    rest_value = get_restfreq(cube) if restfreq is None else restfreq
    return cube.with_spectral_unit(vel_units, velocity_convention='radio',
                                   rest_value=rest_value).spectral_axis

def freq_axis(cube: SpectralCube,
              freq_units: u.Unit = u.GHz) -> u.Quantity:
    """Get the spectral axis in frequency units.

    Args:
      cube: spectral cube to extract the spectral axis.
      freq_units: optional; frequency axis units.
    """
    # Current spectral axis
    spaxis = cube.spectral_axis

    # Cases if input is already in correct units
    if spaxis.unit.is_equivalent(freq_units):
        return spaxis.to(freq_units)
    else:
        # Convert velocity to frequency
        rest_value = get_restfreq(cube)
        return cube.with_spectral_unit(freq_units, velocity_convention='radio',
                                       rest_value=rest_value).spectral_axis

def to_common_beam(cube: SpectralCube, log: Callable = print):
    """Convolve input cube to common beam."""
    log(f'Convolving to common beam:')
    try:
        smallest = cube.beams.smallest_beam()
        sm_maj = smallest.major.to(u.arcsec).value
        sm_min = smallest.minor.to(u.arcsec).value
        largest = cube.beams.largest_beam()
        la_maj = largest.major.to(u.arcsec).value
        la_min = largest.minor.to(u.arcsec).value
        log(('Beam extrema: '
             f"{sm_maj:.4f}'' x {sm_min:.4f}'' -- "
             f"{la_maj:.4f}'' x {la_min:.4f}''"))
        common_beam = cube.beams.common_beam(auto_increase_epsilon=True,
                                             tolerance=5e-5)
        cb_maj = common_beam.major.to(u.arcsec).value
        cb_min = common_beam.minor.to(u.arcsec).value
        log(f"Common beam: {cb_maj:.4f}'' x {cb_min:.4f}''")
        cube.allow_huge_operations = True
        cube.meta['comment'] = 'Convolved to common beam'
        return cube.convolve_to(common_beam, allow_huge=True)
    except AttributeError:
        log('Already with single beam')
        return cube

def moments01(cube: SpectralCube,
              low: Optional[u.Quantity] = None,
              up: Optional[u.Quantity] = None,
              filenames: Sequence[Union[str, pathlib.Path]] = ()
              )-> List[Projection]:
    """Calculate zeroth and first moment maps from cube.

    Args:
      cube: spectral cube.
      low: optional; lower flux limit.
      up: optional; upper flux limit.

    Returns:
      A moment 0 and a moment 1 maps
    """
    # Moment 0
    mom0 = cube.moment(order=0)

    # Moment1
    if low is not None:
        mask = cube.mask & (cube >= low)
    else:
        mask = cube.mask
    if up is not None:
        mask = mask & (cube <= up)
    subcube = cube.with_mask(mask)
    mom1 = subcube.moment(order=1)

    if len(filenames) == 2:
        mom0.write(filenames[0], overwrite=True)
        mom1.write(filenames[1], overwrite=True)
    elif filenames:
        for i, m in enumerate([mom0,mom1]):
            m.write(filenames[0] % i, overwrite=True)
    else:
        pass

    return mom0, mom1

def limits_to_chan_range(cube: SpectralCube,
                         freq_range: Optional[u.Quantity] = None,
                         vel_range: Optional[u.Quantity] = None,
                         chan_range: Optional[List[int]] = None,
                         win_halfwidth: Optional[int] = None,
                         vlsr: Optional[u.Quantity] = None,
                         linefreq: Optional[u.Quantity] = None,
                         log: Callable = print,
                         ) -> Tuple[int, int]:
    """Convert specral limits to channel range.

    Args:
      cube: spectral cube.
      freq_range: frequency range.
      vel_range: velocity range.
      chan_range: channel range.
      win_halfwidth: number of channel in half the spectral range.
      vlsr: required by win_halfwidth; vlsr velocity.
      linefreq: required by win_halfwidth; line frequency.
      log: optional; logging function.
    """
    if chan_range is not None or win_halfwidth is not None:
        return get_spectral_limits(cube, chan_range=chan_range,
                                   win_halfwidth=win_halfwidth, vlsr=vlsr,
                                   linefreq=linefreq, log=log)
    elif freq_range is not None or vel_range is not None:
        spectral = cube.spectral_axis
        equiv = u.doppler_radio(get_restfreq(cube))
        if freq_range is not None:
            if not spectral.unit.is_equivalent(freq_range.unit):
                spectral = spectral.to(freq_range.unit, equivalencies=equiv)
            mask = ((spectral > np.min(freq_range)) &
                    (spectral < np.max(freq_range)))
        else:
            if not spectral.unit.is_equivalent(vel_range.unit):
                if linefreq is not None:
                    equiv = u.doppler_radio(linefreq)
                spectral = spectral.to(freq_range.unit, equivalencies=equiv)
            mask = ((spectral > np.min(vel_range)) &
                    (spectral < np.max(vel_range)))
        chans = np.indices(spectral.shape)[mask]
        ch1, ch2 = np.min(chans), np.max(chans)
        log(f'Spectral range equivalent to channels: {ch1} {ch2}')
        return ch1, ch2
    else:
        raise NotImplementedError('Cannot determine channel range')

def get_spectral_limits(cube: SpectralCube,
                        freq_range: Optional[u.Quantity] = None,
                        vel_range: Optional[u.Quantity] = None,
                        chan_range: Optional[List[int]] = None,
                        win_halfwidth: Optional[int] = None,
                        vlsr: Optional[u.Quantity] = None,
                        linefreq: Optional[u.Quantity] = None,
                        log: Callable = print,
                        ) -> Union[u.Quantity, list]:
    """Convert input ranges to ranges that can be applied to a spectral cube.

    For frequency and velocity ranges, it converts the values to the spectral
    units of the cube. Channel range is sorted and returned. Channel half width
    requires a line frequency and a vlsr to determine the range from the
    spectral axis around a specific line. If velocity range is given
    together with a line frequency, then this frequency is used as the rest
    frequency to convert from frequency to velocity.

    Args:
      cube: spectral cube.
      freq_range: frequency range.
      vel_range: velocity range.
      chan_range: channel range.
      win_halfwidth: number of channel in half the spectral range.
      vlsr: required by win_halfwidth; vlsr velocity.
      linefreq: required by win_halfwidth; line frequency.
      log: optional; logging function.
    """
    # Turn everything into frequency or velocity
    spec_unit = cube.spectral_axis.unit
    rangefmt = '{0.value[0]:.4f} {0.value[1]:.4f} {0.unit}'
    if freq_range is not None:
        range_str = rangefmt.format(freq_range.to(u.GHz))
        log(f'Using frequency range = {range_str}')
        if spec_unit.is_equivalent(freq_range.unit):
            return freq_range[0], freq_range[1]
        else:
            # Convert to spectral units
            restfreq = get_restfreq(cube)
            freq_to_vel = u.doppler_radio(restfreq)
            lim = freq_range.to(spec_unit, equivalencies=freq_to_vel)
            return np.min(lim), np.max(lim)
    elif vel_range is not None:
        range_str = rangefmt.format(vel_range.to(u.km/u.s))
        log(f'Using velocity range = {range_str}')
        if spec_unit.is_equivalent(vel_range.unit):
            return vel_range[0], vel_range[1]
        else:
            # Convert to spectral units
            if linefreq is not None:
                restfreq = linefreq
            else:
                restfreq = get_restfreq(cube)
            freq_to_vel = u.doppler_radio(restfreq)
            lim = vel_range.to(spec_unit, equivalencies=freq_to_vel)
            return np.min(lim), np.max(lim)
    elif chan_range is not None:
        log(f'Channel range = {chan_range[0]} {chan_range[1]}')
        return sorted(chan_range)
    elif win_halfwidth is not None and linefreq and vlsr:
        # Get spectral axis
        spaxis = cube.spectral_axis
        restfreq = get_restfreq(cube)

        # Convert to velocity from line
        vel_to_freq = u.doppler_radio(restfreq)
        spaxis = spaxis.to(linefreq.unit, equivalencies=vel_to_freq)
        freq_to_vel = u.doppler_radio(linefreq)
        spaxis = spaxis.to(vlsr.unit, equivalencies=freq_to_vel)
        spaxis = spaxis - vlsr

        # Closest value to vlsr
        ind = np.nanargmin(np.abs(spaxis.value))
        chmin = ind - win_halfwidth
        chmax = ind + win_halfwidth
        chmin = max(chmin, 0)
        chmax = min(chmax, len(spaxis)-1)
        log(f'Using channel range = {chmin} {chmax}')
        return chmin, chmax
    else:
        log('No spectral limit')
        return None, None

def get_subcube(cube: SpectralCube,
                freq_range: Optional[u.Quantity] = None,
                vel_range: Optional[u.Quantity] = None,
                chan_range: Optional[List[int]] = None,
                win_halfwidth: Optional[int] = None,
                vlsr: Optional[u.Quantity] = None,
                linefreq: Optional[u.Quantity] = None,
                blc_trc: Optional[List[int]] = None,
                xy_ranges: Optional[List[int]] = None,
                put_rms: bool = False,
                common_beam: bool = False,
                filenamebase: Optional[Union[str, pathlib.Path]] = None,
                log: Callable = print) -> SpectralCube:
    """Extract a sub-spectral cube in any of the axes.

    The spectral cube can be saved by providing a filenamebase. The output
    filename will replace the suffix (extension) of filenamebase with
    '.subcube.fits'.

    Args:
      cube: spectral cube.
      freq_range: frequency range.
      vel_range: velocity range.
      chan_range: channel range.
      win_halfwidth: number of channel in half the spectral range.
      vlsr: required by win_halfwidth; vlsr velocity.
      linefreq: required by win_halfwidth; line frequency.
      blc_trc: optional; positions of the bottom left and top right corners.
      xy_ranges: optional; ranges in x and y axes.
      put_rms: optional; put the rms in the header?
      common_beam: optional; calculate and convolve sub-cube to common beam?
      filenamebase: optional; base of the output filename.
      log: optional; logging function.
    """
    # Extract the spectral slab
    low, high = get_spectral_limits(cube,
                                    freq_range=freq_range,
                                    vel_range=vel_range,
                                    chan_range=chan_range,
                                    win_halfwidth=win_halfwidth,
                                    vlsr=vlsr,
                                    linefreq=linefreq,
                                    log=log)
    if hasattr(low, 'unit'):
        log('Applying spectral slab')
        subcube = cube.spectral_slab(low, high)
    elif low is None:
        subcube = cube
    else:
        log('Applying channel range')
        subcube = cube[low:high+1, :, :]

    # Extract spatial box
    if blc_trc:
        xmin, ymin, xmax, ymax = blc_trc
    elif xy_ranges:
        xmin, xmax, ymin, ymax = xy_ranges
    else:
        xmin = ymin = xmax = ymax = None
    if xmin is not None:
        log(f'Selecting spatial box: {xmin}:{xmax}, {ymin}:{ymax}')
        subcube = subcube[:, ymin:ymax+1, xmin:xmax+1]

    # Minimal subcube
    subcube = subcube.minimal_subcube()
    log(f'New sub cube shape: {subcube.shape}')

    # Common beam
    if common_beam:
        subcube = to_common_beam(subcube, log=log)

    # Copy RMS
    if 'RMS' in cube.header:
        subcube.meta['RMS'] = cube.header['RMS']
        log("RMS in cube header: {cube.header['RMS']} {cube.unit}")
    elif put_rms:
        # Use original cube to measure rms
        rms = get_cube_rms(cube, log=log)
        if hasattr(rms, 'unit'):
            log(f'Cube rms = {rms.value} {rms.unit}')
            rms = rms.value
        else:
            log(f'Cube rms = {rms}')
        subcube.meta['RMS'] = rms
    else:
        pass

    # Save
    if filenamebase is not None:
        filename = pathlib.Path(filenamebase)
        filename = filename.expanduser().resolve().with_suffix('.subcube.fits')
        log(f'Saving sub-cube: {filename}')
        subcube.write(filename, overwrite=True)

    return subcube

def moment_from_config(cube: SpectralCube,
                       mom: int,
                       config: ConfigParserAdv,
                       vlsr: Optional[u.Quantity] = None,
                       filenamebase: Optional[Path] = None,
                       filename: Optional[Path] = None) -> Projection:
    """Calculate the moments with parameters from config file.

    The parameters are passed to the `get_moment` function.

    Args:
      cube: spectral cube.
      mom: moment to calculate.
      config: advanced configuration proxy.
      vlsr: optional; vlsr velocity.
      filenamebase: optional; base of the output file name.
      filename: optional; output moment file name.
    """
    # Rest frequency
    try:
        linefreq = config.getquantity('freq')
    except KeyError:
        linefreq = config.getquantity('restfreq', fallback=None)

    # Read arguments for get_moment
    kwargs = {'linefreq': linefreq,
              'filenamebase': filenamebase,
              'filename': filename,
              'lower_limit': config.getquantity('lower_limit', fallback=None),
              'upper_limit': config.getquantity('upper_limit', fallback=None),
              'nsigma': config.getfloat('nsigma', fallback=5.),
              'rms': config.getquantity('rms', fallback=None),
              'auto_rms': 'nsigma' in config}

    # For subcube
    subcube_keys = ['freq_range', 'vel_range', 'chan_range', 'win_halfwidth',
                    'blc_trc', 'xy_ranges']
    for key in filter(lambda opt: opt in config, subcube_keys):
        if key in ['freq_range', 'vel_range']:
            kwargs[key] = config.getquantity(key)
        elif key == 'win_halfwidth':
            kwargs[key] = config.getint(key)
            kwargs['vlsr'] = vlsr or config.getquantity('vlsr')
        else:
            kwargs[key] = config.getintlist(key)

    return get_moment(cube, mom, **kwargs)

def get_moment(cube: SpectralCube,
               mom: int,
               linefreq: Optional[u.Quantity] = None,
               filenamebase: Optional[Path] = None,
               filename: Optional[Path] = None,
               lower_limit: Optional[u.Quantity] = None,
               upper_limit: Optional[u.Quantity] = None,
               nsigma: float = 5.,
               rms: Optional[u.Quantity] = None,
               auto_rms: bool = False,
               skip_beam_error: bool = False,
               log: Callable = print,
               **kwargs) -> Projection:
    """ Calculate a moment map.

    If filenamebase is given, the final file name will be filenamebase with
    suffix (extension) replaced by '.subcube.fits'.

    Note that if the linefreq is not given, the rest frequency value in the
    cube header will be used, so the velocity values will be calculated with
    respect to that value.

    If beams differ more than 1% `spectral_cube` will raise `ValueError`. The
    moment calculation can be skipped in this case by setting `skip_beam_error`
    to `True`.

    Args:
      cube: spectral cube.
      mom: moment to calculate.
      linefreq: optional; line frequency.
      filenamebase: optional; base of the output file name.
      filename: optional; output moment file name.
      lower_limit: optional; intensity lower limit.
      upper_limit: optional; intensity upper limit.
      nsigma: determine lower_limit from number of sigma (rms) values.
      rms: optional; cube rms.
      auto_rms: optional; calculate the cube rms.
      skip_beam_error: optional; raise exception if beams differ more than 1%?
      log: optional; logging function.
      kwargs: additional parameters for get_subcube function.
    """
    # Get subcube if needed
    if filenamebase:
        subcube = pathlib.Path(filenamebase).expanduser().resolve()
        subcube = subcube.with_suffix('.subcube.fits')
    elif filename:
        subcube = pathlib.Path(filename).expanduser().resolve()
    else:
        subcube = None
    if len(kwargs) != 0:
        if subcube and subcube.is_file():
            log(f'Reading sub-cube: {subcube}')
            subcube = SpectralCube.read(subcube)
        else:
            log('Obtaining sub-cube')
            subcube = get_subcube(cube, filenamebase=filenamebase,
                                  common_beam=True, **kwargs)
        if filenamebase:
            filenamebase = Path(filenamebase).expanduser().resolve()
            filenamebase = filenamebase.with_suffix('.subcube.fits')
    else:
        subcube = cube.minimal_subcube()
        try:
            common_beam = subcube.beams.common_beam()
            subcube = subcube.convolve_to(common_beam)
        except AttributeError:
            pass

    # Convert to velocity
    if linefreq is None:
        if mom == 1:
            log('Moment 1 centered around cube 0 vel!')
        linefreq = get_restfreq(cube)
    subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                         rest_value=linefreq)

    # Flux mask
    comment1 = None
    if mom > 0:
        if lower_limit:
            log('Using lower flux limit: %s',
                f'{lower_limit.value:.3f} {lower_limit.unit}')
            mask = subcube >= lower_limit
            comment1 = f'Flux lower limit {lower_limit}'
        elif (rms is not None or 'RMS' in subcube.header or
              'RMS' in subcube.meta or auto_rms):
            if rms is not None:
                log(f'Using input rms: {rms.value:.3e} {rms.unit}')
            elif 'RMS' in subcube.header:
                rms = float(subcube.header['RMS']) * subcube.unit
                log(f'Using header rms: {rms.value:.3e} {rms.unit}')
            elif 'RMS' in subcube.meta:
                rms = float(subcube.meta['RMS']) * subcube.unit
                log(f'Using header rms: {rms.value:.3e} {rms.unit}')
            else:
                # Get rms from original cube
                rms = get_cube_rms(cube)
                log(f'Using cube rms: {rms.value:.3e} {rms.unit}')
            low = nsigma * rms
            log(f'{nsigma}sigma lower flux limit: {low.value:.3f} {low.unit}')
            comment1 = f'Flux lower limit {nsigma} x {rms}'
            mask = subcube >= low
        else:
            mask = subcube.mask
        if upper_limit:
            mask = mask & (subcube <= upper_limit)
        subcube = subcube.with_mask(mask)

    # Moment
    try:
        if mom == 2:
            mmnt = subcube.linewidth_fwhm()
        else:
            mmnt = subcube.moment(order=mom)
    except ValueError as exc:
        if skip_beam_error:
            log('Beams differ more than 1%, skipping moment calculation')
            return None
        else:
            raise ValueError from exc

    # Write comments
    if comment1 is not None:
        mmnt.header['comment'] = comment0
    mmnt.header['comment'] = ('Spectral axis range: '
                              f'{subcube.spectral_axis[0]} -- '
                              f'{subcube.spectral_axis[-1]}')

    # RMS of moment 0
    if mom == 0:
        rms = quick_rms(mmnt.hdu.data)
        mmnt.header['RMS'] = rms

    # Save
    if filenamebase:
        filename = pathlib.Path(filenamebase).expanduser().resolve()
        filename = filename.with_suffix(f'.moment{mom}.fits')
        log(f'Saving moment: {filename}')
        mmnt.write(filename)
    elif filename:
        log(f'Saving moment: {filename}')
        mmnt.write(pathlib.Path(filename).expanduser().resolve())

    return mmnt

def spectrum_at_position(cube: SpectralCube,
                         position: Position,
                         spectral_axis_unit: Optional[u.Unit] = None,
                         restfreq: Optional[u.Quantity] = None,
                         vlsr: Optional[u.Quantity] = None,
                         radius: Optional[u.Quantity] = None,
                         size: Optional[u.Quantity] = None,
                         area_pix: Optional[u.Quantity] = None,
                         filename: Optional[Union[pathlib.Path, str]] = None,
                         log: Optional[Callable] = print,
                         ) -> Tuple[u.Quantity]:
    """Extract spectrum at position.

    Args:
      cube: spectral cube.
      position: coordiante or pixel where to extract the spectrum from.
      spectral_axis_unit: optional; unit of the spectral axis.
      restfreq: optional; rest frequency.
      vlsr: optional; LSR velocity.
      radius: optional; average pixels around this radius.
      size: optional; use ellipse major and minor axes to get radius.
      area_pix: optional; source area in pixels.
      filename: optional; output filename.
      log: optional; logging function.
    Returns:
      xaxis: an array with the spectral axis.
      spec: the spectrum.
    """
    # Restfreq
    if restfreq is None:
        restfreq = get_restfreq(cube)

    # Spectral axis unit
    if spectral_axis_unit is None:
        spectral_axis_unit = cube.spectral_axis.unit

    # Spectral axis
    if spectral_axis_unit.is_equivalent(u.km/u.s):
        aux_cube = cube.with_spectral_unit(spectral_axis_unit,
                                           velocity_convention='radio',
                                           rest_value=restfreq)
    elif (spectral_axis_unit.is_equivalent(u.Hz) and  vlsr is not None and
          restfreq is not None):
        # vlsr to freq
        freq_to_vel = u.doppler_radio(restfreq)
        flsr = vlsr.to(spectral_axis_unit, equivalencies=freq_to_vel)

        # Convert to velocity
        aux_cube = cube.with_spectral_unit(u.km/u.s,
                                           velocity_convention='radio',
                                           rest_value=flsr)

        # Convert back
        aux_cube = aux_cube.with_spectral_unit(spectral_axis_unit,
                                               velocity_convention='radio',
                                               rest_value=restfreq)
    elif spectral_axis_unit.is_equivalent(u.Hz):
        # Convert to requested unit
        aux_cube = cube.with_spectral_unit(spectral_axis_unit,
                                           velocity_convention='radio',
                                           rest_value=restfreq)
    else:
        aux_cube = cube

    # Spectra position
    try:
        wcs = cube.wcs.sub(['longitude', 'latitude'])
        #x, y = wcs.all_world2pix(
        #    [[position.ra.degree, position.dec.degree]], 0)[0]
        #x, y = int(x), int(y)
        x, y = aputils.skycoord_to_pixel(position, wcs)
        x, y = int(x), int(y)
    except (TypeError, AttributeError):
        x, y = tuple(map(int, position))

    # Spectrum
    xaxis = aux_cube.spectral_axis
    if radius is not None or size is not None or area_pix is not None:
        if size is not None:
            radius = np.sqrt(size[0] * size[1]) / 2
        if radius is not None:
            wcs = cube.wcs.sub(['longitude', 'latitude'])
            pixsize = np.sqrt(np.abs(np.linalg.det(wcs.pixel_scale_matrix)))
            pixsize = pixsize * u.Unit(wcs.world_axis_units[0])
            radius_pix = radius.to(pixsize.unit) / pixsize
            radius_pix = radius_pix.decompose().value
            log(f'Source radius: {radius.value:.1f} {radius.unit}')
            log(f'Pixel size: {pixsize.value} {pixsize.unit}')
        else:
            radius_pix = np.sqrt(area_pix / np.pi)
        log(f'Source radius: {radius_pix} pixels')
        shape = aux_cube.shape[-2:]
        yy, xx = np.indices(shape)
        mask = ((yy - y)**2 + (xx - x)**2)**0.5
        mask = mask <= radius_pix
        masked_cube = aux_cube.with_mask(mask)
        masked_cube.allow_huge_operations = True
        spec = masked_cube.mean(axis=(1, 2))
    else:
        spec = aux_cube[:, y, x]

    # Shift velocity
    if xaxis.unit.is_equivalent(u.km/u.s) and vlsr is not None:
        xaxis = xaxis - vlsr

    if filename is not None:
        filename = pathlib.Path(filename).expanduser().resolve()
        with filename.open('w', encoding='utf-8') as out:
            out.write('#v\tF\n')
            out.write(f'#{xaxis.unit}\t{spec.unit}\n')
            for dt in zip(xaxis, spec[:]):
                out.write(f'{dt[0].value:f}\t{dt[1].value:f}\n')

    return xaxis, spec
