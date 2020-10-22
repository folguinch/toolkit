import os

import astropy.units as u
import numpy as np

from ..logger import get_logger
from ..math import quick_rms

LOG = get_logger(__name__)

def get_restfreq(cube):
    try:
        restfreq = cube.header['RESTFRQ'] * u.Hz
    except KeyError:
        restfreq = cube.header['RESTFREQ'] * u.Hz
    return restfreq

def get_cube_rms(cube):
    try:
        rms = quick_rms(cube.unmasked_data)
    except TypeError:
        rms = quick_rms(cube.unmasked_data[:])
    return rms

def vel_axis(cube, restfreq=None, vel_units=u.km/u.s):
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

def freq_axis(cube, freq_units=u.GHz):
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

def moments01(cube, low=None, up=None, filenames=[]):
    # Moment 0
    mom0 = cube.moment(order=0)
    
    # Moment1
    if 'low' is not None:
        mask = cube.mask & (subcube>=low)
    else:
        mask = cube.mask
    if 'up' is not None:
        mask = mask & (subcube<=up)
    subcube = cube.with_mask(mask)
    mom1 = subcube.moment(order=1)

    if len(filenames)==2:
        mom0.write(filenames[0], overwrite=True)
        mom1.write(filenames[1], overwrite=True)
    elif basefilename:
        for i, m in enumerate([mom0,mom1]):
            m.write(basefilename % i, overwrite=True)
    else:
        pass

    return mom0, mom1

def get_spectral_limits(cube, freq_range=None, vel_range=None, chan_range=None,
        chan_halfwidth=None, vlsr=None, linefreq=None):
    # Turn everything into frequency or velocity
    spec_unit = cube.spectral_axis.unit
    rangefmt = '{0.value[0]:.4f} {0.value[1]:.4f} {0.unit}'
    if freq_range is not None:
        LOG.info('Using frequency range = %s', 
                rangefmt.fmt(freq_range.to(u.GHz)))
        if spec_unit.is_equivalent(freq_range.unit):
            return freq_range[0], freq_range[1]
        else:
            # Convert to spectral units
            restfreq = get_restfreq(cube)
            freq_to_vel = u.doppler_radio(restfreq)
            lim = freq_range.to(spec_unit, equivalencies=freq_to_vel)
            return np.min(lim), np.max(lim)
    elif vel_range is not None:
        LOG.info('Using velocity range = %s', 
                rangefmt.fmt(vel_range.to(u.km/u.s)))
        if spec_unit.is_equivalent(vel_range.unit):
            return vel_range[0], vel_range[1]
        else:
            # Convert to spectral units
            restfreq = get_restfreq(cube)
            freq_to_vel = u.doppler_radio(restfreq)
            lim = vel_range.to(spec_unit, equivalencies=freq_to_vel)
            return np.min(lim), np.max(lim)
    elif chan_range is not None:
        LOG.info('Channel range = %i %i', *chan_range)
        return sorted(chan_range)
    elif chan_halfwidth is not None and linefreq and vlsr:
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
        chmin = ind-chan_halfwidth
        chmax = ind + chan_halfwidth
        if chmin < 0:
            chmin = 0
        if chmax >= len(spaxis):
            chmax = len(spaxis)-1
        LOG.info('Using channel range = %i %i', chmin, chmax)
        return chmin, chmax
    else:
        LOG.info('No spectral limit')
        return None, None

def get_subcube(cube, freq_range=None, vel_range=None, chan_range=None,
        chan_halfwidth=None, blc_trc=None, xy_ranges=None, vlsr=None, linefreq=None, 
        put_rms=False, filenamebase=None):
    # Extract the spectral slab
    low, high = get_spectral_limits(cube, freq_range=freq_range, vel_range=vel_range,
            chan_range=chan_range, chan_halfwidth=chan_halfwidth, vlsr=vlsr,
            linefreq=linefreq)
    if hasattr(low, 'unit'):
        subcube = cube.spectral_slab(low, high)
    elif low is None:
        subcube = cube
    else:
        subcube = cube[low:high+1,:,:]

    # Extract spatial box
    if blc_trc:
        xmin, ymin, xmax, ymax = blc_trc
    elif xy_ranges:
        xmin, xmax, ymin, ymax = xy_range
    else:
        xmin = ymin = xmax = ymax = None
    if xmin is not None:
        subcube = subcube[:,ymin:ymax+1,xmin:xmax+1]
    
    # Copy RMS
    if 'RMS' in cube.header:
        subcube.header['RMS'] = cube.header['RMS']
    elif put_rms:
        # Use original cube to measure rms
        rms = get_cube_rms(cube)
        if hasattr(rms, 'unit'):
            rms = rms.value
        subcube.header['RMS'] = rms
    else:
        pass

    # Save
    if filenamebase:
        filename = filenamebase+'.subcube.fits'
        LOG.info('Saving sub-cube: %s', filename)
        subcube.write(filename, overwrite=True)

    return subcube

def moment_from_config(cube, mom, config, vlsr=None, filenamebase=None,
        filename=None):
    # Rest frequency
    try:
        linefreq = config.getquantity('freq')
    except:
        linefreq = config.getquantity('restfreq', fallback=None)

    # Read arguments for get_moment
    kwargs = {'linefreq':linefreq, 'filenamebase': filenamebase,
            'filename':filename,
            'lower_limit':config.getquantity('lower_limit', fallback=None), 
            'upper_limit':config.getquantity('upper_limit', fallback=None),
            'nsigma':config.getfloat('nsigma', fallback=5.),
            'rms':config.getquantity('rms', fallback=None),
            'auto_rms':'nsigma' in config}

    # For subcube
    subcube_keys = ['freq_range', 'vel_range', 'chan_range', 'chan_halfwidth',
            'blc_trc', 'xy_ranges']
    aux = [ key for key in subcube_keys if key in config ]
    if len(aux) >= 1:
        for key in aux:
            if key in ['freq_range', 'vel_range']:
                kwargs[key] = config.getquantity(key)
            elif key == 'chan_halfwidth':
                kwargs[key] = config.getint(key)
                kwargs['vlsr'] = vlsr or config.getquantity('vlsr')
            else:
                kwargs[key] = config.getintlist(key)

    return get_moment(cube, mom, **kwargs)

def get_moment(cube, mom, linefreq=None, filenamebase=None, filename=None,
        lower_limit=None, upper_limit=None, nsigma=5., rms=None, 
        auto_rms=False, **kwargs):
    # Get subcube if needed
    if filenamebase:
        subcube = os.path.expanduser(filenamebase) + '.subcube.fits'
    else:
        subcube = None
    if len(kwargs)!=0:
        if subcube and os.path.exists(subcube):
            LOG.info('Reading sub-cube: %s', subcube)
            subcube = SpectralCube.read(subcube)
        else:
            LOG.info('Obtaining sub-cube')
            subcube = get_subcube(cube, filenamebase=filenamebase, 
                    **kwargs)
        if filenamebase:
            filenamebase += '.subcube'
    else:
        subcube = cube

    # Convert to velocity
    if linefreq is None:
        if mom == 1:
            LOG.warn('Moment 1 centered around cube 0 vel')
        linefreq = get_restfreq(cube)
    subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
            rest_value=linefreq)

    # Flux mask
    if mom>0:
        if lower_limit:
            LOG.info('Using lower flux limit: %s', 
                    '{0.value:.3f} {0.unit}'.format(lower_limit))
            mask = subcube >= lower_limit
        elif rms or 'RMS' in subcube.header or auto_rms:
            if rms:
                LOG.info('Using input rms: %s',
                        '{0.value:.3e} {0.unit}'.format(rms))
            elif 'RMS' in subcube.header:
                rms = subcube.header['RMS'] * subcube.unit
                LOG.info('Using header rms: %s',
                        '{0.value:.3e} {0.unit}'.format(rms))
            else:
                # Get rms from original cube
                rms = get_cube_rms(cube)
                LOG.info('Using cube rms: %s',
                        '{0.value:.3e} {0.unit}'.format(rms))
            low = nsigma * rms
            LOG.info('%isigma lower flux limit: %s', nsigma, 
                    '{0.value:.3f} {0.unit}'.format(low))
            mask = subcube >= low
        else:
            mask = subcube.mask
        if upper_limit:
            mask = mask & (subcube<=upper_limit)
        subcube = subcube.with_mask(mask)

    # Moment
    if mom==2:
        mmnt = subcube.linewidth_fwhm()
    else:
        mmnt = subcube.moment(order=mom)

    # RMS of moment 0
    if mom==0:
        rms = quick_rms(mmnt.hdu.data)
        mmnt.header['RMS'] = rms

    # Save
    if filenamebase:
        filename = os.path.expanduser(filenamebase) + '.moment%i.fits' % mom
        LOG.info('Saving moment: %s', filename)
        mmnt.write(filename)
    elif filename:
        LOG.info('Saving moment: %s', filename)
        mmnt.write(os.path.expanduser(filename))

    return mmnt

