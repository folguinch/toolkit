import astropy.units as u

def get_restfreq(cube):
    return cube.header['RESTFREQ'] * u.Hz

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

