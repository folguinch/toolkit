import numpy as np

def cart_to_sph(x, y, z, pos0=[0.,0.,0.]):
    """Convert cartesian coordinates to spherical.

    Parameters:
        x, y, z: cartesian coordinates.
    """
    r = np.sqrt((x-pos0[0])**2 + (y-pos0[1])**2 + (z-pos0[2])**2)
    th = np.arctan2(np.sqrt((x-pos0[0])**2+(y-pos0[1])**2), z-pos0[2])
    phi = np.arctan2(y-pos0[1], x-pos0[0])

    return r, th, phi

def vel_sph_to_cart(vr, vth, vphi, th, phi):
    """Convert velocity from spherical to cartesian.

    Parameters:
        vr, vth, vphi: velocity components.
        th, phi: angles at each point.
    """
    vx = vr*np.sin(th)*np.cos(phi) + vth*np.cos(th)*np.cos(phi) - \
            vphi*np.sin(phi)
    vy = vr*np.sin(th)*np.sin(phi) + vth*np.cos(th)*np.sin(phi) + \
            vphi*np.cos(phi)
    vz = vr*np.cos(th) - vth*np.sin(th)

    return vx, vy, vz
