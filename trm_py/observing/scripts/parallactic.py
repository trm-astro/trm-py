#!/usr/bin/env python

from __future__ import print_function

usage = \
"""Computes the parallactic angle of a target at a particular time.

Example:

 parallactic.py "05 23 12.1 +30 40 22.2" WHT 40

Returns the parallactic angle of the supplied position
40 minutes ahead of the current time.
"""
import math

import argparse

import numpy as np

from astropy import time, coordinates as coord, units as u
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, FK5

from trm_py import observing
from trm_py.observing import SITES

def parallactic():

    # arguments
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # positional
    parser.add_argument(
        'position', help='position, RA, Dec, sexagesimal'
    )

    parser.add_argument(
        'telescope', help='Telescope name, e.g. WHT, NTT, VLT, TNT'
    )

    parser.add_argument(
        'offset', help='time offset from present time in minutes', type=float,
    )

    # parse them
    args = parser.parse_args()

    # a few checks
    assert(args.telescope in SITES)

    # Get time
    when = time.Time.now() + time.TimeDelta(60*args.offset, format='sec')
    print('\nUT  =',when)

    # Get location
    info = SITES[args.telescope]
    site = coord.EarthLocation.from_geodetic(
        info['long'], info['lat'], info['height']
    )
    print('Lat = {:0.2f} degrees'.format(site.lat.value))
    print('Lon = {:0.2f} degrees'.format(site.lon.value))

    # Get position
    position = SkyCoord(args.position, unit=(u.hourangle, u.deg))

    lst = when.sidereal_time(kind='apparent', longitude=site.lon)
    print('LST =',lst)
    FK5_when = FK5(equinox=when)
    position_when = position.transform_to(FK5_when)
    ha = lst - position_when.ra.to(u.hourangle)
    print('HA  =',ha)

    # now convert to alt/az
    altazframe = AltAz(obstime=when, location=site)
    altaz = position_when.transform_to(altazframe)

    print('Az  = {:0.2f} degrees'.format(altaz.az.value))
    print('Alt = {:0.2f} degrees'.format(altaz.alt.value))

    # Compute latitude, alt, az, dec and ha, all in radians
    lat = math.radians(site.lat.value)
    alt = math.radians(altaz.alt.value)
    az  = math.radians(altaz.az.value)
    delt = math.radians(position_when.dec.value)
    ha = math.radians(15*ha.value)

    # trig
    caz = math.cos(az)

    cdelt = math.cos(delt)
    sdelt = math.sin(delt)

    cha = math.cos(ha)
    sha = math.sin(ha)

    clat = math.cos(lat)
    slat = math.sin(lat)

    calt = math.cos(alt)
    salt = math.sin(alt)

    fac = calt*slat-salt*clat*caz
    dy = cdelt*calt-fac*(slat*cdelt-clat*sdelt*cha)
    dx = clat*sha*fac

    pa = math.degrees(math.atan2(dy, dx))
    if pa < 0: pa += 180
    if pa > 180: pa -= 180
    print(
        '\nParallactic angle  = {:0.2f} or {:0.2f} or {:0.2f} degrees (N-->E)\n'.format(pa-180,pa,pa+180)
    )


if __name__ == '__main__':
    import sys
    parallactic()
    sys.exit(0)
