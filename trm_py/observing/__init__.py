#!/usr/bin/env python3

"""
Module to support observing-related routines such as eclipse predicting.
Provides some classes for reading in ephemeris and positional files and
a routine specific to the Thai telescope
"""

import math as m
import numpy as np
from trm import subs

from astropy import time, coordinates as coord, units as u
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz, SkyCoord
from astropy.time import TimeISO

# Longitude, latitude, height, zenith hole (deg) keyed by name
SITES = {

    'WHT' : {
        'long' : '-17 52 53.9', 'lat' : '+28 45 38.3',
        'height' : 2332., 'zhole' : 3.
    },

    'GTC' : {
        'long' : '-17 53 31', 'lat' : '+28 45 24',
        'height' : 2300., 'zhole' : 8.
    },

    'NTT' : {
        'long' : '-70 44 00.', 'lat' : '-29 15 00.', 
        'height' : 2400., 'zhole' : 5.
    },

    'VLT' : {
        'long' : '-70 24 9.9', 'lat' : '-24 37 30.3', 
        'height' : 2635., 'zhole' : 4.
    },

    'Keck' : {
        'long' : '-155 28 30.04', 'lat' : '+19 49 34.9', 
        'height' : 4145., 'zhole' : 4.
    },

    'TNT' : {
        'long' : '+98 28 00', 'lat' : '+18 34 00',
        'height' : 2457., 'zhole' : 2.
    },

    'A2.3' : {
        'long' : '+22 11 46', 'lat' : '+37 59 04',
        'height' : 2340., 'zhole' : 2.
    },

    'SAAO' : {
        'long' : '+20 48 42', 'lat' : '-32 23 14',
        'height' : 1798., 'zhole' : -1.
    },

    'Magellan' : {
        'long' : '-70 42 05.9', 'lat' : '-29 00 35.8',
        'height' : 2275., 'zhole' : 4.
    },
}

class Sdata (object):
    """
    Stores positional, ephemeris and line weight data on a star.
    """
    def __init__(self, position, ephem, lweight):
        self.position = position
        self.eph = ephem
        self.lw = lweight

class Ephemeris (object):
    """
    Stores an ephemeris.
    """
    def __init__(self, string):
        """
        Constructs an Ephemeris when passed a string having the form:

          Time Poly T0 eT0 P eP [Q eQ]

        Time    -- a time type: 'BMJD', 'HJD', or 'HMJD'
        Poly    -- 'linear', 'quadratic'
        T0, eT0 -- constant term and uncertainty
        P, eP   -- linear coefficient and uncertainty
        Q, eQ   -- quadtraic term and uncertainty
        """

        subv = string.split()

        self.time   = subv[0]
        if self.time != 'HJD' and self.time != 'BJD' and \
           self.time != 'HMJD' and self.time != 'BMJD':
            raise Exception('Ephem: unrecognised time type = ' + self.time)

        self.poly   = subv[1]
        if self.poly == 'linear':
            if len(subv) != 6:
                raise Exception('Ephem: linear ephemerides require 4 numbers')
        elif self.poly == 'quadratic':
            if len(subv) != 8:
                raise Exception('Ephem: quadratic ephemerides require 6 numbers')
        else:
            raise Exception("Ephem: only 'linear' or 'quadratic' recognised ephemeris types")

        self.coeff  = [float(s) for s in subv[2::2]]
        self.ecoeff = [float(s) for s in subv[3::2]]
        if self.time == 'HJD' or self.time == 'BJD':
            if self.coeff[0] < 2400000.:
                raise Exception('Ephem: for HJD or BJD expect T0 > 2400000')
        elif self.time == 'HMJD' or self.time == 'BMJD':
            if self.coeff[0] > 70000.:
                raise Exception('Ephem: for HMJD or BMJD expect T0 < 70000')
        else:
            raise Exception('Ephem: recognised times are HJD, BJD, HMJD or BMJD')

    def __repr__(self):
        return 'Ephemeris(time={!r}, poly={!r}, coeff={!r}, ecoeff={!r})'.format(
            self.time,self.poly,self.coeff,self.ecoeff)

    def phase(self, time):
        """
        Computes phase corresponding to a given time
        """
        pnew = (time - self.coeff[0]) / self.coeff[1]
        if self.poly == 'quadratic':
            pold = pnew - 1
            while np.abs(pnew-pold) > 1.e-8:
                pold = pnew
                pnew = (time - self.coeff[0] - self.coeff[1]*pold**2) / \
                       self.coeff[1]
        return pnew

    def etime(self, cycle):
        """
        Computes uncertainty in time of ephemeris at a given phase
        """
        esum = 0
        fac  = 1.
        for ecf in self.ecoeff:
            esum += fac*ecf**2
            fac  *= cycle**2
        return m.sqrt(esum)

def load_pos_eph(fname):
    """
    Loads positional / ephemeris data from a file, returns in the form of a
    dictionary keyed by the target names. File has following format::

    ES Cet
    02:00:52.17 -09:24:31.7
    null

    Gaia14aae | 2
    16:11:33.97 +63:08:31.81
    BMJD linear 56980.0557197 0.0000013 0.034519487 0.000000016

    The position can be given in any form accepted by
    astropy.coordinates.SkyCoord wbut with the RA and Dec specified in
    sexagesimal form. The "| 2" on the second one is an option to increase the
    line weight

    """
    peinfo = {}
    count = 0
    nline = 0
    name  = None
    with open(fname) as fin:
        for line in fin:
            nline += 1
            try:
                if not line.startswith('#') and not line.isspace():
                    count += 1
                    if count == 1:
                        arr = line.split('|')
                        name = arr[0].strip()
                        lw = 1 if len(arr) == 1 else int(arr[1])
                    elif count == 2:
                        position = SkyCoord(line.strip(), unit=(u.hourangle, u.deg))
                    elif count == 3:
                        try:
                            eph = Ephemeris(line)
                            peinfo[name] = Sdata(position, eph, lw)
                        except:
                            print('No valid ephemeris data found for',name)
                            peinfo[name] = Sdata(position, None, lw)
                        count = 0
            except Exception as err:
                print(err)
                print('Line number',nline)
                print('Line =',line.strip())
                if name:
                    print('Name = ' + name)
                else:
                    print('Name undefined')
                    print('Program aborted.')
                exit(1)

    print('Data on',len(peinfo),'stars loaded.')
    return peinfo

class Switch (object):
    """
    Stores switch target data.
    """
    def __init__(self, line):
        name, ut, delta = line.split('|')
        self.name = name.strip()
        utv = [int(i) for i in ut.split(':')]
        utc = float(utv[0])
        if len(utv) > 1:
            utc += float(utv[1])/60.
        if len(utv) > 2:
            utc += float(utv[2])/3600.
        self.utc   = utc
        self.delta = float(delta)/60.

def load_switches(fname, peinfo):
    swinfo = []
    if fname is not None:
        first = True
        with open(fname) as fin:
            for line in fin:
                if not line.startswith('#') and not line.isspace():
                    swinfo.append(Switch(line))
                    if swinfo[-1].name != 'None' and swinfo[-1].name not in peinfo:
                        raise Exception('switch star: ' + swinfo[-1].name + \
                                        ' not found in position/ephemeris file.')
                    if first:
                        utold = swinfo[-1].utc
                    else:
                        utnew = swinfo[-1].utc
                        if utnew < utold:
                            raise Exception('switch: times not increasing.')
                        utold = utnew
        print(len(swinfo),'target switches loaded.')
    else:
        print('No target switches loaded.')

    return swinfo

class Prange(object):
    """
    Stores phase or time range data
    Very little to this class.
    """
    def __init__(self, name):
        self.name   = name
        self.prange = []

    def add(self, line):
        """
        Given a line of information with 4 components,
        a phase or time range, a pgpplot colour index and a line
        width, this stores the parameters in a an internal list
        prange. Distinguish phase from time (JD) by > or < 1000.
        """
        p1, p2, col, lw = line.split()
        p1  = float(p1)
        p2  = float(p2)
        if p1 < 1000.:
            p2 = p2 - m.floor(p2-p1)
            p_or_t = 'Phase'
        else:
            p_or_t = 'Time'

        col = int(col)
        lw  = int(lw)
        self.prange.append([p1, p2, col, lw, p_or_t])

def load_ptranges(fname, peinfo):
    # Load phase / time ranges
    count = 0
    prinfo = {}
    with open(fname) as fin:
        for line in fin:
            try:
                if line.startswith('#') or line.isspace():
                    if count:
                        if name in peinfo:
                            prinfo[name] = pr
                        else:
                            print(name,'not found in position/ephemeris file and will be skipped.')
                        count = 0
                else:
                    count += 1
                    if count == 1:
                        name = line.strip()
                        pr = Prange(name)
                    elif count > 1:
                        pr.add(line)

            except ValueError:
                print('Could not interpret',line.strip(),'as a phase or time range')
                exit(1)

    print('Data on',len(prinfo),'phase ranges loaded.')
    return prinfo

def tnt_alert(alt, az):
    """
    TNT horizon is complicated by the TV mast. This warns
    that one is being obscured.

    Arguments::

      alt : altitude, degrees
      az  : azimuth, degrees

    Returns with True if the mast is in the way.
    """
    # three critical points defining mast are in (alt,az) (21.0,25.5) (lower
    # left-hand point), (73.5,33.5) (apex), (21.0,50.) (lower right-hand
    # point). Assume that these define two great circles, and find the axial
    # vectors of these great circles.
    #ALT  = np.radians(np.array([21.0,73.5,21.0]))
    #AZ   = np.radians(np.array([25.5,33.5,50.0]))

    # these ones have been expanded 2 degrees because we have had
    # cases where targets have been affected by the mast but eplanner
    # didn't say so i also moved the mast a bit to match better what
    # happened on one of Steven's targets
    ALT  = np.radians(np.array([21.0,76,21.0]))
    AZ   = np.radians(np.array([25.,35.5,54.5]))
    calt, salt = np.cos(ALT), np.sin(ALT)
    caz,  saz  = np.cos(AZ), np.sin(AZ)
    v1  = subs.Vec3(saz[0]*calt[0], caz[0]*calt[0], salt[0])
    v2  = subs.Vec3(saz[1]*calt[1], caz[1]*calt[1], salt[1])
    v3  = subs.Vec3(saz[2]*calt[2], caz[2]*calt[2], salt[2])

    # a1, a2 defined to be downward pointing axial vectors corresponding to
    # great cirles representing each extreme of the mast. Inside mast if
    # actual vector gives positive dot product with both axial vectors.
    a1  = subs.cross(v1,v2)
    a2  = subs.cross(v2,v3)

    ralt, raz = m.radians(alt), m.radians(az)
    v = subs.Vec3(m.sin(raz)*m.cos(ralt), m.cos(raz)*m.cos(ralt), m.sin(ralt))

    return subs.dot(a1,v) > 0 and subs.dot(a2,v) > 0

def load_prefixes(fname, peinfo):
    prefixes = {}
    if fname is not None:
        first = True
        with open(fname) as fin:
            for line in fin:
                if not line.startswith('#') and not line.isspace():
                    if first:
                        name = line.strip()
                        if name not in peinfo:
                            print('Prefix target =',name,'not found in peinfo')
                        first = False
                    else:
                        prefixes[name] = line.strip()
                        first = True
        print(len(prefixes),'target prefixes loaded.')
    else:
        print('No target prefixes loaded.')
    return prefixes


def sun_at_alt(time1, time2, site, alt, tol=0.05):
    """Given two astropy.time.Time values that bracket the times when the Sun is
    at altitude "alt" degrees, this returns the astropy.time.Time when it
    equals alt to with "tol" degrees. If the times do not bracket alt, the
    routine raises an Exception. The times should be close enough that there
    is only one crossing point, but this will not be checked.

    Arguments::

        time1 : (astropy.time.Time)
             first time when the Sun's altitude is either > or < alt

        time2 : (astropy.time.Time)
             second time when the Sun's altitude is either < or > alt

        site : (astropy.coordinates.EarthLocation)
             the position of the observing site.

        alt : (float)
             the altitude of interest in degrees

        tol : (float)
             tolerance of the altitude in degrees
    """

    altazframe1 = AltAz(obstime=time1, location=site)
    sunaltaz1 = get_sun(time1).transform_to(altazframe1)

    altazframe2 = AltAz(obstime=time2, location=site)
    sun = get_sun(time2)
    sunaltaz2 = sun.transform_to(altazframe2)

    if (sunaltaz1.alt.value < alt and sunaltaz2.alt.value < alt) or \
            (sunaltaz1.alt.value > alt and sunaltaz2.alt.value > alt):
        raise Exception('sun_at_alt: times do no bracket the Sun crossing altitude = ' + str(alt))

    if sunaltaz1.alt.value < alt:
        rising = True
    else:
        rising  = False

    # check if we are close enough already (unlikely)
    if abs(sunaltaz2.alt.value-sunaltaz1.alt.value) < tol:
        return time1

    # otherwise, binary chop
    while abs(sunaltaz2.alt.value-sunaltaz1.alt.value) > tol:
        tmid = time.Time((time1.mjd+time2.mjd)/2.,format='mjd')
        altazframe = AltAz(obstime=tmid, location=site)
        sun = get_sun(tmid)
        sunaltaz = sun.transform_to(altazframe)

        if (sunaltaz.alt.value < alt and rising) or (sunaltaz.alt.value > alt and not rising):
            sunaltaz1 = sunaltaz
            time1 = tmid
        else:
            sunaltaz2 = sunaltaz
            time2 = tmid

    return tmid

