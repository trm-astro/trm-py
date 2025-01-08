"""
Roche geometry routines

This module provides some easy-to-use routines for computations
to do with Roche geometry. All routines work in coordinates scaled
by the binary separation, and time units such that the binary angular
frequency = 1. The one exception to the rule are the stream routines
which have a step parameter scaled by the distance to the inner
Lagrangian point.


Functions
=========

bsphases   -- computes ingress & egress phases of bright-spot
face       -- computes position, orientation and gravity of element of specified Roche potential
fblink     -- computes whether a point is eclipsed or not
findi      -- computes inclination given mass ratio and deltaphi
findq      -- computes mass ratio given deltaphi and inclination
findphi    -- computes deltaphi given mass ratio and inclination
ineg       -- calculate the ingress/egress phase of a given point.
lobe1      -- the primary star's Roche lobe
lobe2      -- the secondary star's Roche lobe
pvstream   -- produce position and velocity arrays of a stream
qirbs      -- computes mass ratio, inclination and radius from bright-spot phases
ref_sphere -- computes reference radius and potential given a filling factor
rcirc      -- Verbunt & Rappaport circularisation radius formula
rpot       -- computes Roche potential of a point
rpot1      -- computes asynchronous Roche potential for star 1
rpot2      -- computes asynchronous Roche potential for star 2
shadow     -- returns arrays representing the eclipsed region at a given phase
stream     -- the gas stream's path in position space
streamr    -- the gas stream's path in position space to a given radius
strmnx     -- position and velocity of n-th turning point of gas stream
vlobe1     -- the primary star's Roche lobe, velocity space
vlobe2     -- the secondary star's Roche lobe, velocity space
vstream    -- gas stream in velocity coordinates
wdphases   -- computes white dwarf third and fourth contact phases
wdradius   -- computes scaled white dwarf radius give mass ratio, inclination, egress phase width
xl1        -- L1 position
xl2        -- L2 position
xl3        -- L3 position
xl11       -- L1 position, asynchronous primary
xl12       -- L1 position, asynchronous secondary

Classes
=======

RocheError -- exception class

"""

from ._roche import *
import math as m
import numpy as np
import trm.subs as subs

def bsphases(q, iangle, rbs):
    """
    (pbsi,pbse) = bsphases(q, iangle, rbs) -- computes bright-spot ingress and egress phases.

    q       -- mass rtaio = M2/M1
    iangle  -- orbital inclination, degrees
    rbs     -- bright-spot radius, units of separation

    Returns ingress and egress phases of the bright-spot
    """

    (x,y,vx,vy) = bspot(q, rbs)
    (pbi,pbe)   = ineg(q,iangle,x,y,0.)
    return (pbi-1.,pbe-1.)

def qirbs(deltaphi, pbi, pbe, ilo=78., ns=200):
    """
    (q,i,rbs) = qirbs(deltaphi, pbi, pbe, ilo=78., ihi=90., rlo=0.1) -- computes mass
    ratio, inclination and the radius of the bright-spot given the phase width of the
    white dwarf's eclipse and the ingress and egress phases of the bright-spot.

    deltaph -- phase width of white dwarf's eclipse
    pbi     -- ingress phase of bright-spot
    pbe     -- egress phase of bright-spot
    ilo     -- lowest value of inclination to consider (degrees)
    ns      -- number of points along the gas stream

    Return q = mass ratio = M2/M1, i = orbital inclination in degrees and
    rbs the bright-spot radius in units of the separation.

    The function works by guessing an inclination angle then computing a corresponding
    mass ratio to match the deltaphi value. It then computes the path of the gas stream
    which it then converts to ingress/egress phase. It then calculates the point of closest
    approach by working out the distance of the bright-spot in ingress/egress space from
    each line segment joining any two points along the path. It then binary chops on the
    orbital inclination angle until the right combination is found that gives the correct
    bright-spot phases. It is possible that no closest approach point is found in which case
    the routine will throw a RocheError.
    """

    ihi = 90.
    while ihi > ilo + 0.0001:
        iangle = (ilo+ihi)/2.

        # Compute the mass ratio and the stream path
        q    = findq(iangle, deltaphi)
        (x,y,vx1,vy1,vx2,vy2) = strmnx(q)
        rmin = m.sqrt(x**2+y**2)
        (xs,ys) = stream(q, 1.01*rmin, ns)

        # Convert to ingress/egress phase
        pi = []
        pe = []
        for (x,y) in zip(xs,ys):
            try:
                (bi,be) = ineg(q,iangle,x,y,0.)
                pi.append(bi-1)
                pe.append(be-1)
            except ValueError:
                break
        pi = np.array(pi)
        pe = np.array(pe)

        # Path is modelled as a + lambda * v where a and v are vectors
        # and lambda is a multiplier that must end up between 0 and 1
        # to be valid at the point of closest approach.
        dsqmin = 1.e30
        found  = False
        for i in range(len(pi)-1):
            x0 = pi[i]
            y0 = pe[i]
            vx = pi[i+1] - pi[i]
            vy = pe[i+1] - pe[i]
            ambsq   = (x0-pbi)**2+(y0-pbe)**2
            vsq     = vx**2+vy**2
            vdotamb = vx*(x0-pbi)+vy*(y0-pbe)
            dsq     = ambsq-vdotamb**2/vsq
            if vdotamb <= 0 and vdotamb > -vsq and dsq < dsqmin:
                dsqmin = dsq
                imin   = i
                check  = -vy*(x0-pbi)+vx*(y0-pbe)
                lam    = vdotamb/vsq
                found  = True

        if not found:
            raise RocheError('trm.roche.qirbs: stream had no closest approach to bs phases for q,i = ' \
                                 + str(q) + ',' + str(iangle))
        if check < 0:
            ilo = iangle
        else:
            ihi = iangle

    # Compute the radius of the bright-spot and return
    rbs = m.sqrt((xs[imin]+lam*(xs[imin+1]-xs[imin]))**2+(ys[imin]+lam*(ys[imin+1]-ys[imin]))**2)
    return (q,iangle,rbs)

def wdphases(q, iangle, r1, r2=-1, ntheta=200):
    """
    phi3, phi4 = wdphases(q, iangle, r1, ntheta=200)

    Returns the third and fourth contact phases of the white dwarf.

    q      -- mass ratio = M2/M1
    iangle -- orbital inclination, degrees
    r1     -- scaled white dwarf radius = R1/a
    r2     -- scaled secondary radius, < 0 for Roche lobe filling
    ntheta -- number of angles to compute at the limb of white dwarf.
              (used over quadrants)

    The routine searches points equally-spaced at quadrants of the limb
    of the white dwarf to determine the contact phases. It will fail if
    there is no eclipse at all by raising a RocheError. For partial eclipses
    there will be a valid 'fourth' contact (marking the end of eclipse still)
    but the third contact will be set = -1.
    """

    if r2 <= 0.:
        ffac = 1.
    else:
        ffac = r2/(1.-xl1(q))

    def xyv(iangle, r1, phase):
        """
        Returns x, y vectors which define the projected limb of white dwarf
        when viewed at orbital inclination = iangle and orbital phase = phase
        """
        cosp = m.cos(2.*m.pi*phase)
        sinp = m.sin(2.*m.pi*phase)
        x    = subs.Vec3(-r1*sinp,r1*cosp,0.)
        cosi = m.cos(m.radians(iangle))
        sini = m.sin(m.radians(iangle))
        y    = subs.Vec3(-r1*cosi*cosp,-r1*cosi*sinp,r1*sini)
        return (x,y)

    def uneclipsed3(q, iangle, phase, r1, ffac, ntheta):
        """
        Says whether any of the upper-left quadrant of the WD is uneclipsed at
        phase = phase 'any' means all of the ntheta points computed uniformly
        around the quadrant. This can be used to define the 3rd contact
        """
        x,y = xyv(iangle, r1, phase)
        for i in range(ntheta):
            theta = (m.pi/2.)*i/float(ntheta-1)
            v = -x*m.cos(theta) + y*m.sin(theta)
            if not fblink(q, iangle, phase, v, ffac, 1.e-5):
                return True
        return False

    def eclipsed4(q, iangle, phase, r1, ffac, ntheta):
        """
        Says whether any of lower-right quadrant of the WD is eclipsed at
        phase = phase 'Any' means any of ntheta points computed uniformly
        around the quadrant. This can be used to define the 4th contact
        """
        x,y = xyv(iangle, r1, phase)
        for i in range(ntheta):
            theta = (m.pi/2.)*i/float(ntheta-1)
            v = x*m.cos(theta) - y*m.sin(theta)
            if fblink(q, iangle, phase, v, ffac, 1.e-5):
                return True
        return False

    # fourth contact
    phi4lo   = 0.00
    phi4hi   = 0.25
    if not eclipsed4(q, iangle, phi4lo, r1, ffac, ntheta):
        raise RocheError('roche.wdphases: no eclipse at all for q,i,r1 = ' +
                         str(q) + ',' + str(iangle) + ',' + str(r1))

    while phi4hi - phi4lo > r1/ntheta/10.:
        phi4 = (phi4lo+phi4hi)/2.
        if eclipsed4(q, iangle, phi4, r1, ffac, ntheta):
            phi4lo = phi4
        else:
            phi4hi = phi4

    # third contact
    phi3lo   = 0.00
    phi3hi   = 0.25
    if uneclipsed3(q, iangle, phi3lo, r1, ffac, ntheta):
        return (-1.,phi4)

    while phi3hi - phi3lo > r1/ntheta/10.:
        phi3 = (phi3lo+phi3hi)/2.
        if uneclipsed3(q, iangle, phi3, r1, ffac, ntheta):
            phi3hi = phi3
        else:
            phi3lo = phi3

    return (phi3,phi4)


def wdradius(q, iangle, dpwd, ntheta=100, dr=1.e-5, rmax=0.1):
    """
    Computes scaled radius of white dwarf, r1 = R1/a given
    the mass ratio, orbital inclination and phase width of the
    white dwarf's egress (or ingress).

    q      -- mass ratio = M2/M1
    iangle -- orbital inclination
    dpwd   -- phase width of white dwarf ingress/egress feature
    ntheta -- number of points on limb of white dwarf when using
              wdphases during this routine.
    dr     -- tolerance on scaled radius
    rmax   -- maximum scaled radius value to assume
    """

    r1lo = 0.0
    r1hi = rmax
    while r1hi - r1lo > dr:
        r1 = (r1lo+r1hi)/2.
        phi3, phi4 = wdphases(q, iangle, r1, ntheta=ntheta)
        dp = phi4-phi3
        if dp > dpwd:
            r1hi = r1
        else:
            r1lo = r1
    return (r1lo+r1hi)/2.

def jacobi(q, r, v):
    """
    Computes Jacobi constant, more or less total energy per
    unit mass.

    q : mass ratio = M2/M1
    r : position vector
    v : velocity vector
    """
    f1  = 1/(1+q)
    f2  = f1*q
    sec = subs.Vec3(1,0,0)
    return (v.sqnorm()-r.y**2-(r.x-f2)**2)/2.-f1/r.norm()-f2/(r-sec).norm()

def rcirc(q):
    """
    Returns circularisation radius from Verbunt & Rappaport (as fraction of binary
    separation)

    q : mass ratio = M2/M1
    """
    lq = np.log10(q)
    return 0.0883+lq*(-0.04858+lq*(0.11489+0.020475*lq))

# Exception class
class RocheError(Exception):
    """For throwing exceptions from the roche module"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

