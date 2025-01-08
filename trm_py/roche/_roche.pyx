#cython: embedsignature=True
# distutils: language = c++
cimport numpy as np
import numpy as np
from cython.operator cimport dereference as deref
from libcpp cimport bool

from croche cimport *
from csubs cimport Vec3

class RocheError(Exception):
    pass

# helper functions
cdef Vec3 rconv(obj):
    """Create a Subs::Vec3 from a subs.Vec3"""
    cdef Vec3 p
    p.set(obj.x, obj.y, obj.z)
    return p

cdef STAR star_enum(star):
    cdef STAR rstar
    rstar = PRIMARY if star==1 else SECONDARY
    return rstar

def bspot(q, rad, accel=1.0e-7):
    """returns position and stream velocity on stream at radius rad

    (x,y,vx,vy) = bspot(q, rad, acc=1.e-7)

    Args:
        q (float): mass ratio = M2/M1
        rad (float): radius to aim for
        acc (float[optional]): computationa accuracy parameter
    """
    cdef Vec3 r, v, rs, vs
    roche_strinit(q,r,v)
    try:
        roche_stradv(q,r,v,rad,accel,1.0e-3)
        return r.x(), r.y(), v.x(), v.y()
    except:
        raise RocheError("roche.stradv: never achieved desired radius")

def face(q,spin,dirn,rref,pref,star=2,acc=1.e-5):
    """
    p,d,r,g = face(q, spin, dirn, rref, pref, star=2, acc=1.e-5), returns position and direction of element of specific Roche potential.

     q      -- mass ratio = M2/M1
     spin   -- ratio spin/orbital frequency
     dirn   -- direction (a Vec3) to take from centre of mass of star in question.
     rref   -- reference radius greater than any radius of potential in question.
     pref   -- the potential to aim for.
     star   -- 1 or 2 for primary or secondary star. acc    -- accuracy in terms of separation of location.

    Returns p = position, d = direction perpendicular to face, r = radius from centre of mass, g = gravity.
    """

    assert q>0, "roche.face: q <= 0"
    assert rref > 0, "roche.face: rref must be > 0"
    assert (acc > 0) or (acc < 0.1), "roche.face: acc <= 0 or acc > 0.1"
    assert star in [1,2], "roche.face: star must be either 1 or 2"

    # compute roche lobe
    cdef Vec3 pvec, dvec
    cdef Vec3 cdirn = rconv(dirn)
    cdef double r, g
    #rconv(dirn, &cdirn)

    rstar = star_enum(star)
    try:
        roche_face(q,rstar,spin,cdirn,rref,pref,acc,pvec,dvec,r,g)
    except Exception as err:
        raise RocheError(err)
    return ((pvec.x(),pvec.y(),pvec.z()), (dvec.x(),dvec.y(),dvec.z()), r, g)

def ref_sphere(q,spin,ffac,star=2):
    """
    (rref,pref) = ref_sphere(q, spin, ffac, star=2), returns reference radius and potential needed for face.

    q      -- mass ratio = M2/M1
    spin   -- ratio spin/orbital frequency
    ffac   -- linear filling factor of star in question, defined as the radius of the star along the line of
              centres towards its companion divided by the Roche lobe radius in that direction. For spin = 1
              the latter is simply the distance to the L1 point, but otherwise you need to use modified L1
              radii as returned by xl11 or xl12.
    star   -- 1 or 2 for primary or secondary star.
    """
    cdef double rref, pref
    rstar = star_enum(star)
    roche_ref_sphere(q,rstar,spin,ffac,rref,pref)
    return rref, pref

def findi(q, deltaphi, acc=1.e-4, di=1.0e-5):
    """
    findi(q, deltaphi, acc=1.e-4, di=1.e-5), computes inclination for a given mass ratio and phase width
    """
    assert q>0, "roche.findi: q must be > 0"
    assert deltaphi > 0 and deltaphi <= 0.25, "roche.findi: pwidth out of range 0 to 0.25"
    assert acc > 0 and acc <= 0.1, "roche.findi: acc <= 0 or acc > 0.1"
    assert di > 0 and di <= 10, "roche.findi: di <= 0 or di > 10."

    cdef double ilo=65., ihi=90.
    cdef double phi = deltaphi/2.
    cdef Vec3 earth1 = roche_set_earth(ilo,phi)
    cdef Vec3 earth2 = roche_set_earth(ihi,phi)
    cdef Vec3 r
    cdef bint elo = roche_fblink(q, SECONDARY, 1.0, 1.0, acc, earth1, r)
    cdef bint ehi = roche_fblink(q, SECONDARY, 1.0, 1.0, acc, earth2, r)

    cdef double iangle = 0
    if elo and ehi:
        iangle = -2
    elif (not elo) and (not ehi):
        iangle = -1
    else:
        while ihi - ilo > di:
            iangle = (ilo+ihi)/2.
            earth1 = roche_set_earth(iangle,phi)
            if (roche_fblink(q,SECONDARY,1.0,1.0,acc,earth1,r)):
                ihi = iangle
            else:
                ilo = iangle
        iangle = (ilo+ihi)/2
    return iangle

def findq(i, deltaphi, acc=1.e-4, dq=1.0e-5, qlo=0.001, qhi=2.):
    """
    findq(i, deltaphi, acc=1.e-4, dq=1.e-5, qlo=0.001, qhi=2.), computes mass ratio q for a given phase width and inclination
    """
    assert i>0 and i<=90, "roche.findq: iangle out of range 0 to 90"
    assert deltaphi > 0 and deltaphi <= 0.25, "roche.findq: pwidth out of range 0 to 0.25"
    assert acc > 0 and acc <= 0.1, "roche.findq: acc <= 0 or acc > 0.1"
    assert dq > 0 and dq <= 0.1, "roche.findq: dq <= 0 or di > 0.1"


    cdef double phi = deltaphi/2.
    cdef Vec3 r
    cdef Vec3 earth = roche_set_earth(i,phi)

    cdef bint elo = roche_fblink(qlo, SECONDARY, 1.0, 1.0, acc, earth, r)
    cdef bint ehi = roche_fblink(qhi, SECONDARY, 1.0, 1.0, acc, earth, r)

    cdef double q=0
    if elo and ehi:
        q = -2
    elif (not elo) and (not ehi):
        q = -1
    else:
        while qhi - qlo > dq:
            q = (qlo+qhi)/2.
            if (roche_fblink(q,SECONDARY,1.0,1.0,acc,earth,r)):
                qhi = q
            else:
                qlo = q
        q = (qlo+qhi)/2
    return q

def findphi(q,iangle,delta=1.0e-6):
    """
    findphi(q, i, delta=1.e-6), computes deltaphi for a given mass ratio and inclination
    """
    assert q>0, "roche.findphi: q <= 0"
    assert iangle>0 and iangle <= 90, "roche.findphi: iangle out of range 0 to 90"
    assert 0 < delta <= 0.001, "roche.findphi: delta <= 0 or delta > 0.001"

    cdef Vec3 r = Vec3(0,0,0)
    cdef double ingress, egress
    cdef bint status = roche_ingress_egress(q,SECONDARY,1.0,1.0,iangle,delta,r,ingress,egress)
    if not status:
        raise RocheError("roche.findphi: the centre of mass of the white dwarf is not eclipsed")
    return egress-ingress

def fblink(q,iangle,phi,r,ffac=1.,acc=1.0e-4,star=2,spin=1):
    """fblink(q, i, phi, r, ffac=1., acc=1.e-4, star=2, spin=1), computes whether a point is eclipsed or not

    Args:
        r (subs.Vec3): position to check
    """
    assert q>0, "roche.fblink: q <= 0"
    assert 0 < iangle <= 90, "roche.fblink: iangle out of range 0 to 90"
    assert 0 < ffac <= 1, "roche.fblink: ffac out of range 0 to 1"
    assert 0 < acc < 0.1, "roche.fblink: acc <= 0 or acc > 0.1"
    assert star==1 or star==2, "roche.fblink: star must be either 1 or 2"

    cdef STAR rstar = star_enum(star)
    cdef Vec3 earth = roche_set_earth(iangle,phi)
    cdef Vec3 rvec = rconv(r)
    if roche_fblink(q,rstar,spin,ffac,acc,earth,rvec):
        return True
    else:
        return False

def ineg(q,iangle,x,y,z=0,ffac=1.,delta=1.0e-7,star=2,spin=1):
    """ in,out = ineg(q, i, x, y, z=0, ffac=1., delta=1.e-7, star=2, spin=1), computes ingress and egress phase of a point"""
    assert 0 < iangle <= 90, "roche.ineg: iangle out of range 0 to 90"
    assert q>0, "roche.ineg: q must be > 0"
    assert 0 < ffac <= 1, "roche.ineg: ffac out of range 0 to 1"
    assert delta > 0, "roche.ineg: delta <= 0"
    assert star==1 or star==2, "roche.ineg: star must be either 1 or 2"

    cdef Vec3 r = Vec3(x,y,z)
    cdef double ingress, egress
    cdef STAR rstar = star_enum(star)
    cdef bint status = roche_ingress_egress(q,rstar,spin,ffac,iangle,delta,r,ingress,egress)
    if not status:
        raise RocheError("roche.ineg: point is not eclipsed")
    return ingress, egress

def lobe1(q,n=200):
    """
    x,y = lobe1(q, n=200), q = M2/M1. Returns arrays of primary star's Roche lobe.
    """
    assert q>0, "roche.lobe1: q must be > 0"
    assert n>=2, "roche.lobe1: n<2"

    # create output array containing x and y
    x = np.zeros(n, dtype=np.dtype("f"))
    y = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] x_view = x
    cdef float[:] y_view = y

    # compute roche lobe
    roche_lobe1(q, &x_view[0], &y_view[0], n)
    return x, y

def lobe2(q,n=200):
    """
    x,y = lobe2(q, n=200), q = M2/M1. Returns arrays of secondary star's Roche lobe.
    """
    assert q>0, "roche.lobe1: q must be > 0"
    assert n>=2, "roche.lobe1: n<2"

    # create output array containing x and y
    x = np.zeros(n, dtype=np.dtype("f"))
    y = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] x_view = x
    cdef float[:] y_view = y

    # compute roche lobe
    roche_lobe2(q, &x_view[0], &y_view[0], n)
    return x, y

def rpot(q, r):
    """
    rp = rpot(q, r), q = M2/M1. Returns Roche potential at position r.

    Args:
        q (float): mass ratio M2/M1
        r (subs.Vec3): position to calculate potential
    Returns:
        (float): roche potential
    """
    assert q>0, "roche.rpot: q must be > 0"
    return roche_rpot(q,rconv(r))

def rpot1(q, spin, r):
    """
    rp = rpot1(q, spin, r), q = M2/M1. Returns asynchronous Roche potential for star 1

    Args:
        q (float): mass ratio M2/M1
        spin (float): spin/orbital ratio
        r (subs.Vec3): position to calculate potential
    Returns:
        (float): roche potential
    """
    assert q>0, "roche.rpot1: q must be > 0"
    return roche_rpot1(q,spin,rconv(r))

def rpot2(q, spin, r):
    """
    rp = rpot2(q, spin, r), q = M2/M1. Returns asynchronous Roche potential for star 2

    Args:
        q (float): mass ratio M2/M1
        spin (float): spin/orbital ratio
        r (subs.Vec3): position to calculate potential
    Returns:
        (float): roche potential
    """
    assert q>0, "roche.rpot2: q must be > 0"
    return roche_rpot2(q,spin,rconv(r))

def drpot(q, r):
    """
    dx,dy,dz = drpot(q, r)  Returns partial derivs of Roche potential at position r
    Args:
        q (float): mass ratio M2/M1
        r (subs.Vec3): position to calculate potential gradient
    Returns:
        (tuple): x,y,z tuple of potential gradient
    """
    assert q>0, "roche.drpot: q must be > 0"
    cdef Vec3 drp = roche_drpot(q,rconv(r))
    return (drp.x(), drp.y(), drp.z())

def drpot1(q, spin, r):
    """
    dx,dy,dz = drpot1(q, spin, r)  Returns partial derivs of Roche potential for star 1 at position r
    Args:
        q (float): mass ratio M2/M1
        spin (float): spin/orbital ratio
        r (subs.Vec3): position to calculate potential gradient
    Returns:
        (tuple): x,y,z tuple of potential gradient
    """
    assert q>0, "roche.drpot1: q must be > 0"
    cdef Vec3 drp = roche_drpot1(q,spin,rconv(r))
    return (drp.x(), drp.y(), drp.z())

def drpot2(q, spin, r):
    """
    dx,dy,dz = drpot2(q, spin, r)  Returns partial derivs of Roche potential for star 2 at position r
    Args:
        q (float): mass ratio M2/M1
        spin (float): spin/orbital ratio
        r (subs.Vec3): position to calculate potential gradient
    Returns:
        (tuple): x,y,z tuple of potential gradient
    """
    assert q>0, "roche.drpot2: q must be > 0"
    cdef Vec3 drp = roche_drpot2(q,spin,rconv(r))
    return (drp.x(), drp.y(), drp.z())

def shadow(q,iangle,phi,n=200,dist=5.,acc=1.e-4):
    """
    Compute roche shadow region in equatorial plane
    x,y,s = shadow(q, iangle, phi, n=200, dist=5., acc=1.e-4),

    Args:
        q (float):  M2/M1.
        iangle (float): inclination
        phi (float): orbital phase
        n (int): number of points in output arrays
        dist (float): maximum distance to search for shadow measured
        acc (float): numerical accuracy parameter
    Returns:
        (np.ndarray[float]): x array of roche lobe shadow
        (np.ndarray[float]): y array of roche lobe shadow
        (np.ndarray[bool]): true/false if genuine shade or not.
                            The array goes all the way round and when not in shade it
                            will be glued to the red star. This array allows you to see if this is the case or not.
    """
    assert q>0, "roche.shadow: q <= 0"
    assert 0 < iangle <= 90, "roche.shadow: iangle out of range 0 to 90"
    assert 0 < acc < 0.1, "roche.shadow: acc <= 0 or acc > 0.1"
    assert n>2, "roche.shadow: n < 2"
    assert dist > 0, "roche.shadow: dist <= 0"

    # create output arrays
    x = np.zeros(n, dtype=np.dtype("f"))
    y = np.zeros(n, dtype=np.dtype("f"))
    # boolean arrays just need one bit, so use char and retype on return
    s = np.zeros(n, dtype=np.uint8)

    # memory views
    cdef float[:] x_view = x
    cdef float[:] y_view = y
    cdef char[:] s_view = s

    # compute roche shadow
    roche_shadow(q,iangle,phi,dist,acc,&x_view[0],&y_view[0],<bool*> &s_view[0],n)

    return x, y, s.view(np.bool)

def streamr(q, rad, n=200):
    """
    x,y = streamr(q, rad, n=200), returns arrays of the gas stream. q = M2/M1, rad = minimum radius to aim for.
    """
    assert q>0, "roche.streamr: q <= 0"
    assert 0 < rad < 1, "roche.streamr: rad < 0 or rad > 1"
    assert n >= 2, "roche.streamr: n < 2"

    # create output arrays
    x = np.zeros(n, dtype=np.dtype("f"))
    y = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] x_view = x
    cdef float[:] y_view = y

    try:
        roche_streamr(q,rad,&x_view[0],&y_view[0],n)
    except Exception as err:
        raise RocheError("roche.streamr: " + err)
    return x, y

def stream(q, step, n=200):
    """
    x,y = stream(q, step, n=200), returns arrays of the gas stream. q = M2/M1, step=distance between adjacent points.
    """
    assert q>0, "roche.stream: q <= 0"
    assert 0 < step < 1, "roche.stream: step <= 0 or rad > 1"
    assert n >= 2, "roche.stream: n < 2"

    # create output arrays
    x = np.zeros(n, dtype=np.dtype("f"))
    y = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] x_view = x
    cdef float[:] y_view = y

    try:
        roche_stream(q,step,&x_view[0],&y_view[0],n)
    except Exception as err:
        raise RocheError("roche.stream: " + err)
    return x, y

def astream(q, type, r0, v0, step, n=200, acc=1.0e-9):
    """
    returns arrays of the gas stream given arbitrary initial conditions

    Args:
        q (float): M2/M1
        type (int): 0, 1, 2 or 3.
            0 gives x,y positions, 1,2,3 give
            different velocities -- see vstream for more detail.
        r0 (subs.Vec3): starting position
        v0 (subs.Vec3): starting velocity
        step (float): step size for integration
        n (int): number of points in integration
        acc (float): numerical accuracy parameter
    Returns:
        (np.ndarray) x position of stream
        (np.ndarray) y position of stream
    """
    assert q>0, "roche.astream: q <= 0"
    assert 0 < step < 1, "roche.astream: step <= 0 or rad > 1"
    assert n >= 2, "roche.astream: n < 2"
    assert 0 <= type <= 3 , "roche.astream: type out of range 0 to 3"
    assert 0 < acc < 0.1, "roche.astream: acc <= 0 or acc > 0.1"

    # create output arrays
    x = np.zeros(n, dtype=np.dtype("f"))
    y = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] x_view = x
    cdef float[:] y_view = y

    # setup some variables
    cdef double xold, yold, apx, apy
    cdef double time, dist, tdid, tnext, frac, ttry
    cdef double vel, smax
    cdef int lp = 0
    smax = min(1.0e-3,step/2.)

    # convert vectors
    cdef Vec3 r = rconv(r0)
    cdef Vec3 v = rconv(v0)
    cdef int num = n
    cdef float cstep = step
    cdef int ctype = type
    cdef float cq = q
    # carry out the stream integration
    try:
        xold = r.x()
        yold = r.y()
        if type == 0:
            x_view[0] = xold
            y_view[0] = yold
        else:
            # calculate stream velocity
            roche_vtrans(cq,type,r.x(),r.y(),v.x(),v.y(),apx,apy)
            x_view[0] = apx
            y_view[0] = apy

        vel = (v.x()**2 + v.y()**2)**.5
        ttry = smax/max(1.0e20,vel)

        # loop and integrate
        while(lp < num-1):

            roche_gsint(cq, r, v, ttry, tdid, tnext, time, acc)
            dist = ( (r.x()-xold)**2 + (r.y()-yold)**2 )**0.5
            if dist > cstep:
                frac = cstep / dist

                if ctype == 0:
                    x_view[lp+1] = x_view[lp] + (r.x() - x_view[lp])*frac
                    y_view[lp+1] = y_view[lp] + (r.y() - y_view[lp])*frac
                else:
                    roche_vtrans(q,type,r.x(),r.y(),v.x(),v.y(),apx,apy)
                    x_view[lp+1] = x_view[lp] + (apx - x_view[lp])*frac
                    y_view[lp+1] = y_view[lp] + (apy - y_view[lp])*frac

                xold = r.x()
                yold = r.y()
                lp = lp + 1

            vel = (v.x()**2 + v.y()**2)**0.5
            ttry = min(smax/vel, tnext)

    except Exception as err:
        raise RocheError("roche.astream: " + repr(err))

    # return the values
    return x, y

def strmnx(q,n=1,acc=1.0e-7):
    """
    Calculates position & velocity of n-th turning point of stream.
    x,y,vx1,vy1,vx2,vy2 = strmnx(q, n=1, acc=1.e-7), q = M2/M1.
    Two sets of velocities are reported, the first for the pure stream, the second for the disk at that point.
    """
    assert q>0, "roche.strmnx: q <= 0"
    assert acc > 0, "roche.strmnx: acc <= 0 "
    assert n >= 1, "roche.stream: n < 1"
    cdef Vec3 r, v

    roche_strinit(q,r,v)
    for i in range(n):
        roche_strmnx(q, r, v, acc)
    cdef double tvx1, tvy1, tvx2, tvy2
    roche_vtrans(q, 1, r.x(), r.y(), v.x(), v.y(), tvx1, tvy1)
    roche_vtrans(q, 2, r.x(), r.y(), v.x(), v.y(), tvx2, tvy2)
    return (r.x(), r.y(), tvx1, tvy1, tvx2, tvy2)

def vlobe1(q,n=200):
    """
    (vx,vy) = vlobe1(q, n=200), q = M2/M1. Returns arrays of primary star's Roche lobe in velocity space.
    """
    assert q>0, "roche.vlobe1: q <= 0"
    assert n >= 2, "roche.vlobe1: n < 2"

    # create output arrays
    vx = np.zeros(n, dtype=np.dtype("f"))
    vy = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] vx_view = vx
    cdef float[:] vy_view = vy

    roche_vlobe1(q, &vx_view[0], &vy_view[0], n)
    return vx, vy

def vlobe2(q,n=200):
    """
    (vx,vy) = vlobe2(q, n=200), q = M2/M1. Returns arrays of secondary star's Roche lobe in velocity space.
    """
    assert q>0, "roche.vlobe2: q <= 0"
    assert n >= 2, "roche.vlobe2: n < 2"

    # create output arrays
    vx = np.zeros(n, dtype=np.dtype("f"))
    vy = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] vx_view = vx
    cdef float[:] vy_view = vy

    roche_vlobe2(q, &vx_view[0], &vy_view[0], n)
    return vx, vy

def vstream(q, step=0.01, vtype=1, n=60):
    """
    Returns arrays of positions of the gas stream in velocity space.

    vx,vy = vstream(q, step=0.01, type=1, n=60)

    Args:
        q (float):  M2/M1
        step (float): step is measured as a fraction of the distance to the inner Lagrangian point from the primary star.
        vtype (int): vtype=1 is the straight velocity of the gas stream while vtype=2 is the velocity of the disc along the stream.
        n (int): number of points in output arrays
    Returns:
        (np.ndarray): array of x-velocities
        (np.ndarray): array of y-velocities
    """
    assert q>0, "roche.vstream: q <= 0"
    assert 0 < step < 1, "roche.vstream: step <= 0 or step > 1"
    assert n >= 2, "roche.vstream: n < 2"
    assert vtype == 1 or vtype == 2, "roche.vstream: vtype must be 1 or 2"

    # create output arrays
    vx = np.zeros(n, dtype=np.dtype("f"))
    vy = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] vx_view = vx
    cdef float[:] vy_view = vy

    try:
        roche_vstrreg(q, step, &vx_view[0], &vy_view[0], n, vtype)
    except Exception as err:
        return RocheError("roche.vstream: " + repr(err))

    return vx, vy

def pvstream(q, step=0.01, vtype=1, n=60):
    """
    Returns arrays of positions, velocities time and jacobi constant along THE GAS STREAM

    x, y, vx, vy, t, jac = pvstream(q, step=0.01, type=1, n=60)

    Args:
        q (float): M2/M1.
        step (float): is measured as a fraction of the distance to the inner Lagrangian point.
        vtype (int): vtype=1 is the straight velocity of the gas stream while vtype=2 is the velocity of the disc along the stream.
        n (int): number of points in output arrays
    Returns:
        (np.ndarray): array of x-positions
        (np.ndarray): array of y-positions
        (np.ndarray): array of x-velocities
        (np.ndarray): array of y-velocities
        (np.ndarray): time
        (np.ndarray): jacobi constant
    """
    assert q>0, "roche.pvstream: q <= 0"
    assert 0 < step < 1, "roche.pvstream: step <= 0 or step > 1"
    assert n >= 2, "roche.pvstream: n < 2"
    assert vtype == 1 or vtype == 2, "roche.pvstream: vtype must be 1 or 2"

    # create output arrays
    vx_arr = np.zeros(n, dtype=np.dtype("f"))
    vy_arr = np.zeros(n, dtype=np.dtype("f"))
    x_arr  = np.zeros(n, dtype=np.dtype("f"))
    y_arr  = np.zeros(n, dtype=np.dtype("f"))
    t_arr  = np.zeros(n, dtype=np.dtype("f"))
    jc_arr  = np.zeros(n, dtype=np.dtype("f"))

    # memory views
    cdef float[:] vx = vx_arr
    cdef float[:] vy = vy_arr
    cdef float[:]  x = x_arr
    cdef float[:]  y = y_arr
    cdef float[:]  t = t_arr
    cdef float[:]  jc = jc_arr

    cdef double TLOC = 1.0e-8
    cdef double RLOC = 1.0e-8
    cdef int i, decr
    cdef double dt, tvx, tvy, rend, rnext
    cdef Vec3 r, v, rm, vm

    cdef double rl1 = roche_xl1(q)
    try:
        # store L1 as first point
        roche_vtrans(q, vtype, rl1, 0., 0., 0., tvx, tvy)
        x[0]  = rl1
        y[0]  = 0.
        vx[0] = tvx
        vy[0] = tvy
        t[0]  = 0.
        i    = 1
        rnext = rl1*(1.0-step)
        decr  = 1

        # initialise stream
        roche_strinit(q, r, v)
        jc[0] = roche_jacobi(q, r, v)

        while (i < n):

            # advance one step
            dt = roche_stradv(q, r, v, rnext, RLOC, 1.0e-3)
            roche_vtrans(q, vtype, r.x(), r.y(), v.x(), v.y(), tvx, tvy)
            x[i] = r.x()
            y[i] = r.y()
            vx[i] = tvx
            vy[i] = tvy
            t[i] = t[i-1] + dt
            jc[i] = roche_jacobi(q, r, v)
            i = i + 1
            rnext = rnext - rl1*step if decr else rnext + rl1*step

            # locate and store next turning point
            rm = r
            vm = v
            roche_strmnx(q, rm, vm, TLOC)
            rend = rm.length()

            # loop over all radii wanted before next turning point
            while (i < n) and (decr and (rnext > rend)) or \
                (not decr and (rnext < rend)):
                    dt = roche_stradv(q, r, v, rnext, RLOC, 1.0e-3)
                    roche_vtrans(q, vtype, r.x(), r.y(), v.x(), v.y(), tvx, tvy)
                    x[i] = r.x()
                    y[i] = r.y()
                    vx[i] = tvx
                    vy[i] = tvy
                    t[i] = t[i-1] + dt
                    jc[i] = roche_jacobi(q, r, v)
                    i = i + 1
                    rnext = rnext - rl1*step if decr else rnext + rl1*step

            # change direction of search and move it to start at turning
            # point
            rnext = rnext + rl1*step if decr else rnext - rl1*step
            r = rm
            v = vm
            decr = not decr
    except Exception as err:
        raise RocheError("roche.pvstream: " + repr(err))

    return x_arr, y_arr, vx_arr, vy_arr, t_arr, jc_arr

def xl1(q):
    """Calculate the inner Lagrangian point distance.

    Args:
        q (float): mass ratio = M2/M1
    """
    assert q > 0, "roche.xl1: q <= 0"
    return roche_xl1(q)

def xl2(q):
    """Calculate the L2 point distance.

    Args:
        q (float): mass ratio = M2/M1
    """
    assert q > 0, "roche.xl2: q <= 0"
    return roche_xl2(q)

def xl3(q):
    """Calculate the L3 point distance.

    Args:
        q (float): mass ratio = M2/M1
    """
    assert q > 0, "roche.xl3: q <= 0"
    return roche_xl3(q)

def xl11(q, spin):
    """Calculate the L1 point distance if primary is asynchronous.

    Args:
        q (float): mass ratio = M2/M1
        spin (float): ratio of spin/orbital of primary
    """
    assert q > 0, "roche.xl11: q <= 0"
    assert 0 < spin <= 1, "roche.xl11: spin <= 0 or spin > 1"
    return roche_xl11(q,spin)

def xl12(q, spin):
    """Calculate the L1 point distance if secondary is asynchronous.

    Args:
        q (float): mass ratio = M2/M1
        spin (float): ratio of spin/orbital of secondary
    """
    assert q > 0, "roche.xl12: q <= 0"
    assert 0 < spin <= 1, "roche.xl12: spin <= 0 or spin > 1"
    return roche_xl12(q,spin)