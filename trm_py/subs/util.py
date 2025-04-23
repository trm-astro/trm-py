import os
import re
import copy
import math as m
import numpy as np

from py_subs import cpp_subs

# Constants 
ME     = 9.1093897e-31     # Mass of the electron, SI
MP     = 1.672621637e-27   # Mass of the proton, SI
E      = 1.602176487e-19   # Charge on the electron, SI
H      = 6.6262e-34        # Planck's constant SI
C      = 2.99792458e8      # Speed of light, SI
K      = 1.3806e-23        # Boltzmann's constant, SI
SIGMA  = 5.6703e-8         # Stefan-Boltzmann constant, SI
RSUN   = 6.9599E8          # Radius of Sun, SI
PARSEC = 3.085678e16       # 1 parsec, SI
AU     = 1.4959787061e11   # 1 AU, SI, matches Gauss gravitational constant
G      = 6.673e-11         # Gravitational constant, SI
MSUN   = 1.989e30          # Mass of the Sun, SI
GMSUN  = 1.32712442099e20  # G*MSUN, SI (m^3 s^-2)
KGAUSS = 0.01720209895     # Gauss' gravitational constant G*MSUN, AU^(3/2) day^-1
MJ     = 1.8986e27         # Mass of Jupiter, SI
DAY    = 86400.            # Length of a day, SI
DYEAR  = 365.25            # Length of a Julian year, in days.
YEAR   = DYEAR*DAY         # Length of a Julian year, SI


# Properly wrapping the C++ functions in the cpp module
def voigt(a, x, eps=1.e-8):
    """voigt(a, x) -- returns the Voigt function for a given a and x
    Inputs:
    a -- the Lorentzian width (double)
    x -- the frequency offset (double or numpy array)
    eps -- the precision of the calculation (double)
    Outputs:
    returns the Voigt function (double or numpy array)    
    """
    return cpp_subs.voigt(a, x)

def gammq(a, x):
    """gammq(a, x) -- returns the incomplete gamma function
    Inputs
    a -- the gamma parameter (double)
    x -- the upper limit of the integral (double or numpy array)
    Outputs
    returns the incomplete gamma function (double or numpy array)
    """
    return cpp_subs.gammq(a, x)


def planck(wave, temp, units='FLAM'):
    """
    Returns Planck spectrum in W/m^2/steradian/m or W/m^2/steradian/Hz.

    Arguments:

    wave     -- wave, a number or a numpy array. Units of nm.
    temp     -- temperature, a number or numpy array. Kelvin
    units    -- 'FLAM' gives flambda per m, 'FNU' gives fnu per Hz.

    NB only one of wave or temp can be an array, not both.
    """

    if isinstance(wave, np.ndarray) and isinstance(temp, np.ndarray):
        raise SubsError('trm.subs.planck: either wave or temp can be an array, not both')

    efac = H*C/K/(1.e-9*wave)/temp
    if units == 'FLAM':
        flux = 2.*H*C*C/pow(1.e-9*wave,5)
    elif units == 'FNU':
        flux = 2.*H*C/pow(1.e-9*wave,3)
    else:
        raise ValueError('planck: units = ' + units + ' not recognised.')
    if isinstance(temp, np.ndarray):
        flux = flux*np.ones_like(temp)

    # Need to take a little care evaluating the exponential
    if isinstance(wave, np.ndarray) or isinstance(temp, np.ndarray):
        ok = efac < 50.
        flux[ok] /= np.exp(efac[ok])-1.
        ok = efac >= 50.
        flux[ok] *= np.exp(-efac[ok])
        return flux
    elif efac < 50.:
        return flux/(np.exp(efac)-1.)
    else:
        return flux*np.exp(-efac)

def slitloss(width, seeing):
    """
    Computes fraction of gaussian seeing profile getting through a slit
    (numerically integrates across slit)

    Arguments:

    width  -- slit width
    seeing -- FWHM seeing.

    """
    import scipy.integrate

    def sfunc(x, fwhm):
        sigma = fwhm/2.3548
        return m.exp(-m.pow(x/sigma,2)/2)/m.sqrt(2.*m.pi)/sigma

    (sloss,err) = scipy.integrate.quad(sfunc, -width/2., width/2., args=(seeing))
    return sloss

def splfit(x, y, ye, nspline, thresh, slow=True):
    """
    This uses scipy.interpolate to fit a spline to data with rejection.
    The spline knots are regularly spaced in x and rejection cycles are used
    to kick out points more than rms sigma away from the fit, one at a time.

    x       -- the X values
    y       -- the Y values
    ye      -- uncertainties in Y, used to weight fit. Must be > 0.
    nspline -- the number of splines
    thresh  -- the rejection threshold
    slow    -- True for worst-point-at-a-time rejection, False for lots.

    Returns (sfit,nrej,rms,ok) where sfit = the final fit, nrej = the number of points rejected,
    rms = the RMS of the sigma normalised residuals (i.e. = 1 in ideal case where scatter matches
    errors), and ok = a mask to show which points were used in the final fit.
    """

    from scipy import interpolate

    if ye.min() <= 0.:
        raise SubsError('Y errors have at least one value <= 0')

    ok   = ye > 0.
    nrej = 1

    while nrej:

        xg  = x[ok]
        yg  = y[ok]
        yeg = ye[ok]

        # space knots (which must be interior knots)
        knots   = np.arange(xg[0], xg[-1], (xg[-1]-xg[0])/nspline)[1:]

        # Fit the spline
        sp = interpolate.splrep(xg, yg, 1./yeg, task=-1, t=knots)

        # Calculate the fit (to all points)
        sfit = interpolate.splev(x, sp)

        resid = np.abs((y - sfit)/ye)
        rms   = np.sqrt((resid[ok]**2).sum()/(len(yg)-nspline))

        if slow:
            worst = resid[ok].max()
            if worst > rms*thresh:                
                ok = ok & (resid != worst)
                nrej = 1
            else:
                nrej = 0
        else:
            nok = len(yg)
            ok  = ok & (resid < rms*thresh)
            nnok = len(y[ok])
            nrej = nok - nnok

    nrej = len(x) - len(x[ok])
    return (sfit,nrej,rms,ok)

def sinfit(x, y, ye, period, pfree=False, mask=None):
    """
    Fits a sinusoid to data.

    Model y = c + a*sin(2*pi*x/period) + b*cos(2*pi*x/period)

    x       -- the X values. Take care with e.g. JD when there is an enormous
               offset as it can slow the routine down. Worth offsetting to make
               the mean value close to 0.
    y       -- the Y values
    ye      -- uncertainties in Y, used to weight fit. Must be > 0.
    period  -- period, same units as X
    pfree   -- True/False tow let the period be free or not
    mask    -- specifies which points to be used for the fit. Fit is
               still evlauted at all points, but the coefficients are
               based only upon the points set to True in this array

    Returns (sfit,c,a,b) where sfit = the final fit, c,a,b are the fit coefficients, or
    (sfit,c,a,b,p) if pfree
    """

    if ye.min() <= 0.:
        raise SubsError('Y errors have at least one value <= 0')
    if mask is None:
        mask = ye > 0.

    if pfree:
        from scipy import optimize

        def func(p,x,y,ye):
            c = p[0]
            a = p[1]
            b = p[2]
            period = p[3]
            phi = 2.*np.pi*x/period
            return (y-c-a*np.cos(phi)-b*np.sin(phi))/ye

        def dfunc(p,x,y,ye):
            c = p[0]
            a = p[1]
            b = p[2]
            period = p[3]
            phi = 2.*np.pi*x/period
            wcp = np.cos(phi)/ye
            wsp = np.sin(phi)/ye
            return [ -1./ye, -wcp, -wsp, (-a*wsp+b*wcp)*phi/period]

        (sfit,c,a,b) = sinfit(x,y,ye,period,False,mask)
        pinit  = np.array((c,a,b,period))
        result = optimize.leastsq(func, pinit, (x[mask],y[mask],ye[mask]), dfunc, col_deriv=1)
        c,a,b,period = result[0]
        phi = 2.*np.pi*x/period
        return (c+a*np.cos(phi)+b*np.sin(phi),c,a,b,period)

    else:
        from scipy import linalg

        phi = 2.*np.pi*x/period
        a   = np.matrix(np.column_stack((1./ye, np.cos(phi)/ye, np.sin(phi)/ye)))
        b   = np.matrix(np.reshape(y/ye,(len(y),1)))
        (c,resid,rank,s) = linalg.lstsq(a[mask,:],b[mask,:])

        sfit = np.array(ye*np.ravel(a*c))

        return (sfit,float(c[0]),float(c[1]),float(c[2]))

def vac2air(wvac):
    """wair = vac2air(wvac): convert vacuum wavelengths (in Angstroms) to air wavelengths

    I use a formula quoted in http://www.sdss.org/dr7/products/spectra/vacwavelength.html
    which I have not however checked independently. It quotes Morton (1991, ApJS, 77, 119)
    where it appears in a less convenient form.
    """
    return wvac/(1.0002735182 + (131.4182+2.76249e8/wvac**2)/wvac**2)

def air2vac(wair):
    """wvac = air2vac(wair): convert air wavelengths (in Angstroms) to vacuum wavelengths

    This uses equation 3 quoted in Morton (1991, ApJS, 77, 119)
    """

    sigmasq = (1.e4/wair)**2

    return wair*(1.000064328 + 2.94981e-2/(146.-sigmasq) + 2.5540e-4/(41.-sigmasq))

def m2min(mf, m1):
    """
    Returns the minimum possible mass(es) of a companion star given a mass function 
    (or a numpy array of them) and a single primary mass. Works by iteration.

    mf -- mass function, i.e. K1**2*(P/2*Pi*G).
    m1 -- mass of primary
    """

    m2 = mf
    for i in range(3):
        m2 = (mf*(m1+m2)**2)**(1./3.)
    for i in range(3):
        m2 -= (m2**3-mf*(m1+m2)**2)/(3*m2**2-2*mf*(m1+m2))
    return m2

def linterp(x, y, xnew):

    """
    Returns linear interpolation of array y onto a new array xnew

    Values off ends are set to the end values.

    Arguments:

    x    -- array on which values y are defined
    y    -- values defined for every x
    xnew -- new values of x for which y values will be found by linear interpolation
    """

    if len(x) != len(y):
        raise Exception('subs.linterp: x and y have different lengths')
    if x.ndim != 1 or y.ndim != 1:
        raise Exception('subs.linterp: x and/or y do not have correct number of dimensions')

    # array such that for all j
    # x[iarr[j]-1] < xnew[j] < x[iarr[j]]

    if isinstance(xnew, np.ndarray):

        iarr = np.digitize(xnew, x)

        # linearly interpolate onto same scale as data
        ok   = (iarr > 0) & (iarr < len(x))

        x1   = xnew[ok] - x[iarr[ok]-1]
        x2   = x[iarr[ok]] - xnew[ok]
        nvec = np.empty_like(xnew)
        nvec[ok] = (x1*y[iarr[ok]]+x2*y[iarr[ok]-1])/(x1+x2)

    # deal with ends
        nvec[iarr == 0]      = y[0]
        nvec[iarr == len(x)] = y[-1]
        return nvec

    else:

        iarr = np.digitize(np.array([xnew]), x)
        if iarr[0] == 0:
            return y[0]
        elif iarr[0] == len(x):
            return y[-1]
        else:
            return ((xnew - x[iarr[0]-1])*y[iarr[0]]+(x[iarr[0]]-xnew)*y[iarr[0]-1])/(x[iarr[0]]-x[iarr[0]-1])

def centroid(xpos, fwhm, y, emission=True, e=None):
    """(xcen,{xerr}) = centroid(xpos, fwhm, y, e=None) -- computes weighted centroid of a feature in a 1D array.

    This uses the classic Young & Schneider cross-correlation with a gaussian method. Note that this routine will
    use all of the arrays passed to it. If it runs slowly you may be be able to speed it up by passing sub-sections
    of the arrays.

    Arguments:

    xpos     -- initial position, with centre of the first pixel defined as 0.0. This must initially be
                within the range 1 to len(y)-2, i.e. more than 1 from either end. The routine fails if it cannot
                locate an extremum of the cross-correlation only shifts
                by 3*fwhm from this position at the most.
    fwhm     -- the full width at half maximum in pixels of the gaussian which must be >=2.
    emission -- True for emission lines, False for absorption
    y        -- the array of values containing the feature. e are 1-sigma uncertainties. Must have at least 3 pixels.
    e        -- 1-sigma uncertainties if available. These are only used to estimate the uncertainty in the answer

    Returns:

    xcen     -- centroid
    xerr     -- 1-sigma uncertainty in centroid if e is not None

    Raises SubsError exceptions if there are problems.
    """

    from scipy.optimize import brentq

    def comperr(xd, y, e, sigma):
        """Computes uncertainty in final position from median"""
        xdsq = np.power(xd,2)
        vc   = np.sum(xdsq*np.exp(-xdsq)*np.power(e,2))
        ddc  = np.sum(y*(xdsq-1)*np.exp(-xdsq/2))
        return sigma*m.sqrt(vc)/abs(ddc)

    def dcorr(xd, y):
        """Returns derivative of correlation with gaussian

        xd -- xd offsets from centre of gaussian, scaled by RMS
        y  -- the corresponding y values
        """
        return np.sum(y*xd*np.exp(-np.power(xd,2)/2))

    # check inputs
    if fwhm < 2:
        raise SubsError('centroid: fwhm = ' + str(fwhm) + ' < 2')
    if xpos < 1 or xpos > len(y)-2:
        raise SubsError('centroid: xpos = ' + str(xpos) + ' is out of allowed range 1 to ' + str(len(y)-2))
    if len(y) < 3:
        raise SubsError('median: length of y array < 3')
    if e is not None and len(e) != len(y):
        raise SubsError('median: length of e array does not match the y array')

    sigma = fwhm/2.3548
    x     = np.arange(len(y),dtype=float)
    xd    = (x - xpos)/sigma
    dc1   = dcorr(xd, y)
    if dc1 == 0:
        if e is not None:
            return (xpos,comperr(xd, y, e, sigma))
        return xpos

    # Try to bracket the peak by looking for closest switch in sign
    # of the derivative of the x-corr.

    x1  = xpos
    x2  = xpos
    dc2 = dc1
    found_switch = False

    # move in direction defined by sign of gradient.
    if (dc1 > 0.0 and emission) or (dc1 < 0.0 and not emission):
        cshift =  0.25
    else:
        cshift = -0.25

    shift  = cshift
    while abs(shift) < 3.1:
        x2  = xpos + fwhm*shift
        xd  = (x - x2)/sigma
        dc2 = dcorr(xd, y)
        if dc2 == 0:
            if e != None:
                return (x2,comperr(xd, y, e, sigma))
            else:
                return x2
        if (dc1 > 0 and dc2 < 0) or (dc1 < 0 and dc2 > 0):
            found_switch = True
            break
        shift += cshift

    if not found_switch:
        raise SubsError('centroid: could not find peak within 3*fwhm of initial position')

    # reorder if necessary
    if x1 > x2:
        (x1,x2)   = (x2,x1)
        (dc1,dc2) = (dc2,dc1)

    if x1 < 0.0 or x2 < 1.0 or x1 > len(y)-2 or x2 > len(y)-1:
        raise SubsError('centroid: bracketting limits = ' + str(x1) + ', ' + str(x2) + ' out of range.')

    # now find root using Brent's method

    def bfunc(xpos, x, y, sigma):
        """Function for the root finder brentq"""
        return dcorr((x - xpos)/sigma, y)

    (xm,r) = brentq(bfunc, x1, x2, args=(x,y,sigma), full_output=True)
    if not r.converged:
        raise SubsError('centroid: brentq did not converge')

    if e is not None:
        xd = (x - xm)/sigma
        return (xm, comperr(xd, y, e, sigma))
    else:
        return xm

class iStr(str):
    """
    Case insensitive strings class.
    Performs like str except comparisons are case insensitive.
    """

    def __init__(self, strMe):
        str.__init__(self, strMe)
        self.__lowerCaseMe = strMe.lower()

    def __repr__(self):
        return "iStr(%s)" % str.__repr__(self)

    def __eq__(self, other):
        return self.__lowerCaseMe == other.lower()

    def __lt__(self, other):
        return self.__lowerCaseMe < other.lower()

    def __le__(self, other):
        return self.__lowerCaseMe <= other.lower()

    def __gt__(self, other):
        return self.__lowerCaseMe > other.lower()

    def __ne__(self, other):
        return self.__lowerCaseMe != other.lower()

    def __ge__(self, other):
        return self.__lowerCaseMe >= other.lower()

    def __cmp__(self, other):
        return cmp(self.__lowerCaseMe, other.lower())

    def __hash__(self):
        return hash(self.__lowerCaseMe)

    def __contains__(self, other):
        return other.lower() in self.__lowerCaseMe

    def count(self, other, *args):
        return str.count(self.__lowerCaseMe, other.lower(), *args)

    def endswith(self, other, *args):
        return str.endswith(self.__lowerCaseMe, other.lower(), *args)

    def find(self, other, *args):
        return str.find(self.__lowerCaseMe, other.lower(), *args)

    def index(self, other, *args):
        return str.index(self.__lowerCaseMe, other.lower(), *args)

    def lower(self):   # Courtesy Duncan Booth
        return self.__lowerCaseMe

    def rfind(self, other, *args):
        return str.rfind(self.__lowerCaseMe, other.lower(), *args)

    def rindex(self, other, *args):
        return str.rindex(self.__lowerCaseMe, other.lower(), *args)

    def startswith(self, other, *args):
        return str.startswith(self.__lowerCaseMe, other.lower(), *args)

class Odict(dict):
    """
    A dictionary which stores a key order which it uses when printing
    """

    ninset  = 0
    nincr   = 5
    nlength = 20

    def __init__(self, dct = None):
        if dct == None:
            self._keys = []
            dict.__init__(self, {})
        else:
            self._keys = dct.keys()
            dict.__init__(self, dct)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        dict.__setitem__(self, key, item)
        # 'try' needed to avoid error with pickling with protocol = 2 
        try:
            if key not in self._keys: self._keys.append(key)
        except AttributeError:
            self._keys = [key]

    def __str__(self):
        st    = ''
        inset = ' '*Odict.ninset
        Odict.ninset += Odict.nincr
        for (key,val) in self.iteritems():
            if isinstance(val, Odict):
                st += ('%s%-' + str(Odict.nlength) + 's\n') % (inset,key)
            else:
                st += ('%s%-' + str(Odict.nlength) + 's ') % (inset,key)
            st += str(val) + '\n'
            
        Odict.ninset -= Odict.nincr
        return st

    def clear(self):
        dict.clear(self)
        self._keys = []

    def copy(self):
        newInstance = Odict()
        newInstance.update(self)
        return newInstance

    def items(self):
        return zip(self._keys, self.values())

    def keys(self):
        return self._keys

    def popitem(self):
        try:
            key = self._keys[-1]
        except IndexError:
            raise KeyError('dictionary is empty')

        val = self[key]
        del self[key]

        return (key, val)

    def setdefault(self, key, failobj = None):
        if key not in self._keys: self._keys.append(key)
        return dict.setdefault(self, key, failobj)

    def update(self, dct):
        for (key,val) in dct.items():
            self.__setitem__(key,val)

    def values(self):
        return map(self.get, self._keys)

    def iteritems(self):
        return _odict_iteritem(self)

class _odict_iteritem:
    """Iterator class for Odict"""

    def __init__(self, od):
        self.last = 0  
        self.odct = od

    def __iter__(self):
        return self

    def next(self):
        if self.last == len(self.odct._keys):
            raise StopIteration
        key = self.odct._keys[self.last]
        val = self.odct[key]
        self.last += 1
        return (key,val)

def str2radec(position, crude=False):
    """
    ra,dec,system = str2radec(position, crude) -- translates an astronomical
    coordinate string to double precision RA and Dec.

    'ra' is the RA in decimal hours; 'dec' is the declination in degrees;
    'system' is one of 'ICRS', 'B1950', 'J2000'. Entering coordinates is an
    error-prone and often irritating chore.  This routine is designed to make
    it as safe as possible while supporting a couple of common formats.

    Here are example input formats, both good and bad:

    12 34 56.1 -23 12 12.1     -- OK. 'system' will default to ICRS
    234.5 34.2  B1950          -- OK. RA, Dec entered in decimal degrees, B1950 to override default ICRS
    11 02 12.1 -23.2 J2000     -- NOK. Cannot mix HMS/DMS with decimals
    11 02 12.1 -23 01 12 J4000 -- NOK. Only 'ICRS', 'B1950', 'J2000' allowed at end.
    1.02323e2 -32.5            -- NOK. Scientific notation not supported.
    11 02 12.1 23 01 12        -- NOK. In HMS mode, the sign of the declination must be supplied
    11 02 12.1 -00 00 02       -- OK. - sign will be picked up
    25 01 61.2 +90 61 78       -- NOK. various items out of range.
    12:32:02.4 -12:11 10.2     -- OK. Colon separators allowed.

    If crude=True low precision coordinates such as "03 50.0 +13 25" will be
    recognised as well

    A SubsError is raised on failure
    """
    
    # Try out three types of match
    m = re.search(r'^\s*(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)\s+([\+-])(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)(?:\s+(\w+))?\s*$', position)
    if m:
        (rah,ram,ras,decsign,decd,decm,decs,system) = m.groups()
        rah  = int(rah)
        ram  = int(ram)
        ras  = float(ras)
        decd = int(decd)
        decm = int(decm)
        decs = float(decs)
        if (rah > 23 and ram > 0 and ras > 0.) or ram > 60 or ras > 60. \
                or (decd > 89 and decm > 0 and decs > 0.) or decm > 60 or decs > 60.:
            raise SubsError('trm.subs.str2radec: one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

        if not system: system = 'ICRS'
        if system != 'ICRS' and system != 'J2000' and system != 'B1950':
            raise SubsError('trm.subs.str2radec: astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + ' is not recognised.')

        ra  = rah + ram/60. + ras/3600.
        dec = decd + decm/60. + decs/3600.
        if decsign == '-': dec = -dec
        return (ra,dec,system)

    # No arcseconds of dec as sometimes is the case with coords from simbad
    m = re.search(r'^\s*(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)\s+([\+-])(\d{1,2})(?:\:|\s+)(\d{1,2}\.\d*)(?:\s+(\w+))?\s*$', position)
    if m:
        (rah,ram,ras,decsign,decd,decm,system) = m.groups()
        rah  = int(rah)
        ram  = int(ram)
        ras  = float(ras)
        decd = int(decd)
        decm = float(decm)
        if (rah > 23 and ram > 0 and ras > 0.) or ram > 60 or ras > 60. \
                or (decd > 89 and decm > 0.) or decm > 60.:
            raise SubsError('trm.subs.str2radec: one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

        if not system: system = 'ICRS'
        if system != 'ICRS' and system != 'J2000' and system != 'B1950':
            raise SubsError('trm.subs.str2radec: astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + ' is not recognised.')

        ra  = rah + ram/60. + ras/3600.
        dec = decd + decm/60.
        if decsign == '-': dec = -dec
        return (ra,dec,system)

    # Decimal entries
    m = re.search(r'^\s*(\d{1,3}(?:\.\d*)?)\s+([+-]?\d{1,2}(?:\.\d*)?)\s+(\w+)?\s*$', position)
    if m:
        print ('matched decimal entries')
        (rad,dec,system) = m.groups()
        ra   = float(rad)/15.
        dec  = float(dec)
        if ra >= 24. or dec < -90. or dec > 90.:
            raise SubsError('trm.subs.str2radec: one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

        if not system: system = 'ICRS'
        if system != 'ICRS' and system != 'J2000' and system != 'B1950':
            raise SubsError('trm.subs.str2radec: astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + ' is not recognised.')

        return (ra,dec,system)

    if crude:
        # crude coords of form "12 34.5 +12 34"
        m = re.search(r'^\s*(\d{1,2})(?:\:|\s+)(\d{1,2}(?:\.\d*)?)\s+([\+-])(\d{1,2})(?:\:|\s+)(\d{1,2})(?:\s+(\w+))?\s*$', position)
        if m:
            (rah,ram,decsign,decd,decm,system) = m.groups()
            rah  = int(rah)
            ram  = float(ram)
            decd = int(decd)
            decm = int(decm)
            if (rah > 23 and ram > 0.) or ram > 60. or (decd > 89 and decm > 0) or decm > 60:
                raise SubsError('trm.subs.str2radec: one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

            if not system: system = 'ICRS'
            if system != 'ICRS' and system != 'J2000' and system != 'B1950':
                raise SubsError('trm.subs.str2radec: astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + \
                                    ' is not recognised.')

            ra  = rah + ram/60.
            dec = decd + decm/60.
            if decsign == '-': dec = -dec
            return (ra,dec,system)

        # even cruder: "12 34 +12.4"
        m = re.search(r'^\s*(\d{1,2})(?:\:|\s+)(\d{1,2})\s+([\+-])(\d{1,2}(?:\.\d*)?)(?:\s+(\w+))?\s*$', position)
        if m:
            rah,ram,decsign,decd,system = m.groups()
            rah  = int(rah)
            ram  = int(ram)
            decd = float(decd)
            if (rah > 23 and ram > 0) or ram > 60 or decd > 90.:
                raise SubsError('trm.subs.str2radec: one or more of the entries in the astronomical coordinates "' + position + '" is out of range')

            if not system: system = 'ICRS'
            if system != 'ICRS' and system != 'J2000' and system != 'B1950':
                raise SubsError('trm.subs.str2radec: astronomical coordinate system must be one of ICRS, B1950, J2000; ' + system + \
                                    ' is not recognised.')

            ra  = rah + ram/60.
            dec = decd
            if decsign == '-': dec = -dec
            return (ra,dec,system)

    raise SubsError('trm.subs.str2radec: could not interpret "' + position + '" as astronomical coordinates')

def observatory(telescope=None):
    """
    (tel,obs,longitude,latitude,height) = observatory(telescope=None)

    telescope = None implies interactive selection of observatory data. The
    user will be prompted to choose preloaded observatory data.

    Returns:

    tel        -- telescope name (which can be used to force selection in future)
    obs        -- observatory
    longitude  -- East-positive longitude in degrees
    latitude   -- latitude in degrees
    height     -- height above sea-level, metres
    """
    
    OBS = [ \
            ['WHT','La Palma','17 52 53.9','W','28 45 38.3','N',2332.], \
            ['INT','La Palma','17 52 40.', 'W','28 45 43.','N',2336.], \
            ['LT','La Palma','17 52 45.', 'W','28 45 44.','N',2336.], \
            ['SAAO','Sutherland','20 48 38.5','E','32 22 46.','S', 1798.], \
            ['NTT', 'La Scilla','70 44 00.','W','29 15 00.','S', 2400.], \
            ['VLT', 'Paranal','70 24 9.9','W','24 37 30.3','S',2635.], \
            ['SRT', 'University of Warwick','01 28 0.0','W','52 25 00','N',20.], \
            ['Greenwich', 'London','00 00 0.0','W','51 29 00','N',20.], \
            ['SDSS', 'Apache Point','105 49 13','W','32 46 49','N',2788.], \
            ['SOAR', 'Cerro Pachon','70 44 01.11','W','30 14 16.41','S',2713.], \
            ['TNT','Doi Inthanon','98 28 56','E','18 34 26','N',2457.], \
            ['ATCA','Narrabri','149 32 56.327','E','30 18 52.048', 'S', 209.3], \
            ['200-in','Palomar','116 51 54','W','33 21 23', 'N', 1713.], \
            ]

    if telescope is not None:
        for obs in OBS:
            if obs[0] == telescope:
                break
        else:
            raise SubsError('observatory: did not recognize telescope = ' + telescope)
        select = obs
    else:
        print ('Observatories: ')
        for n in range(len(OBS)):
            (tel,obs,lng,eorw,lat,nors,height) = OBS[n]
            if nors == 'S':
                lats = '-'
            else:
                lats = '+'
            print ('%2d) telescope = %-5s  observatory = %-12s longitude = %-11s%s  latitude = %s%-11s  height = %6.1f' % (n+1,tel,obs,lng,eorw,lats,lat,height))

        try:
            n = int(raw_input('Select number of observatory: '))
        except:
            raise SubsError('Failed to interpret selection as an integer')
        if n < 1 or n > len(OBS):
            raise SubsError('Observatory number = ' + str(n) + ' is out of range.')
        select = OBS[n-1]

    # format for output
    (tel,obs,lng,eorw,lat,nors,height) = select
    (ld,lm,ls) = lng.split()
    longitude = float(ld) + float(lm)/60. + float(ls)/3600.
    if eorw == 'W':
        longitude = -longitude
    elif eorw != 'E':
        raise SubsError('Pre-stored observatory East-or-West string = ' + eorw + ' is not one of E or W')
    
    (ld,lm,ls) = lat.split()
    latitude = float(ld) + float(lm)/60. + float(ls)/3600.
    if nors == 'S':
        latitude = -latitude
    elif nors != 'N':
        raise SubsError('Pre-stored observatory North-or-South string = ' + nors + ' is not one of N or S')

    return (tel,obs,longitude,latitude,float(height))

def d2hms(hour, precision=2, sep=':', dp=3, sign=False, nearest=False):
    """hms = d2hms(hour, precision=2) -- returns HH MM SS.SS time string from a decimal hour

    hour      -- decimal hour
    precision -- 0 gives hh:mm, 1 gives hh:mm:ss, 2 gives hh:mm:ss.sss
    sep       -- separator
    dp        -- number of decimal places in the case of precision = 2
    sign      -- True if you want to make sure a '+' is appended at the start
    nearest   -- If you want only limited precision this controls whether one
                 truncates or takes the closest value
    """
    if nearest:
        add = 0.5
    else:
        add = 0.00001
    h = int(abs(hour))
    if precision == 0:
        m = int(60.*(abs(hour)-h)+add)
        st = '%02.2d%s%02.2d' % (h,sep,m)
    else:
        m = int(60.*(abs(hour)-h))
    if precision == 1:
        s  = int(3600.*(abs(hour)-h-m/60.)+add)
        st = '%02.2d%s%02.2d%s%02.2d' % (h,sep,m,sep,s)
    else:
        s = 3600.*(abs(hour)-h-m/60.)
        if dp:
            st = ('%02.2d%s%02.2d%s%0' + str(dp+3) + '.' + str(dp) +'f') % (h,sep,m,sep,s)
        else:
            st = ('%02.2d%s%02.2d%s%02d') % (h,sep,m,sep,s)

    if sign:
        if hour < 0:
            st = '-' + st
        else:
            st = '+' + st
    return st

def hms2d(utc):
    """hour = hms2d(utc) -- returns decimal hour from and HH MM SS.SS time string. Can also cope
    with initial + or - sign.

    utc    -- time in HH MM SS.SSS format. HH MM and SS must be separated by exactly 1 character.
    """
    if utc[:1] == '-':
        sign = -1.
        off = 1
    elif utc[:1] == '+':
        sign = +1.
        off = 1
    else:
        sign = +1.
        off = 0
    h = int(utc[off:2+off])
    m = int(utc[3+off:5+off])
    s = float(utc[6+off:])
    return sign*(h + m/60. + s/3600.)

# Exception class
class SubsError(Exception):
    """For throwing exceptions from the subs module"""
    def __init__(self, value):
        self.value = value
            
    def __str__(self):
        return repr(self.value)

class Poly(np.ndarray):
    """
    Class for storing a polynomial which is an numpy.ndarray of
    coefficients 'a', an offset 'x0' and scale 'xscale' so that
    y = sum_{i=0}^{n-1} a[i]*((x-x0)/xscale)**i
    """

    def __new__(subtype, n, x0, xscale):
        """
        Need to supply the number of coefficients, the offset and scale
        """
        # Stuff below cribbed off web without any understanding on my part

        # Make sure we are working with an array, and copy the data if requested
        subarr = np.array((n), np.float)
    
        # Transform 'subarr' from an ndarray to our new subclass.
        subarr = subarr.view(subtype)
   
        subarr.x0     = x0
        subarr.xscale = xscale
        return subarr

    def __array_finalize__(self,obj):
           # We use the getattr method to set a default if 'obj' doesn't have the 'info' attribute
           self.info = getattr(obj, 'info', {})
           # We could have checked first whether self.info was already defined:
           #if not hasattr(self, 'info'):
           #    self.info = getattr(obj, 'info', {})
  
    def __repr__(self):
        desc="""
array(data=
         %(data)s,
  x0    =%(x0)s,
  xscale=%(xscale)s)
  """
        return desc % {'data': str(self), 'x0': self.x0, 'xscale': self.xscale }

    def __str__(self):
        return 'coeffs=\n%s\nx0=%f, xscale=%f' % (np.ndarray(self).__str__(), self.x0, self.xscale)

class Vec3(object):
    """
    A simple 3D vector class.

    Examples:

    r = Vec3(1,2,3) # sets r = (1,2,3)
    r = a + b       # adds two Vec3s to make another.
    x = r[0]        # gets the x ordinate
    r = 2.*r        # multiply all ordinates by 2
    print r.norm()  # Print length of a Vec3

    See also dot and cross.

    NB. This does not have many checks, and you should try to get the
    code right.
    """

    def __init__(self, *args):
        """
        If no arguments supplied, the vector is set to (0,0,0).

        If 3 arguments are supplied, they specify (x,y,z)

        """
        if len(args) == 0:
            self.x = self.y = self.z = 0
        elif len(args) == 3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        else:
            raise SubsError('Vec3(): can only have 0 or 3 arguments')

    def __repr__(self):
        return '(%f,%f,%f)' % (self.x,self.y,self.z)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __add__(self, other):
        temp = copy.copy(self)
        temp += other
        return temp

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __sub__(self, other):
        temp = copy.copy(self)
        temp -= other
        return temp

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        self.z *= other
        return self

    def __mul__(self, other):
        temp = copy.copy(self)
        temp *= other
        return temp

    def __rmul__(self, other):
        temp = copy.copy(self)
        temp *= other
        return temp

    def __itruediv__(self, other):
        self.x /= other
        self.y /= other
        self.z /= other
        return self

    def __truediv__(self, other):
        temp = copy.copy(self)
        temp /= other
        return temp

    def __neg__(self):
        temp = copy.copy(self)
        temp *= -1
        return temp

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        else:
            raise SubsError('Index out of range 0:3')

    def sqnorm(self):
        """
        Returns Euclidean length squared of the vector
        """
        return self.x**2 + self.y**2 + self.z**2

    def norm(self):
        """
        Returns Euclidean length of the vector
        """
        return m.sqrt(self.sqnorm())

    def unit(self):
        """
        Returns vector as a unit vector in same direction
        """
        leng = m.sqrt(self.sqnorm())
        if leng == 0.:
            raise SubsError("Vec3.unitv: zero length vector")
        return self / leng

    def dot(self, other):
        """
        Returns the scalar or dot product of self with other, a 3-vector
        """
        return self.x*other.x+self.y*other.y+self.z*other.z

    def cross(self, other):
        """
        Computes the vector or cross product of self with other, a 3-vector
        """
        return Vec3(self.y*other.z-self.z*other.y,
                    self.z*other.x-self.x*other.z,
                    self.x*other.y-self.y*other.x)


def dot(a: Vec3, b: Vec3) -> float:
    """
    Computes the scalar or dot product of two 3-vectors
    """
    assert isinstance(a, Vec3) and isinstance(b, Vec3),\
        "Arguments must be of type Vec3"
    return a.dot(b)


def cross(a: Vec3, b: Vec3) -> Vec3:
    """
    Computes the vector or cross product of two 3-vectors
    """
    assert isinstance(a, Vec3) and isinstance(b, Vec3),\
        "Arguments must be of type Vec3"
    return a.cross(b)


def mr_wd_eggleton(mass):
    """
    Returns radius of a white dwarf in solar radii given its mass in solar masses
    This is according to a formula of Eggleton's quoted in Verbunt and Rappaport (1988)
    """

    if isinstance(mass, np.ndarray):
        if (mass <= 0).any() or (mass >= 1.44).any():
            raise SubsError("mr_wd_eggleton: mass array contains at least element outside range 0 to 1.44")
    elif mass <= 0 or mass >= 1.44:
        raise SubsError("mr_wd_eggleton: mass = " + str(mass) + " out of range.");

    fac1 = np.power(mass/1.44,2./3.)
    fac2 = 0.00057/mass;
    return 0.0114*np.sqrt(1./fac1-fac1)*np.power(1.+3.5*np.power(fac2,2./3.)+fac2,-2./3.)


def orbital_separation(m1, m2, period):
    """
    Returns orbital semi-major axis in solar radii given masses in solar masses
    and an orbital period in days, i.e. Kepler's third law. Works with scalar
    or numpy array inputs.
    """

    return (G*MSUN*(m1+m2)*(DAY*period/(2.*m.pi))**2)**(1./3.)/RSUN


def orbital_period(m1, m2, a):
    """
    Returns orbital period in days given masses in solar masses
    and the orbital semi-major axis in solar radii, i.e. Kepler's third law
    """

    return 2.*m.pi/m.sqrt(G*MSUN*(m1+m2)/m.pow(RSUN*a,3))/DAY


def jdotgr(m1, m2, a):
    """
    Returns Jdot/J due to GR. m1, m2, a in solar units, result
    has SI units
    """
    return -32./5.*(G*MSUN/C/RSUN)**3/C**2/RSUN*m1*m2*(m1+m2)/a**4


def jorb(m1, m2, a):
    """
    Returns orbital angular momentum, (SI), given m1, m2, a in solar units
    """
    return MSUN*m1*m2*m.sqrt(G*MSUN*RSUN*a/(m1+m2))


def jdotrvj(m2, r2, gamma, porb):
    """
    Returns Rappaport, Verbut, Joss ang mom loss rate (SI)
    m2, r2 solar, porb days.
    """
    return -3.8e-26*m2*MSUN*RSUN**4*r2**gamma*(2.*m.pi/porb/DAY)**3


def rlobe_eggleton(q):
    """
    Returns scaled Roche lobe radius of star 2 according to Eggleton's (1983) formula given
    q = M2/M1. Works for scalar or array values.
    """
    return 0.49*q**(2./3.)/(0.6*q**(2./3.)+np.log(1.+ q**(1./3.)))


def date2dmy(dstring):
    """Returns (year,month,day) as numbers given input dates such as '11 Apr 2004' or 'Apr 11 2004'"""

    rec = re.compile('^\s*(\d{1,2})\s+([a-zA-Z]{3})\s+(\d\d\d\d)\s*$')
    m   = rec.match(dstring)

    if m:
        day   = int(m.group(1))
        year  = int(m.group(3))
        month = month2int(m.group(2))
    else:
        rec = re.compile('^\s*([a-zA-Z]{3})\s+(\d\d)\s*(\d\d\d\d)\s*$')
        m   = rec.match(dstring)
        if m:
            day   = int(m.group(2))
            year  = int(m.group(3))
            month = month2int(m.group(1))
        else:
            raise SubsError('Date = ' + dstring + ' not in a recognised format.')
    return (year, month, day)


def month2int(name):
    """
    i = month2int(name) -- returns integer 1 to 12 given first three
    characters of a month name, case insensitively
    """

    names = {'jan' : 1, 'feb' : 2, 'mar' : 3, 'apr' : 4, 'may' : 5, 'jun' : 6, 'jul' : 7, 'aug' : 8, 'sep' : 9, 'oct' : 10, 'nov' : 11, 'dec' : 12}
    lname = name.lower()[:4]
    if lname in names:
        month = names[lname]
    else:
        raise SubsError('subs.month2int: month = ' + str(lname) + ' not recognised')
    return month


def int2month(month):
    """name = int2month(month) -- returns 3 letter month name given integer 1 to 12"""
    names = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
    if month < 1 or month > 12:
        raise SubsError('subs.int2month: month = ' + str(month) + ' out of range 1 to 12')
    return names[month-1]


def inpoly(poly, x, y):
    """Determines whether a point (x,y) is inside a polygon defined by
    the numpy 2D arrays poly which lists the vertices of the polygon
    such that poly[0,0] is the first x coordinate, poly[0,1] is the
    first y coordinate, poly[1,0] is the second x coordinate etc. It
    returns True if the point is inside the poly.

    See this web page for lengthy description of the method which I
    translated from the short piece of C-code listed there:

    http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

    I have yet to test whether or Python does the same as C in not
    bothering with the second test in an 'and' conditional if the
    first fails.

    See also: minpoly
    """
    #TODO: Likely this is better bound into some C for speed, or use shapley

    outside = True
    j = poly.shape[0] - 1
    for i in range(poly.shape[0]):
        if ( (poly[i,1] > y) != (poly[j,1] > y) ) and \
                (x < (poly[j,0]-poly[i,0]) * (y-poly[i,1]) / (poly[j,1]-poly[i,1]) + poly[i,0]):
            outside = not outside
        j = i
    return not outside


def minpoly(poly, x, y):
    """Determines whether points in numpy arrays (x,y) are inside a
    polygon defined by the numpy 2D arrays poly which lists the
    vertices of the polygon such that poly[0,0] is the first x
    coordinate, poly[0,1] is the first y coordinate, poly[1,0] is the
    second x coordinate etc. It returns True if the point is inside
    the poly.

    See this web page for lengthy description of the method which I 
    translated from the short piece of C-code listed there:

    http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

    I have yet to test whether or Python does the same as C in not
    bothering with the second test in an 'and' conditional if the
    first fails.

    See also: inpoly
    """
    #TODO: Likely this is better bound into some C for speed, or use shapley


    outside = np.ones_like(x,dtype=bool)
    j = poly.shape[0] - 1
    for i in range(poly.shape[0]):
        switch =  ((poly[i,1] > y) != (poly[j,1] > y)) & \
                  (x < (poly[j,0]-poly[i,0]) * (y-poly[i,1]) / (poly[j,1]-poly[i,1]) + poly[i,0])
        outside[switch] = ~outside[switch]
        j = i
    return ~outside


def sigma_reject(data, thresh, onebyone):
    """
    This routine computes the raw mean and RMS and the mean and RMS after
    rejection of points more than thresh*RMS from the mean. The rejection 
    can be carried out one bad pixel value at a time (onebyone=True) or all bad pixels
    at a time. The first method (which can involve mutiple pixels of the sam evalue) 
    is slow but robust whereas the second is fast but potentially dangerous.

    Arguments:

    data     -- a numpy array
    thresh   -- threshold number of RMS e.g. 3.
    onebyone -- True for slow single worst pixel reject methos, False for fast, 
                all bad pixels reject method.

    Returns (rawmean,rawrms,clipmean,cliprms,nreject,ncycles)
    """
    
    if not isinstance(data, np.ndarray):
        raise SubsError("sigma_reject: first input is not a numpy.ndarray")

    ncycles = 0
    nrej = 1
    mask = (data == data)
    nrold = 0

    while nrej:
        mean = data[mask].mean()
        rms  = data[mask].std()
    
        if ncycles == 0:
            rawmean = mean
            rawrms  = rms

        # Update mask
        diff = np.abs(data - mean)
        if onebyone:
            if diff[mask].max() > thresh*rms:
                mask = mask & (diff < diff[mask].max()) 
        else:
            mask = mask & (diff <= thresh*rms)

        # Check on number of rejects
        nrnew = np.size(data) - np.size(data[mask])
        nrej  = nrnew - nrold
        nrold = nrnew 
        ncycles += 1

    return (rawmean,rawrms,mean,rms,nrnew,ncycles)


def find_exec(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or raises a 
    SubsError if not found.
    """
    if path is None:
        path = os.environ['PATH']

    paths = path.split(os.pathsep)

    for p in paths:
        f = os.path.join(p, executable)
        if os.path.isfile(f):
            return f

    raise SubsError('Failed to find command = ' + executable + ' in path = ' + path)

def zeta_wd_eggleton(mass):
    """
    Returns d log(R)/d log(M) of a white dwarf according to a formula of Eggleton's quoted in Verbunt and Rappaport (1988)
    """

    if isinstance(mass, np.ndarray):
        if (mass <= 0).any() or (mass >= 1.44).any():
            raise SubsError("zeta_wd_eggleton: mass array contains at least element outside range 0 to 1.44")
    elif mass <= 0 or mass >= 1.44:
        raise SubsError("zeta_wd_eggleton: mass = " + str(mass) + " out of range.");

    fac1 = np.power(mass/1.44,4./3.)
    fac2 = 0.00057/mass;
    fac3 = np.power(fac2,2./3.);
    return -(1.+fac1)/(1.-fac1)/3.+ 2.*(7.*fac3/3.+fac2)/(1.+3.5*fac3+fac2)/3.


def pts2cont(x, y, x1, x2, nx, y1, y2, ny, xb=0., yb=0., frac=None):
    """
    Converts a series of points into a binned image and optionally
    blurring them. Can also compute contour levels marking particular 
    fractions of points, e.g. 0.68, 0.95
    
    Arguments:

    x, y       -- arrays / lists of x and y values of points
    x1, x2, nx -- range and number of pixels for the X axis. 
                  Leftmost pixel starts at x1, rightmost ends at x2.
    y1, y2, ny -- range and number of pixels for the Y axis
    xb, yb     -- amount of blurring in X and Y (RMS in pixels).
    frac       -- Array / list of fractions to determine equivalent contour levels.

    Returns: (shist, cont) where shist is a smoothed histogrammed version of
    the points and comt is an array of contour levels corresponding to frac.
    """

    import scipy.signal

    # Define edges, accumulate histogram and transpose it
    xe = np.linspace(x1, x2, nx+1)
    ye = np.linspace(y1, y2, ny+1)
    hist,xe,ye = np.histogram2d(x,y,(xe,ye))
    hist = hist.transpose()

    if xb > 0. or yb > 0.:
        # compute blurring array
        nxwid = 2*int(3.*xb)+1
        nywid = 2*int(3.*yb)+1
        (x,y) = np.mgrid[0:nywid,0:nxwid]
        if xb > 0. and yb > 0.:
            fuzz = np.exp(-((x-nxwid//2)/xb)**2/2.-((y-nywid//2)/yb)**2/2.)
        elif xb > 0.:
            fuzz = np.exp(-((x-nxwid//2)/xb)**2/2.)
        else:
            fuzz = np.exp(-((y-nywid//2)/yb)**2/2.)

        # apply blurr
        shist = scipy.signal.convolve2d(hist, fuzz, mode='same')
    else:
        nxwid = 1
        nywid = 1
        shist = hist

    def _levdet(shist, hist, clev):
        """
        Determine level which captures a fraction clev of the data.
        Supply smoothed and raw versions of the histogram. The level
        returned applies to the smoothed version.
        """
        lower = 0.
        upper = 1.1*shist.max()
        nlo   = 0
        nhi   = hist.sum()
        naim  = int(nhi*(1.-clev)+0.5)

        nlold = nlo + 1
        nhold = nhi + 1
        while nlo != nlold or nhi != nhold:
            nlold  = nlo
            nhold  = nhi
            guess  = (lower+upper)/2.
            mask   = shist < guess
            nbelow = hist[mask].sum()
            if nbelow > naim:
                nhi   = nbelow
                upper = guess
            else:
                nlo   = nbelow
                lower = guess
        return guess

    if frac is not None:
        cont = []
        for f in frac:
            cont.append(_levdet(shist, hist, f))
    else:
        cont = None
    return (shist, cont)

#def fastpure(x, y, e, nelem, fmax, nfreq, freq, pgram):
#    """
#    Implementation of Press & Rybicki using scipy rather than
#    home-cooked FFT routines.
#    """
#
#    import numpy as np
#    from scipy import fftpack
#
#    xmin = x[e > 0].min()
#    xmax = x[e > 0].max()
#    ofac = nfreq/(xdif*fmax)
#
#    MACC    = 4
#    nfreqt  = 2*MACC*nfreq
#
#    # First power of 2 > nfreqt
#    ndim = 2**int(np.ceil(np.log(nfreqt)/np.log(2.)))
#
#    fac   = ndim*(fmax/nfreq)
#
#    def spread(y, yy, n, x, m):
#        """
#        Spreads an array 

def scargle(x, y, e, f1, f2, nf):
    """
    Compute weighted floating mean Lomb-Scargle periodogram. This computes
    it directly. You should be considering the Press-Rybicki alternative if
    you have lots of points and are not miles away from regular spacing.

    x   -- input X values
    y   -- input Y values
    e   -- input Y uncertainies (only those > 0 will be used)
    f1  -- initial frequency, cycles per unit x
    f2  -- initial frequency, cycles per unit x
    nf  -- number of frequencies (>1)

    Returns (f,p) arrays of frequencies and powers

    Raises SubsError s in nf < 2
    """

    if nf < 2:
        raise SubsError("scargle: nf < 2")

    # Use Clenshaw's recurrence to speed evaluation of sine and cosines
    # initialise them first

    phi   = 2.*np.pi*x[e > 0]
    theta = (f2-f1)/(nf-1)*phi
    alpha = -2.*np.sin(0.5*theta)**2
    beta  =  np.sin(theta)
    theta = f1*phi
    czero = np.cos(theta)
    szero = np.sin(theta)
    
    w  = 1./e[e>0]**2

    f = np.linspace(f1,f2,nf)
    p = np.empty((nf,))
    for n in xrange(nf):

        cc = (w*czero**2).sum()
        cs = (w*czero*szero).sum()
        ss = (w*szero**2).sum()
        cy = (w*czero*y[e>0]).sum()        
        sy = (w*szero*y[e>0]).sum()
        
        p[n] = (ss*cy*cy-2.*cs*sy*cy+cc*sy*sy)/(cc*ss-cs*cs)/2.
        
        # recurrence update
        tmp    = np.copy(czero)
        czero += alpha*czero-beta*szero
        szero += alpha*szero+beta*tmp
        
    return (f,p)


def convex_hull(points):
    """
    Returns the points on the convex hull of a 2D set of points in 
    counter-clockwise order. "points" is an Nx2 list of 2D 
    points, which can be obtained by zip-ing an x and y array 
    together. An Mx2 list od 2D points will be returned by this 
    routine
    """

    TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

    def _turn(p, q, r):
        """Returns -1, 0, 1 if p,q,r forms a right, straight, or left turn."""
        return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)

    def _dist(p, q):
        """Returns the squared Euclidean distance between p and q."""
        dx, dy = q[0] - p[0], q[1] - p[1]
        return dx * dx + dy * dy

    def _next_hull_pt(points, p):
        """Returns the next point on the convex hull in CCW from p."""
        q = p
        for r in points:
            t = _turn(p, q, r)
            if t == TURN_RIGHT or t == TURN_NONE and _dist(p, r) > _dist(p, q):
                q = r
        return q

    hull = [min(points)]
    for p in hull:
        q = _next_hull_pt(points, p)
        if q != hull[0]:
            hull.append(q)
    return hull

def expint2(x):
    """
    Compute E_2(x), exponential integral. x must be > 0
    """
    A = np.array((-0.57721566, 0.99999193, -0.24991055, 0.05519968,\
                       -0.00976004, 0.00107857, 8.5733287401, 18.0590169730, 8.6347608925,\
                       0.2677737343, 9.5733223454, 25.6329561486, 21.0996530827, 3.9584969228))


    ex  = np.exp(-x)
    if isinstance(x, float) or isinstance(x, int):
        if x == 0.:
            return 1.
        if x <= 1.:
            e1 = A[0]+x*(A[1]+x*(A[2]+x*(A[3]+x*(A[4]+x*A[5]))))-np.log(x)
        else:
            e1 = ex*(A[9]+x*(A[8]+x*(A[7]+x*(A[6]+x))))/(x*(A[13]+x*(A[12]+x*(A[11]+x*(A[10]+x))))) 
    else:
        e1          = np.empty_like(x)
        e1[x == 0.] = 0.

        ok = (x > 0.) & (x <= 1.)
        s  = x[ok]
        e1[ok] = A[0]+s*(A[1]+s*(A[2]+s*(A[3]+s*(A[4]+s*A[5]))))-np.log(s)

        ok = x > 1.
        s  = x[ok]
        e1[x > 1.] = ex[ok]*(A[9]+s*(A[8]+s*(A[7]+s*(A[6]+s))))/(s*(A[13]+s*(A[12]+s*(A[11]+s*(A[10]+s))))) 

    return ex-x*e1
