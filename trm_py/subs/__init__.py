"""
a module of general utility routines, defined mostly in util

subs is intended to provide generically useful classes, functions and
data. See also the sub-packages below.

Functions
=========
air2vac            -- convert from air to vacuum wavelengths
centroid           -- Schneider & Young-style single-gaussian cross-correlation
convex_hull        -- returns convex hull in 2D
date2dmy           -- returns (year,month,day) from a date string
d2hms              -- produce an hh:mm:ss.ss time string from a decimal hour.
find_exec          -- searches for an executable
gammq              -- incomplete gamma function (pybind function from C++)
hms2d              -- produce a decimal hour from an hh:mm:ss.ss time string
inpoly             -- determines whether a point is inside or outside a polygon
int2month          -- return 3-letter name of a month from an integer 1 to 12
linterp            -- linear interpolation of an array onto a different sampling
minpoly            -- determines whether points are inside or outside a polygon
m2min              -- computes minimum mass of a companion star given m1 and a mass function.
month2int          -- return an integer 1 to 12 from name of a month
mr_wd_eggleton     -- radius of white dwarf according to a formula of Eggleton's
observatory        -- interactive selection of observatory information
orbital_separation -- orbital separation given masses and period
orbital_period     -- orbital period given masses and the separation
planck             -- compute Planck function
pts2cont           -- converts x,y points into contourable image
rlobe_eggleton     -- Eggleton's relative Roche lobe radius formula
sigma_reject       -- computes clipped mean of a numpy array
sinfit             -- fits a sinusoid plus a constant to some data
slitloss           -- compute slitloss for gaussian seeing
splfit             -- fits a spline to a 1D array
str2radec          -- convert from a string to numerical RA, Dec
vac2air            -- convert from vacuum to air wavelengths
voigt              -- voigt function (pybind function from C++)
zeta_wd_eggleton   -- logarithmic radius vs mass derivative

Classes
=======
iStr        -- case insensitive string class
Odict       -- ordered dictionary class     
Fname       -- handles filenames with standard extensions
SubsError   -- exception class for the package
Poly        -- Polynomial class
Vec3        -- Cartesian 3-vectors

Data
====
ME     -- mass of the electron, SI
MP     -- mass of the proton, SI
E      -- charge on the electron, SI
H      -- Planck's constant SI
C      -- Speed of light, SI
K      -- Boltzmann's constant, SI
SIGMA  -- Stefan-Boltzmann, SI
RSUN   -- Radius of Sun, SI
PARSEC -- 1 parsec, SI
AU     -- 1 AU, SI
G      -- Gravitational constant, SI
MSUN   -- Mass of the Sun, SI
GMSUN  -- Gravitational parameter of the Sun, SI
KGAUSS -- Gauss' gravitational constant G*MSUN, AU^(3/2) day^-1
MJ     -- Mass of Jupiter, SI
DAY    -- Length of a day, SI
YEAR   -- Length of a year, SI

Sub-packages
============
cpp    -- some C++ helper routines
dvect  -- data vectors [data, errors, label]
input  -- parameter input with default storage 
smtp   -- provides one function useful for handling smtp-based e-mail

Withdrawn functions
===================

fasper  -- please see the package trm.pgram instead.
medfilt -- use scipy.signal.medfilt instead
"""

from . import cpp
from . import dvect
from . import input
from . import smtp
from . import util
# import above from .util

# Functions
from .util import air2vac
from .util import centroid
from .util import convex_hull
from .util import date2dmy
from .util import d2hms
from .util import find_exec
from .util import gammq
from .util import hms2d
from .util import inpoly
from .util import int2month
from .util import linterp
from .util import minpoly
from .util import m2min
from .util import month2int
from .util import mr_wd_eggleton
from .util import observatory
from .util import orbital_separation
from .util import orbital_period
from .util import planck
from .util import pts2cont
from .util import rlobe_eggleton
from .util import sigma_reject
from .util import sinfit
from .util import slitloss
from .util import splfit
from .util import str2radec
from .util import vac2air
from .util import voigt
from .util import zeta_wd_eggleton

# Classes
from .util import iStr
from .util import Odict
from .util import SubsError
from .util import Poly
from .util import Vec3

# Data
from .util import ME
from .util import MP
from .util import E
from .util import H
from .util import C
from .util import K
from .util import SIGMA
from .util import RSUN
from .util import PARSEC
from .util import AU
from .util import G
from .util import MSUN
from .util import GMSUN
from .util import KGAUSS
from .util import MJ
from .util import DAY
from .util import YEAR


# Expose Fname from input (for backward compatibility)
class Fname(input.Fname):
    # simply inherit from input.Fname
    # add a warning that this should be imported from input instead
    def __init__(*args, **kwargs):
        import warnings
        warnings.warn(
            "Fname should be imported from subs.input instead of subs",
            DeprecationWarning
        )
        super().__init__(*args, **kwargs)


__all__ = [
    # sub-packages
    'cpp', 'dvect', 'input', 'smtp',

    # Functions
    'air2vac', 'centroid', 'convex_hull', 'date2dmy', 'd2hms', 'find_exec', 
    'gammq', 'hms2d', 'inpoly', 'int2month', 'linterp', 'minpoly', 'm2min', 
    'month2int', 'mr_wd_eggleton', 'observatory', 'orbital_separation', 
    'orbital_period', 'planck', 'pts2cont', 'rlobe_eggleton', 'sigma_reject', 
    'sinfit', 'slitloss', 'splfit', 'str2radec', 'vac2air', 'voigt', 
    'zeta_wd_eggleton',

    # Classes
    'iStr', 'Odict', 'Fname', 'SubsError', 'Poly', 'Vec3',

    # Data
    'ME', 'MP', 'E', 'H', 'C', 'K', 'SIGMA', 'RSUN', 'PARSEC', 'AU', 'G', 
    'MSUN', 'GMSUN', 'KGAUSS', 'MJ', 'DAY', 'YEAR'

]


# # Withdrawn functions

def fasper(*args, **kwargs):
    # Raise an import error here
    raise ImportError("fasper no longer supported \
                      please see the package trm.pgram instead")


def medfilt(*args, **kwargs):
    # Raise an import error here
    raise ImportError("mdefilt no longer supported \
                      use scipy.signal.medfilt instead")
