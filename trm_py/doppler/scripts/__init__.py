"""Scripts sub-module of trm.doppler contains all the commands used
from the terminal.

"""

from .comdat import comdat
from .comdef import comdef
from .makedata import makedata
from .makemap import makemap
from .memit import memit
from .optgam import optgam
from .optscl import optscl


__all__ = ['makedata', 'makemap', 'memit', 'optgam', 'optscl',
           'comdat', 'comdef']
