# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Doppler tomography

This is an implementation of Doppler tomography in Python. It
allows flexible configuration of both images and input data.
"""

from .._cpp import _cpp_doppler

from .core import afits, acfg
from .data import Spectra, Data
from .map import Image, Map
from .grid import Grid
from .derived import genmat, genvec, svd

# scripts is a sub-package
from . import scripts

__all__ = ['Image', 'Map', 'Spectra', 'Data', 'Grid',
           'genmat', 'genvec', 'svd', 'afits', 'acfg',
           ]
