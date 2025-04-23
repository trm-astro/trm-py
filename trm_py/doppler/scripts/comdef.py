#!/usr/bin/env python
import copy
import argparse
from .. import (
    Map,
    afits,
    cpp_doppler as doppler,
)

__all__ = ['comdef',]


def comdef(args=None):
    """
    comdef computes the default equivalent to an image.
    """

    parser = argparse.ArgumentParser(description=comdef.__doc__)

    # positional
    parser.add_argument('map',   help='name of the input map')
    parser.add_argument('dout',  help='default output file')

    # OK, done with arguments.
    args = parser.parse_args()

    # load map
    dmap  = Map.rfits(afits(args.map))

    # copy the map to compute the entropy
    mcopy = copy.deepcopy(dmap)

    # compute default
    # This is from the C++ code
    doppler.comdef(dmap)

    # write the result to a FITS file
    dmap.wfits(afits(args.dout))
