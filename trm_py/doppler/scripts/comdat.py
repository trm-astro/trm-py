import copy
import argparse
import numpy as np
from .. import (
    Map,
    Data,
    cpp_doppler as doppler
)


__all__ = ['comdat',]


def comdat(args=None):
    """comdat computes the data corresponding to a Doppler map, with the option
    of adding noise. Use '-h' to see help.

    """

    parser = argparse.ArgumentParser(description=comdat.__doc__)

    # positional
    parser.add_argument('map',   help='name of the input map')
    parser.add_argument('dtemp', help='data template file')
    parser.add_argument('dout',  help='data output file')

    # optional
    parser.add_argument(
        '-n', dest='noise', action='store_true',
        help='add noise according to uncertainty array in template'
    )

    # OK, done with arguments.
    args = parser.parse_args()

    # load map and data
    dmap  = Map.rfits(doppler.afits(args.map))
    dtemp = Data.rfits(doppler.afits(args.dtemp))

    flux = dtemp.data[0].flux
    ferr = dtemp.data[0].ferr

    # compute data
    dcopy = copy.deepcopy(dtemp)
    
    # This is from the C++ code
    doppler.comdat(dmap, dtemp)

    # optionally add noise
    if args.noise:
        for spectra in dtemp.data:
            spectra.flux = np.random.normal(spectra.flux, np.abs(spectra.ferr))
    else:
        chisq, ndata = doppler.chisquared(dcopy, dtemp)
        print('Chi**2 = ',chisq,', chi**2/N =',chisq/ndata,', N =',ndata)

    # Write to a fits file
    dtemp.wfits(doppler.afits(args.dout))
    print('Written computed data to',args.dout)

