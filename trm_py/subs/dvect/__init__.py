"""Dvect

:author: Tom Marsh
"""

__author__ = "Tom Marsh"

#import .dvect
from .dvect import *

__all__ = ['dvect']
__all__.extend(dvect.__all__)

if __name__ == "__main__":

    import numpy as np
    from numpy.ma import MaskedArray

    try:
        import matplotlib.pyplot as plt
        x  = np.arange(100,dtype=np.float)
        y  = x*x
        ey = 100.*np.ones_like(x)
        mx = (x > 20) & (x < 40)
        my = (x > 60) & (x < 65)
        dx = Dvect(x, mask=mx)
        dy = Dvect(y, err=ey, mask=mx)
        plt.plot(dx, dy)
        errorbar(plt, dx, dy)
        plt.show()
    except ImportError:
        raise ImportError("matplotlib is needed for testing plotting")
