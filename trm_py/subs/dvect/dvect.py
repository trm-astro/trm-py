"""
trm.subs.dvect : a package to represent data plus errors with masking of bad data.

It is common to have data in the form of a series of values plus errors.  Very
often some of the data is bad but you may not want to remove it in order to be
able to keep it in step with some corresponding bit of other data.  You thus
want to mask it. Representing such data vectors is the purpose of this
package. It also makes some limited attempts to combine errors of multiple
data vectors in a consistent manner.

Example
=======

from trm.subs.dvect.Dvect as Dvect 
import numpy as np

x  = np.array([1,2,3,4])
e  = np.array([0.1,0.2,0.1,0.15])
m  = np.array([0,1,0,0])
d  = Dvect(x, err=e, mask=mx,name='Wavelength',units='nm')

a  = 2*d

Note that the operations one performs on a Dvect may often invalidate 
the name and units but no attempt is made to do anything more than pass
these along. 

Warning on errors
=================

One must beware blindly assuming that the routines will return the right uncertainties.
For instance comsider:

a = b*b

vs

a = b**2

They will give the same data values but potentially different values for the errors if
it is not realised that the two inputs in the first case are in fact the same and therefore 
100% correlated. The second is safer, although in this particular simple example, a check
that the two inputs are identical will be applied to obtain the correct result. In general
however, such checks cannot be done, e.g. the two arguments of the * operation in

a = (b+c)*c

are correlated but not identical and it would be hard to get the routine to
spot this automatically. In such cases you will have to cook up a bespoke
routine to get the right errors.
"""

from __future__ import with_statement, print_function 

__author__ = "Tom Marsh"

import numpy as np
from numpy import array, ndarray
from numpy.ma import MaskedArray, masked_array, nomask, mask_or, getdata, masked, getmaskarray

__all__ = ['Dvect', 'DvectError', 'errorbar']

# used to indicate that there are no errors available.
noerr   = None

# default error value for filling bad data
fillerr = -1.

class _binarywrap(object):
    """
    Wrapper class for two-argument maths functions (add, mul, div, etc) of 
    the form z = f(x, y). When called, returns a new Dvect object with the
    result of the specified operation applied to the two inputs.

    The errors are usually combined assuming that the errors of the two inputs
    are independent except for the special case of identical inputs. The
    constructor expects three arguments: f, fx and fy where the last two are
    partial derivatives wrt x and y respectively. The name and units are taken
    from the first Dvect argument encountered, with no checks for
    incompatibilities.

    Parameters:

    f   -- string of corresponding binary operation to apply, e.g. '__add__'
    fx  -- function corresponding to partial derivative wrt first argument
    fy  -- function corresponding to partial derivative wrt second argument
    inp -- carry out operation in place

    If either fx or fy is None, no errors will be propagated.
    """
    def __init__ (self, f, fx=None, fy=None, inp=False):
        self.f   = f
        self.fx  = fx
        self.fy  = fy
        self.inp = inp
        self.obj = None

    def __get__(self, obj, objtype=None):
        "Gets the calling object."
        self.obj = obj
        return self

    def __call__ (self, two):
        "Executes the call behavior."

        # first argument
        one = self.obj

        # carry out basic operation, transfer name and units
        func = getattr(super(Dvect, one), self.f)
        if self.inp:
            func(two)
            result = one
        else:
            result = MaskedArray(func(two), subok=True).view(type(one))
            result._name  = one._name
            result._units = one._units

        # handle the errors. They are worked out in a linear approximation
        # which requires the user to supply two partial derivative functions
        # in addition to the basic function.
        if self.fx is not None and self.fy is not None:
            if isinstance(two, Dvect):
                if one._err is noerr:
                    result._err = getdata(np.abs(self.fy(one.dat, two.dat))*two._err)
                elif two._err is noerr:
                    result._err = getdata(np.abs(self.fx(one.dat, two.dat))*one._err)
                else:
                    if one is two:
                        # Two inputs are identical and not independent
                        if result.mask is nomask:
                            result._err = np.hypot(self.fx(one.dat, two.dat), self.fy(one.dat, two.dat))*one._err
                        else:
                            result._err = np.where(result.mask, fillerr, 
                                                   getdata(np.hypot(self.fx(one.dat, two.dat), self.fy(one.dat, two.dat)*two._err)))
                        pass
                    else:
                        # Two inputs are assumed to be independent.
                        if result.mask is nomask:
                            result._err = getdata(np.hypot(self.fx(one.dat, two.dat)*one._err,
                                                           self.fy(one.dat, two.dat)*two._err))
                        else:
                            result._err = np.where(result.mask, fillerr, 
                                                   getdata(np.hypot(self.fx(one.dat, two.dat)*one._err,
                                                                    self.fy(one.dat, two.dat)*two._err)))
            else:
                result._err = getdata(np.abs(self.fx(one.dat, two))*one._err)
        else:
            result._err = noerr

        return result

class Dvect(MaskedArray):

    """
    Class for storing data based upon numpy arrays. It allows for masking by 
    inheriting from numpy masked arrays, with the additional attributes of errors and
    name and unit strings.

    Construction:
        x = Dvect(data, err=None, name='', units='', mask=np.ma.nomask, dtype=None, edtype=None)

    Parameters
    ----------
    data : {array_like}
         Input data, array object
    err : {None, array_like}, optional
         Errors. Shape must match the shape of input_array or a DvectError will be raised.
    name : {string}, optional
         name of the data, e.g. 'Wavelength'
    units : {string}, optional
         units for the data, e.g.'nm'
    mask : {sequence}, optional
         data mask. True elements will be masked.
    dtype : {dtype}, optional   
         Data type of data
    edtype : {dtype}, optional
         Data type of the errors
    """

    def __new__(cls, data, err=noerr, name='', units='', mask=nomask, dtype=None, edtype=None):
        """
        Construct a Dvect from an input array which we try to make into a masked array as far as possible
        """

        # Try to make a masked array out of the input
        #obj = masked_array(data, dtype=dtype, mask=mask, subok=True).view(cls)
        obj = MaskedArray.__new__(cls, data, dtype=dtype, mask=mask)
        if not isinstance(obj, Dvect):
            obj = obj.view(cls)

        # Add attributes
        if err is not noerr:
            obj._err = array(err, dtype=edtype)
            if obj._err.shape != obj.shape:
                raise DvectError('Dvect.__new__: errors and data are incompatible')
        obj._name  = name
        obj._units = units

        return obj

    def __array_finalize__(self,obj):
        """
        This function is needed for numpy subclassing because of the ways numpy arrays
        can be constructed.
        """
        if obj is None: return
        self._err   = getattr(obj, 'err', noerr)
        self._name  = getattr(obj, 'name', '')
        self._units = getattr(obj, 'units', '')
        MaskedArray.__array_finalize__(self, obj)

    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]
        
        Returns the item described by indx. Not a copy.
        """

        _data = ndarray.__getattribute__(self, '_data')
        _mask = ndarray.__getattribute__(self, '_mask')
        _err  = ndarray.__getattribute__(self, '_err')

        output = _data.__getitem__(indx)
        if not getattr(output, 'ndim', False):
            if _mask is not nomask and _mask[indx]:
                return masked
            if _err is not noerr:
                return (output, _err[indx])
            else:
                return (output, None)

        newdvect = _data.__getitem__(indx).view(type(self))
        newdvect._mask = _mask[indx]
        if _err is noerr:
            newdvect._err = noerr
        else:
            newdvect._err = _err[indx]
        newdvect._name  = self._name
        newdvect._units = self._units
        return newdvect

    def __setitem__(self, indx, value):
        """x.__setitem__(i, y) <==> x[i]=y

        Sets item described by indx. If value is masked, masks those locations.
        Errors are also replaced if there are errors in both source and destination.
        """

        MaskedArray.__setitem__(self, indx, value)
        if isinstance(value, Dvect) and self._err is not noerr and value._err is not noerr:
            print (len(indx), len(self._err), len(value._err))
            self._err[indx] = value._err

    def __repr__(self):
        """Literal string representation"""
        clss   = 'Dvect'
        bclss  = ' '*len(clss)
        strg   = clss +  '(data   = ' + str(self.dat) + '\n'
        if self.err is noerr:
            strg  += bclss + ' errors = ' + str(self.err) + '\n'
        else:
            strg  += bclss + ' errors = ' + str(MaskedArray(self.err, mask=self.mask, copy=False)) + '\n'
        strg  += bclss + ' mask   = ' + str(self.mask) + '\n'
        strg  += bclss + ' name   = "' + self._name + '"\n'
        strg  += bclss + ' units  = "' + self._units + '")\n'        
        return strg

    def __str__(self):
        """String representation. Data then errors, appropriately masked"""
        strg   = self.label + ', ' + str(self.dat)
        if self.err is not noerr:
            strg  += ', ' + str(MaskedArray(self.err, mask=self.mask, copy=False))
        return strg

    # Access to the errors
    def __seterr__(self, err):
        """Set the errors"""

        if err is noerr:
            self._err = noerr
            return
        
        # Try to make an array from the input
        err = array(err, copy=False)
        if not err.ndim:
            err = array([err.item()*len(self)])
        elif isinstance(err, np.ndarray):
            if err.shape != self.shape:
                raise DvectError('dvect.Dvect.__seterr_: errors array must be same length as data')
            err = array(err)
        self._err = err
        return

    def _geterr(self):
        """Gets the errors"""
        if self._err is noerr:
            return noerr
        else:
            return MaskedArray(self._err, mask=self.mask)

    err = property(_geterr, __seterr__, doc='The errors of a Dvect')

    # Access to the data
    def _getdata(self):
        """Returns the data as a pure masked array"""
        return self.view(MaskedArray)

    dat = property(_getdata, doc='Copy of the data of the Dvect as a MaskedArray')

    # Protected access to the name
    def _getname(self):
        """Returns name of the Dvect"""
        return self._name

    def _setname(self, name):
        """Sets name of the Dvect"""
        self._name = str(name)
        return

    name = property(_getname, _setname, doc='The name of a Dvect')

    # Protected access to the units
    def _getunits(self):
        """Returns units of the Dvect"""
        return self._units

    def _setunits(self, units):
        """Sets units of the Dvect"""
        self._units = str(units)
        return

    units = property(_getunits, _setunits, doc='The units of a Dvect')


    def _getlabel(self):
        """Gets label of a Dvect"""
        return self._name + ' [' + self._units + ']'
    label = property(_getlabel, doc='The label of a Dvect')

    # Binary functions
    # To get the reversed versions one must transpose the last two arguments
    # and swap x for y and y for x in the body of the function. In commutative 
    # cases, this has no overall effect
    __add__       = _binarywrap('__add__', lambda x, y: 1.0, lambda x, y: 1.0)
    __radd__      = __add__
    __mul__       = _binarywrap('__mul__', lambda x, y: y, lambda x, y: x)
    __rmul__      = __mul__
    __sub__       = _binarywrap('__sub__', lambda x, y: 1.0, lambda x, y: -1.0)
    __rsub__      = _binarywrap('__rsub__', lambda x, y: -1.0, lambda x, y: 1.0)
    __div__       = _binarywrap('__div__', lambda x, y: 1./y, lambda x, y: -x/y**2)
    __rdiv__      = _binarywrap('__rdiv__', lambda x, y: -y/x**2, lambda x, y: 1./x)
    __pow__       = _binarywrap('__pow__', lambda x, y: y*pow(x,y-1), lambda x, y: log(x)*pow(x,y))
    __truediv__   = _binarywrap('__truediv__', lambda x, y: np.reciprocal(y), lambda x, y: -np.divide(x,np.square(y**2)))
    __rtruediv__  = _binarywrap('__rtruediv__', lambda x, y: -np.divide(y,np.square(x**2)), lambda x, y: np.reciprocal(x))
    __floordiv__  = _binarywrap('__floordiv__')
    __rfloordiv__ = _binarywrap('__rfloordiv__')
    __iadd__      = _binarywrap('__iadd__', lambda x, y: 1.0, lambda x, y: 1.0, True)
    __imul__      = _binarywrap('__imul__', lambda x, y: y, lambda x, y: x, True)
    __isub__      = _binarywrap('__isub__', lambda x, y: 1.0, lambda x, y: -1.0, True)
    __idiv__      = _binarywrap('__idiv__', lambda x, y: 1./y, lambda x, y: -x/y**2, True)
    __lt__        = _binarywrap('__lt__')
    __gt__        = _binarywrap('__gt__')
    __le__        = _binarywrap('__le__')
    __ge__        = _binarywrap('__ge__')
    __eq__        = _binarywrap('__eq__')

# Exception class
class DvectError(Exception):
    """Class for dvect related errors"""
    pass

# Related functions
def errorbar(ax, x, y, xerr=True, yerr=True, fmt='-', ecolor=None, elinewidth=None, capsize=3, 
             barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, **kwargs):
    """Plots two Dvects against each other. Requires matplotlib.
    
    See matplotlib.axes.Axes.errorbar for an explanation of most fields. The only differences are:

    ax   -- the Axes you are plotting in, needed to allow a call of the form ax.errorbar(...) to the
            standard matplotlib errorbar plotter
    x    -- x Dvect
    y    -- y Dvect

    Quite often you may want to set fmt=.' to get points. The lines are for compatibility with
    the matplotlib errorbar routine.
    """
    import warnings
    
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise DvectError("dvect.errorbar: matplotlib needed for plotting") 

    # Identify joint OK part
    ok = ~(getmaskarray(x) | getmaskarray(y))

    # Catch xerr and yerr which are re-defined cf maplotlib errorbar since 
    # Dvects carry their own errors
    xerr = x.err.data[ok] if xerr and x.err is not None and isinstance(x,Dvect) else None
    yerr = y.err.data[ok] if yerr and y.err is not None and isinstance(y,Dvect) else None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax.errorbar(x.data[ok], y.data[ok], xerr=xerr, yerr=yerr, fmt=fmt, ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, 
                    barsabove=barsabove, lolims=lolims, uplims=uplims, xlolims=xlolims, xuplims=xuplims, **kwargs)
