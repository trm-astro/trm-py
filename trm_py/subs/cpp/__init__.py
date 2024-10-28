import os
import struct
import numpy as npy
import trm.subs as subs

"""
Various routines to help interface with data written by my C++
programs like ULTRACAM and Ponto
"""

def read_string(fobj, endian=''):
    """
    Reads a string written in binary format by my C++ code

    fobj   -- file object opened for binary input
    endian -- '>' for big-endian, '' for little-endian.
    """
    (nchar,)  = struct.unpack(endian + 'i', fobj.read(4))
    (string,) = struct.unpack(endian + str(nchar) + 's', fobj.read(nchar))
    return string

def open_ponto(fname):
    """
    Opens a Ponto file for reading and returns a file object

    Returns (fobj,endian) where fobj is the file object and endian is a string to be passed
    to later routines indicating endian-ness
    """
    MAGIC = 7856871
    fobj = open(fname, 'rb')

    # read the format code
    fbytes = fobj.read(4)
    (fcode,) = struct.unpack('i',fbytes)
    if fcode != MAGIC:
        (fcode,) = struct.unpack('>i',fbytes)
        if fcode != MAGIC:
            fobj.close()
            raise CppError('open_ponto: could not recognise first 4 bytes of ' + fname + ' as a Ponto file')
        endian = '>'
    else:
        endian = ''
    return (fobj,endian)

def write_string(fobj, strng):
    """
    Writes a string in binary format for my C++ code

    fobj         -- file object opened for binary output
    strng        -- string to file object opened for binary output
    """
    nchar = len(strng)
    fobj.write(struct.pack('i' + str(nchar) + 's',nchar, strng)) 

def read_units(fobj, endian=''):
    """
    Reads Units object written by my C++ routines from disk file
    positioned at start of header.

    fobj   -- file object
    endian -- '>' for big-endian, '' for little-endian.
    """

    (nmass,dmass,nlength,dlength,ntime,dtime,ncharge,dcharge,nsolid,dsolid,scale) = struct.unpack(endian + '10id',fobj.read(48))  
    pgname = read_string(fobj, endian)
    (pgset,)  = struct.unpack(endian + 'B',fobj.read(1))
    return (pgname,scale,nmass,dmass,nlength,dlength,ntime,dtime,ncharge,dcharge,nsolid,dsolid)

def read_buffer1d(fobj, dtype, endian=''):
    """
    Reads buffer1d object written by my C++ routines from disk file
    positioned at start of it

    fobj   -- file object
    dtype  -- 'float' or 'double'
    endian -- '>' for big-endian, '' for little-endian.
    """

    (npix,) = struct.unpack(endian + 'i', fobj.read(4))
    if dtype == 'float':
        arr = npy.fromfile(file=fobj, dtype=npy.float32, count=npix)
    elif dtype == 'double':
        arr = npy.fromfile(file=fobj, dtype=npy.float64, count=npix)
    else:
        raise CppError('read_buffer1d: do not recogniise dtype = ' + str(dtype))
    return arr
    
def read_header(fobj, endian=''):
    """
    Reads header written by my C++ routines from disk file
    positioned at start of header. The header is returned in
    an trm.subs.Odict which loses the comments and so this
    needs upgrading but it will do for now.

    fobj   -- file object
    endian -- '>' for big-endian, '' for little-endian.
    """    

    # read the header
    lstr = fobj.read(4)
    if lstr == '':
        raise EOFError('read_header: EOF encountered at start of header read')
    (lmap,) = struct.unpack(endian + 'i', lstr)
    
    head = subs.Odict()
    for i in xrange(lmap):
        name  = read_string(fobj, endian)
        (itype,) = struct.unpack(endian + 'i', fobj.read(4))
        comment = read_string(fobj, endian)
        
        if itype == 0: # double
            (value,) = struct.unpack(endian + 'd', fobj.read(8))
        elif itype == 1: # char
            raise CppError('read_header: char not enabled')
        elif itype == 2: # int
            (value,) = struct.unpack(endian + 'i', fobj.read(4))
        elif itype == 3: # uint
            raise CppError('read_header: uint not enabled')
        elif itype == 4: # lint
            raise CppError('read_header: linit not enabled')
        elif itype == 5: # ulint
            raise CppError('read_header: ulint not enabled')
        elif itype == 6: # float
            (value,) = struct.unpack(endian + 'f', fobj.read(4))
        elif itype == 7: # string
            value = read_string(fobj, endian)
        elif itype == 8: # bool
            (value,) = struct.unpack(endian + 'B', fobj.read(1))
        elif itype == 9: # directory
            value = subs.Odict()
        elif itype == 10: # date
            raise CppError('read_header: date not enabled')
        elif itype == 11: # time
            (mjd,)  = struct.unpack(endian + 'i', fobj.read(4))
            (hour,) = struct.unpack(endian + 'd', fobj.read(8))
            value   = (mjd, hour)
        elif itype == 12: # position
            value = subs.Odict()
            (value['RA'],)       = struct.unpack(endian + 'd', fobj.read(8))
            (value['Dec'],)      = struct.unpack(endian + 'd', fobj.read(8))
            value['System']      = 'ICRS'
            (value['Epoch'],)    = struct.unpack(endian + 'd', fobj.read(8))
            (value['PmRA'],)     = struct.unpack(endian + 'f', fobj.read(4))
            (value['PmDec'],)    = struct.unpack(endian + 'f', fobj.read(4))
            (value['Parallax'],) = struct.unpack(endian + 'f', fobj.read(4))
            (value['RV'],)       = struct.unpack(endian + 'f', fobj.read(4))
        elif itype == 13: # dvector
            raise CppError('read_header: dvector not enabled')
        elif itype == 14: # uchar
            (value,) = struct.unpack(endian + 'c', fobj.read(1))
        elif itype == 15: # telescope
            tname = read_string(fobj, endian)
            sname = read_string(fobj, endian)
            (longitude,) = struct.unpack(endian + 'd', fobj.read(8))
            (latitude,)  = struct.unpack(endian + 'd', fobj.read(8))
            (height,)    = struct.unpack(endian + 'f', fobj.read(4))
            value = subs.Odict()
            value['Name']        = tname
            value['Observatory'] = sname
            value['Longitude']   = longitude
            value['Latitude']    = latitude
            value['Height']      = height
        else:
            raise CppError('read_header: itype = ' + str(itype) + ' not recognised.')

        clist = name.split('.')
        head_set(head, clist, value)
    
    return head

def head_set(head, clist, value):
    """
    Recursively sets a header item from a list of string 'clist' in a series of 'Odict' objects
    e.g. if clist = ['Town','Street'] and value = 'High Street', then this sets
    head['Town']['Street'] = 'High Street'
    """

    if len(clist) > 1:
        if clist[0] not in head:
            raise CppError('head_set: key = ' + str(clist[0]) + ' is not defined.')
        head_set(head[clist[0]], clist[1:], value)
    elif len(clist) == 1:
        head[clist[0]] = value
    else:
        raise CppError('head_set: cannot call with zero length list length')

def read_ponto_axis(fobj, endian):
    """
    Reads a Ponto Axis.

    fobj   -- file object
    endian -- endian character
    
    Returns (name,units,data,errors). 'errors' = None if there are none.
    """

    axtype = read_string(fobj, endian)

    # skip 4 byte repeat of npix
    (npix,) = struct.unpack(endian + 'i', fobj.read(4))

    axname = read_string(fobj, endian)

    axunits = read_units(fobj, endian)

    # ignore wrap string and level
    read_string(fobj, endian)
    fobj.seek(4, os.SEEK_CUR)

    # read in data and errors
    if axtype == 'DOUBLE':
        data   = read_buffer1d(fobj, 'double', endian)
        errors = None
    elif axtype == 'DOUBLEERROR':
        data   = read_buffer1d(fobj, 'double', endian)
        errors = read_buffer1d(fobj, 'float', endian)
    elif axtype == 'FLOAT':
        data   = read_buffer1d(fobj, 'float', endian)
        errors = None
    elif axtype == 'FLOATERROR':
        data   = read_buffer1d(fobj, 'float', endian)
        errors = read_buffer1d(fobj, 'float', endian)
    else:
        return CppError('rponto: sorry have not implemented reading of Axis type = ' + axtype + ' yet.')
    
    return (axname, axunits[0], data, errors)
        

class CppError(Exception):
    pass


