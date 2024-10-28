"""
smtp related routines
"""

def ttime(tstr):
    """
    ttime translates the time string signifying arrival time of an email
    into a float so that the messages can be ordered by arrival time
    
    Argument:

    tstr  -- date/time string from the result of applying parsestr
             of an email.Parser.Parser object to a message fetched with
             .fetch method of an imaplib.IMAP4_SSL object.

    Return 

    """

    days = (0,31,59,90,120,151,181,212,243,273,304,334)
    (wday,mday,month,year,utc,corr,ttyp) = tstr.split()
    (hour,min,sec) = utc.split(':')
    iy = int(year)

    # start by adding up days since 2000 Jan 0
    if iy < 2000:
        raise SubsError('ttime: year < 2000')
    day = 0.
    for y in range(2000, iy):
        if y % 4 == 0:
            day += 366.
        else:
            day += 365.

    if iy % 4 == 0:
        lday = 1
    else:
        lday = 0
        
    if month == 'Jan':
        im = 0
    elif month == 'Feb':
        im = 1
    elif month == 'Mar':
        im = 2
    elif month == 'Apr':
        im = 3
    elif month == 'May':
        im = 4
    elif month == 'Jun':
        im = 5
    elif month == 'Jul':
        im = 6
    elif month == 'Aug':
        im = 7
    elif month == 'Sep':
        im = 8
    elif month == 'Oct':
        im = 9
    elif month == 'Nov':
        im = 10
    elif month == 'Dec':
        im = 11

    day = float(365*(iy - 2000)) + float(days[im]) + float(mday) - 1 + (int(hour) + int(min)/60. + int(sec)/3600.)/24.
    # correct for leap day of the current year
    if im > 1:
        day += 1.
    return day

