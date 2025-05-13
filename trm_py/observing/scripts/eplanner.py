#!/usr/bin/env python3

usage = \
""" Plots out observing schedules, with a focus on periodic events. You
supply a file of target positions and ephemerides, another defining phase
ranges of interest and optionally a third specifying when to switch targets,
and this will make a graphical representation of the results.

The tracks for the stars start and stop when they reach a user-defined
airmass, set to 2.2 by default, or when the Sun is 10 degrees below
the horizon, whichever is more constraining. Sun at -10 is optimistic
for most targets.

If you see horizontal black error bars at the far right, these
indicate +/- 1 sigma uncertainties on ephemerides. As quoted
uncertainties can barely ever be trusted, they are indicative only. If
they are big, be careful: your eclipse may not appear when you expect.

Grey boxes represent zenith holes for Alt/Az telescopes which can't
track very close to the zenith. These are not always hard limits, but
you may not be able to observe during these intervals. For the TNT
there is also a special indicator of close approach to the TV mast
(only approximate).

A curved, red dashed line represents the elevation of the Moon. An
indication of its illuminated percentage is written at the top and its
minimum separation from targets is indicated if it is below a
user-defined minimum.

Black dots indicate phase 0 for any targets with
ephemerides. Similar-sized open circles mark phase 0.5, but only if
coverage of that phase has been specified.  """
import argparse
import datetime
import pathlib

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter

from astropy import time, coordinates as coord, units as u
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz
from astropy.utils import iers

from trm_py import observing
from trm_py.observing import SITES


def utc_formatter(x, pos):
    """
    Make sure UTC labels are always in range 1-24
    """
    if int(x) > 24:
        x -= 24
    return '{:d}'.format(int(x))


def iso_datehm(t):
    """
    Get round stricter astropy constraints on out_subfmt
    """
    try:
        # works with astropy > 4.0.1
        return t.to_value('iso', out_subfmt='date_hm')
    except (AttributeError, TypeError):
        # fallback for older versions
        out = t.copy(format='iso')
        out.out_subfmt = 'date_hm'
        return out.iso


def sexagesimal(ftime):
    """
    Convert time in hours to sexagesimal string
    """
    if ftime >= 24.: ftime -= 24
    hours = int(ftime)
    minutes = int(60*(ftime-hours))
    seconds = 3600*(ftime-hours-minutes/60)
    return f'{hours:02d}:{minutes:02d}:{seconds:04.1f}'


def eplanner():

    # arguments
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # positional
    parser.add_argument(
        'stardata',
        type=lambda p: pathlib.Path(p).absolute(),
        help='file containing positions and ephemerides. It should have the same format of my C++ ephemeris file.')

    parser.add_argument(
        'telescope', help='Telescope name, e.g. WHT, NTT, VLT, TNT')

    parser.add_argument(
        'pranges',
        type=lambda p: pathlib.Path(p).absolute(),
        help='file of stars to include and phase ranges of interest. Each entry should start with the name of the star on its own on a line, and then follow it with the phase ranges in the form:\np1 p2 col lw\nwhere p1 p2 is the range of phase to indicate, col is the pgplot colour index (2 = red, etc), lw is the pgplot line width. It is OK to have p1 > p2; the programme will map it into a positive range.')

    parser.add_argument(
        'date', help='date at start of night, format: "YYYY-MM-DD"')

    # optional
    parser.add_argument(
        '-a', dest='airmass', type=float, default=2.2,
        help='airmass limit.')

    parser.add_argument(
        '-c', dest='csize', type=float, default=1.0,
        help='default character size for star names')

    parser.add_argument(
        '-d', dest='divide', type=float, default=0.18,
        help='divide point between names on left and plot on right as fraction of total width')

    parser.add_argument(
        '-f', dest='hcopy', default=None,
        help='Name for hard copy file, e.g. plot.pdf, plot.png there will be no interactive plot in this case)')

    parser.add_argument(
        '-v', dest='verbose', action='store_true',
        help='Causes verbose output (at the moment just lists UTs corresponding to phase ranges)'
    )

    parser.add_argument(
        '-i', dest='iers', default='',
        help='URL of IERS table. should normally set one by default but if the server is down, you want another such as http://toshi.nofs.navy.mil/ser7/finals2000A.all'
    )

    parser.add_argument(
        '-m', dest='mdist', type=float, default=25,
        help='separation below which to highlight that the Moon gets close [degrees]')

    parser.add_argument(
        '-o', dest='offset', type=float, default=0.3,
        help='offset as fraction of plot width at which to print duplicate names.')

    parser.add_argument(
        '-p', dest='prefix', default=None,
        help='target name prefixes. These are placed before the left-hand list of target names to allow one to add e.g. priorities.')

    parser.add_argument(
        '-r', dest='reduce', type=float, default=None,
        help='reduce plots as far as possible. Suppress any targets with phase ranges provides where the indicated phases do not crop up on the night, and also targets closer than the optional number of degrees from the Moon (default=0)')

    parser.add_argument(
        '-s', dest='switch', default=None,
        help='switch targets data file name. Each line should have the format:\nstar name | start switch time (e.g. 12:34) | dead time\nwhere the dead time is the time taken (in minutes) to make the switch. The line will be split on the pipe | characters which must therefore not appear in the star name. Use a star name = None as the final switch to terminate early. The times should increase monotonically, so use times > 24 if necessary.')

    parser.add_argument(
        '-t', dest='twilight', type=float, default=-15,
        help='altitude of Sun for twilight [degrees]')

    parser.add_argument(
        '-x', dest='width', type=float, default=11.69,
        help='plot width, inches')

    parser.add_argument(
        '-y', dest='height', type=float, default=8.27,
        help='plot height, inches'
    )

    parser.add_argument(
        '--hfrac', type=float, default=0.84,
        help='plot height as fraction of total')

    parser.add_argument(
        '--bfrac', type=float, default=0.08,
        help='offset of bottom of plot as fraction of total')

    parser.add_argument(
        '--xmajor', type=float, default=2,
        help='spacing of labelled major ticks in X [hours]')

    parser.add_argument(
        '--xminor', type=float, default=1,
        help='spacing of minor ticks in X [hours]')

    # parse them
    args = parser.parse_args()

    # a few checks
    assert(args.airmass > 1 and args.airmass < 6)
    assert(args.width > 0 and args.height > 0)
    assert(args.telescope in SITES)

    if args.iers != '':
        iers.conf.iers_auto_url = args.iers

    # Interpret date.
    date = time.Time(args.date)

    # Location
    info = SITES[args.telescope]
    site = coord.EarthLocation.from_geodetic(
        info['long'], info['lat'], info['height']
    )

    # Load position and ephemeris data
    peinfo = observing.load_pos_eph(args.stardata)

    # Load phase / time ranges
    prinfo = observing.load_ptranges(args.pranges, peinfo)

    # Load switches
    swinfo = observing.load_switches(args.switch, peinfo)

    # Load prefixes
    prefixes = observing.load_prefixes(args.prefix, peinfo)

    # Compute maximum length of the names
    if len(prefixes):
        mpre = max([len(n) for n in prefixes.values()]) + 1
    else:
        mpre = 0

    left = max([len(n) for n in peinfo.keys()]) + mpre

    # Rather primitive times to bracket Sun down and up; won't work properly
    # in far north or south or around the dateline, but should be OK for most
    # sites.  'date' is set at the UTC of the start of the entered date,
    # i.e. 0 hours. To get to the day time before the night of interest, we
    # advance by 0.5 days and apply an offset in longitude. This should be the
    # local mid-day, then half day steps are used to locate midnight and the
    # following mid day. Armed with these we can narrow down on the sunset /
    # rise etc times.
    toffset = time.TimeDelta(0.5 - site.lon.degree/360., format='jd')

    day1  = date + toffset
    night = day1 + time.TimeDelta(0.5, format='jd')
    day2  = night + time.TimeDelta(0.5, format='jd')

    # sunset & sunrise times
    sunset = observing.sun_at_alt(day1, night, site, -0.25)
    sunrise = observing.sun_at_alt(night, day2, site, -0.25)

    # integer offset
    isun = int(sunset.mjd)
    utc1 = 24.*(sunset.mjd-isun)
    utc2 = 24.*(sunrise.mjd-isun)
    print('Sun is at alt=-0.25 at {0:s} and {1:s}'.format(
        iso_datehm(sunset), iso_datehm(sunrise))
    )

    # the lines in the plot will terminate at the point the Sun is at -10 as
    # the absolute limit for observation.
    rset = observing.sun_at_alt(day1, night, site, -10.)
    rrise = observing.sun_at_alt(night, day2, site, -10.)
    utc5 = 24.*(rset.mjd-isun)
    utc6 = 24.*(rrise.mjd-isun)
    print('Sun is at alt=-10.0 at {0:s} and {1:s}'.format(
        iso_datehm(rset), iso_datehm(rrise)))

    # user-defined twilight -- will be indicated on the plot with dashed lines
    twiset = observing.sun_at_alt(day1, night, site, args.twilight)
    twirise = observing.sun_at_alt(night, day2, site, args.twilight)
    utc3 = 24.*(twiset.mjd-isun)
    utc4 = 24.*(twirise.mjd-isun)
    print('Sun is at alt={0:5.1f} at {1:s} and {2:s}'.format(
        args.twilight, iso_datehm(twiset), iso_datehm(twirise)))

    # simulate the old PGPLOT colours
    cols = {
        2 : (0.7,0,0), # red
        3 : (0,0.7,0), # green
        4 : (0,0,0.7), # blue
        5 : (0.7,0.7,0.7), # grey
        9 : (0.8,0.8,0.8),
    }

    # more memorable names
    COLS = {
        'red' : (0.7,0,0), # red
        'green' : (0,0.7,0), # green
        'blue' : (0,0,0.7), # blue
        'grey1' : (0.7,0.7,0.7), # grey
        'grey2' : (0.8,0.8,0.8),
        'highlight' : (1,0.8,0.5),
    }

    # Start the plot
    fig = plt.figure(figsize=(args.width,args.height),facecolor='white')

    # left and right axes. Left for labels, right for the main stuff
    edge = 0.04
    axl = plt.axes([edge, args.bfrac, args.divide, args.hfrac])
    axr = plt.axes([args.divide+edge, args.bfrac, 1-2*edge-args.divide,
                    args.hfrac], frameon=False)

    # correct switch times to lie within range
    for sw in swinfo:
        if sw.utc < utc1 and sw.utc+24 < utc2:
            sw.utc += 24.

    # set up general scale, draw vertical dashed lines every minor
    # tick spacing
    utstart = utc1 - 0.5
    utend = utc2 + 0.5

    n1 = int(np.ceil(utc1/args.xminor))
    n2 = int(np.ceil(utc2/args.xminor))

    for n in range(n1,n2):
        axr.plot([n*args.xminor,n*args.xminor],[0,1],'--',color=COLS['grey2'])

    # mark sunset / twiset / twirise / sunrise

    # first with vertical lines
    kwargs = {'color' : COLS['red']}
    axr.plot([utc1,utc1],[0,1],'--',**kwargs)
    axr.plot([utc2,utc2],[0,1],'--',**kwargs)
    axr.plot([utc3,utc3],[0,1],'--',**kwargs)
    axr.plot([utc4,utc4],[0,1],'--',**kwargs)

    # then labels
    kwargs = {'color' : COLS['red'], 'ha' : 'center', 'va' : 'center',
              'size' : 10}
    axr.text(utc1, 1.02, 'sunset', **kwargs)
    axr.text(utc2, 1.02, 'sunrise', **kwargs)
    axr.text(utc3, 1.02, str(args.twilight), **kwargs)
    axr.text(utc4, 1.02, str(args.twilight), **kwargs)

    # 600 points from start to end defined by the Sun at -10, i.e.
    # roughly 1 per minute.
    utcs = np.linspace(utc5,utc6,600)
    mjds = isun + utcs/24.
    mjds = time.Time(mjds, format='mjd')

    # Calculate MJD in middle of visibility period and the position of the
    # Moon at this time in order to calculate a representative illumination
    # for the Moon which is added at the top of the plot.
    mjd_mid = time.Time(isun + (utc5+utc6)/48., format='mjd')

    moon = get_moon(mjd_mid, location=site)
    sun = get_sun(mjd_mid)
    elong = sun.separation(moon)
    mphase = np.arctan2(sun.distance*np.sin(elong),
                        moon.distance - sun.distance*np.cos(elong))
    illum = (1 + np.cos(mphase))/2.0
    axr.text((utc5+utc6)/2., 1.02,
             'Moon: {:3d}%'.format(int(round(100*illum.value))), **kwargs)

    # Loop through the stars listed in the ephemeris file, computing their hour
    # angles at midnight to establish the plot order
    has = []
    lst = mjd_mid.sidereal_time(kind='mean', longitude=site.lon)

    for key, star in peinfo.items():
        ha = (lst - star.position.ra).value
        if ha > 12.: ha -= 24
        if ha < -12.: ha += 24
        has.append((ha, key))

    # first go at plot order. might be reduced below
    keys = [item[1] for item in sorted(has, key=lambda ha: ha[0], reverse=True)]

    # compute altitude of Moon through the night. Add as red dashed line
    # scaled so that 90 = top of plot.
    moon = get_moon(mjds, location=site)
    altazframe = AltAz(obstime=mjds, location=site)
    altaz = moon.transform_to(altazframe)
    alts = altaz.alt.value
    plt.plot(utcs,alts/90.,'--',color=COLS['red'])

    # Loop through all targets
    stored_info = {}
    skipped_keys = []
    for key in keys:
        star = peinfo[key]

        # Compute airmasses of the star for all times (altazframe
        # encodes them all))
        altaz = star.position.transform_to(altazframe)
        airmasses = altaz.secz.value
        alts = altaz.alt.value
        azs = altaz.az.value
        ok = alts > 90.-np.degrees(np.arccos(1./args.airmass))
        if len(ok[ok]) < 2:
            # must have some points
            print(f'reduce: skipping {key} as it is never observable (airmass < {args.airmass}) during the night')
            skipped_keys.append(key)
            continue

        # Compute minimum distance to the Moon during period target is
        # above airmass limit
        seps = moon.separation(star.position).degree[ok]
        if len(seps):
            sepmin = seps.min()
            moon_close = sepmin < args.mdist
            if moon_close:
                col_moon = COLS['red']
            else:
                col_moon = 'k'
        else:
            moon_close = False

        if args.reduce:
            # section designed to reduce clutter. see whether we can
            # skip any targets according to whether they are too close
            # to the Moon or they have indicated phase ranges that do
            # not appear.

            if args.reduce > sepmin:
                # suppress targets too close to the Moon
                print(f'reduce: skipping {key} as minimum distance from Moon = {sepmin} < {args.reduce}')
                skipped_keys.append(key)
                continue

            if key in prinfo:
                # target has phase or time range info
                pr = prinfo[key]
                pranges = pr.prange

                # determine start and stop times of visibility
                # period. a bit complex because can have targets that
                # are visible at start and end but not the middle
                first = True
                for n, (flag, utc, mjd) in enumerate(zip(ok, utcs, mjds)):
                    if flag:
                        if first:
                            # start a visibility interval
                            utc_start = utc_end = utc
                            mjd_start = mjd_end = mjd
                            first = False
                        else:
                            # extend a visibility interval
                            utc_end = utc
                            mjd_end = mjd

                    if (not flag or n == len(ok)-1) and not first:
                        # end a visibility interval
                        first = True

                        # now test whether any phase or time ranges overlap
                        pstart = None
                        for pt1, pt2, col, lw, p_or_t in pranges:

                            if p_or_t == 'Time':
                                utc1, utc2 = 24.*(pt1-isun), 24.*(pt2-isun)
                                if utc1 < utc_end and utc2 > utc_start:
                                    ut1  = max(utc1, utc_start)
                                    ut2  = min(utc2, utc_end)
                                    if ut1 < ut2:
                                        # yes, we have overlap
                                        break

                            elif p_or_t == 'Phase':
                                # Now the phase info
                                if pstart is None:
                                    # must compute start and stop phases
                                    eph = star.eph
                                    times = time.Time((mjd_start,mjd_end))
                                    if eph.time.startswith('H'):
                                        times += times.light_travel_time(star.position, 'heliocentric', location=site)
                                    elif eph.time.startswith('B'):
                                        times += times.light_travel_time(star.position, location=site)
                                    else:
                                        raise Exception('Unrecognised type of time = ' + eph.time)

                                    if eph.time == 'HJD' or eph.time == 'BJD':
                                        pstart = eph.phase(times[0].jd)
                                        pend   = eph.phase(times[1].jd)
                                    elif eph.time == 'HMJD' or eph.time == 'BMJD':
                                        pstart = eph.phase(times[0].mjd)
                                        pend   = eph.phase(times[1].mjd)
                                    else:
                                        raise Exception('Unrecognised type of time = ' + eph.time)

                                d1 = pstart + (pt1 - pstart) % 1 - 1
                                d2 = pend + (pt1 - pend) % 1
                                nphs = int(np.ceil(d2 - d1))
                                for n in range(nphs):
                                    ut1 = utc_start + (utc_end-utc_start)*(d1 + n - pstart)/(pend-pstart)
                                    ut2  = ut1 + (utc_end-utc_start)/(pend-pstart)*(pt2-pt1)
                                    ut1  = max(ut1, utc_start)
                                    ut2  = min(ut2, utc_end)
                                    if ut1 < ut2:
                                        break
                                else:
                                    # no phase ranges have overlapped, so just continue
                                    continue

                                # only get here from the break statement a few lines back
                                # so break again to escape outer loop
                                break
                        else:
                            # if we get here, there has been no
                            # overlap of any kind so skip the target
                            skipped_keys.append(key)
                            print(f'reduce: skipping {key} as it has no observable phase or time ranges this night')

                            # continue to next target
                            continue

        # Store computed data for later retrieval
        stored_info[key] = {
            'airmasses' : airmasses,
            'alts' : alts,
            'azs' : azs,
            'ok' : ok,
            'moon_close' : moon_close,
            'sepmin' : sepmin if moon_close else None
        }

    # Remove skipped targets
    for key in skipped_keys:
        print(f'Removing {key} from list of stars to plot')
        keys.remove(key)

    # Now compute plot positions
    ys = {}
    for i,key in enumerate(keys):
        ys[key] = (len(keys)-i)/float(len(keys)+1)

    # Plot a line to represent which target to observe
    if args.switch is not None:
        kwargs = {'color': COLS['grey1'], 'lw': 5}

        first = True
        for sw in swinfo:
            if first:
                xstart, ystart = sw.utc, ys[sw.name]
                first = False
            else:
                if sw.name == 'None':
                    axr.plot([xstart,sw.utc],[ystart,ystart],**kwargs)
                    break
                else:
                    xend, yend = sw.utc+sw.delta, ys[sw.name]
                    axr.plot([xstart, sw.utc, xend],[ystart, ystart, yend],**kwargs)
                    xstart, ystart = xend, yend

    # Finally, loop through the stars listed in the phase ranges.
    lbar = min(0.01, 1./max(1,len(keys))/3.)
    for key in keys:
        star = peinfo[key]
        y = ys[key]

        # Retrieve values saved earlier
        sinfo = stored_info[key]
        airmasses = sinfo['airmasses']
        alts = sinfo['alts']
        azs = sinfo['azs']
        ok = sinfo['ok']
        moon_close = sinfo['moon_close']
        sepmin = sinfo['sepmin']

        # Compute minimum distance to the Moon during period target
        # is above airmass limit
        seps = moon.separation(star.position).degree[ok]
        if len(seps):
            sepmin = seps.min()
            moon_close = sepmin < args.mdist
            if moon_close:
                col_moon = COLS['red']
            else:
                col_moon = 'k'
        else:
            moon_close = False

        # now start plotting stuff associated with individual targets
        first = True
        afirst = True
        for n, (flag, utc, mjd) in enumerate(zip(ok, utcs, mjds)):

            if flag:
                if first:
                    utc_start = utc_end = utc
                    n_start = n
                    first = False
                    if afirst:
                        mjd_first = mjd_last = mjd
                        utc_first = utc_last = utc
                        afirst = False
                else:
                    utc_end = utc_last = utc
                    mjd_last = mjd

            if (not flag or n == len(ok)-1) and not first:
                first = True
                if key not in prinfo:
                    # highlight anytime objects
                    axr.plot([utc_start,utc_end],[y,y],color=COLS['highlight'],lw=6,zorder=-10)

                # plot visibility period, highlighted if close to the Moon
                axr.plot([utc_start,utc_end],[y,y],'--',color=col_moon)
                axr.plot([utc_start,utc_start],[y-lbar,y+lbar],color=col_moon)
                axr.plot([utc_end,utc_end],[y-lbar,y+lbar],color=col_moon)

                if flag:
                    n_end = n+1
                else:
                    n_end = n

                if utc_start > utstart + args.offset*(utend-utstart):
                    # repeat target name if the start is delayed to make it
                    # easier to line up names and tracks
                    kwargs = {'ha' : 'right', 'va' : 'center',
                              'size' : 9*args.csize}
                    axr.text(utc_start-0.2, y, key, **kwargs)

        if afirst:
            # never found any ok bit; move on ...
            continue

        if moon_close:
            axr.text(
                utc_last+0.07, y,
                '${:d}^\circ$'.format(int(round(sepmin))),
                ha='left', va='center', size=9*args.csize,
                color=COLS['red']
            )

        # zenith holes
        start = True
        for alt, utc in zip(alts[ok], utcs[ok]):
            if start and alt > 90.-info['zhole']:
                air_start = utc
                air_end = utc
                start = False
            elif not start:
                if alt > 90.-info['zhole']:
                    air_end = utc
                else:
                    break

        if not start:
            plt.fill(
                [air_start,air_end,air_end,air_start],
                [y-lbar,y-lbar,y+lbar,y+lbar],
                color=COLS['grey2'], zorder=10, alpha=0.65
            )
            plt.plot(
                [air_start,air_end,air_end,air_start,air_start],
                [y-lbar,y-lbar,y+lbar,y+lbar,y-lbar],
                color='k',lw=0.8, zorder=11
            )

        # TNT has a stupid TV aerial
        if args.telescope == 'TNT':
            start = True
            aerial = []
            for alt,az,utc in zip(alts[ok],azs[ok],utcs[ok]):
                if az < 0.: az += 360.

                if observing.tnt_alert(alt, az):
                    if start:
                        # start bad period
                        air_start = utc
                        air_end   = utc
                        start     = False
                    else:
                        # update end time of bad period
                        air_end = utc
                else:
                    if not start:
                        # We are out of it. Record bad period.
                        aerial.append((air_start,air_end))
                        start = True

            if not start:
                # We were still in a bad period at the end. Record it.
                aerial.append((air_start,air_end))

            # Plot the bad periods
            for air_start, air_end in aerial:
                axr.fill([air_start,air_end,air_end,air_start],
                         [y-lbar,y-lbar,y+lbar,y+lbar], color=COLS['grey2'], zorder=10, alpha=0.65)
                axr.plot([air_start,air_end,air_end,air_start,air_start],
                         [y-lbar,y-lbar,y+lbar,y+lbar,y-lbar], color='k',lw=0.8, zorder=11)

        if key in prinfo:
            # handles the time ranges only
            pr = prinfo[key]

            pranges = pr.prange
            for t1, t2, col, lw, p_or_t in pranges:
                if p_or_t == 'Time':
                    utc1, utc2 = 24.*(t1-isun), 24.*(t2-isun)
                    if utc1 < utc_last and utc2 > utc_first:
                        ut1  = max(utc1, utc_first)
                        ut2  = min(utc2, utc_last)
                        if ut1 < ut2:
                            plt.plot([ut1,ut2],[y,y],color=cols[col],lw=lw)

        if star.eph:
            # Now the phase info

            eph = star.eph
            times = time.Time((mjd_first,mjd_last))
            if eph.time.startswith('H'):
                times += times.light_travel_time(star.position, 'heliocentric', location=site)
            elif eph.time.startswith('B'):
                times += times.light_travel_time(star.position, location=site)
            else:
                raise Exception('Unrecognised type of time = ' + eph.time)

            if eph.time == 'HJD' or eph.time == 'BJD':
                pstart = eph.phase(times[0].jd)
                pend   = eph.phase(times[1].jd)
            elif eph.time == 'HMJD' or eph.time == 'BMJD':
                pstart = eph.phase(times[0].mjd)
                pend   = eph.phase(times[1].mjd)
            else:
                raise Exception('Unrecognised type of time = ' + eph.time)

            plothalf = False
            if key in prinfo:
                pr   = prinfo[key]

                # Draw phase ranges of interest
                pranges = pr.prange
                for p1, p2, col, lw, p_or_t in pranges:

                    if (p1 < 0.5 and p2 > 0.5) or (p1 < -0.5 and p2 > -0.5) or (p1 < 1.5 and p2 > 1.5):
                        # determine whether or not to draw circles at phase 0.5
                        plothalf = True

                    if p_or_t == 'Phase':
                        d1 = pstart + (p1 - pstart) % 1 - 1
                        d2 = pend + (p1 - pend) % 1
                        nphs = int(np.ceil(d2 - d1))
                        for n in range(nphs):
                            st1 = ut1 = utc_first + (utc_last-utc_first)*(d1 + n - pstart)/(pend-pstart)
                            st2 = ut2  = ut1 + (utc_last-utc_first)/(pend-pstart)*(p2-p1)
                            ut1  = max(ut1, utc_first)
                            ut2  = min(ut2, utc_last)
                            if ut1 < ut2:
                                plt.plot([ut1,ut2],[y,y],color=cols[col],lw=lw)
                                if args.verbose:
                                    delta = 1440.*min(eph.etime((pstart+pend)/2.), eph.coeff[1]/2.)
                                    print(
                                        f'Target {key} is at phase range {p1} to {p2} at'
                                        f' UT {sexagesimal(st1)} to {sexagesimal(st2)} (+/- {delta:.1f} mins), '
                                        f'[visible range: {sexagesimal(utc_first)} to {sexagesimal(utc_last)}]'
                                    )

            # Compute uncertainty in predictions
            delta = 24.*min(eph.etime((pstart+pend)/2.), eph.coeff[1]/2.)
            if delta > 0.03:
                plt.plot([utend-2.*delta,utend],[y,y],'k',lw=2)
                plt.plot([utend-2.*delta,utend-2*delta],[y-0.005,y+0.005],'k',lw=2)
                plt.plot([utend,utend],[y-0.005,y+0.005],'k',lw=2)

            # draws dots at phase zero
            d1 = np.ceil(pstart)
            d2 = np.floor(pend)
            nphs = int(np.ceil(d2 - d1))+1
            for n in range(nphs):
                ut = utc_first + (utc_last-utc_first)*(d1 + n - pstart)/(pend-pstart)
                plt.plot(ut,y,'ok',ms=4)

            # draws circles at phase 0.5
            if plothalf:
                d1 = np.ceil(pstart-0.5)
                d2 = np.floor(pend-0.5)
                nphs = int(np.ceil(d2 - d1))+1
                for n in range(nphs):
                    ut = utc_first + (utc_last-utc_first)*(d1 + n - pstart + 0.5)/(pend-pstart)
                    plt.plot(ut,y,'o',ms=4,mfc='w',mec='k')

        # draw vertical bar at meridian crossing based on what the
        # azimuth does. bit horrible because of odd way azimuth can
        # behave.
        if ok.any():
            az_start, az_end = azs[n_start], azs[n_end-1]
            az_min, az_max = azs[n_start:n_end].min(), azs[n_start:n_end].max()

            if az_start < 180 and az_end > 180:
                # OK, we do cross the meridian
                razs = azs[n_start:n_end]
                rutcs = utcs[n_start:n_end]
                az0 = razs.min()
                az180 = np.abs(razs-180).min()
                if az0 < az180:
                    # we cross to the North of zenith
                    # chop down to actual moment of crossing
                    corr = razs > 180
                    razs[corr] -= 360
                    nmin = np.argmin(razs[razs > 0])
                    nmax = min(len(razs), nmin+2)
                    if nmax == nmin+2:
                        razs = razs[nmin:nmax]
                        rutcs = rutcs[nmin:nmax]
                        utc_mer = np.interp(0.,razs[::-1],rutcs[::-1])
                        plt.plot([utc_mer,utc_mer],[y-1.3*lbar,y+1.3*lbar],'k',lw=2,zorder=20)
                else:
                    # we cross to the South of zenith
                    # chop down to actual moment of crossing
                    nmin = np.argmax(razs[razs < 180])
                    nmax = min(len(razs), nmin+2)
                    if nmax == nmin+2:
                        razs = razs[nmin:nmin+2]
                        rutcs = rutcs[nmin:nmin+2]
                        utc_mer = np.interp(180.,razs,rutcs)
                        plt.plot([utc_mer,utc_mer],[y-1.3*lbar,y+1.3*lbar],'k',lw=2,zorder=20)

    # finish off
    axr.set_xlabel('UTC')
    date.out_subfmt='date'
    axr.set_title('{0!s} ({1:s}, airmass < {2:3.1f})'.format(
        date, args.telescope, args.airmass))

    axr.set_xlim(utstart, utend)
    axr.set_ylim(0,1.05)

    axr.xaxis.set_major_locator(MultipleLocator(args.xmajor))
    axr.xaxis.set_minor_locator(MultipleLocator(args.xminor))
    axr.xaxis.set_major_formatter(FuncFormatter(utc_formatter))
    axr.get_xaxis().tick_bottom()
    axr.axes.get_yaxis().set_visible(False)
    fig.canvas.draw()

    # add target names at left
    for key in keys:
        y = ys[key]
        if key in prefixes:
            name = prefixes[key] + ' ' + key
        elif len(prefixes):
            name = mpre*' ' + key
        else:
            name = key
        axl.text(0.95,y,name,ha='right',va='center',size=9*args.csize)
    axl.set_xlim(0,1)
    axl.set_ylim(0,1.05)
    axl.set_axis_off()

    xmin, xmax = axr.get_xaxis().get_view_interval()
    ymin, ymax = axr.get_yaxis().get_view_interval()
    axr.add_artist(
        Line2D((xmin, xmax), (ymin, ymin), color='black',
               linewidth=2))

    if args.hcopy:
        plt.savefig(args.hcopy)
    else:
        plt.show()


if __name__ == '__main__':
    import sys
    # run the main function
    eplanner()
    sys.exit(0)
