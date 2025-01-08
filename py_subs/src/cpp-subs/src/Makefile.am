## Process this file with automake to generate Makefile.in
##
## This is the file that must be edited if you are changing anything in the source directory

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS   = weekday tcorr gap test

weekday_SOURCES = weekday.cc
tcorr_SOURCES   = tcorr.cc
gap_SOURCES     = gap.cc
test_SOURCES    = test.cc 

AM_CPPFLAGS = -I../include -I../.

LDADD = libsubs.la

DATE = $(shell date)

install-data-hook:
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION

## Library

lib_LTLIBRARIES = libsubs.la

#libsubs_la_SOURCES = gauss.cc

libsubs_la_SOURCES = gauss.cc ran.cc ephem.cc gammln.cc gammp.cc gammq.cc \
gser.cc hitem.cc header.cc date.cc time.cc position.cc star.cc binary.cc \
telescope.cc svdfit.cc filnam.cc tchi.cc input.cc string_funcs.cc \
centroid.cc bsstep.cc mmid.cc pzextr.cc gcf.cc lud.cc gaussj.cc jacob.cc \
sleep.cc sigma_reject.cc what_colour.cc planck.cc fft.cc llsqr.cc fasper.cc \
complex.cc plot.cc formula.cc fraction.cc units.cc format.cc poly.cc rebin.cc \
amoeba.cc genetic.cc rtsafe.cc brent.cc dbrent.cc mnbrak.cc powell.cc \
safunc.cc poisson.cc extinct.cc byte_swap.cc endian.cc boxcar.cc numdiff.cc \
factln.cc runge_kutta.cc voigt.cc stoerm.cc

libsubs_la_LDFLAGS = -version-info 1:0:0

## special rule for this one
##getcomm.o: getcomm.c
##	$(COMPILE) -DPOSIX -c getcomm.c
