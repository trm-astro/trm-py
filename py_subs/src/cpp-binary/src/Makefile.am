## Process this file with automake to generate Makefile.in
##
## This is the file that must be edited if you are changing anything in the source directory

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS   = binfo edisc wdsolve

binfo_SOURCES   = binfo.cc
edisc_SOURCES   = edisc.cc
wdsolve_SOURCES = wdsolve.cc

AM_CPPFLAGS = -I../include -I../.

LDADD = libbinary.la

## Library

lib_LTLIBRARIES = libbinary.la 

libbinary_la_SOURCES = eclipses.cc  orbital_ang_mom.cc rcirc.cc jdotgr.cc  mr_main_sequence.cc  \
orbital_omega.cc rmin.cc kepler_omega.cc mr_wd_eggleton.cc orbital_period.cc zeta_wd_eggleton.cc \
mr_wd_nauenberg.cc orbital_separation.cc zeta_main_sequence.cc zeta_wd_nauenberg.cc gyrate.cc vellipse.cc

libbinary_la_LDFLAGS = -version-info 1:0:0

DATE = $(shell date)

install-data-hook: 
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION
