## Process this file with automake to generate Makefile.in
##
## This is the file that must be edited if you are changing anything in the source directory

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS    = algol impact lagrange pdwd cphase pstream pcont pphi equalrv vspace pshade speed_test

algol_SOURCES    = algol.cc 
cphase_SOURCES   = cphase.cc
equalrv_SOURCES  = equalrv.cc
impact_SOURCES   = impact.cc
lagrange_SOURCES = lagrange.cc
pcont_SOURCES    = pcont.cc
pdwd_SOURCES     = pdwd.cc 
pphi_SOURCES     = pphi.cc
pshade_SOURCES   = pshade.cc
pstream_SOURCES  = pstream.cc
vspace_SOURCES   = vspace.cc
speed_test_SOURCES = speed_test.cc

AM_CPPFLAGS = -I../include

LDADD = libroche.la

DATE = $(shell date)

install-data-hook:
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION

## Library

lib_LTLIBRARIES = libroche.la 

libroche_la_SOURCES = blink.cc hits.cc rpot.cc gsint.cc rdot.cc strinit.cc rlobe_eggleton.cc irrad.cc rocacc.cc \
xl1.cc xl2.cc xl3.cc lobe1.cc lobe2.cc flobe1.cc streamr.cc stradv.cc strmnx.cc drpot.cc vtrans.cc vlobe1.cc vlobe2.cc \
zeta_rlobe_eggleton.cc contact.cc vstream.cc sphere_eclipse.cc ref_sphere.cc pot_min.cc rpot_val.cc rpot_grad.cc \
rpot_val_grad.cc fblink.cc ingress_egress.cc face.cc set_earth.cc roche_shadow.cc disc_eclipse.cc stream.cc vstrreg.cc \
jacobi.cc xl11.cc xl12.cc



