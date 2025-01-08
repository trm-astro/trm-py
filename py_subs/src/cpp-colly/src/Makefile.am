## Process this file with automake to generate Makefile.in
##

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS   = correct grab group reorder mlist

correct_SOURCES = correct.cc
grab_SOURCES    = grab.cc
group_SOURCES   = group.cc
mlist_SOURCES   = mlist.cc
reorder_SOURCES = reorder.cc

AM_CPPFLAGS = -I../include -I../.

LDADD = libcolly.la

## Library

lib_LTLIBRARIES = libcolly.la 

libcolly_la_SOURCES = write_molly_head.cc read_molly_head.cc find_groups.cc rejig_molly_head.cc skip_molly_data.cc skip_bytes.cc


DATE = $(shell date)

install-data-hook: 
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION

