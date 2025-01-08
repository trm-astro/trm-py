dnl Process this file with autoconf to produce a configure script.

AC_INIT([colly], [1.1])

AC_CONFIG_SRCDIR([src/grab.cc])

AM_INIT_AUTOMAKE

dnl Get rid of all the -D stuff
AM_CONFIG_HEADER(config.h)

AC_PROG_MAKE_SET

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC

dnl Make C++ the default
AC_LANG_CPLUSPLUS
AC_REQUIRE_CPP

dnl Initialise libtool
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_LIBTOOL

dnl headers

AC_CHECK_HEADERS([trm/subs.h trm/position.h trm/constants.h trm/array1d.h],[],
		 [AC_MSG_ERROR([cannot find 'subs' headers])])

AC_CHECK_HEADERS([stdlib.h float.h math.h sstream string fstream iostream], [],
                 [AC_MSG_ERROR(missing header; please fix)])

AC_CHECK_HEADERS([pcrecpp.h sofa.h], [],
                 [AC_MSG_ERROR(missing header; please fix)])

dnl libraries

AC_CHECK_LIB([pcrecpp], [main], [],
             [AC_MSG_ERROR(cannot link to the pcrecpp library)])

AC_CHECK_LIB([sofa_c], [main], [],
             [AC_MSG_ERROR(cannot find the sla C library)])

AC_CHECK_LIB([subs], [main], [],
             [AC_MSG_ERROR([cannot link to 'subs' library])])

dnl PGPLOT has its own macro 'cos its a pain
TRM_LIB_PGPLOT

dnl Installation program
AC_PROG_INSTALL

dnl 
AC_CONFIG_SUBDIRS([subs])

dnl The Makefiles to create
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile doc/Makefile]) 

AC_OUTPUT
