dnl Process this file with autoconf to produce a configure script.

AC_INIT([roche], [1.0.12])

AC_CONFIG_SRCDIR([src/gsint.cc])

dnl Get rid of all the -D stuff
AM_CONFIG_HEADER(config.h)

dnl automake stuff
AM_INIT_AUTOMAKE

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
                 [AC_MSG_ERROR(cannot find headers associated with subs)])

AC_CHECK_HEADERS([trm/binary.h],[],
                 [AC_MSG_ERROR(cannot find the binary header)])

AC_CHECK_HEADERS([stdlib.h float.h math.h sstream string fstream iostream], [],
                 [AC_MSG_ERROR(missing header; please fix)])

AC_CHECK_HEADERS([pcrecpp.h sofa.h], [],
                 [AC_MSG_ERROR(missing header; please fix)])

dnl libraries

AC_CHECK_LIB([pcrecpp], [main], [],
             [AC_MSG_ERROR(cannot link to the pcrecpp library)])

AC_CHECK_LIB([sofa_c], [main], [],
             [AC_MSG_ERROR(cannot find the sofa C library)])

AC_CHECK_LIB([subs], [main], [],
             [AC_MSG_ERROR(cannot link to the subs library)])

AC_CHECK_LIB([binary], [main], [],
             [AC_MSG_ERROR(cannot link to the binary library)])

dnl PGPLOT has its own macro 'cos its a pain
TRM_LIB_PGPLOT

dnl Installation program
AC_PROG_INSTALL

dnl The Makefiles to create

AC_OUTPUT([Makefile include/Makefile src/Makefile doc/Makefile]) 
