dnl Process this file with autoconf to produce a configure script.

dnl Initialise
AC_INIT([subs],[2.1])

AC_CONFIG_SRCDIR([src/gauss.cc])

AM_INIT_AUTOMAKE

dnl Avoid all the -D stuff
AC_CONFIG_HEADERS([config.h])

AC_PROG_MAKE_SET

dnl Compiler checks
AC_PROG_CXX
AC_PROG_CC

dnl keep correct libtool macros in tree
AC_CONFIG_MACRO_DIR([m4])

dnl Initialise libtool
AC_PROG_LIBTOOL

dnl Make C++ the default
AC_LANG([C++])
AC_REQUIRE_CPP

dnl headers
AC_CHECK_HEADERS([pcrecpp.h sofa.h], [], [AC_MSG_ERROR(missing header; please fix)])

AC_CHECK_HEADERS([stdlib.h float.h math.h sstream string fstream iostream], [], [AC_MSG_ERROR(missing header; please fix)])

AC_CHECK_HEADERS([pcrecpp.h sofa.h], [], [AC_MSG_ERROR(missing header; please fix)])


dnl libraries.
AC_CHECK_LIB([sofa_c], [main], [], [AC_MSG_ERROR(cannot link to the sofa C library)])

AC_CHECK_LIB([pcrecpp], [main], [], [AC_MSG_ERROR(cannot link to the pcrecpp library)])

dnl PGPLOT has its own macro 'cos its a pain
TRM_LIB_PGPLOT

dnl Installation program
AC_PROG_INSTALL

dnl The Makefiles to create

dnl AC_CONFIG_FILES([Makefile src/Makefile include/Makefile doc/Makefile])
AC_CONFIG_FILES([Makefile src/Makefile include/Makefile])

AC_OUTPUT
