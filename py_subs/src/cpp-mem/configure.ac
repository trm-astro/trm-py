dnl Process this file with autoconf to produce a configure script.

AC_INIT([mem], [1.0.3])

AC_CONFIG_SRCDIR([src/memsys.cc])

AM_INIT_AUTOMAKE

dnl Get rid of all the -D stuff
AM_CONFIG_HEADER([config.h])

AC_PROG_MAKE_SET

dnl Compiler checks
AC_PROG_CXX
AC_PROG_CC

dnl keep correct libtool macros in tree
AC_CONFIG_MACRO_DIR([m4])

dnl Initialise libtool
AC_PROG_LIBTOOL

dnl Make C++ the default
AC_LANG_CPLUSPLUS
AC_REQUIRE_CPP

dnl Check for some header files.
AC_CHECK_HEADERS([stdlib.h float.h math.h sstream string fstream iostream], [], [AC_MSG_ERROR(missing header; please fix)])

dnl Installation program
AC_PROG_INSTALL

dnl The Makefiles to create
AC_CONFIG_FILES([Makefile src/Makefile include/Makefile])

AC_OUTPUT
