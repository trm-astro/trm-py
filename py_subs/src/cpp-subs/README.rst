
.. |subs| replace:: ``subs``

cpp-subs: astronomically-focussed utility routines
==================================================

|subs| is a library of general utility software used by many of my
other packages. After setting up a few third-party packages it will
generally be the first of my packages that you want to install as it
depends upon no others, but is used by many of them. It is therefore
likely to be the first where you run into installation problems, so
this pages includes a section at the end on what you are likely to
encounter and how to fix it.

Installation
------------

You should first make sure of the following third-party packages:
'pcre' (Perl Compatible Regular Expressions), 'slalib' (C version) and
'pgplot'. Make sure these are in place first. 'slalib' is
unfortunately not open source, but its author Pat Wallace allows me to
distribute an obfusticated version, but you will need to ask me for
it. Apologies for this. You will also need a suite of installation
routines called 'autotools'. These will almost certainly be available
on any managed system, but on your own laptop you may need to install
them. Search in whatever package manager you have (yum, yast, zypper,
etc) for keywords such as `autoconf`, `automake` and `libtool`; you
will need all of them. Also if you run your own system, note that you
will the "development" version of pcre which will typically be called
pcre-dev or pcre-devel because you are compiling and linking with it,
not simply using a program already linked to it.

If installing for first time, run './bootstrap' then follow the
instructions in INSTALL. Since these are generic to the autotools,
they are a bit forbidding, so below is some more trouble-shooting
advice which you might want to refer back to as you install others of
my packages, although generally once you have |subs| installed, you
will probably find the others relatively straightforward.

Brief explanation of 'configure'
--------------------------------

If you look in the INSTALL doc you will see that at some point you
will have to invoke the command ``configure`` which will run a script
built by the autotools according to your system. The idea is to try to
cope with all the variations between platforms that makes software
installation painful. The key aspect you need to be aware of is the
configure tries to track down where the header (XXXX.h) and library
(XXXX.a, XXXX.so) files needed to compile and link the software are
located and it performs simple tests to ensure that they work. As it
does so you get lots of 'yes' and 'no' statements printed to the
terminal. For something widely used like ``pcre``, it will probably
locate it automatically because it is located in standard directories
where header or library files are located. For something you might
have installed yourself in a non-standard directory however, configure
has no way of knowing what to do. This is the source of perhaps the
most common problem when installing |subs| which is that although you
have just installed pgplot, it says "PGPLOT header not found; please
fix." or "could not find or work out libraries for PGPLOT".

If you have encountered this problem, then please look up the meaning
and use of the environment variables ``CPPFLAGS`` and ``LDFLAGS`` because
their precise purpose is to direct the compiler and linker to the
right directories. I set these for myself in my .bashrc /
.bash_profile scripts. They depend upon *your* installation and
therefore cannot be set by the installation software. You will only be
able to side-step them if you install in a standard location as root,
but even on my own laptop, I never do this because I like to separate
such "user-software" from the proper stuff.

A typical setting for ``CPPFLAGS`` is::

  export CPPFLAGS="-I/path/to/dir/with/headers1 -I/path/to/dir/with/headers2"

(bash) or equivalently::

  setenv CPPFLAGS "-I/path/to/dir/with/headers1 -I/path/to/dir/with/headers2"

(csh). The paths often end with "include". ``CPPFLAGS`` therefore can be
used to add non-standard locations for header files, e.g. "cpgplot.h"
for compiling PGPLOT. ``LDFLAGS`` is similar but is invoked at linking time
to find libraries. Typical invocations for this are::

  export LDFLAGS="-L/path/to/dir/with/libraries1 -I/path/to/dir/with/libraries2"

or::

  export LDFLAGS=-L/path/to/dir/with/libraries

(bash) or equivalently::

  setenv LDFLAGS "-L/path/to/dir/with/libraries1 -I/path/to/dir/with/libraries2"

or::

  setenv LDFLAGS -L/path/to/dir/with/libraries

(csh). These paths often end with "lib", and should contain files like
libpgplot.a, libpgplot.so, libcpgplot.a.

Thus when configure fails, try to work out whether it has failed to
find a header or a library and act accordingly. In such cases it can
often help to load up 'config.log' in an editor and search for the
package that failed as you will see the explicit test which failed
detailed there.

Good luck!

Tom Marsh
Warwick
