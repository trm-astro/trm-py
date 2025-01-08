## Process this file with automake to generate Makefile.in

EXTRA_DIST = bootstrap

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src include

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; ln -sf $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz"

## File of aliases

ALIASES = Subs

.PHONY : $(ALIASES)

DATE    = $(shell date)

progdir = @bindir@/@PACKAGE@

prog_DATA = $(ALIASES)

$(ALIASES):
	@echo '#' > $(ALIASES)
	@echo '# This file defines aliases for the $(PACKAGE) v$(VERSION) software.' >> $(ALIASES)
	@echo '# Source it to define them.' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo '# this to allow both bash and csh to work' >> $(ALIASES)
	@echo 'test "$$?BASH_VERSION" = "0" || eval '\''alias() { command alias "$$1=$$2"; }'\' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo 'alias tcorr     $(progdir)/tcorr'       >> $(ALIASES)
	@echo 'alias weekday   $(progdir)/weekday'     >> $(ALIASES)
	@echo 'alias gap       $(progdir)/gap'         >> $(ALIASES)
	@echo ' ' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Welcome to $(PACKAGE) version $(VERSION), built $(DATE)."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Commands available are: gap, interp, tcorr, weekday."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Point a browser at ${prefix}/html/$(PACKAGE)/index.html for help."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)

clean-local:
	-rm -f $(ALIASES)
	-rm -rf autom4te.cache
	-rm -f compile config.guess config.h config.log config.status libtool config.h.in config.h.in~ config.sub configure configure~ configure.scan autoscan.log Makefile Makefile.in include/Makefile.in include/Makefile src/Makefile.in src/Makefile missing COPYING ltmain.sh INSTALL install-sh depcomp aclocal.m4 m4/libtool.m4 m4/lt~obsolete.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 stamp-h1



