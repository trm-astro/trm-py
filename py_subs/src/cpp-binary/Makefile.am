## Process this file with automake to generate Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src include doc

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; ln -sf $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz"

ALIASES = Binary

.PHONY : $(ALIASES)

DATE   = $(shell date)

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
	@echo 'alias binfo     $(progdir)/binfo'      >> $(ALIASES)
	@echo 'alias edisc     $(progdir)/edisc'      >> $(ALIASES)
	@echo 'alias wdsolve   $(progdir)/wdsolve'    >> $(ALIASES)
	@echo ' ' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Welcome to $(PACKAGE) version $(VERSION), built $(DATE)"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Commands available are: binfo, edisc, wdsolve"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Point a browser at ${prefix}/html/$(PACKAGE)/index.html for help."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)

clean-local:
	-rm -f $(ALIASES)
