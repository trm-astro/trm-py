## Process this file with automake to generate Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src doc include

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; \rm -f $(PACKAGE).tar.gz ; ln -s $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz; \rm -fr $(PACKAGE); mkdir $(PACKAGE)"
	scp -r doc/html $(WEB_SERVE):$(WEB_PATH)/software/$(PACKAGE)/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- roche version --><td>.*</td>%<\!-- roche version --><td>$(VERSION)</td>%' index.html > junk.html; mv junk.html index.html"
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- roche date --><td>.*</td>%<\!-- roche date --><td>$(DATE)</td>%' index.html > junk.html; mv junk.html index.html"

ALIASES = Roche

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
	@echo 'alias algol     $(progdir)/algol'    >> $(ALIASES)
	@echo 'alias cphase    $(progdir)/cphase'   >> $(ALIASES)
	@echo 'alias equalrv   $(progdir)/equalrv'  >> $(ALIASES)
	@echo 'alias impact    $(progdir)/impact'   >> $(ALIASES)
	@echo 'alias lagrange  $(progdir)/lagrange' >> $(ALIASES)
	@echo 'alias pdwd      $(progdir)/pdwd'     >> $(ALIASES)
	@echo 'alias pcont     $(progdir)/pcont'    >> $(ALIASES)
	@echo 'alias pphi      $(progdir)/pphi'     >> $(ALIASES)
	@echo 'alias pshade    $(progdir)/pshade'   >> $(ALIASES)
	@echo 'alias pstream   $(progdir)/pstream'  >> $(ALIASES)
	@echo 'alias vspace    $(progdir)/vspace'   >> $(ALIASES)
	@echo ' ' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Welcome to $(PACKAGE) version $(VERSION), built $(DATE)"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Commands available are: algol, cphase, equalrv, impact, lagrange, pdwd, pcont, pshade, pstream, vspace"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Point a browser at ${prefix}/html/$(PACKAGE)/index.html for help."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
