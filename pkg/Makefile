# Makefile
# Makefile for package checking building, installing, uninstalling, etc.

VERSION=1.5
VERSIONSBML=0.1

DEFAULT:
	make check
	make install

check:
	R CMD check smfsb
	R CMD check smfsbSBML

check-cran:
	make build
	R CMD check --as-cran smfsb_$(VERSION).tar.gz
	R CMD check --as-cran smfsbSBML_$(VERSIONSBML).tar.gz

check-all:
	R CMD check --as-cran --run-donttest smfsb

build:
	R CMD build smfsb
	R CMD build smfsbSBML

install:
	make build
	R CMD INSTALL smfsb_$(VERSION).tar.gz
	R CMD INSTALL smfsbSBML_$(VERSIONSBML).tar.gz

web:
	make build
	scp smfsb_$(VERSION).tar.gz smfsbSBML_$(VERSIONSBML).tar.gz @unix.ncl.ac.uk:public_html/smfsb/3e/

install-pkg:
	R CMD INSTALL smfsb
	R CMD INSTALL smfsbSBML

remove:
	R CMD REMOVE smfsb
	R CMD REMOVE smfsbSBML

clean:
	rm -rf *~ smfsb*.tar.gz *.Rcheck


update:
	svn update
	svn log|less

commit:
	svn commit
	make update



# eof

