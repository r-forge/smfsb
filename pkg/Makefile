# Makefile
# Makefile for package checking building, installing, uninstalling, etc.

VERSION=1.3
VERSIONSBML=0.1

DEFAULT:
	make check
	make install

check:
	R CMD check smfsb

check-cran:
	R CMD check --as-cran smfsb

build:
	R CMD build smfsb

install:
	make build
	R CMD INSTALL smfsb_$(VERSION).tar.gz

install-pkg:
	R CMD INSTALL smfsb

remove:
	R CMD REMOVE smfsb

clean:
	rm -rf *~ smfsb_*.tar.gz smfsb.Rcheck


update:
	svn update
	svn log|less

commit:
	svn commit
	make update



# eof

