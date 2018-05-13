# Makefile
# Makefile for package checking building, installing, uninstalling, etc.

VERSION=1.2

DEFAULT:
	make check
	make install

check:
	R CMD check smfsb

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
	rm -rf *~ smfsb_*.tar.gz pkg.Rcheck


update:
	svn update
	svn log|less

commit:
	svn commit
	make update



# eof

