# Makefile
# Makefile for package checking building, installing, uninstalling, etc.


check:
	R CMD check pkg

build:
	R CMD build pkg

clean:
	rm -rf *~ smfsb_*.tar.gz pkg.Rcheck

update:
	svn update
	svn log|less

commit:
	svn commit
	make update



# eof

