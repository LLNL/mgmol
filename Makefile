# Makefile
# $Id: Makefile,v 1.10 2009/11/19 17:36:07 jeanluc Exp $
mgmol:
	( cd src; $(MAKE) opt )

tests:
	( cd testing; $(MAKE) tests )

mgmol_debug:
	( cd src; $(MAKE) debug )

mgmol_serial:
	( cd src; $(MAKE) mgmol_serial )

mgmol_serial_debug:
	( cd src; $(MAKE) mgmol_serial_dbg )

comparewf:
	( cd src; $(MAKE) comparewf )

comparewf_mpi:
	( cd src; $(MAKE) comparewf_mpi )

#------------------------------------------------------------------------------
clean:
	( cd src; $(MAKE) clean)
	( cd testing; $(MAKE) clean)

depend:
	( cd src; $(MAKE) depend)

test_mpi:
	( cd testing; $(MAKE) test_mpi)
test_clean:
	( cd testing; $(MAKE) test_clean)

ctags :
	/usr/bin/ctags src/*.cc src/*.C include/*.h
