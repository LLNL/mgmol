# Makefile
mgmol:
	( cd src; $(MAKE) opt )

debug:
	( cd src; $(MAKE) debug )

#------------------------------------------------------------------------------
clean:
	( cd src; $(MAKE) clean)

depend:
	( cd src; $(MAKE) depend)

