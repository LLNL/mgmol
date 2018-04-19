# Makefile
# $Id:$
CURRENT_DIR = $(PWD)
INSTALL_DIR = $(CURRENT_DIR)/../install

include $(INSTALL_DIR)/Makefile.in


readhdf5 = ../bin/read_hdf5_$(PLAT)

READHDF_OBJECTS = read_hdf5.o

LD = $(CXX)
CPPFLAGS = $(MACHINE_FLAGS) $(DEFS) -I$(HDF_DIR_SERIAL)/include
 
$(readhdf5): $(READHDF_OBJECTS) 
	$(LD) $(LDFLAGS) $(READHDF_OBJECTS) \
		$(HDF_LIB_SERIAL) $(COMMONLIBS)\
		-o $(readhdf5)
