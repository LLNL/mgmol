MAIN = ./main
GCCDIR=/usr/local/tools/gcc-4.6.1
LD_LIBRARY_PATH = -L$(GCCDIR)/lib
LD = $(CXX)


OBJECTS = \
 spline.o \
 splint.o \
 log2lingrid.o
 
%.o: %.cc
	$(CXX) -g $(CPPFLAGS) -c $<


$(MAIN): $(OBJECTS)
	$(LD) $(LD_LIBRARY_PATH) $(LDFLAGS) $(OBJECTS) -lm \
                   -o $(MAIN)
 
clean:
	rm -f *.o 
	rm -f main
