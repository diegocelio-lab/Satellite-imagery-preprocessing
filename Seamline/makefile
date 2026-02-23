# COMPATIBLE WITH 64BIT 

# C source code
CSRC	= 

# C++ source code
CXXSRC	= main.cpp seamline.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= seamline

default	: $(BIN)

# use DEBUG mode with `make DEBUG=1`
ifdef DEBUG
  # C optimization flags
  COPT = -g
else
  # C optimization flags
  COPT	= -O3 -ftree-vectorize -funroll-loops
endif


# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT)  -std=c99 -m64
# C++ compilation flags
# option -DNO_XML needed to compile fast bilateral filter 
# implementation by Sylvain Paris and Frédo Durand
CXXFLAGS = $(CXXOPT) -std=c++17 -DNO_XML -m64 $(shell gdal-config --cflags)
# link flags
LDFLAGS	= -lpng -ltiff -m64 $(shell gdal-config --libs)


CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas


# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS) 
# partial compilation of C++ source code
%.o: %.cpp
	$(CXX) -c -o $@  $<  $(CXXFLAGS) 

# link all the object code
$(BIN): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
distclean	: clean
	$(RM) $(BIN) *.o *~ 

