# C++ source code
CXXSRC  = library/parser.cpp simplest_colour_balance.cpp
# C source code
CSRC    = library/io_png.c library/auxiliary_lib.c

# All objects
CXXOBJ  = $(CXXSRC:.cpp=.o)
COBJ    = $(CSRC:.c=.o)
OBJ     = $(CXXOBJ) $(COBJ)

# Binary target
BIN     = simplest_colour_balance

# C++ compilation flags
CXXFLAGS= -O3 -Wall -Wextra -Wno-write-strings -std=c++17 -I/opt/local/include/ -I/usr/local/include/ $(shell gdal-config --cflags)
# C compilation flags
CFLAGS  = -O3 -Wall -Wextra -Wno-write-strings -std=c99 -I/opt/local/include/ -I/usr/local/include/ $(shell gdal-config --cflags)

# Link flags
LDFLAGS = -lfftw3 -lpng -ltiff -lm -L/opt/local/lib/ -L/usr/local/lib/ $(shell gdal-config --libs)

# Compile C++ source
%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# Compile C source
%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

# Link all objects
$(BIN): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: clean distclean
clean:
	$(RM) $(OBJ)
distclean: clean
	$(RM) $(BIN)
