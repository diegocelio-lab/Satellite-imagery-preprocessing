# Makefile for LLCC
#
# Original Copyright 2011 Jose-Luis Lisani, joseluis.lisani@uib.es
# Modifications 2026 by Diego Celio, diego.celio@me.com

# 
# This updated version has been modified to support:
#   - C++17 compilation
#   - GDAL and TIFF image support
#   - Qt5 XML (QDom) support
#   - PNG support
#   - WSL / modern Ubuntu environments
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the original copyright notice and this notice are preserved.
# This file is offered as-is, without any warranty.


# C source code
CSRC = library/io_png.c \
       library/MCMMondelli/FDS_MCM.c \
       library/gaussianGetreuer/gaussian_conv_sii.c \
       library/gaussianGetreuer/gaussian_conv.c \
       
# C++ source code
CXXSRC = llcc.cpp

# All source code
SRC = $(CSRC) $(CXXSRC)

# C objects
COBJ = $(CSRC:.c=.o)
# C++ objects
CXXOBJ = $(CXXSRC:.cpp=.o)
# All objects
OBJ = $(COBJ) $(CXXOBJ)

# Binary target
BIN = llcc

# Default target
default: $(BIN)

# Debug mode
ifdef DEBUG
  COPT = -g
else
  COPT = -O3 -ftree-vectorize -funroll-loops
endif

# Compiler flags
CFLAGS = $(COPT) -std=c99 -m64 -Wno-unknown-pragmas \
         $(shell gdal-config --cflags) \
         $(shell pkg-config --cflags Qt5Xml)

CXXFLAGS = $(COPT) -std=c++17 -DNO_XML -m64 -Wno-unknown-pragmas \
           $(shell gdal-config --cflags) \
           $(shell pkg-config --cflags Qt5Xml)

# Linker flags
LDFLAGS = -lpng -ltiff $(shell gdal-config --libs) $(shell pkg-config --libs Qt5Xml) -m64

# Compile C files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

# Compile C++ files
%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# Link all objects
$(BIN): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Housekeeping
.PHONY: clean distclean
clean:
	$(RM) $(OBJ)
distclean: clean
	$(RM) $(BIN) *.o *~