# Copyright 2009, 2010 IPOL Image Processing On Line http://www.ipol.im/
# Author: Thierry Guillemot <thierry.guillemot.work@gmail.com> and Julie Delon <julie.delon@parisdescarte.fr>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.

# -----------------------------------------
# Makefile for Midway_img_equ (PNG + TIFF support)
# AddressSanitizer enabled
# -----------------------------------------

# Compiler
CXX = g++
CC  = gcc

# Directories
DIR       = bin obj
SRCUTILS  = src/Utilities
SRCMAIN   = src/Midway

# Compiler options
CXXOPT    = -O3 -funroll-loops -fomit-frame-pointer
SANITIZE_FLAGS = -g -fsanitize=address -fno-omit-frame-pointer

# Include paths
CXXFLAGS  = $(CXXOPT) -Wall -Wextra -I $(SRCUTILS) $(shell gdal-config --cflags) $(SANITIZE_FLAGS)
CFLAGS    = $(CXXOPT) -Wall -Wextra -I $(SRCUTILS) $(shell gdal-config --cflags) $(SANITIZE_FLAGS)

# Linker flags
LDFLAGS   = -lpng -ltiff $(shell gdal-config --libs) $(SANITIZE_FLAGS)

# Targets
all: $(DIR) bin/midway
default: all
.PHONY: clean distclean

# --- Create directories ---
$(DIR):
	mkdir -p $(DIR)

# --- C files ---
obj/mt19937ar.o	: $(SRCUTILS)/mt19937ar.c $(SRCUTILS)/mt19937ar.h
	$(CC) -c -o $@ $(CFLAGS) $<

obj/io_png.o : $(SRCUTILS)/io_png.c $(SRCUTILS)/io_png.h
	$(CC) -c -o $@ $(CFLAGS) $<

# --- CPP files ---
obj/LibImages.o	: $(SRCUTILS)/LibImages.cpp $(SRCUTILS)/LibImages.h
	$(CXX) -c -o $@ $(CXXFLAGS) $<

obj/midway.o : $(SRCMAIN)/midway.cpp $(SRCMAIN)/midway.h $(SRCUTILS)/utils.h
	$(CXX) -c -o $@ $(CXXFLAGS) $<

obj/main.o	: src/main.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<

# --- Link binaries ---
bin/midway : obj/mt19937ar.o obj/LibImages.o obj/io_png.o obj/midway.o obj/main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

# --- Clean rules ---
clean :
	$(RM) -r obj
distclean : clean
	$(RM) -r bin


