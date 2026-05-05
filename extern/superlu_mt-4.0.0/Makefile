############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   August 15, 1997
#
#  Modified:   	    September 1, 1999 version 1.0
#                   March 20, 2013  version 2.1
#
############################################################################

include make.inc

all: install lib 

lib: superlulib 

clean: cleanlib

install:
	( cd INSTALL; $(MAKE) )
#	( cd INSTALL; cp lsame.c ../SRC/; \
#	  cp dlamch.c ../SRC/; cp slamch.c ../SRC/ )

blaslib:
	 ## built blas lib routes to "lib"...error thrown if directory does NOT exist
	mkdir -p lib
	( cd CBLAS; $(MAKE) )

superlulib:
	( cd SRC; $(MAKE) )

cleanlib:
	( cd SRC; $(MAKE) clean )
	( cd CBLAS; $(MAKE) clean )
	( cd lib; rm -f *.a )
	( rm -f *.a )

