# Makefile created by JK, 22/02/09

.SUFFIXES: .cxx

CFLAGS = -O2 # For Cygwin

GLLIBS =

libs = ${GLLIBS} # Cygwin or Mingw

EXEEXT=

all:	sprngTest

sprngTest:	../geometry.o	main.o
	g++ $(CFLAGS) geometry.o main.o $(libs) -o sprngTest
	cp sprngTest$(EXEEXT) ~/bin

.cxx.o:
	g++ -c $(CFLAGS) $*.cxx

clean:
	/bin/rm -f *.o
	/bin/rm -f sprngTest$(EXEEXT)
	/bin/rm -f ~/bin/sprngTest$(EXEEXT)
