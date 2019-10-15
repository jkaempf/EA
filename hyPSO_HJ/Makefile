######################################################################
#
# Makefile to compile hyPSO_HJ on Linux machines
# Jerome.Kaempf@epfl.ch					    2008-06-24
#
######################################################################

PROG     = hyPSO_HJ
CPPEXT	 = cxx
OBJEXT	 = o
BINDIR   = ~/bin
LIBDIR   = /usr/lib

CPPDIR	 = /usr/bin
CPP	 = $(CPPDIR)/g++

# flags

CPPFLAGS = -DNTHREADS=16 -O2 -march=native -c
LDFLAGS =  -o $@
LIBS    =  -L$(LIBDIR) -lpthread

# sources, objects and so on

SRCS =	../geometry.$(CPPEXT) ../problem.$(CPPEXT) heuristic.$(CPPEXT) evalpthread.$(CPPEXT) main.$(CPPEXT)

OBJS =	geometry.$(OBJEXT) problem.$(OBJEXT) heuristic.$(OBJEXT) evalpthread.$(OBJEXT) main.$(OBJEXT)

###################################################
all: $(PROG)
	cp $(PROG) $(BINDIR)/$(PROG)


$(PROG): $(OBJS) 
	$(CPP) $(LDFLAGS) $(OBJS) $(LIBS) 


.$(CPPEXT).$(OBJEXT):
	$(CPP) $(CPPFLAGS) $<


clean:
	rm -f $(PROG) *.o $(BINDIR)/$(PROG)


.SUFFIXES : .$(OBJEXT) .$(CPPEXT)

###################################################

main.$(OBJEXT): geometry.$(OBJEXT) problem.$(OBJEXT) evalpthread.$(OBJEXT) heuristic.$(OBJEXT)

geometry.$(OBJEXT): ../geometry.h ../geometry.$(CPPEXT)
	$(CPP) $(CPPFLAGS) ../geometry.$(CPPEXT) -o geometry.$(OBJEXT) 	

problem.$(OBJEXT): ../problem.h ../problem.$(CPPEXT)
	$(CPP) $(CPPFLAGS) ../problem.$(CPPEXT) -o problem.$(OBJEXT)

evalpthread.$(OBJEXT): ../evalpthread.h ../evalpthread.$(CPPEXT)
	$(CPP) $(CPPFLAGS) ../evalpthread.$(CPPEXT) -o evalpthread.$(OBJEXT)

heuristic.$(OBJEXT): ../heuristic.h ../heuristic.$(CPPEXT) evalpthread.$(OBJEXT)
	$(CPP) $(CPPFLAGS) ../heuristic.$(CPPEXT) -o heuristic.$(OBJEXT)

