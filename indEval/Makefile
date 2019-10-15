######################################################################
#
# Makefile to compile IndEval on Linux machines
# Jerome.Kaempf@epfl.ch					    2008-06-24
#
######################################################################

PROG     = indEval
CPPEXT	 = cxx
OBJEXT	 = o
BINDIR   = ~/bin
LIBDIR   = /usr/lib

CPPDIR	 = /usr/bin
CPP	 = $(CPPDIR)/g++

# flags

CPPFLAGS = -O2 -march=native -c
LDFLAGS =  -o $@
LIBS    =  -L$(LIBDIR)

# sources, objects and so on

SRCS =	../geometry.$(CPPEXT) ../problem.$(CPPEXT) main.$(CPPEXT)

OBJS =	geometry.$(OBJEXT) problem.$(OBJEXT) main.$(OBJEXT)

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

main.$(OBJEXT): geometry.$(OBJEXT) problem.$(OBJEXT)

geometry.$(OBJEXT): ../geometry.h ../geometry.$(CPPEXT)
	$(CPP) $(CPPFLAGS) ../geometry.$(CPPEXT) -o geometry.$(OBJEXT) 	

problem.$(OBJEXT): ../problem.h ../problem.$(CPPEXT)
	$(CPP) $(CPPFLAGS) ../problem.$(CPPEXT) -o problem.$(OBJEXT)
