OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
        `root-config --cflags` \
        -O3 -fno-common -c \
        -I$(INCDIR) -I$(HEEDDIR)

# Debug flags
# CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm
# LDFLAGS += -g

all: thgem analysis

thgem: thgem.C 
	$(CXX) $(CFLAGS) thgem.C
	$(CXX) -o thgem thgem.o $(LDFLAGS)
	rm thgem.o

analysis: analysis.c
	$(CXX) $(CFLAGS) analysis.C
	$(CXX) -o analysis analysis.o $(LDFLAGS)
	rm analysis.o

clean:
	rm thgem
