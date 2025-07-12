# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = bands extract noise perturbation spec2tab

# Installation directory...
DESTDIR ?= ../bin

# Library directories...
LIBDIR += -L ../libs/build/lib

# Include directories...
INCDIR += -I ../libs/build/include

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags)

# Compiler...
CC ?= gcc
MPICC ?= mpicc.openmpi

# CFLAGS...
CFLAGS ?= $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums $(OPT) -g -fopenmp

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lcoda -lm

# LDFLAGS for retrieval...
LDFLAGS_RET ?= $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lm

# Optimization information...
ifeq ($(INFO),1)
  CFLAGS += -fopt-info
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Static compilation...
ifeq ($(STATIC),1)
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doxygen indent install lizard mkdocs strip uninstall

all: $(EXC) retrieval
	rm -f *~

$(EXC): %: %.c libiasi.o jurassic.o
	$(CC) $(CFLAGS) -o $@ $< libiasi.o jurassic.o $(LDFLAGS) -static

libiasi.o: libiasi.c libiasi.h Makefile
	$(CC) $(CFLAGS) -c -o libiasi.o libiasi.c 

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

retrieval: retrieval.c jurassic.o
	$(MPICC) $(CFLAGS) -o retrieval retrieval.c jurassic.o $(LDFLAGS_RET)

clean:
	rm -f $(EXC) retrieval *.o *~

cppcheck:
	cppcheck --enable=all ./

dist:
	zip iasi_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC) retrieval

uninstall:
	cd $(DESTDIR) && rm $(EXC)
