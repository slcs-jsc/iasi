# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = bands extract noise perturbation spec2tab

# Library directories...
LIBDIR = -L ../lib/build/lib

# Include directories...
INCDIR = -I ../lib/build/include

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc
MPICC = mpicc

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED $(DEFS) -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp -g

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lcoda -lm

# LDFLAGS for retrieval...
LDFLAGS_RET = $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lm

# Profiling...
ifdef PROF
  CFLAGS += -O2 -pg 
else
  CFLAGS += -O3
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

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

dist:
	zip iasi_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC) retrieval
