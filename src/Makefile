
CC   = gcc
COPT = -g -Wall -O3 -pg #-O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

INCLUDEDIR = -I../../libcpropep/include -I../../libcompat/include -I../../libnum/ -I../include/

DEF = -DGCC #-DTRUE_ARRAY

LIBNAME  = libthermo.a

LIBOBJS  = load.o thermo.o

.SUFFIXES: .c

all: $(LIBNAME)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(LIBNAME): $(LIBOBJS)
	ar -r $@ $(LIBOBJS)
	ranlib $@
	
clean:
	rm -f *.o *~

deep-clean: clean
	rm -f $(LIBNAME)
