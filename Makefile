
CC   = gcc
COPT = -g -Wall -pg -O6\
 -mpentium -ffast-math -funroll-loops -fnonnull-objects\
 -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
 -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

LIB  =  -lcpropep -lm -lnum -lcruft -lf2c
LIBDIR = -L../libnum/ -L.
INCLUDEDIR = -I../libnum/ -I.

DEF = -DGCC -DTRUE_ARRAY -D_GNU

LIBNAME = libcpropep.a

PROG = cpropep

LIBOBJS = equilibrium.o load.o print.o performance.o derivative.o

OBJS = cpropep.o getopt.o 

.SUFFIXES: .c

all: $(LIBNAME) $(PROG)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(LIBNAME): $(LIBOBJS)
	ar -r $@ $(LIBOBJS)
	ranlib $@

$(PROG): $(OBJS) $(LIBNAME)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@

clean:
	rm -f *.o *~
	make -C cgi clean

deep-clean: clean
	rm -f $(PROG) $(LIBNAME)
	make -C cgi deep-clean
