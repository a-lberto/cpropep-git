


CC   = gcc
COPT = -g -Wall
LIB  =  -lcpropep -lm -lnum -lcruft -lf2c
LIBDIR = -L../libnum/ -L.
INCLUDEDIR = -I../libnum/ -I.

LIBNAME = libcpropep.a

PROG = cpropep

LIBOBJS = equilibrium.o load.o print.o performance.o derivative.o

OBJS = cpropep.o getopt.o 

.SUFFIXES: .c

all: $(LIBNAME) $(PROG)

.c.o:
	$(CC) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

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
