
CC   = gcc
COPT = -g -Wall -O3 #-pg -O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

LIB  =  -lcpropep -lthermo -lm -lnum #-lcruft -lf2c
LIBDIR = -L../../libnum/ -L../lib/
INCLUDEDIR = -I../../libnum/ -I../lib/

DEF = -DGCC #-DTRUE_ARRAY

PROG = cpropep

OBJS = cpropep.o getopt.o 

.SUFFIXES: .c

all: $(PROG)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(PROG): $(OBJS) $(LIBNAME)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@

clean:
	rm -f *.o *~

deep-clean: clean
	rm -f $(PROG)
