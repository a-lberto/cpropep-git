
CC   = gcc
COPT = -g -Wall -O3 #-pg -O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

LIB  =  -lthermo #-lcpropep -lm -lnum #-lcruft -lf2c
LIBDIR = -L../../libnum/ -L../lib/
INCLUDEDIR = -I../../libnum/ -I../lib/

DEF = -DGCC #-DTRUE_ARRAY

PROG1 = prop
PROG2 = plot

OBJS1 = prop.o 
OBJS2 = plot.o

.SUFFIXES: .c

all: $(PROG1) $(PROG2)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(PROG1): $(OBJS1) $(LIBNAME)
	$(CC) $(COPT) $(OBJS1) $(LIBDIR) $(LIB) -o $@

$(PROG2): $(OBJS2) $(LIBNAME)
	$(CC) $(COPT) $(OBJS2) $(LIBDIR) $(LIB) -o $@

clean:
	rm -f *.o *~ *.m

deep-clean: clean
	rm -f $(PROG1) $(PROG2)
