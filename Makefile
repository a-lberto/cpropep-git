


CC   = gcc
COPT = -g -Wall
LIB  =  -lm -lnum 
LIBDIR = -L../libnum/
INCLUDEDIR = -I../libnum/

PROG = cpropep

OBJS = equilibrium.o cpropep.o load.o getopt.o 

.SUFFIXES: .c

.c.o:
	$(CC) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm -f *.o *~
