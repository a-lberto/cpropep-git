


CC   = gcc
COPT = -g -Wall
LIB  = -lm -lnum
LIBDIR = -L../libnum/
INCLUDEDIR = -I../libnum/

PROG = cpropep

OBJS = cpropep.o load.o

.SUFFIXES: .c

.c.o:
	$(CC) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm *.o *~