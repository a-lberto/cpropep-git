CC   = gcc
COPT = -g -Wall
LIB  =  -lcpropep -lm -lnum -lcgi
LIBDIR = -L../../libnum/ -L../
INCLUDEDIR = -I../../libnum/ -I../

PROG = cgitest

OBJS = cgitest.o


.SUFFIXES: .c

all: $(PROG)

.c.o:
	$(CC) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o


$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm -f *.o *~