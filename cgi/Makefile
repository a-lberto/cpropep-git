CC   = gcc
COPT = -g -Wall
LIB  =  -lcpropep -lm -lnum -lcgi #-lcruft -lf2c
LIBDIR = -L../../libnum/ -L../ -L./
INCLUDEDIR = -I../../libnum/ -I../ -I./

PROG = cgitest

OBJS = cgitest.o


.SUFFIXES: .c

all: $(PROG)

.c.o:
	$(CC) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o


$(PROG): $(OBJS)
	make -C ../ all
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm -f *.o *~

deep-clean: clean
	rm -f cgitest
