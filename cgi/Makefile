CC   = gcc
COPT = -g -Wall
LIB  =  -lcpropep -lthermo -lm -lnum -lcgi #-lcruft -lf2c
LIBDIR = -L../../libnum/ -L../lib/ -L./
INCLUDEDIR = -I../../libnum/ -I../lib/ -I./

DEF = -DGCC -DSOURCEFORGE 

PROG = cgitest

OBJS = cgitest.o


.SUFFIXES: .c

all: $(PROG)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o


$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm -f *.o *~

deep-clean: clean
	rm -f cgitest

upload:
	scp $(PROG) antoine@rocketworkbench.sourceforge.net:/home/groups/rocketworkbench/cgi-bin/
