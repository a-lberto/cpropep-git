


CC   = gcc
COPT = -g -Wall
LIB  = -lcruft -lf2c


PROG = num

OBJS = num.o libnum.o
LIBNUM = libnum.a

.SUFFIXES: .c

all: $(PROG) $(LIBNUM)

.c.o:
	$(CC) $(COPT) -c $*.c -o $*.o

$(LIBNUM): libnum.o
	ar -r $@ libnum.o
	ranlib $@

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIB) -o $@


clean:
	rm -f *.o *~
