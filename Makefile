


CC   = gcc
COPT = -g -Wall
LIB  = -lm


PROG = cpropep

OBJS = cpropep.o

.SUFFIXES: .c

.c.o:
	$(CC) $(COPT) -c $*.c -o $*.o

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIB) -o $@


clean:
	rm *.o *~