


CC   = gcc
COPT = -g -Wall
LIB  = 


PROG = num

OBJS = num.o libnum.o

.SUFFIXES: .c

.c.o:
	$(CC) $(COPT) -c $*.c -o $*.o

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIB) -o $@


clean:
	rm -f *.o *~
