
CC   = gcc
COPT = -g -Wall -O3 -pg #-O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

LIB  = -lnum -lm # -lcruft -lf2c
LIBDIR = -L.

DEF = -DLINUX #-DTRUE_ARRAY
#DEF += -DLAPACK
PROG = test

OBJS = test.o
LIBOBJS = lu.o rk4.o rkf.o general.o print.o sec.o newton.o ptfix.o
LIBNUM = libnum.a

.SUFFIXES: .c

all: $(LIBNUM) $(PROG)

.c.o:
	$(CC) $(DEF) $(COPT) -c $*.c -o $*.o

$(LIBNUM): $(LIBOBJS)
	ar -r $@ $(LIBOBJS) 
	ranlib $@

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm -f *.o *~

deep-clean: clean
	rm -f $(LIBNUM) $(PROG) 
