
CC   = gcc
COPT = -g -Wall -O3 -pg #-O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

INCLUDEDIR = -I../../libcompat/include -I../../libthermo/include -I../../libnum/ -I../include/

DEF = -DGCC #-DTRUE_ARRAY

LIBNAME = libcpropep.a

LIBOBJS = equilibrium.o print.o performance.o derivative.o

.SUFFIXES: .c

all: $(LIBNAME)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(LIBNAME): $(LIBOBJS)
	ar -r $@ $(LIBOBJS)
	ranlib $@

clean:
	rm -f *.o *~

deep-clean: clean
	rm -f $(LIBNAME) 
