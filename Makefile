
CC   = gcc
COPT = -g -Wall -O3 #-pg -O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

DEF = -DGCC #-DTRUE_ARRAY

all:
	make -C lib all
	make -C cpropep all
	make -C cgi all
	make -C prop all


clean:
	make -C lib clean
	make -C cpropep clean
	make -C cgi clean
	make -C prop clean

deep-clean: clean
	make -C lib deep-clean
	make -C cpropep deep-clean
	make -C cgi deep-clean
	make -C prop deep-clean
