
all:
	make -C src all
	
clean:
	make -C src clean

deep-clean: clean
	make -C src deep-clean
