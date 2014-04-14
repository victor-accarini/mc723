all: mips1

mips1:
	acsim mips1.ac -abi
	make -f Makefile.archc

run:
	./mips1.x --load=$(P)
