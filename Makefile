# The simplest makefile you've ever seen...
all:
	make -j4 -C src
	make -j4 -C app

clean:
	make clean -C src
	make clean -C app
