all: init

init:
	git submodule update --init
	cd sonLib && make
	cd jobTree && make
	cd hal && make
