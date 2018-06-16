include ../../../makefile.setup
include ../../../makefile.sys

.PHONY: all clean run run_paradisimp

all: Home_imp.py

Home_imp.py: ../../../implicit/Include/Home.h
	@echo create Home_imp.py based on Home.h
	python2 ../../../python/ctypesgen/ctypesgen.py --cpp="gcc -E -D_IMPLICIT -D_RETROCOLLISIONS -I ../../../include" $^ -o $@

pyshdis: pyshdis.py Home_imp.py frank_read_src_imp.ctrl frank_read_src_imp.data
	python pyshdis.py frank_read_src_imp
clean: 
	rm -f *.pyc  Home_imp.py

