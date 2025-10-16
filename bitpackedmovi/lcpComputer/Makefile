SHELL := /bin/bash

all: lcpComputer

.PHONY: all clean test

numIterRand=10

testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

optIncDir = ../include
incdir = ../thirdparty/include
libdir = ../thirdparty/lib
ropebwt3dir = ../../ropebwt3

objs = $(ropebwt3dir)/bre.o\
	   $(ropebwt3dir)/dawg.o\
	   $(ropebwt3dir)/fm-index.o\
	   $(ropebwt3dir)/io.o\
	   $(ropebwt3dir)/kalloc.o\
	   $(ropebwt3dir)/kthread.o\
	   $(ropebwt3dir)/libsais16.o\
	   $(ropebwt3dir)/libsais16x64.o\
	   $(ropebwt3dir)/misc.o\
	   $(ropebwt3dir)/mrope.o\
	   $(ropebwt3dir)/rld0.o\
	   $(ropebwt3dir)/rle.o\
	   $(ropebwt3dir)/rope.o\
	   $(ropebwt3dir)/ssa.o

lcpHeads = $(optIncDir)/optbwtrl/util.h\
           $(optIncDir)/optbwtrl/lcpComputer.h\
           $(optIncDir)/optbwtrl/moveStructure.h\

libs = -lsdsl -lpthread -lz -lm -latomic
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(ropebwt3dir)

developmentflags+= -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+=$(developmentflags) -O3 -DNDEBUG
#cppflags+=$(developmentflags) -Wall -g

ifneq ($(asan),)
	cppflags+=-fsanitize=address
	libs+=-fsanitize=address -ldl
endif

ifneq ($(omp),0)
	cppflags+=-DLIBSAIS_OPENMP
	cppflags+=-fopenmp
	libs+=-fopenmp
endif

clean:
	-rm lcpComputer

test: lcpComputer
	set -e && \
	for testind in $(testinds); do \
		./lcpComputer $$testind $$testind; \
	done

lcpComputer: $(lcpHeads) $(objs) lcpComputer.cpp
	$(CXX) $(cppflags) $(incs) -L$(libdir) $^ -o $@ $(libs)
