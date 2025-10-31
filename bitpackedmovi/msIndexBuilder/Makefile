SHELL := /bin/bash

all: msIndexBuilder

.PHONY: all clean test lcpComputerTest

numIterRand=10

testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

optIncDir = ../include
incdir = ../thirdparty/include
libdir = ../thirdparty/lib
ropebwt3dir = ../../ropebwt3

ropebwt3objs = $(shell grep -oPz "OBJS=.*\\\\\n.*\n" $(ropebwt3dir)/Makefile | tail -c +6 | tr -d "\\\\")

objs = $(foreach obj, $(ropebwt3objs), $(ropebwt3dir)/$(obj))

msHeads = $(optIncDir)/util/util.h\
		  $(optIncDir)/lcpComputer/lcpComputer.h\
		  $(optIncDir)/matchingStatistics/msIndex.h\
		  $(optIncDir)/moveStructure/moveStructure.h

libs = -lsdsl -lpthread -lz -lm -latomic
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(ropebwt3dir)

developmentflags+= -std=c++17 -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
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
	-rm msIndexBuilder

lcpComputerTest:
	make -C ../lcpComputer/ testIndexOnly

test: msIndexBuilder lcpComputerTest
	set -e && \
	for testind in $(testinds); do \
		./msIndexBuilder -i ../lcpComputer/$$testind.lcp_index -o $$testind; \
	done

msIndexBuilder: msIndexBuilder.cpp $(msHeads) $(objs)
	$(CXX) $(cppflags) $(incs) -L$(libdir) $(objs) $< -o $@ $(libs)
