SHELL := /bin/bash

all: TeraLCP

.PHONY: all clean test testIndexOnly

numIterRand=10

testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

optIncDir = ../include
incdir = ../thirdparty/include
libdir = ../thirdparty/lib
ropebwt3dir = ../../ropebwt3

ropebwt3objs = $(shell grep -oPz "OBJS=.*\\\\\n.*\n" $(ropebwt3dir)/Makefile | tail -c +6 | tr -d "\\\\")

objs = $(foreach obj, $(ropebwt3objs), $(ropebwt3dir)/$(obj))

lcpHeads = $(optIncDir)/util/util.h\
           $(optIncDir)/TeraLCP/TeraLCP.h\
           $(optIncDir)/moveStructure/moveStructure.h\

libs = -lsdsl -lpthread -lz -lm -latomic
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(ropebwt3dir)

developmentflags+= -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+=$(developmentflags) -O3 -DNDEBUG -std=c++11
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
	-rm TeraLCP

testIndexOnly: TeraLCP
	set -e && \
	for testind in $(testinds); do \
		./TeraLCP -f fmd -i $$testind -oindex $$testind -t testTEMPFILE; \
        rm testTEMPFILE; \
	done

test: TeraLCP
	set -e && \
	for testind in $(testinds); do \
		./TeraLCP -f fmd -i $$testind -orlcp test -t testTEMPFILE; \
        rm testTEMPFILE; \
	done

TeraLCP: TeraLCP.cpp $(lcpHeads) $(objs)
	$(CXX) $(cppflags) $(incs) -L$(libdir) $(objs) $< -o $@ $(libs)
