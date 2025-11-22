SHELL := /bin/bash 

all: TeraMS

.PHONY: all clean test

#numIterRand=10

#testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

optIncDir = ../include
incdir = ../thirdparty/include
libdir = ../thirdparty/lib
kseqdir = ../thirdparty/kseq
ropebwt3dir = ../../ropebwt3

ropebwt3objs = $(shell grep -oPz "OBJS=.*\\\\\n.*\n" $(ropebwt3dir)/Makefile | tail -c +6 | tr -d "\\\\")

objs = $(foreach obj, $(ropebwt3objs), $(ropebwt3dir)/$(obj))

msHeads = $(optIncDir)/util/util.h\
           $(optIncDir)/TeraLCP/TeraLCP.h\
           $(optIncDir)/TeraIndex/TeraIndex.h\
           $(optIncDir)/moveStructure/moveStructure.h\

libs = -lsdsl -lpthread -lz -lm -latomic
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(ropebwt3dir) -isystem$(kseqdir)

developmentflags+= -std=c++17 -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+=$(developmentflags) -O3 -DNDEBUG
#cppflags+=$(developmentflags) -Wall -g -fsanitize=thread

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
	-rm TeraMS

#test: TeraMS
#	set -e && \
#	for testind in $(testinds); do \
#		./TeraLCP -f fmd -i $$testind -orlcp test.rlcp -t testTEMPFILE; \
#	done

TeraMS: TeraMS.cpp $(msHeads) $(objs)
	$(CXX) $(cppflags) $(incs) -L$(libdir) $(objs) $< -o $@ $(libs)
