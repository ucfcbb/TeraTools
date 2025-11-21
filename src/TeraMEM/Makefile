all: MEM

.PHONY: all clean test msIndexBuilderTest

testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

LEN = 2000

optIncDir = ../include
incdir = ../thirdparty/include
libdir = ../thirdparty/lib
ropebwt3dir = ../../ropebwt3

ropebwt3objs = $(shell grep -oPz "OBJS=.*\\\\\n.*\n" $(ropebwt3dir)/Makefile | tail -c +6 | tr -d "\\\\")

objs = $(foreach obj, $(ropebwt3objs), $(ropebwt3dir)/$(obj))

MEMHeads = $(optIncDir)/util/util.h\
           $(optIncDir)/matchingStatistics/msIndex.h\
           $(optIncDir)/moveStructure/moveStructure.h\

libs = -lsdsl -lpthread -lz -lm
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(ropebwt3dir)

developmentflags+= -std=c++17 -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+=$(developmentflags) -O3 -DNDEBUG
#cppflags+=-Wall -g

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
	-rm MEM

msIndexBuilderTest:
	make -C ../msIndexBuilder/ test

test: MEM msIndexBuilderTest
	set -e && \
	for testind in $(testinds); do \
		./MEM -MEM SM -i $$testind.ms_index -o $$testind.SM.out; \
		./MEM -MEM LM -i $$testind.ms_index -o $$testind.LM.out -L $(LEN); \
	done

MEM: MEM.cpp $(MEMHeads) $(objs)
	$(CXX) $(cppflags) $(incs) -L$(libdir) $(objs) $< -o $@ $(libs)
