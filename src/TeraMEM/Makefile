all: TeraMEM

.PHONY: all clean test TeraIndexTest

testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

LEN = 2000

optIncDir = ../include
incdir = ../thirdparty/include
libdir = ../thirdparty/lib
ropebwt3dir = ../../ropebwt3

ropebwt3objs = $(shell if [ ! -d $(ropebwt3dir) ]; then git clone https://github.com/lh3/ropebwt3 $(ropebwt3dir); fi &&\
	grep -oPz "OBJS=.*\\\\\n.*\n" $(ropebwt3dir)/Makefile | tail -c +6 | tr -d "\\\\")

objs = $(foreach obj, $(ropebwt3objs), $(ropebwt3dir)/$(obj))

TeraMEMHeads = $(optIncDir)/util/util.h\
           $(optIncDir)/TeraIndex/TeraIndex.h\
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
	-rm TeraMEM

TeraIndexTest:
	make -C ../TeraIndex/ test

test: TeraMEM TeraIndexTest
	set -e && \
	for testind in $(testinds); do \
		./TeraMEM -MEM SM -i $$testind.ms_index -o $$testind.SM.out; \
		./TeraMEM -MEM LM -i $$testind.ms_index -o $$testind.LM.out -L $(LEN); \
	done

$(objs): 
	make -C $(ropebwt3dir)

$(libdir) $(incdir):
	make -C ../thirdparty/

TeraMEM: TeraMEM.cpp $(TeraMEMHeads) $(objs) $(libdir) $(incdir)
	$(CXX) $(cppflags) $(incs) -L$(libdir) $(objs) $< -o $@ $(libs)
