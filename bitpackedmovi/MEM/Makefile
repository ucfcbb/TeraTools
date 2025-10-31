all: MEM

.PHONY: all clean test msIndexBuilderTest

testinds = ../../data/testing/idx.fmd ../../data/testing/mtb152.fmd

LEN = 2000

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

MEMHeads = $(optIncDir)/util/util.h\
           $(optIncDir)/matchingStatistics/msIndex.h\
           $(optIncDir)/moveStructure/moveStructure.h\

libs = -lsdsl -lpthread -lz -lm
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(ropebwt3dir)

cppflags+=-Wall -O3 -DNDEBUG
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
