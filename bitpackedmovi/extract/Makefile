all: extract extractMultiple

.PHONY: all clean test

testinds = "../../data/testing/idx.len ../../data/testing/idx.fmd.optbwtrl" "../../data/testing/mtb152.fmd.len ../../data/testing/mtb152.fmd.optbwtrl"

sdsl_optbwt_incdir = ../include
sdsl_optbwt_libdir = ../lib
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

libs = -lsdsl -lpthread -lz -lm
incs = -I$(sdsl_optbwt_incdir) -I$(ropebwt3dir)

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
	-rm extract extractMultiple

test: extract
	set -e && \
	for testind in $(testinds); do \
		./extract 0 $$testind; \
	done

extract: ../include/optbwtrl/optbwtrl.h ../include/optbwtrl/util.h $(objs) extract.cpp
	$(CXX) $(cppflags) $(incs) -L$(sdsl_optbwt_libdir) $^ -o $@ $(libs)

extractMultiple: ../include/optbwtrl/optbwtrl.h ../include/optbwtrl/util.h $(objs) extractMultiple.cpp
	$(CXX) $(cppflags) $(incs) -L$(sdsl_optbwt_libdir) $^ -o $@ $(libs)
