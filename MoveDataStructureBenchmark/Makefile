all: bench

.PHONY: all clean test

developmentflags+= -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+= -std=c++20
cppflags+=$(developmentflags) -O3 -DNDEBUG
#cppflags+=$(developmentflags) -Wall -g

optIncDir = ../bitpackedmovi/include
thDir = ../bitpackedmovi/thirdparty
incdir = $(thDir)/include
libdir = $(thDir)/lib
move-rdir = Move-r/include
move-rincs = malloc_count concurrentqueue gtl/include ips4o/include
move-rlibs = 

#ltbb and latomic for Move-r and its dependences (ips40)
libs = -lsdsl -lpthread -lz -lm -ltbb -latomic
incs = -I$(optIncDir) -isystem$(incdir) -isystem$(move-rdir) $(foreach d, $(move-rincs), -isystemMove-r/external/$(d))

ifneq ($(asan),)
	cppflags+=-fsanitize=address
	libs+=-fsanitize=address -ldl
endif

#openmp required for Move-r (through ips40::parallel::sort)
ifneq ($(omp),0)
	cppflags+=-DLIBSAIS_OPENMP
	cppflags+=-fopenmp
	libs+=-fopenmp
endif

benchHeads = $(optIncDir)/optbwtrl/util.h\
             $(optIncDir)/optbwtrl/moveStructure.h
clean cleanours:
	-rm bench

bench: $(incdir) $(libdir) $(benchHeads) $(move-rdir)

bench: bench.cpp
	$(CXX) $(cppflags) $(incs) -L$(libdir) bench.cpp -o $@ $(libs)

$(move-rdir): Move-r
	cd Move-r && mkdir build && cd build && cmake .. && cp -r ../patched-files/* .. && make

Move-r:
	git clone --recurse-submodules https://github.com/asanaullah2015/Move-r.git


$(incdir) $(libdir) &:
	$(MAKE) -C $(thDir)
