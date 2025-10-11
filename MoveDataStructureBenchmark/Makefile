all: bench

.PHONY: all clean test

developmentflags+= -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+=$(developmentflags) -O3 -DNDEBUG
#cppflags+=$(developmentflags) -Wall -g

optIncDir = ../bitpackedmovi/include
thDir = ../bitpackedmovi/thirdparty
incdir = $(thDir)/include
libdir = $(thDir)/lib

libs = -lsdsl -lpthread -lz -lm
incs = -I$(optIncDir) -isystem$(incdir)

ifneq ($(asan),)
	cppflags+=-fsanitize=address
	libs+=-fsanitize=address -ldl
endif

benchHeads = $(optIncDir)/optbwtrl/util.h\
             $(optIncDir)/optbwtrl/moveStructure.h
clean:
	-rm bench

bench: $(incdir) $(libdir) $(benchHeads) 

bench: bench.cpp
	$(CXX) $(cppflags) $(incs) -L$(libdir) bench.cpp -o $@ $(libs)

$(incdir) $(libdir) &:
	$(MAKE) -C $(thDir)
