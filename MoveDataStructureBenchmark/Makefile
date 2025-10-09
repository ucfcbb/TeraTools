all: bench

.PHONY: all clean test

developmentflags+= -Wall -Wpedantic -Wextra -pedantic-errors -Wold-style-cast
cppflags+=$(developmentflags) -O3 -DNDEBUG
#cppflags+=$(developmentflags) -Wall -g

optIncDir = ../bitpackedmovi/include
incdir = ../bitpackedmovi/thirdparty/include
libdir = ../bitpackedmovi/thirdparty/lib

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

bench: $(benchHeads) bench.cpp
	$(CXX) $(cppflags) $(incs) -L$(libdir) $^ -o $@ $(libs)
