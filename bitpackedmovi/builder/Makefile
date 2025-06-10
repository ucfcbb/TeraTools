all: builder

.PHONY: all clean

sdslincdir = ../include
sdsllibdir = ../lib

libs = -lsdsl

clean:
	-rm builder

builder: builder.cpp
	$(CXX) -I$(sdslincdir) -L$(sdsllibdir) $^ -o $@ $(libs)
