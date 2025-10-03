all: kactl_impl

.PHONY: all clean 

clean: 
	-rm kactl_impl

kactl_impl: algo.cpp
	$(CXX) $^ -o $@
