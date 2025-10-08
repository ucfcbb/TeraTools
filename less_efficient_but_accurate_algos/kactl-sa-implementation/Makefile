all: kactl_impl minLCP

.PHONY: all clean 

clean: 
	-rm kactl_impl minLCP

kactl_impl: algo.cpp
	$(CXX) $^ -o $@

minLCP: minLCP.cpp
	$(CXX) $^ -o $@
