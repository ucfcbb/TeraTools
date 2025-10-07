all: brute_force_repeats

.PHONY: all clean

clean:
	-rm brute_force_repeats

brute_force_repeats: brute_force_repeats.cpp
	$(CXX) -std=c++17 $^ -o $@
