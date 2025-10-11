all: syntheticgen

.PHONY: all clean

clean:
	-rm syntheticgen

syntheticgen: syntheticgen.cpp
	$(CXX) $^ -o $@
