all: builder/

.PHONY: all clean builder/

clean:
	-rm -r sdsl-lite lib include
	make -C builder/ clean

sdsl-lite/:
	git clone https://github.com/simongog/sdsl-lite

lib/ include/ &: sdsl-lite/
	cd sdsl-lite && ./install.sh ..

builder/: lib/ include/
	make -C builder/
