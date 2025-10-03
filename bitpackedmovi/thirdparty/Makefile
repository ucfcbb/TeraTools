all: lib/ include/

.PHONY: all clean

clean: 
	-rm -r sdsl-lite lib include

sdsl-lite/:
	git clone https://github.com/simongog/sdsl-lite

lib/ include/: sdsl-lite/
	cd sdsl-lite && ./install.sh ..
