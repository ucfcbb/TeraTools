all: thirdparty/ builder/ loader/

.PHONY: all clean thirdparty/ builder/ loader/

clean:
	make -C thirdparty/ clean
	make -C builder/ clean
	make -C loader/ clean

thirdparty/:
	make -C thirdparty/

builder/: thirdparty/
	make -C builder/

loader/: thirdparty/
	make -C loader/
