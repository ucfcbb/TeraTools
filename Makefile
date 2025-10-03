all: bitpackedmovi less_efficient_but_accurate_algos

.PHONY: all clean bitpackedmovi less_efficient_but_accurate_algos

clean:
	-rm -r ropebwt3
	make -C bitpackedmovi/ clean
	make -C less_efficient_but_accurate_algos/ clean

bitpackedmovi: ropebwt3
	make -C bitpackedmovi/

ropebwt3:
	git clone https://github.com/lh3/ropebwt3
	make -C ropebwt3/

less_efficient_but_accurate_algos:
	make -C less_efficient_but_accurate_algos/
