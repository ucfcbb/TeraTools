all: bitpackedmovi

.PHONY: all clean bitpackedmovi

clean:
	-rm -r ropebwt3
	make -C bitpackedmovi/ clean

bitpackedmovi: ropebwt3
	make -C bitpackedmovi/

ropebwt3:
	git clone https://github.com/lh3/ropebwt3
	make -C ropebwt3/
