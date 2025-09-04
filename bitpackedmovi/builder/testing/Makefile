all: test

.PHONY: all clean ropebwt3 builder kactl-sa-implmentation generation

builder: 
	make -C ../

loader:
	make -C ../../loader

ropebwt3: 
	make -C ../../../ropebwt3/

generation:
	make -C ../../../data/generation/

kactl-sa-implementation:
	make -C ../../../kactl-sa-implementation

test: ropebwt3 builder kactl-sa-implementation generation loader
	./tester.sh ../../../ropebwt3/ropebwt3 ../builder ../../../kactl-sa-implementation/kactl_impl ../../../data/generation/syntheticgen ../../loader/loader
