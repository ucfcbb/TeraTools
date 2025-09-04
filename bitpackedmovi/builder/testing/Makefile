all: test

.PHONY: all clean ropebwt3 builder kactl-sa-implmentation generation

numIter = $(shell seq 1 10)

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
	for i in $(numIter); do \
		./tester.sh ../../../ropebwt3/ropebwt3 ../builder ../../../kactl-sa-implementation/kactl_impl ../../../data/generation/syntheticgen ../../loader/loader ; \
	done
