SHELL := /bin/bash

all: test

.PHONY: all clean ropebwt3 builder kactl-sa-implmentation generation

numIter=10

builder: 
	make -C ../

loader:
	make -C ../../loader

ropebwt3: 
	make -C ../../../ropebwt3/

generation:
	make -C ../../../data/generation/

kactl-sa-implementation:
	make -C ../../../less_efficient_but_accurate_algos/kactl-sa-implementation

test: ropebwt3 builder kactl-sa-implementation generation loader
	test=0;\
	while ./tester.sh ../../../ropebwt3/ropebwt3 ../builder ../../../less_efficient_but_accurate_algos/kactl-sa-implementation/kactl_impl ../../../data/generation/syntheticgen ../../loader/loader; do \
		echo passed test $$((++test)); \
		if [[ $$test == $(numIter) ]]; then \
			echo passed all $(numIter) tests; \
			break; \
		fi \
	done

