SHELL := /bin/bash

all: test

syngenparam = "10 10 10 10"

ropebwt3 = ../../../ropebwt3/
lcpComputer = ../
minLCP = ../../../less_efficient_but_accurate_algos/kactl-sa-implementation/
syntheticgen = ../../../data/generation/

numIter=10

.PHONY: all test ropebwt3 builder kactl_impl lcpComputer minLCP syntheticgen loader

ropebwt3 builder kactl_impl lcpComputer minLCP syntheticgen loader:
	make -C $($@) $@

test: ropebwt3 lcpComputer minLCP syntheticgen
	test=0;\
	while ./tester.sh $(foreach prog, $^, $($(prog))$(prog)) $(syngenparam) ; do \
		echo passed test $$((++test)); \
		if [[ $$test == $(numIter) ]]; then \
			echo passed all $(numIter) tests; \
			break; \
		fi \
	done
