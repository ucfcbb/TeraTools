SHELL := /bin/bash

all: test

syngenparam = "10 50 20 80"
modes = SM LM
LEN = 10
numIter=10

ropebwt3 = ../../../ropebwt3/
TeraLCP = ../../TeraLCP/
msIndexBuilder = ../../msIndexBuilder/
MEM = ../
brute_force_repeats = ../../../less_efficient_but_accurate_algos/Repeats/
syntheticgen = ../../../data/generation/

allprogs = ropebwt3 TeraLCP msIndexBuilder MEM brute_force_repeats syntheticgen

.PHONY: all test ropebwt3 TeraLCP msIndexBuilder MEM brute_force_repeats syntheticgen

$(allprogs):
	make -C $($@) $@

test: $(allprogs)
	for mode in $(modes); do \
		test=0;\
		while ./tester.sh $(foreach prog, $^, $($(prog))$(prog)) $(syngenparam) $$mode $(LEN); do \
			echo passed $$mode test $$((++test)); \
			if [[ $$test == $(numIter) ]]; then \
				echo passed all $(numIter) $$mode tests; \
				break; \
			fi \
		done; \
        if [ $$test -ne $(numIter) ]; then echo FAILED $$mode test $$((++test)); break; fi \
	done
