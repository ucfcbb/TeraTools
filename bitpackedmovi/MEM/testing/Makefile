SHELL := /bin/bash

all: test

syngenparam = "10 10 10 10"
modes = SM LM
LEN = 2
numIter=10

ropebwt3 = ../../../ropebwt3/
lcpComputer = ../../lcpComputer/
msIndexBuilder = ../../msIndexBuilder/
MEM = ../
brute_force_repeats = ../../../less_efficient_but_accurate_algos/Repeats/
syntheticgen = ../../../data/generation/

allprogs = ropebwt3 lcpComputer msIndexBuilder MEM brute_force_repeats syntheticgen

.PHONY: all test ropebwt3 lcpComputer msIndexBuilder MEM brute_force_repeats syntheticgen

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
	done
