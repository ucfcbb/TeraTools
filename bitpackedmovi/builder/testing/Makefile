all: test

.PHONY: all clean builder kactl-sa-implmentation generation

builder: 
	make -C ../

generation:
	make -C ../../../data/generation/

kactl-sa-implementation:
	make -C ../../../kactl-sa-implementation

test: builder kactl-sa-implementation generation
	./tester.sh ../builder ../../../kactl-sa-implementation/ ../../../data/generation
