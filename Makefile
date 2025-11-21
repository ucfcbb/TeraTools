subdirs = src data less_efficient_but_accurate_algos #MoveDataStructureBenchmark

all: $(subdirs)

.PHONY: all clean $(subdirs)

clean:
	-rm -r ropebwt3
	for dir in $(subdirs); do \
        $(MAKE) -C $$dir clean; \
    done

cleanours:
	for dir in $(subdirs); do \
		$(MAKE) -C $$dir cleanours; \
	done

src: ropebwt3

ropebwt3:
	git clone https://github.com/lh3/ropebwt3
	$(MAKE) -C ropebwt3/

$(subdirs): 
	$(MAKE) -C $@
