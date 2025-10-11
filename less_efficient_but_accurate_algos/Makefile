subdirs = kactl-sa-implementation Repeats
all: kactl-sa-implementation Repeats

.PHONY: all clean $(subdirs)

clean cleanours:
	for dir in $(subdirs); do \
        $(MAKE) -C $$dir clean; \
    done

$(subdirs): 
	$(MAKE) -C $@
