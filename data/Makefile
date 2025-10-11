subdirs = generation
all: $(subdirs)

.PHONY: all clean $(subdirs)

clean cleanours:
	for dir in $(subdirs); do \
        $(MAKE) -C $$dir clean; \
    done

$(subdirs):
	$(MAKE) -C $@
