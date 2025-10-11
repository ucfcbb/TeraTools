ourdirs = builder extract lcpComputer loader msComputer msIndexBuilder repeats stats 
subdirs = $(ourdirs) thirdparty

all: $(subdirs)

.PHONY: all clean $(subdirs)

clean:
	for dir in $(subdirs); do \
        $(MAKE) -C $$dir clean; \
    done

cleanours:
	for dir in $(ourdirs); do \
        $(MAKE) -C $$dir clean; \
    done

$(ourdirs): thirdparty ../ropebwt3

../ropebwt3:
	$(MAKE) -C .. ropebwt3

$(subdirs):
	$(MAKE) -C $@
