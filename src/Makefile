QUACK_ROOT=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))/..))
export QUACK_ROOT

SDIR=$(QUACK_ROOT)/src
ALL_DIRS=$(filter-out modules,$(patsubst %/,%,$(wildcard */)))

# Rules for Modules
###################

MOD_DIRS=numgrid


# Rules for Libraries
#####################

LDIR =$(QUACK_ROOT)/lib
LIB_DIRS=$(filter-out $(MAIN_DIRS), $(ALL_DIRS))

FORCE:

$(SDIR)/%/Makefile: FORCE
	@(echo QUACK_ROOT=$(QUACK_ROOT) > $@)
	@(echo TARGET=$*.a >> $@)
	@(echo 'include $$(QUACK_ROOT)/src/Makefile.include' >> $@)

$(LDIR)/%.a: FORCE $(SDIR)/%/Makefile $(MOD_DIRS)
	$(MAKE) -C $(SDIR)/$* lib



# Rules for executables
#######################

MAIN_DIRS=QuAcK eDFT

BDIR=$(QUACK_ROOT)/bin
ALL_EXEC=$(patsubst %, $(BDIR)/%, $(MAIN_DIRS) )

$(BDIR)/%: $(patsubst %,$(LDIR)/%.a,$(LIB_DIRS)) $(SDIR)/%/Makefile
	$(MAKE) -C $* $@




# Rules for both
################

default: $(ALL_EXEC) $(wildcard $(SDIR)/*/Makefile) 

debug: 
	DEBUG=1 $(MAKE) default

.DEFAULT_GOAL := default
.PHONY: default debug
.PRECIOUS: $(LDIR)/%.a %/Makefile 

clean: 
	rm -f -- $(patsubst %/,$(LDIR)/%.a,$(wildcard */)) ; \
	rm -f -- $(ALL_EXEC) ; \
	for i in $(MAIN_DIRS) $(LIB_DIRS) ; do rm -f -- $$i/obj/* ; done

