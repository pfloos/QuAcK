ALL_DIRS=$(filter-out modules,$(patsubst %/,%,$(wildcard */)))

# Rules for Modules
###################

MOD_DIRS=numgrid


# Rules for Libraries
#####################

LDIR =../lib
LIB_DIRS=$(filter-out $(MAIN_DIRS), $(ALL_DIRS))
MAKEFILES=$(patsubst %,%/Makefile, $(ALL_DIRS))

FORCE:

%/Makefile:
	@(echo TARGET=$*.a > $*/Makefile && echo include ../Makefile.include >> $@)

$(LDIR)/%.a: FORCE %/Makefile $(MOD_DIRS)
	$(MAKE) -C $* lib



# Rules for executables
#######################

MAIN_DIRS=QuAcK eDFT

BDIR=../bin
ALL_EXEC=$(patsubst %, $(BDIR)/%, $(MAIN_DIRS) )

$(BDIR)/%: $(patsubst %,$(LDIR)/%.a,$(LIB_DIRS)) %/Makefile
	$(MAKE) -C $* ../$@




# Rules for both
################

default: $(ALL_EXEC) 

.DEFAULT_GOAL := default
.PHONY: default
.PRECIOUS: $(LDIR)/%.a %/Makefile 

clean: 
	rm -f -- $(patsubst %/,$(LDIR)/%.a,$(wildcard */)) ; \
	rm -f -- $(ALL_EXEC) ; \
	for i in $(MAIN_DIRS) $(LIB_DIRS) ; do rm -f -- $$i/obj/* ; done

