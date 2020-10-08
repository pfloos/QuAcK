# Frontend to dune.
.NOTPARALLEL:

TARGETS=quack_input \
        quack_integrals

.PHONY: default build install uninstall test clean                                            

%.exe: build
	dune build $@

$(TARGETS): $(patsubst %, %.exe, $(TARGETS))
	for i in $* ; do ln _build/default/$$i $${i%.exe} ; done

default: $(TARGETS)

build:
	dune build

clean:
	dune clean


