# Frontend to dune.
.NOTPARALLEL:

TARGETS=quack_input \
        quack_integrals \
        GoDuck

.PHONY: default build install uninstall test clean

default: $(TARGETS)
	mv GoDuck ../

%.exe: build
	dune build $@

$(TARGETS): $(patsubst %, %.exe, $(TARGETS))
	for i in $(TARGETS) ; do cp _build/default/$${i}.exe $$i ; done

build:
	dune build

clean:
	dune clean


