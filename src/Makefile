.SECONDARY:

QUACK_ROOT=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))/..))
export QUACK_ROOT

default: build.ninja
	python3 make_ninja.py
	ninja

clean: build.ninja
	python3 make_ninja.py
	ninja -t clean
	for i in */ ; do ninja -t clean ; done

