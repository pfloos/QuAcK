.SECONDARY:

#QUACK_ROOT=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))/..))
#export QUACK_ROOT

default:
	./utils/gitversion.sh
	mv gitver.f90 ./utils/
	python3 make_ninja.py
	ninja

debug:
	./utils/gitversion.sh
	mv gitver.f90 ./utils/
	python3 make_ninja.py -d
	ninja

clean: 
	python3 make_ninja.py
	for i in */ ; do cd $$i ; ninja -t clean ; cd .. ; done
