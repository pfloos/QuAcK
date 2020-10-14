#!/bin/bash
set -e

: ${LIBXC_VERSION:=5.0.0}
: ${LXC:=libxc-${LIBXC_VERSION}}
: ${TARGZ:=${LXC}.tar.gz}

if [[ $1 == "install" ]] ; then
  make -C $LXC install
  exit 0
fi

set -x
[[ -f $TARGZ ]] || wget https://gitlab.com/libxc/libxc/-/archive/${LIBXC_VERSION}/${TARGZ}
tar -xzf $TARGZ
cd $LXC
autoreconf -i
./configure --prefix=$QUACK_ROOT
make
#make check
make install
set +x
