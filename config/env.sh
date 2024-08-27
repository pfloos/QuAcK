#!/bin/env

CURRENT_HOSTNAME=$(hostname)
echo "Current machine: $CURRENT_HOSTNAME"

export QUACK_ROOT=/users/p18005/ammar/tmpdir/QuAcK

case $CURRENT_HOSTNAME in
        *turpan*)
        ;;
        *olympe*)
		module purge
		module load python/3.9.5
		module load intel/18.2.199
		export PATH=$PATH:/users/p18005/ammar/qp2/bin
		export LD_LIBRARY_PATH=/usr/local/miniconda/4.9.2/envs/python-3.9.5/lib:$LD_LIBRARY_PATH
		export LD_LIBRARY_PATH=$QUACK_ROOT/lib:$LD_LIBRARY_PATH
		export PYTHONPATH=$QUACK_ROOT/src/pyscf:$PYTHONPATH
        ;;
        *)
                echo "Unknown hostname: $CURRENT_HOSTNAME. No modules loaded."
        ;;
esac

