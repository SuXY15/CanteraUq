#!/bin/bash
# .bashrc
source /HOME/pp336/scripts/cn-module.sh

module load Python/3.6.3
export PYTHON_CMD=/WORK/app/Python/3.6.3/bin/python
export PYTHONPATH=$HOME/.local/lib/python3.6/site-packages
export PATH=$PATH:$HOME/.local/bin:PYTHON_CMD
export LD_LIBRARY_PATH=$HOME/.local/lib/sundials-2.7.0:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/.local:$HOME/.local/include/ncurses:$LD_LIBRARY_PATH

yhrun -N 1 -n 24 python3 src/activeSubspace.py DME sampling


# yhbatch -N 1 -p paratera -J DMEsk42 run.sh
