#!/bin/bash
tar -xzf python37.tar.gz
export PATH=$PWD/python/bin:$PATH
export PYTHONPATH=$PWD
export HOME=$PWD
cd /software/mrblackburn/pymol/
./pymol -c $_CONDOR_SCRATCH_DIR/SASAquatch.py $1
