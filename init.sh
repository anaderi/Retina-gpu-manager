#!/bin/bash

export LD_LIBRARY_PATH=:/usr/local/cuda/lib64
export PATH=$PATH:/usr/local/cuda/bin
source /cvmfs/lhcb.cern.ch/lib/LbLogin.sh
export User_release_area=`pwd`
SetupProject Brunel v47r2p1
