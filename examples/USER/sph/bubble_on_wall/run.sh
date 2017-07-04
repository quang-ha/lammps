#! /bin/bash

set -e
set -u

lmp=../../../../src/lmp_mpi
mpirun=mpirun
dname=data

${lmp} -in bubble.lmp -var dname ${dname} 
