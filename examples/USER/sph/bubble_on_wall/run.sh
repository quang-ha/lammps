#! /bin/bash

set -e
set -u

lmp=../../../../src/lmp_mpi
mpirun=mpirun
dname=data

mpirun -np 8 ${lmp} -in bubble.lmp -var dname ${dname} 
