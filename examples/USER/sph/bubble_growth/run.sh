#! /bin/bash

set -e
set -u

ndim=3
nx=20
lmp=../../../../src/lmp_mpi

dname=data-ndim${ndim}-nx${nx}
mpirun -np 8 ${lmp} -var nx ${nx} -var ndim ${ndim} -var dname ${dname} -in bubble.lmp
    
    
