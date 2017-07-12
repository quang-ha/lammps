#! /bin/bash

dname=data
lmp=../../../../src/lmp_mpi
${lmp} -in infslab.lmp -var dname ${dname}

