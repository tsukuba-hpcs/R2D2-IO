#!/bin/bash
#------ pjsub option -------#
#PJM -N "d093"
##PJM -L "rscgrp=fx-debug"
##PJM -L "rscgrp=fx-large"
#PJM -L "rscgrp=fx-small"
#PJM -L "node=1"
#PJM --mpi "proc=8"
#PJM -L "elapse=00:15:00"
##PJM -L "elapse=010:00:00"
#PJM -j
#PJM -o "log.txt"
#PJM --step
#------ Program execution -------#
export PARALLEL=8
export OMP_NUM_THREADS=8
#rm -rf prof
#rm -rf data
mpiexec ./a.out
#fapp -C -d prof -I mpi -I hwm mpiexec ./a.out
