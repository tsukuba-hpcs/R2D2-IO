#!/bin/bash
#PJM -L "node=1"
#PJM -L "rscunit=rscunit_ft01"
#PJM -L "rscgrp=eap-int"
#PJM -L "elapse=00:10:00"
#PJM -o log.txt
#PJM --mpi "max-proc-per-node=4"
#PJM -s

export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}

# execute job
#mpiexec -n 4 ./a.out

rm -rf prof rep* *csv
rm -rf data
fapp -C -d ./prof -Icpupa mpiexec ./a.out
