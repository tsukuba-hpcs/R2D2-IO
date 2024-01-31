#!/bin/bash
#------ pjsub option -------#
#PJM -N "d018"
#PJM -L "rscgrp=fx-debug"
##PJM -L "rscgrp=fx-large"
##PJM -L "rscgrp=fx-extra"
##PJM -L "rscgrp=fx-middle"
##PJM -L "rscgrp=fx-small"
#PJM -L "node=1x2x2"
#PJM --mpi "proc=16"
#PJM -L "elapse=01:00:00"
#PJM -j
#PJM -o "log.txt"
#PJM --step
#------ Program execution -------#
#export PARALLEL=12
export OMP_NUM_THREADS=12
#rm -rf prof
#rm -rf data
#mpiexec ./a.out
#fapp -C -d prof -I mpi -I hwm mpiexec ./a.out

mpiexec ./a.out

#rm -rf prof rep* *csv
#rm -rf data
#echo prof start
#fapp -C -d ./prof -Icpupa mpiexec ./a.out

#for i in `seq 1 11`
#do
#    rm -rf data
#    echo $i start
#    fapp -C -d ./rep$i -Hevent=pa$i mpiexec ./a.out
#done

