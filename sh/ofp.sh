#!/bin/sh -l

#-------- pjsub option -------#
#PJM -N "d001"
#PJM -L rscgrp=regular-flat
##PJM -L rscgrp=debug-flat
#PJM -L node=16
#PJM --mpi proc=128
#PJM --omp thread=16
#PJM -L elapse=06:00:00
##PJM -L elapse=12:00:00
##PJM -L elapse=12:00:00
#PJM -g hp120286
#PJM -j
#PJM --restart
#PJM --step

module load fftw

export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST="0,1,68,69,136,137,204,205"
#export I_MPI_PIN_DOMAIN=256
#export I_MPI_PERHOST=1

# 2 thread/mpi is expected

## omp thread = 1
#export I_MPI_PIN_DOMAIN=2
#export I_MPI_PERHOST=128

## omp thread = 2
#export I_MPI_PIN_DOMAIN=4
#export I_MPI_PERHOST=64

## omp thread = 4
#export I_MPI_PIN_DOMAIN=8
#export I_MPI_PERHOST=32

## omp thread=8
#export I_MPI_PIN_DOMAIN=16
#export I_MPI_PERHOST=16

## omp thread 16
export I_MPI_PIN_DOMAIN=32
export I_MPI_PERHOST=8

#export I_MPI_DEBUG=5
export KMP_HW_SUBSET=2T
export KMP_AFFINITY=verbose

#-------- Program execution -------#
#mpiexec.hydra -n ${PJM_MPI_PROC} ./a.out > log.txt
mpiexec.hydra -n ${PJM_MPI_PROC} numactl --preferred=1 ./a.out > log.txt
#mpiexec.hydra -n ${PJM_MPI_PROC} numactl --membind=1 ./a.out > log.txt

#mpiexec.hydra -n ${PJM_MPI_PROC} numactl --membind=1 ./a.out > log0.txt
#module load advisor
#mpiexec.hydra -n ${PJM_MPI_PROC} advixe-cl -collect survey -project-dir result ./a.out > log.txt
#mpiexec.hydra -n ${PJM_MPI_PROC} advixe-cl -collect tripcounts -flop -project-dir result -no-auto-finalize ./a.out > log.txt
