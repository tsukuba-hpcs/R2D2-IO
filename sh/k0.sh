#!/bin/sh -x
#PJM --name "s002"
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=2048"
#PJM --rsc-list node-quota=29G
#PJM --rsc-list "elapse=00:30:00"
#PJM --mpi "assign-online-node"
#PJM --mpi "proc=2048"
#PJM -S
#PJM --spath ana.txt
#PJM --step
#
### STAGING DISCRIPTION
#
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
### STAGE IN ###
#PJM --stgin "rank=* ./a.out %r:./"
#PJM --stgin "rank=* ./input_data/value_cart.dac %r:./input_data/value_cart.dac"
#PJM --stgin "rank=* ./input_data/eos_table_sero.dac %r:./input_data/eos_table_sero.dac"
#PJM --stgin "rank=* ./input_data/eos_table_enro.dac %r:./input_data/eos_table_enro.dac"
#PJM --stgin "rank=* ./input_data/params.txt %r:./input_data/params.txt"
#
### STAGE OUT ###
#PJM --stgout "rank=* %r:./data/* ./data/"
#PJM --stgout "rank=* %r:./data/time/*  ./data/time/"
#PJM --stgout "rank=* %r:./data/remap/* ./data/remap/"
#PJM --stgout "rank=* %r:./data/param/* ./data/param/"
#PJM --stgout "rank=* %r:./log.txt.%r ./data/log/log.txt.%r"
##PJM --stgout "rank=* %r:./data/qq/qq.dac.e.%08r           ./data/qq/qq.dac.e.%08r"
##PJM --stgout "rank=* %r:./data/qq/qq.dac.o.%08r           ./data/qq/qq.dac.o.%08r"
#PJM --stgout "rank=* %r:./data/qq/*           ./data/qq/"
##PJM --stgout "rank=* %r:./out/* ./out/"
#PJM --stgout "rank=* %r:./out/* ./job%j/out/"
#PJM --stgout "rank=0 %r:../* ./data/remap/"
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#
#mpiexec -of-proc log.txt ./a.out
fapp -C -d out  -Ihwm -L1 mpiexec -of-proc log.txt ./a.out
