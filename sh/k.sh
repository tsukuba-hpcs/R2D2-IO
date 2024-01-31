#!/bin/sh -x
#PJM --name "s156"
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=512"
#PJM --rsc-list node-quota=29G
#PJM --rsc-list "elapse=05:00:00"
#PJM --mpi "assign-online-node"
#PJM --mpi "proc=512"
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
#PJM --stgin "rank=0 data/time/t.dac.e      0:data/time/t.dac.e"
#PJM --stgin "rank=0 data/time/t.dac.o      0:data/time/t.dac.o"
#PJM --stgin "rank=0 data/param/nd.dac           0:data/param/nd.dac"
#PJM --stgin "rank=0 data/param/params.dac       0:data/param/params.dac"
#PJM --stgin "rank=* data/qq/qq.dac.e.%08r           %r:data/qq/qq.dac.e.%08r"
#PJM --stgin "rank=* data/qq/qq.dac.o.%08r           %r:data/qq/qq.dac.o.%08r"
#
### STAGE OUT ###
#PJM --stgout "rank=* %r:./data/* ./data/"
#PJM --stgout "rank=* %r:./data/time/*  ./data/time/"
#PJM --stgout "rank=* %r:./data/remap/* ./data/remap/"
#PJM --stgout "rank=0 0:./data/param/nd.dac ./data/param/nd.dac"
#PJM --stgout "rank=0 0:./data/param/params.dac ./data/param/params.dac"
#PJM --stgout "rank=* %r:./log.txt.%r ./data/log/log.txt.%r"
#PJM --stgout "rank=* %r:./data/qq/* ./data/qq/"
#PJM --stgout "rank=* %r:./out/* ./job%j/out/"
#PJM --stgout "rank=0 %r:../* ./data/remap/"
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
#
mpiexec /work/system/bin/msh "mkdir ./data"
mpiexec /work/system/bin/msh "mkdir ./data/qq"
mpiexec /work/system/bin/msh "mkdir ./data/time"
mpiexec /work/system/bin/msh "mkdir ./data/param"
mpiexec /work/system/bin/msh "mkdir ./data/remap"
mpiexec /work/system/bin/msh "touch ./data/param/nd.dac"
#fapp -C -d out  -Ihwm -L1 mpiexec -of-proc log.txt ./a.out
mpiexec -of-proc log.txt ./a.out