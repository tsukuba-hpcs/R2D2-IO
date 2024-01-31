#!/bin/bash
#PJM -N "d001_star"
### Yin-Yangの場合: node = ix0 x jx0 x kx/2
### Yin-Yangでない場合: node = ix0 x jx0/2 x kx0/2
#PJM -L "node=8x8x8"
#PJM -L "rscunit=rscunit_ft01"
##PJM -L "rscgrp=small"
#PJM -L "rscgrp=large"
#PJM -L "elapse=00:30:00"
#PJM -o out.txt
#PJM -e err.txt
#PJM --mpi "max-proc-per-node=4"
#PJM --mpi "proc=2048"
#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM --step
# LLIOの領域確保
#PJM --llio localtmp-size=20Gi
#PJM --llio sharedtmp-size=20Gi

export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}

# a.outを全てのSIOノードにコピー
llio_transfer ./a.out

# execute job
mpiexec -n 2048 -stdout-proc ./output.%j/%/1000R/%m/stdout -stderr-proc ./output.%j/%/1000R/%m/stderr ./a.out
