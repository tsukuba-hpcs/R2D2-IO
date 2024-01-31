#PBS -N d003
#PBS -l nodes=2
#PBS -l walltime=2:00:00
##PBS -q large-a
#PBS -q bulk-a
cd ${PBS_O_WORKDIR}
export ATP_ENABLED=1

#export OMP_NUM_THREADS=1
#aprun -n 80 -N 40 -d ${OMP_NUM_THREADS} -cc depth ./a.out > log.txt
export OMP_NUM_THREADS=5
aprun -n 16 -N 8 -d ${OMP_NUM_THREADS} -cc depth ./a.out > log.txt

## -n total number of MPI thread
## -N number of MPI thread per node
