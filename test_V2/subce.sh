#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=my_job
#SBATCH --license=SCRATCH
#SBATCH --ntasks-per-node=32
#SBATCH -C haswell

########
#Hybrid MPI/OMP applications
########
## #SBATCH --cpus-per-task=64(total logical cores per node, 32*2)/(MPI tasks per node)
#export OMP_NUM_THREADS=2
#export OMP_PROC_BIND=spread   # new recommendations for hybrid MPI/OpenMP
#export OMP_PLACES=threads

########
#to run with pure MPI
########
#export OMP_NUM_THREADS=1

#######
# When Not All CPUs on a Node are Used
# if #MPI_per_node is not a divisor of 64, -n value times -c value does 
# not equal to 64, meaning not all CPUs on a node are used.
#######
#srun -n 12 -c 10 --cpu_bind=cores check-hybrid.intel.cori  # -c is set to 64/#MPI_per_node

#srun ./calc_eps.cplx.x chi0mat.h5 > ce0.out
#sleep 1
srun ./EpsInv.x chimat.h5 > EpsInv.out
