#!/bin/bash
#SBATCH --job-name=disk_comb_1mil         # Job name    (default: sbatch)
#SBATCH --output=disk_comb_1mil.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_comb_1mil.err         # Error file  (default: slurm-%j.out)
#SBATCH --cpus-per-task=4         # Number of CPUs per task
#SBATCH --nodes=4
#SBATCH --mem-per-cpu=4G          # Memory per CPU
#SBATCH --time=02:00:00           # Wall clock time limit

module load openmpi
module load hdf5

OMP_NUM_THREADS=4
export OMP_NUM_THREADS

# ./main/src/sphexa/sphexa --init '/home/lwatan/data/cloud.h5' --prop std-angmom -s 3 -w 1 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o ../output/run_cloud_3_5.hdf5
mpirun -np 4 build/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5.hdf5' --prop std-angmom -s 1e6 -w 10 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o ../output/run_disk_comb_1mil.hdf5
