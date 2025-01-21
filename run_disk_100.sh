#!/usr/bin/env bash
#SBATCH --job-name=disk_comb_100         # Job name    (default: sbatch)
#SBATCH --output=disk_comb_100.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_comb_100.err         # Error file  (default: slurm-%j.out)
#SBATCH --cpus-per-task=2         # Number of CPUs per task
#SBATCH --ntasks=2                # Number of tasks
#SBATCH --ntasks-per-node=1      # Number of tasks per node§
#SBATCH --mem-per-cpu=4G          # Memory per CPU
#SBATCH --time=01:00:00           # Wall clock time limit

module load openmpi
module load hdf5
pmix_info
srun --mpi=list

OMP_NUM_THREADS=8

srun --mpi=pmix_v5  build-release/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5.hdf5' --prop std-angmom -s 100 -w 10 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o '/home/lwatan/data/SPH-EXA-fork/output/runs/run_disk_comb_100.hdf5' --nthreads=$OMP_NUM_THREADS --mpi=pmix

