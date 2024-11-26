#!/usr/bin/env bash
#SBATCH --job-name=disk_mom_20        # Job name    (default: sbatch)
#SBATCH --output=disk_mom_20.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_mom_20_%j.err         # Error file  (default: slurm-%j.out)
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --ntasks-per-node=1      # Number of tasks per node§
#SBATCH --mem-per-cpu=4G          # Memory per CPU
#SBATCH --time=18:00:00           # Wall clock time limit

module load openmpi
module load hdf5

OMP_NUM_THREADS=16
export OMP_NUM_THREADS

srun build-release/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5.hdf5' --prop std-angmom -s 20000 -w 10 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o '/home/lwatan/data/SPH-EXA-fork/output/run_disk_mom_20.hdf5'

