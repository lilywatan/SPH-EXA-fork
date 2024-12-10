#!/usr/bin/env bash
#SBATCH --job-name=disk_comb_100         # Job name    (default: sbatch)
#SBATCH --output=disk_comb_100.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_comb_100.err         # Error file  (default: slurm-%j.out)
#SBATCH --cpus-per-task=6         # Number of CPUs per task
#SBATCH --ntasks=8              # Number of tasks
#SBATCH --nodes=2   
#SBATCH --ntasks-per-socket=2
#SBATCH --cores-per-socket=12   
#SBATCH --mem-per-cpu=4G          # Memory per CPU
#SBATCH --time=01:00:00           # Wall clock time limit

export OMP_NUM_THREADS=6
export OMP_PLACES=cores
module load stack/.2024-06-silent  gcc/12.2.0
module load openmpi/4.1.6
module load hdf5/1.14.3

srun --cpus-per-task=6 build-release/main/src/sphexa/sphexa --init '../disk5.hdf5' --prop std-planet -s 100 -w 1 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o './output/runs/run_disk_comb_100.hdf5' --nthreads=$OMP_NUM_THREADS 

