#!/usr/bin/env bash
#SBATCH --job-name=disk_subgrid_beta_planet_star0_2         # Job name    (default: sbatch)
#SBATCH --output=disk_subgrid_beta_planet_star0_2.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_subgrid_beta_planet_star0_2.err         # Error file  (default: slurm-%j.err)
#SBATCH --cpus-per-task=32       # Number of CPUs per task
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --ntasks-per-node=1      # Number of tasks per node§
#SBATCH --mem-per-cpu=1G          # Memory per CPU
#SBATCH --time=10:00:00           # Wall clock time limit

module load openmpi/5.0.5
module load hdf5/1.12.2

# Set the library path for MPI
export LD_LIBRARY_PATH=/apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-5.0.5-nwkqils3tmp73uobzl5gaktusuwblarf/lib:$LD_LIBRARY_PATH

# Verify the library path and confirm libmpi.so.40 is accessible
if [ ! -f /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/openmpi-5.0.5-nwkqils3tmp73uobzl5gaktusuwblarf/lib/libmpi.so.40 ]; then
    echo "Error: libmpi.so.40 not found in the specified library path."
    exit 1
fi

echo "Library path set to: $LD_LIBRARY_PATH"
echo "libmpi.so.40 found and accessible."

OMP_NUM_THREADS=32
export OMP_NUM_THREADS
EXEC_PATH="/home/lwatan/data/SPH-EXA-fork/build-subgrid/main/src/sphexa/sphexa"

srun $EXEC_PATH --init '/home/lwatan/data/SPH-EXA-fork/subgrid_init_planet_star0.hdf5' --prop std-subgrid -s 1e6 -w 1000 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o '/home/lwatan/scratch/run_subgrid_beta_planet_star0_2.hdf5'

