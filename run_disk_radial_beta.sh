#!/usr/bin/env bash
#SBATCH --job-name=disk_radial_1e6_beta       # Job name    (default: sbatch)
#SBATCH --output=disk_radial_1e6_beta_%j.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_radial_1e6_beta_%j.err         # Error file  (default: slurm-%j.out)
#SBATCH --cpus-per-task=32         # Number of CPUs per task
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --ntasks-per-node=1      # Number of tasks per node§
#SBATCH --mem-per-cpu=1G          # Memory per CPU
#SBATCH --time=80:00:00           # Wall clock time limit

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

srun build-radial/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5_beta.hdf5' --prop std-planet -s 1e6 -w 1000 -f m,c,x,y,z,rho,p,vx,vy,vz,h -o '/home/lwatan/scratch/run_disk_radial_1e6_beta.hdf5'

