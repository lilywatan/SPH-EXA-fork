#!/usr/bin/env bash
#SBATCH --job-name=disk_calibration_beta_planet         # Job name    (default: sbatch)
#SBATCH --output=disk_calibration_beta_planet.out        # Output file (default: slurm-%j.out)
#SBATCH --error=disk_calibration_beta_planet_%j.err         # Error file  (default: slurm-%j.err)
#SBATCH --cpus-per-task=32       # Number of CPUs per task
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

srun build-no-accretion/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5_beta.hdf5' --prop std-planet -s 1e6 -w 1000 -f m,x,y,z,rho,vx,vy,vz,h,u,temp,alpha,du_m1,x_m1,y_m1,z_m1 -o '/home/lwatan/data/SPH-EXA-fork/output/runs/run_disk_cal_1e6_beta_planet.hdf5'
# srun build-no-accretion/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5_beta.hdf5' --prop std-hydro -s 1e6 -w 1000 -f m,x,y,z,rho,vx,vy,vz,h,u,temp,alpha,du_m1,x_m1,y_m1,z_m1 -o '/home/lwatan/data/SPH-EXA-fork/output/runs/run_disk_cal_1e6_beta_hydro.hdf5'

