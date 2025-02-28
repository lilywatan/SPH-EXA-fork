#!/usr/bin/env bash
#SBATCH --job-name=disk_calibration_beta_planet        # Job name    (default: sbatch)
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

echo "28-02-2025: Running SPH-EXA with 1e6 particles, planet propagator, h = 1.2, star_inner = 0"
srun build-star0/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5_beta.hdf5' --prop std-planet -s 1e6 -w 1000 -f c,m,x,y,z,rho,vx,vy,vz,h,u,temp,alpha,du_m1,x_m1,y_m1,z_m1 -o '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_star0.hdf5'
echo "28-02-2025: Running SPH-EXA with 1e6 particles, planet propagator, h = 1.2, star_inner = 0, no removal limit h in accretion"
srun build-no-hlim/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5_beta.hdf5' --prop std-planet -s 1e6 -w 1000 -f c,m,x,y,z,rho,vx,vy,vz,h,u,temp,alpha,du_m1,x_m1,y_m1,z_m1 -o '/home/lwatan/scratch/run_disk_cal_1e6_beta_no_hlim.hdf5'
echo "28-02-2025: Running SPH-EXA with 1e6 particles, planet propagator, h = 5.0, star_inner = 0"
srun build-star0/main/src/sphexa/sphexa --init '/home/lwatan/data/disk5_beta_noh.hdf5' --prop std-planet -s 1e6 -w 1000 -f c,m,x,y,z,rho,vx,vy,vz,h,u,temp,alpha,du_m1,x_m1,y_m1,z_m1 -o '/home/lwatan/scratch/run_disk_cal_1e6_beta_h5.hdf5'
