#!/bin/bash
module load anaconda3
source activate disk-analysis

echo "set up done – running python scripts"

python subgrid_disk.py
echo "subgrid disk done"

#python mass_accretion_cluster.py
#echo "mass accretion done"

#python vertical_profile.py
#echo "vertical profile done"

#python kin_visc.py
#echo "kinematic viscosity done"

#python surface_density.py
#echo "surface density done"

