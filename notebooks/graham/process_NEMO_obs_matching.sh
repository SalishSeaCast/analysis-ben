#!/bin/bash
#SBATCH --account=rrg-allen
#SBATCH --job-name="process_NEMO_obs_matching"
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-user=bmoorema@eoas.ubc.ca
#SBATCH --mail-type=ALL

module load glost/0.3.1
module load python/3.9.6
module load netcdf/4.7.4
virtualenv --no-download ~/env/
source ~/env/bin/activate
python3 -m pip install --no-index --upgrade pip
python3 -m pip install --no-index numpy bottleneck netCDF4 pandas scipy xarray

echo "Starting run at: `date`"
srun glost_launch process_NEMO_obs_matching.list
echo "Finished run at: `date`"
