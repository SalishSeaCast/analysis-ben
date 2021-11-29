#!/bin/bash
#SBATCH --account=rrg-allen
#SBATCH --job-name="process_NEMO_drifters"
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-user=bmoorema@eoas.ubc.ca
#SBATCH --mail-type=ALL

#module load glost/0.3.1
module load python/3.9.6
module load netcdf/4.7.4
virtualenv --no-download ~/env/
source ~/env/bin/activate
python3 -m pip install --no-index --upgrade pip
python3 -m pip install --no-index numpy bottleneck netCDF4 pandas scipy xarray

echo "Starting run at: `date`"
srun python3 process_NEMO_drifters.py charnock 20160401 20160930 
#srun glost_launch process_NEMO_drifters.list
echo "Finished run at: `date`"
