import timeseries_extras
import xarray as xr

# Load mask and grid
mask = xr.open_dataset(
    '/data/bmoorema/MEOPAR/NEMO-forcing/grid/mesh_mask_downbyone2.nc'
)

# Timerange
timerange = ['2016 Jan 1 00:00', '2016 Dec 31 23:00']

# Loop through requests
dims = ['depth', 'depth', 'y', 'y', 'y']
indices = [0, 20, 450, 520, 680]
variables = ['votemper', 'vosaline', 'vozocrtx', 'vomecrty']
timeseries_extras.iterate_NEMO_timeseries(
    timerange, variables, mask, dims, indices,
)
