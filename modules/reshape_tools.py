import numpy as np
import xarray as xr
from salishsea_tools import viz_tools


def reshape_grid(tmask, deptht, gridY, gridX, index=0, dim=0, spacing=1):
    """Prepare the mask and grid for the selected timeseries slice, and reshape into 1 spatial dimension
    """
    
    tmask = tmask.take(index, axis=dim)[::spacing, ::spacing]
    ngrid = tmask.shape[0] * tmask.shape[1]
    tmask = tmask.reshape(ngrid)
    ngrid_water = tmask.sum()
    gridX = gridX.take(index, axis=dim)[::spacing, ::spacing].reshape(ngrid)[tmask]
    gridY = gridY.take(index, axis=dim)[::spacing, ::spacing].reshape(ngrid)[tmask]
    deptht = deptht.take(index, axis=dim)[::spacing, ::spacing].reshape(ngrid)[tmask]
    
    GRIDMASK = {'tmask': tmask, 'deptht': deptht, 'gridY': gridY,
                'gridX': gridX, 'ngrid': ngrid, 'ngrid_water': ngrid_water}
    
    return GRIDMASK


def reshape_data(data, GRIDMASK, dim, index=0, spacing=1, unstagger_dim=None):
    """
    """
    
    field = data.isel(**{dim: index})
    if unstagger_dim is not None:
        field = viz_tools.unstagger_xarray(field, unstagger_dim)
    field = field.values[:, ::spacing, ::spacing]
    field = field.reshape((-1, GRIDMASK['ngrid']))
    field_trim = np.zeros((field.shape[0], GRIDMASK['ngrid_water']))
    
    return field, field_trim
