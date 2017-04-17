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


def calc_rho(Sal, TempC, P):
    """ Calculate rho: Based on SOG code
    """
    
    # Calculate the square root of the salinities
    sqrSal = np.sqrt(Sal)

    # Calculate the density profile at the grid point depths
    # Pure water density at atmospheric pressure
    # (Bigg P.H., (1967) Br. J. Applied Physics 8 pp 521-537)
    R1 = ((((6.536332e-9 * TempC - 1.120083e-6) * TempC + 1.001685e-4)
           * TempC - 9.095290e-3) * TempC + 6.793952e-2) * TempC - 28.263737
    R2 = (((5.3875e-9 * TempC - 8.2467e-7) * TempC + 7.6438e-5)
          * TempC - 4.0899e-3) * TempC + 8.24493e-1
    R3 = (-1.6546e-6 * TempC + 1.0227e-4) * TempC - 5.72466e-3

    # International one-atmosphere equation of state of seawater
    SIG = (4.8314e-4 * Sal + R3 * sqrSal + R2) * Sal + R1

    # Specific volume at atmospheric pressure
    V350P = 1.0 / 1028.1063
    SVA   = -SIG * V350P / (1028.1063 + SIG)

    # Density anomoly at atmospheric pressure
    rho = 28.106331 - SVA / (V350P * (V350P + SVA)) + 1000
    
    return rho
